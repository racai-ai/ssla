# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

from cnn_wrapper import CNN
import dynet as dy
import math
import numpy as np
import random
import sys

class AttentionNetwork:
    def __init__(self, dataset, config):
        self.clip = 5 * dataset.sample_rate # set to zero for full training
        self.sample_rate = dataset.sample_rate
        self.config = config
        self.dataset = dataset
        
        self.model = dy.Model()
        self.trainer = dy.AdamTrainer(self.model)
        
        #lookups
        self.phoneme_lookup = self.model.add_lookup_parameters((len(dataset.phoneme2int), config.phone_embeddings_size))
        self.context_lookup = [self.model.add_lookup_parameters((len(ctx), config.context_embeddings_size)) for ctx in dataset.context2int]
        
        #encoder
        inp_sz = config.phone_embeddings_size * 5 #+ len(dataset.context2int) * config.context_embeddings_size
        self.encoder_fw = [dy.LSTMBuilder(1, inp_sz, config.encoder_size, self.model)]
        self.encoder_bw = [dy.LSTMBuilder(1, inp_sz, config.encoder_size, self.model)]
        [self.encoder_fw.append(dy.LSTMBuilder(1, config.encoder_size * 2, config.encoder_size, self.model)) for _ in range(config.encoder_layers-1)]
        [self.encoder_bw.append(dy.LSTMBuilder(1, config.encoder_size * 2, config.encoder_size, self.model)) for _ in range(config.encoder_layers-1)]
        
        #receptive network
        receptive_input = config.receptive_input
        self.receptive_w = []
        self.receptive_b = []
        for size in config.receptive_layers:
            receptive_output = size
            self.receptive_w.append(self.model.add_parameters((receptive_output, receptive_input)))
            self.receptive_b.append(self.model.add_parameters((receptive_output)))
            receptive_input = receptive_output
            
        #receptive network
        attention_input = config.receptive_input
        self.attention_w = []
        self.attention_b = []
        for size in config.attention_layers:
            attention_output = size
            self.attention_w.append(self.model.add_parameters((attention_output, attention_input)))
            self.attention_b.append(self.model.add_parameters((attention_output)))
            attention_input = attention_output
            
        #decoder
        self.decoder = dy.GRUBuilder(config.decoder_layers, config.encoder_size * 2 + config.receptive_layers[-1], config.decoder_size, self.model)
        
        #attention
        self.att_w1 = self.model.add_parameters((config.att_proj_size, config.encoder_size * 2))
        self.att_w2 = self.model.add_parameters((config.att_proj_size, config.decoder_size))
        self.att_w3 = self.model.add_parameters((config.att_proj_size, config.att_lsa_filters))
        self.att_w4 = self.model.add_parameters((config.att_proj_size, config.attention_layers[-1]))
        self.att_v = self.model.add_parameters((1, config.att_proj_size))
        self.cnn_attention = CNN(self.model)
        self.cnn_attention.add_layer_conv(config.att_lsa_input_size, 1, 1, 1, config.att_lsa_filters, same=True)
        
        #output
        
        
        presoftmax_input = config.decoder_size + config.receptive_layers[-1] + config.sample_trail_size
        self.presoftmax_w = []
        self.presoftmax_b = []
        for size in config.presoftmax_layers:
            presoftmax_output = size
            self.presoftmax_w.append(self.model.add_parameters((presoftmax_output, presoftmax_input)))
            self.presoftmax_b.append(self.model.add_parameters((presoftmax_output)))
            presoftmax_input = presoftmax_output
        
        self.softmax_w = self.model.add_parameters((257, presoftmax_input))
        self.softmax_b = self.model.add_parameters((257))
        
        
        
    def argmax(self, data):
        max = 0
        for zz in range(1, len(data)):
            if data[zz] > data[max]:
                max = zz
                
        return max
    
    def _attend(self, input_vectors, state, prev_att, prev_att_expr, receptive, compute_attention):
        if compute_attention or prev_att_expr is None:
            w1 = self.att_w1.expr()
            w2 = self.att_w2.expr()
            w3 = self.att_w3.expr()
            w4 = self.att_w4.expr()
            v = self.att_v.expr()
            attention_weights = []
            att_cnn = self.cnn_attention.apply(dy.reshape(prev_att, (len(input_vectors), 1)))
            att_cnn = dy.reshape(att_cnn, (len(input_vectors), self.config.att_lsa_filters))


            w2dt = w2 * state.h()[-1]
            w4dt = w4 * receptive
            for cnn, input_vector in zip (att_cnn, input_vectors):
                attention_weight = v * dy.tanh(w1 * input_vector + w2dt + w3 * cnn + w4dt)
                attention_weights.append(attention_weight)

            attention_weights = dy.softmax(dy.concatenate(attention_weights))
            #print attention_weights.value()
        else:
            attention_weights = prev_att_expr
        
        output_vectors = dy.esum([vector * attention_weight for vector, attention_weight in zip(input_vectors, attention_weights)])
        return output_vectors, attention_weights

    def _step(self, prev_samples, encoder_output, decoder_state, prev_att, prev_att_expr, runtime, compute_attention):
        if prev_att is None:
            prev_att = dy.inputVector([0] * len(encoder_output))
        else:
            prev_att = dy.inputVector(prev_att) #this truncates backpropagation - don't know if it is ok to do that
        
        
        
        #input from receptive network
        while len(prev_samples) < self.config.receptive_input:
            prev_samples = [0] + prev_samples
        input_vect2 = dy.inputVector(prev_samples[-self.config.receptive_input:])
        input_vect3 = dy.inputVector(prev_samples[-self.config.sample_trail_size:])
        for w, b in zip (self.receptive_w, self.receptive_b):
            input_vect2 = dy.rectify(w.expr() * input_vect2 + b.expr())
            if not runtime:
                input_vect2 = dy.dropout(input_vect2, self.config.receptive_dropout)
        #input from encoder
        if compute_attention or prev_att_expr is None:
            att_vect = dy.inputVector(prev_samples[-self.config.receptive_input:])
            for w, b in zip (self.attention_w, self.attention_b):
                att_vect = dy.rectify(w.expr() * att_vect + b.expr())
        else:
            att_vect=None
                    
        input_vect1, prev_att = self._attend(encoder_output, decoder_state, prev_att, prev_att_expr, att_vect, compute_attention)
        
        decoder_state = decoder_state.add_input(dy.concatenate([input_vect1, input_vect2]))
        presoftmax = dy.concatenate([decoder_state.output(), input_vect2, input_vect3])
        for w, b in zip (self.presoftmax_w, self.presoftmax_b):
            presoftmax = dy.rectify(w.expr() * presoftmax + b.expr())
            
        softmax = dy.softmax(self.softmax_w.expr() * presoftmax + self.softmax_b.expr())
        return softmax, decoder_state, prev_att.value(), prev_att
    
    
    def _make_input(self, labels):
        x_list = []
        for iLabel in range (len(labels)):
            ppP = "#"
            pP = "#"
            cP = labels[iLabel].phoneme
            nP = "#"
            nnP = "#"
            if iLabel-1 >= 0:
                pP = labels[iLabel-1].phoneme
            if iLabel-2 >= 0:
                ppP = labels[iLabel-2].phoneme
            if iLabel + 1 < len(labels):
                nP = labels[iLabel + 1].phoneme
            if iLabel + 2 < len(labels):
                nnP = labels[iLabel + 2].phoneme
            x1 = self.phoneme_lookup[self.dataset.phoneme2int[ppP]]
            x2 = self.phoneme_lookup[self.dataset.phoneme2int[pP]]
            x3 = self.phoneme_lookup[self.dataset.phoneme2int[cP]]
            x4 = self.phoneme_lookup[self.dataset.phoneme2int[nP]]
            x5 = self.phoneme_lookup[self.dataset.phoneme2int[nnP]]
            x = dy.concatenate([x1, x2, x3, x4, x5])#, self._get_x_context(labels[iLabel])])
            x_list.append(x)
        return x_list
    
    def _pick_sample(self, probs):
        #print np.sum(probs)
        probs=np.asarray(probs[:-1])
        probs=probs/np.sum(probs)
        return np.random.choice(np.arange(256), p=probs)
    
    def _predict(self, labels, start_sample, stop_sample, samples_cont, prev_att, last_decoder_state, runtime=False):
        #compute encoder
        x_list = self._make_input(labels)
        for fw, bw in zip (self.encoder_fw, self.encoder_bw):
            x_fw = fw.initial_state().transduce(x_list)
            x_bw = list(reversed(bw.initial_state().transduce(reversed(x_list))))
            x_list = [dy.concatenate([x1, x2]) for x1, x2 in zip (x_fw, x_bw)]
        encoder_out = x_list 
        pred_samples = []
        
        if last_decoder_state is None:
            rnn = self.decoder.initial_state().add_input(dy.vecInput(self.config.encoder_size * 2 + self.config.receptive_layers[-1]))
        else:
            s=[dy.inputVector(state) for state in last_decoder_state]
            rnn=self.decoder.initial_state().set_s(s)
        att = prev_att
        att_expr=None#dy.inputVector(prev_att)
        for iSample in range (start_sample, stop_sample):
            prev_samples = samples_cont[:iSample]
            if iSample-start_sample % self.config.compute_attention_samples == 0:
                compute_attention = True
            else:
                compute_attention = False
            new_sample, rnn, att, att_expr = self._step(prev_samples, encoder_out, rnn, att, att_expr, runtime, compute_attention)
            pred_samples.append(new_sample)
            if runtime:
                sample = self._pick_sample(new_sample.value())#self._argmax(new_sample.value())
                samples_cont.append(float(sample) / 127 - 1.0)
            
        return pred_samples, att, samples_cont, rnn
    
    def _get_x_context(self, label):
        zeroVec = dy.vecInput(self.config.context_embeddings_size)
        context_list = []
        for i in range (len(label.context)):
            ctx = label.context[i]
            if ctx in self.dataset.context2int[i]:
                context_list.append(self.context_lookup[i][self.dataset.context2int[i][ctx]])
            else:
                context_list.append(zeroVec)
        return dy.concatenate(context_list)

    def _get_x_phoneme(self, label):
        return self.phoneme_lookup[self.dataset.phoneme2int[label.phoneme]]
    
    def synthesize(self, labels, max_samples):
        num_batches = (max_samples) / self.config.batch_size
        if (max_samples) % self.config.batch_size != 0:
            num_batches += 1
        
        att = None
       
        
        last_proc = 0
        f_cont = []
        f_disc = []
        last_decoder_state=None
        for iBatch in range (num_batches):
            proc = iBatch * 100 / num_batches
            while last_proc + 10 < proc:
                last_proc += 10
                sys.stdout.write(" " + str(last_proc))
                sys.stdout.flush()
            dy.renew_cg()
            start_sample = iBatch * self.config.batch_size
            stop_sample = start_sample + self.config.batch_size
            if stop_sample > max_samples:
                stop_sample = max_samples
            samples_disc, att, samples_cont, last_decoder_state = self._predict(labels, start_sample, stop_sample, f_cont, att, last_decoder_state, True)
            last_decoder_state=[s.value() for s in last_decoder_state.s()]
            for iSample in range (stop_sample -start_sample):
                #sample = self.argmax(pred_samples[iSample].value())
                #f_disc.append(sample)
                #norm_sample = float(sample) / 127-1.0
                norm_sample = samples_cont[iSample]
                f_cont.append(norm_sample)
                
        return f_cont
    
    def learn(self, labels, f_cont, f_disc):
        #trunchiem pauzele:
        f_cont = f_cont[labels[0].stop * self.dataset.sample_rate / 1000-80:labels[-1].start * self.dataset.sample_rate / 1000 + 80]
        f_disc = f_disc[labels[0].stop * self.dataset.sample_rate / 1000-80:labels[-1].start * self.dataset.sample_rate / 1000 + 80]
        
        clipped = False
        if self.clip != 0 and len(f_cont) > self.clip:
            clipped = True
            f_cont = f_cont[:self.clip]
            f_disc = f_disc[:self.clip]
        
        num_batches = (len(f_cont) + 1) / self.config.batch_size
        if (len(f_cont) + 1) % self.config.batch_size != 0:
            num_batches += 1
        
        att = None
        total_loss = 0
       
        
        last_proc = 0
        last_decoder_state=None
        for iBatch in range (num_batches):
            proc = iBatch * 100 / num_batches
            while last_proc + 10 < proc:
                last_proc += 10
                sys.stdout.write(" " + str(last_proc))
                sys.stdout.flush()
            dy.renew_cg()
            start_sample = iBatch * self.config.batch_size
            stop_sample = start_sample + self.config.batch_size
            if stop_sample > len(f_disc) + 1:
                stop_sample = len(f_disc) + 1
            pred_samples, att, samples_cont, rnn = self._predict(labels, start_sample, stop_sample, f_cont, att, last_decoder_state, False)
            last_decoder_state=[s.value() for s in rnn.s()]
            losses = []
            for iSample in range (stop_sample -start_sample):
                if iSample + start_sample != len(f_cont):
                    losses.append(-dy.log(dy.pick(pred_samples[iSample], f_disc[iSample + start_sample])))
                elif not clipped:
                    losses.append(-dy.log(dy.pick(pred_samples[-1], 256))) #special end of sequence
            loss = dy.esum(losses)
            total_loss += loss.value()
            loss.backward()
            self.trainer.update()
            p_one = 1
            if clipped:
                p_one = 0
            
        return total_loss, len(f_disc) + p_one
    
    def store_model(self, path):
        sys.stdout.write("Storing " + path + "\n")
        self.model.save(path)
        
    def load_model(self, path):
        sys.stdout.write("Loading vocoder from " + path + "\n")
        self.model.populate(path)
    
