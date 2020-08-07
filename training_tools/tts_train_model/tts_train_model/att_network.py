# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import dynet as dy
import sys
import numpy as np

class AttentionTTS:
    def __init__(self, ds, config):
        
        self.model = dy.Model()
        self.trainer = dy.AdamTrainer(self.model, alpha=0.002)
        self.config = config
        self.dataset = ds
        self.sample_embedding_size=64
        #wave discreete encoding
        self.sample_lookup=self.model.add_lookup_parameters((256,self.sample_embedding_size))
        
        
        #label transduction and encoding
        self.phoneme_lookup = self.model.add_lookup_parameters((len(ds.phoneme2int), 32))
        self.context_lookup = []
        for zz in range(len(ds.context2int)):
            self.context_lookup.append(self.model.add_lookup_parameters((len(ds.context2int[zz]), 8)))
            
        self.input_size = 32 + len(self.context_lookup) * 8
        
        self.encoder_fw = [dy.LSTMBuilder(1, self.input_size, config.encoder_size, self.model)]
        self.encoder_bw = [dy.LSTMBuilder(1, self.input_size, config.encoder_size, self.model)]
        
        for zz in range(1, config.encoder_layers):
            self.encoder_fw.append(dy.LSTMBuilder(1, config.encoder_size * 2, config.encoder_size, self.model))
            self.encoder_bw.append(dy.LSTMBuilder(1, config.encoder_size * 2, config.encoder_size, self.model))
            
        #state decoding and synthesis
        #self.decoder = [dy.GRUBuilder(1, config.encoder_size * 2, config.decoder_size, self.model)]
        #for zz in range(1, config.decoder_layers):
        #    self.decoder.append(dy.GRUBuilder(1, config.decoder_size, config.decoder_size, self.model))
        self.decoder = dy.LSTMBuilder(config.decoder_layers, config.encoder_size * 2, config.decoder_size, self.model)
        
        self.presoftmaxW = self.model.add_parameters((256, config.decoder_size))
        self.presoftmaxB = self.model.add_parameters((256))
        
        self.softmaxW = self.model.add_parameters((256, 256))
        self.softmaxB = self.model.add_parameters((256))
        
        #easy peasy attention
        self.attention_window=3#
        self.att_w1 = self.model.add_parameters((config.encoder_size * 2, config.encoder_size * 2))
        self.att_w2 = self.model.add_parameters((config.encoder_size * 2, config.decoder_size))
        self.att_v = self.model.add_parameters((1, config.encoder_size * 2))
        
    
    def get_x(self, pi):
        zero_vec = dy.vecInput(8)
        zero_vec.set([0] * 8)
        concat = []
        if pi.phoneme in self.dataset.phoneme2int:
            phoneme_x = self.phoneme_lookup[self.dataset.phoneme2int[pi.phoneme]]
            #print pi.phoneme
        concat.append(phoneme_x)
        for zz in range(len(pi.context)):
            c_x = zero_vec
            if (pi.context[zz] in self.dataset.context2int[zz]):
                c_x = self.context_lookup[zz][self.dataset.context2int[zz][pi.context[zz]]]
            
            concat.append(c_x)
        x = dy.concatenate(concat)
        return x
    
    def _gaussian(self, x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    
    def _attend(self, input_vectors, state):
        w1 = self.att_w1.expr()
        w2 = self.att_w2.expr()
        v = self.att_v.expr()
        attention_weights = []

        w2dt = w2 * state.h()[-1]
        for input_vector in input_vectors:
            attention_weight = v * dy.tanh(w1 * input_vector + w2dt)
            attention_weights.append(attention_weight)
        attention_weights = dy.softmax(dy.concatenate(attention_weights))
        pos=self.argmax(attention_weights.value())
        #print pos
        att_inp=[]
        for x in range(pos-self.attention_window, pos+self.attention_window+1):
            gaussian_value=self._gaussian(x, pos, self.attention_window)
            #print gaussian_value
            if x>=0 and x<len(input_vectors):
                vector=input_vectors[x]
            else:
                vector=dy.vecInput(self.config.encoder_size*2)
            att_inp.append(vector*gaussian_value)
        #output_vectors = dy.esum([vector * attention_weight for vector, attention_weight in zip(input_vectors, attention_weights)])
        output_vectors=dy.esum(att_inp)
        return output_vectors
    
    def argmax(self, data):
        max = 0
        for zz in range(1, len(data)):
            if data[zz] > data[max]:
                max = zz
                
        return max
    
    def _encode(self, label_list):
        x_list = []

        for label in label_list:
            x = self.get_x(label)
            x_list.append(x)
        
        for zz in range(self.config.encoder_layers):
            lstm_fw = self.encoder_fw[zz].initial_state()
            lstm_bw = self.encoder_bw[zz].initial_state()
            l_fw = lstm_fw.transduce(x_list)
            l_bw = list(reversed(lstm_bw.transduce(reversed(x_list))))
            
            x_list = [dy.concatenate ([fw, bw]) for fw, bw in zip(l_fw, l_bw)]
        
        return x_list
            
    def _decode(self, input_list, stored_state, batch_size, last_output, gt_samples=None):
        
        #first layer has to attend input
        state = []
        rnn = self.decoder.initial_state()
        #print stored_state
        if stored_state!=None:
            rnn=rnn.set_s([dy.inputVector(state) for state in stored_state])
        else:
            rnn=rnn.add_input(dy.vecInput(self.config.encoder_size*2))
        
        num_predictions = batch_size * self.config.win_size * self.dataset.sample_rate / 1000
        soft_out = []
        x=None
        #last_out=last_output
        
        
        for zz in range(num_predictions):
            inp = self._attend(input_list, rnn)
            x_inp=inp
            rnn = rnn.add_input(x_inp)
            #x=dy.concatenate([rnn.output(), self.sample_lookup[last_out]])
            x=rnn.output()
            psm=dy.tanh(self.presoftmaxW.expr() * x + self.presoftmaxB.expr())
            soft_out.append(dy.softmax(self.softmaxW.expr() * psm + self.softmaxB.expr()))
#            if gt_samples!=None and zz<len(gt_samples):
#                last_out=gt_samples[zz]
#            else:
#                x_out=self.argmax(soft_out[-1].value())
#                last_out=x_out
#                #print last_out
        return soft_out, rnn
        
    
    def learn_full(self, frames_disc, frames_cont, continous_signal, labels):
        dy.renew_cg()
        t_list = self._encode(labels)
        
        stored_initial_state = None
        num_batches = (len(frames_disc) / self.config.batch_size) + 1
        losses = []
        total_loss = 0
        num_frames_per_batch = self.config.batch_size
        last_proc = 0
        min=len(frames_disc)
        last_out=127
        total_frames=0
        #start pause
        for iBatch in range(num_batches):
            
            proc = iBatch * 100 / num_batches
            if proc % 5 == 0 and proc != last_proc:
                sys.stdout.write(" " + str(proc))
                sys.stdout.flush()
                last_proc = proc
            
            start = iBatch * num_frames_per_batch
            stop = start + num_frames_per_batch
            if stop >= min:
                stop = min
                
            gt_values=[]
            for iFrame in range(start,stop):
                f_disc=frames_disc[iFrame]
                for pp in range(len(f_disc)):
                    gt_values.append(f_disc[pp])
                    
            if len(gt_values)>0:
                soft_out, final_state = self._decode(t_list, stored_initial_state, self.config.batch_size, last_out, gt_values)
                last_out=gt_values[-1]

                index = 0
                for iFrame in range(start, stop):
                    #do not propagate loss for start and end silence
                    if iFrame>labels[0].stop/5 and iFrame<labels[-1].start/5:
                        f_disc = frames_disc[iFrame]
                        for tt in range (len(f_disc)):
                            losses.append(-dy.log(dy.pick(soft_out[index], f_disc[tt])))
                            index += 1
                    else:
                        index+=len(frames_disc[iFrame])

                if len(losses) != 0:
                    total_frames+=len(losses)
                    lss = dy.esum(losses)
                    total_loss += lss.value()
                    lss.backward()
                    self.trainer.update()
                    losses = []
                    #salvam ultima stare
                    #print final_state.h()
                    stored_initial_state = [state.value(recalculate=True) for state in final_state.s()]
                    dy.renew_cg()
                    t_list = self._encode(labels)
        
        return total_loss, total_frames
    
    def predict(self, label_list):
        dy.renew_cg()
        t_list = self._encode(label_list)
        
        #stored_initial_state = [[0] * self.config.decoder_size] * self.config.decoder_layers
        wav_out_network = []
        num_frames=(label_list[-1].stop-label_list[0].start)/5
        soft_out, gru_state = self._decode(t_list, None, num_frames, 127)
        for tt in range (len(soft_out)):
            wav_out_network.append(self.argmax(soft_out[tt].value()))
        
        return wav_out_network