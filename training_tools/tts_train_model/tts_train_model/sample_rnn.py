import dynet as dy
import numpy as np
import sys


class SampleRNN:
    def __init__(self, config):
        self.config = config
        self.model = dy.Model()
        self.trainer = dy.AdamTrainer(self.model)
        self.rnn = []
        self.rnn_upsample_w = []
        self.rnn_linear_w = [None]
        rnn_input_size = config.FS[0]
        first = True
        for ls, fs, ups in zip(config.rnn_layers, config.FS, config.upsample):
            self.rnn.append(dy.GRUBuilder(1, rnn_input_size, ls, self.model))
            self.rnn_upsample_w.append([self.model.add_parameters((ls, ls))] * ups)
            if first:
                first = False
            else:
                self.rnn_linear_w.append(self.model.add_parameters((ls, fs)))
            rnn_input_size = ls

        layer_is = rnn_input_size
        self.mlp_w = []
        self.mlp_b = []
        for layer_os in config.mlp:
            self.mlp_w.append(self.model.add_parameters((layer_os, layer_is)))
            self.mlp_b.append(self.model.add_parameters((layer_os)))
            layer_is = layer_os
        self.mlp_w.append(self.model.add_parameters((256, layer_is)))
        self.mlp_b.append(self.model.add_parameters((256)))

    def learn(self, f_cont, f_disc, batch_size=512):
        num_batches = len(f_cont) / batch_size
        if len(f_cont) / batch_size != 0:
            num_batches += 1

        last_proc = 0
        rnn_states = None
        total_loss = 0
        last_proc=0
        for batch in xrange(num_batches):

            proc = batch * 100 / num_batches
            while last_proc + 10 < proc:
                last_proc += 10
                sys.stdout.write(" " + str(last_proc))
                sys.stdout.flush()

            dy.renew_cg()
            start_sample = batch * batch_size
            stop_sample = (batch + 1) * batch_size
            if stop_sample >= len(f_cont) - 1:
                stop_sample = len(f_cont) - 2
            # for sample in range (start_sample, stop_sample):
            softmax_list, sample_list, rnn_states = self._predict(f_cont=f_cont[start_sample:],
                                                                  f_disc=f_disc[start_sample:],
                                                                  rnn_states=rnn_states, num_samples=batch_size,
                                                                  runtime=False)
            losses = []
            for target, output in zip(f_disc[start_sample + 1:stop_sample + 1], softmax_list):
                losses.append(-dy.log(dy.pick(output, target)))
            loss = dy.esum(losses)
            total_loss += loss.value()
            loss.backward()
            self.trainer.update()

        return total_loss / len(f_disc)

    def _pick_sample(self, probs):
        # print np.sum(probs)
        # probs=np.asarray(probs[:-1])
        probs = probs / np.sum(probs)
        return np.random.choice(np.arange(256), p=probs)

    def _predict(self, f_cont=None, f_disc=None, rnn_states=None, num_samples=512, runtime=True):
        softmax_list = []
        sample_list = []
        if rnn_states is None:
            tiers = [rnn.initial_state() for rnn in self.rnn]
        else:
            tiers = []
            for rnn, state in zip(self.rnn, rnn_states):
                rnn = rnn.initial_state()
                rnn = rnn.set_s([dy.inputVector(state[0])])
                tiers.append(rnn)

        if f_cont is None:
            f_cont = [0]
            f_disc = [127]

        ct_tiers = [None] * (len(self.rnn) + 1)
        for iSample in xrange(num_samples):
            if runtime:
                trail_cont = f_cont
                trail_disc = f_disc
            else:
                trail_cont = f_cont[:iSample]
                trail_disc = f_disc[:iSample]
            if len(trail_cont) < self.config.maxFS:
                while len(trail_cont) < self.config.maxFS:
                    trail_cont.insert(0, 0)
                    trail_disc.insert(0, 127)

            # tiers
            for fs, w, tier_index in zip(self.config.FS, self.rnn_linear_w, xrange(len(self.rnn))):
                tier = tiers[tier_index]
                # make the input
                fs_input = dy.inputVector(trail_cont[-fs:])
                if w is not None:
                    fs_input = w.expr() * fs_input

                ct = ct_tiers[tier_index]
                if ct != None:
                    fs_input += ct[iSample % fs]

                tier = tier.add_input(fs_input)
                tiers[tier_index] = tier
                if iSample % fs == 0:  # update linear projections for next tier
                    ct = [w_ups.expr() * tier.output() for w_ups in self.rnn_upsample_w[tier_index]]
                    ct_tiers[tier_index + 1] = ct
            # postnetwork
            fs_input = dy.inputVector(trail_cont[-self.config.FS[-1]:])
            fs_input = self.rnn_linear_w[-1].expr() * fs_input
            ct_input = ct_tiers[-1][iSample % self.config.FS[-1]]
            input = fs_input + ct_input
            for iLayer in xrange(len(self.mlp_w)):
                w = self.mlp_w[iLayer]
                b = self.mlp_b[iLayer]
                input = w.expr() * input + b.expr()
                if iLayer == len(self.mlp_w) - 1:
                    input = dy.softmax(input)
                else:
                    input = dy.tanh(input)
            softmax_list.append(input)
            if runtime:
                sample_list.append(self._pick_sample(softmax_list[-1]))
                f_disc.append(sample_list[-1])
                f_cont.append(float(sample_list[-1]) / 127.0 - 1.0)

        rnn_states = [[tier.s()[0].value()] for tier in tiers]
        #print rnn_states
        return softmax_list, sample_list, rnn_states

    def synthesize(self, num_frames):
        return 0


class SampleRNNConfig:
    def __init__(self, filename=None):
        self.rnn_layers = [1024, 1024]
        self.FS = [8, 2, 2]
        self.maxFS = 8
        self.upsample = [8, 2]  # de recalculat
        self.mlp = [1024, 1024, 256]
