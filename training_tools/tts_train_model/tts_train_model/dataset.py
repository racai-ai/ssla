# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

if __name__ == "__main__":
    print "Hello World"

import numpy as np
from os import listdir
from os.path import isfile
from os.path import join
import scipy
from scipy.cluster.vq import kmeans
import scipy.io
import scipy.io.wavfile
import sys

class Dataset:
    def __init__(self, train_path, dev_path):
        self.train_list = []
        self.dev_list = []
        sys.stdout.write("scanning " + train_path + "...")
        onlyfiles = [f for f in listdir(train_path) if isfile(join(train_path, f))]
        for f in onlyfiles:
            fn = join (train_path, f)
            if fn.endswith(".wav"):
                self.train_list.append(fn.replace(".wav", ""))
        sr, wav = scipy.io.wavfile.read(self.train_list[0] + ".wav")
        self.sample_rate = sr
        sys.stdout.write(" found " + str(len(self.train_list)) + " training files\n")
        
        sys.stdout.write("scanning " + dev_path + "...")
        onlyfiles = [f for f in listdir(dev_path) if isfile(join(dev_path, f))]
        for f in onlyfiles:
            fn = join (dev_path, f)
            if fn.endswith(".wav"):
                self.dev_list.append(fn.replace(".wav", ""))
        sys.stdout.write(" found " + str(len(self.dev_list)) + " dev files\n")
        
        sys.stdout.write("Building feature maps...")
        self.phoneme2int = {}
        self.phoneme2int["#"] = 0
        self.duration2int = {}
        self.context2int = []
        self.total_states = 0
        for fname in self.train_list:
            pi_list = self.read_lab(fname)
            for pi in pi_list:
                self.total_states += (pi.stop-pi.start) / 5
                phoneme = pi.phoneme
                if phoneme not in self.phoneme2int:
                    self.phoneme2int[phoneme] = len(self.phoneme2int)
                duration = pi.duration
                if duration not in self.duration2int:
                    self.duration2int[duration] = len(self.duration2int)
                if len(self.context2int) == 0:
                    print "Creating context of size " + str(len(pi.context))
                    for zz in range (len(pi.context)):
                        self.context2int.append({})
                for zz in range(len(pi.context)):
                    cont = pi.context[zz]
                    if cont not in self.context2int[zz]:
                        self.context2int[zz][cont] = len(self.context2int[zz])
                
        sys.stdout.write("done\n")
        sys.stdout.write("Statistics:\n")
        
        sys.stdout.write(" \t" + str(len(self.phoneme2int)) + " unique phonemes\n")
        sys.stdout.write(" \t" + str(len(self.duration2int)) + " unique durations\n")
        sys.stdout.write(" \t" + str (len(self.context2int)) + " feature maps\n")
        sys.stdout.write(" \t" + str(self.total_states) + " total states\n")
        self.state2int = None
        #print self.context2int
        #for zz in range(len(self.context2int)):
        #    print len(self.context2int[zz])
        
    def ulaw_encode(self, data):
        out_discreete = []
        out_continous = []
        for zz in range(len(data)):
            f = float(data[zz]) / 32768
            sign = np.sign(f)
            encoded = sign * np.log(1.0 + 255.0 * np.abs(f)) / np.log(1.0 + 255.0)
            encoded_d = int((encoded + 1) * 127)
            out_discreete.append(encoded_d)
            out_continous.append(encoded)
            
        return [out_discreete, out_continous]
    def ulaw_decode(self, data, discreete=True):
        out = []
        for zz in range (len(data)):
            if discreete:
                f = float(data[zz]) / 128-1.0
            else:
                f = data[zz]
            sign = np.sign(f)
            decoded = sign * (1.0 / 255.0) * (pow(1.0 + 255, abs(f))-1.0)
            decoded = int(decoded * 32768)
            out.append(decoded)
        return out
    
    def read_wave(self, filename):
        sr, wav = scipy.io.wavfile.read(filename + ".wav")
        out = self.ulaw_encode(wav)
        #rec=self.ulaw_decode(out)
        #data2 = np.asarray(rec, dtype=np.int16)
        #scipy.io.wavfile.write("test.wav", 16000, data2)
        return out
    
    def read_lab(self, filename):
        out = []
        with open (filename + ".lab") as f:
            lines = f.readlines()
            for line in lines:
                line = line.replace("\n", "")
                parts = line.split(" ")
                start = int(parts[0]) / 10000
                stop = int(parts[1]) / 10000
                pp = parts[2].split(":")
                phon = pp[0]
                context = parts[2][parts[2].find(":") + 2:]
                phon = phon.split("-")[1]
                phon = phon.split("+")[0]
                pi = PhoneInfo(phon, context, start, stop)
                out.append(pi)
        return out
    
    def make_dataframes(self, audio, frame_len_ms):
        frames = []
        frame_len_samp = self.sample_rate / 1000 * frame_len_ms
        num_frames = len(audio) / frame_len_samp
        if len(audio) % frame_len_samp != 0:
            num_frames += 1
        index = 0
       
        for frame in range (num_frames):
            frm = []
            for zz in range (frame_len_samp):
                if index < len(audio):
                    frm += [audio[index]]
                else:
                    frm += [0]
                index += 1
            frames.append(frm)
            
        return frames
    
    def squared_distance(self, p1, p2):
        dist = 0
        for zz in range(len(p1)):
            dist += (p1[zz]-p2[zz]) ** 2
        return np.sqrt(dist)
    
    def spectrum2int(self, spec):
        
        best_index = 0
        best_score = self.squared_distance(spec, self.state2int[0])
        for zz in range(1, self.num_clusters):
            score = self.squared_distance(spec, self.state2int[zz])
            if score < best_score:
                best_score = score
                best_index = zz
        return best_index
            
class PhoneInfo:
    context2int = {}
    def __init__(self, phoneme, context, start, stop):
        self.phoneme = phoneme
        self.context = context.split("/")
        self.start = start
        self.stop = stop
        self.duration = (stop-start)
            