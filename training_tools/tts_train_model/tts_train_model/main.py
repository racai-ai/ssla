# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
#import dynet_config
#dynet_config.set(mem=1024,random_seed=9, weight_decay=1e-6)
#dynet_config.set_gpu()

import dynet_config
import sys

dynet_config.set(mem=4 * 512, random_seed=9, autobatch=True)
dynet_config.set_gpu()
import time

from dataset import Dataset
import dynet as dy
import numpy as np
import scipy.io.wavfile
import sys
import os


from config import Config
from vocoder import WorldVocoder
from tts_network import AttentionNetwork
from sample_rnn import SampleRNN, SampleRNNConfig

def display_help():
    print ("Neural TTS Model Trainer version 0.9 beta.")
    print ("Usage:")
    print ("\t--train <train folder> <dev folder> <model output base> <num itt no improve>")
    
    
def generate_raw_audio(network, ds):
    dy.renew_cg()
    lab_list = ds.read_lab(ds.dev_list[0])
    out = synth(lab_list, network,)#generate 5 seconds of speech
    rec = ds.ulaw_decode(out)
    data2 = np.asarray(rec, dtype=np.int16)
    scipy.io.wavfile.write("test.wav", network.sample_rate, data2)
    
def render_spec(spec, fname):
    from PIL import Image
    img = Image.new('RGB', (spec.shape[0], spec.shape[1]), "black")
    pixels = img.load() # create the pixel map
    min_val = spec[0][0]
    max_val = spec[0][0]
    for i in range(img.size[0]):
        for j in range(img.size[1]):
            if spec[i][j] > max_val:
                max_val = spec[i][j]
            if spec[i][j] < min_val:
                min_val = spec[i][j]

    for i in range(img.size[0]):
        for j in range(img.size[1]):
            color = (spec[i][j]-min_val) / (max_val-min_val + 1)
            c = int(color * 255) & 0xff
            pixels[i, img.size[1]-j-1] = (c, c, c)
    out_file = open(fname, 'wb')
    img.save(out_file, "PNG")
    out_file.flush()
    out_file.close()
    
def recode(vocoder, ds, output_base):
    sys.stdout.write("Recoding development data")
    start = time.time()
    directory = output_base + "/dev-recoded"
    errors = 0
    total = 0
    if not os.path.exists(directory):
        os.makedirs(directory)
    for zz in range(len(ds.dev_list)):
        fname = ds.dev_list[zz]
        [disc, cont] = ds.read_wave(fname)
        #draw power spectrum
        from PIL import Image
        fn = fname[fname.rfind('/'):]
        spec, aper, f0 = vocoder.extract_spectrum(cont)
        render_spec(spec, directory + fn + "_spec.png")
        render_spec(aper, directory + fn + "_aper.png")
        
        for i in range (aper.shape[0]):
            f0[i]=200
            #for j in range (aper.shape[1]):
            #    aper[i][j]=0
        
        wav = vocoder.synthesize(spec, aper, f0)
        wav_decoded = np.asarray(wav * 32767, dtype=np.int16)
        
        scipy.io.wavfile.write(directory + fn + ".wav", vocoder.sample_rate, wav_decoded)
        
        
    stop = time.time()
    total += 1
    sys.stdout.write(" acc=" + str((1.0-float(errors) / total)) + " execution time=" + str(stop-start) + "\n")
    return (1.0-float(errors) / total)




def synthesize(network, ds, output_base):
    sys.stdout.write("===Recoding development data===\n")
    start = time.time()
    directory = output_base + "/dev-synth"
    errors = 0
    total = 0
    if not os.path.exists(directory):
        os.makedirs(directory)
    for zz in range(len(ds.dev_list)):
        
        fname = ds.dev_list[zz]
        sys.stdout.write("\t" + str(zz + 1) + "/" + str(len(ds.dev_list)) + " " + fname)
        fn = fname[fname.rfind('/'):]
        [disc, cont] = ds.read_wave(fname)
        labels = ds.read_lab(fname)
        wav_out=network.synthesize(labels, len(disc)+4000)
        loss=0
        num_frames=len(wav_out)
        sys.stdout.write(" loss=" + str(loss / num_frames) + "\n")
        wav=ds.ulaw_decode(wav_out, discreete=False)
        wav_decoded = np.asarray(wav, dtype=np.int16)
        
        scipy.io.wavfile.write(directory + fn + ".wav", network.sample_rate, wav_decoded)
        
        
    stop = time.time()
    total += 1
    sys.stdout.write(" acc=" + str((1.0-float(errors) / total)) + " execution time=" + str(stop-start) + "\n")
    return (1.0-float(errors) / total)

def train_decoder(train, dev, model_output_base, itt_no_improve, config):
    sys.stdout.write("Testing world vocoder on devset\n");
    ds = Dataset(train, dev)
    
    worldVocoder = WorldVocoder(ds.sample_rate, config.win_size)
    recode(worldVocoder, ds, model_output_base)
        
def train(train, dev, model_output_base, itt_no_improve, config):
    sys.stdout.write("Pretraining powerspectrum decoder on static embeddings\n");
    ds = Dataset(train, dev)
    
    #coder = WorldVocoder(ds.sample_rate, config.win_size)
    network = SampleRNN(SampleRNNConfig()) #vocoder #f0 softmax #vuv
    epoch = 0
    best_acc = 0
    itt = itt_no_improve
    while itt > 0:
        sys.stdout.write("Starting epoch " + str(epoch) + "\n")
        #acc = synthesize(decoder, ds, model_output_base)
        nf = len(ds.train_list)
        cf = 0
        for fname in ds.train_list:
            
            sys.stdout.write(str(cf) + "/" + str(nf) + " " + fname + "...")
            sys.stdout.flush()
            [disc, cont] = ds.read_wave(fname)
            if len(disc)<=10*ds.sample_rate: #de scos dupa
                cf += 1
                start = time.time()
                labels = ds.read_lab(fname)
                loss, num_frames = network.learn(cont, disc)

                stop = time.time()
                sys.stdout.write(" avg loss=" + str(loss/(num_frames)) + " execution time=" + str(stop-start) + "\n")
                if cf % 10 == 0:
                    acc = synthesize(network, ds, model_output_base)
                    if acc > best_acc:
                        best_acc = acc
                        network.store_model(model_output_base + "/att_network.bestAcc")
                    network.store_model(model_output_base + "/att_network.last")
            else:
                sys.stdout.write("too long. skipping\n")
        
        
        network.store_model(model_output_base + "/att_network.last")
        #network.store_model(model_output_base + "/encoder.network")
        acc = synthesize(network, ds, model_output_base)
        if acc > best_acc:
            best_acc = acc
            network.store_model(model_output_base + "/att_network.bestAcc")
        sys.stdout.write("\n\nDevset accuracy is " + str(acc) + "\n\n")
        itt -= 1
        epoch += 1
        
if len(sys.argv) == 1:
    display_help()
else:
    if (sys.argv[1] == "--train" and len(sys.argv) == 6):
        config = Config()
        train(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), config)
    else:
        display_help()