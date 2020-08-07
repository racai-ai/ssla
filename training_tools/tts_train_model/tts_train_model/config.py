# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

class Config:
    def __init__(self):
        #static
        self.encoder_size = 200
        self.encoder_layers = 2
        self.decoder_size = 1000
        self.decoder_layers = 2
        self.batch_size = 512 #32 ms at 16Khz
        
        self.att_proj_size = 100
        self.att_lsa_filters = 10
        self.att_lsa_input_size = 5
        self.receptive_input = 320 #20ms at 16Khz
        self.receptive_layers = [512, 256]
        self.receptive_dropout = 0.5
        self.attention_layers = [512, 256]
        self.phone_embeddings_size = 32
        self.context_embeddings_size = 8
        self.presoftmax_layers = [1024, 1024, 256]
        self.compute_attention_samples = 15 * 16 # 15 ms
        self.sample_trail_size = 2 #append samples to final MLP
        