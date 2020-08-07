/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.baseprocessors;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import ro.racai.tts.WS.NLP.language.BaseProcessor;
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_TAGGER;
import ro.racai.tts.WS.NLP.language.Token;
import ro.racai.tts.WS.NLP.machinelearning.ID3;

/**
 * Implements standard POS tagging functionality. Romanian accuracy is 94.19%
 *
 * @author tibi
 */
public class BasicTagger implements BaseProcessor {

    @Override
    public int getProcessorType() {
        return PROCESSOR_TAGGER;
    }

    @Override
    public void processTokens(Token[] tokens) {
        String pt = "_";
        String ppt = "_";

        for (int i = 0; i < tokens.length; i++) {
            String pw = "_";
            String ppw = "_";
            String nw = "_";
            String nnw = "_";
            String w = tokens[i].word;
            if (i > 0) {
                pw = tokens[i - 1].word;
            }
            if (i > 1) {
                ppw = tokens[i - 2].word;
            }
            if (i < tokens.length - 1) {
                nw = tokens[i + 1].word;
            }
            if (i < tokens.length - 2) {
                nnw = tokens[i + 2].word;
            }

            //trasaturi (le facem de kkt momentan)
            String feats = "";
            //feats += " w:" + w;
            feats += " pt:" + pt;
            feats += " ppt:" + ppt;
            if (w.toLowerCase().equals(w)) {
                feats += " w_caps:lower";
            } else if (w.toUpperCase().equals(w)) {
                feats += " w_caps:upper";
            } else {
                feats += " w_caps:mixed";
            }
            w = w.toLowerCase();
            if (nw.toLowerCase().equals(nw)) {
                feats += " nw_caps:lower";
            } else if (nw.toUpperCase().equals(nw)) {
                feats += " nw_caps:upper";
            } else {
                feats += " nw_caps:mixed";
            }

            w = w.toLowerCase();
            feats += " ws:" + w.length();
            for (int k = 0; k < w.length(); k++) {
                feats += " wh_" + k + ":" + w.substring(k, k + 1);
            }
            /*feats += " w_s0:" + w.substring(0, 1);
             if (w.length() > 1) {
             feats += " w_s1:" + w.substring(1, 2);
             }
             if (w.length() > 2) {
             feats += " w_s2:" + w.substring(2, 3);
             }
             if (w.length() > 3) {
             feats += " w_s3:" + w.substring(3, 4);
             }*/

            feats += " w_e0:" + w.substring(w.length() - 1, w.length());
            if (w.length() > 1) {
                feats += " w_e1:" + w.substring(w.length() - 2, w.length() - 1);
            } else {
                feats += " w_e1:_";
            }
            if (w.length() > 2) {
                feats += " w_e2:" + w.substring(w.length() - 3, w.length() - 2);
            } else {
                feats += " w_e2:_";
            }
            if (w.length() > 3) {
                feats += " w_e3:" + w.substring(w.length() - 3, w.length() - 3);
            } else {
                feats += " w_e3:_";
            }
            feats += " nw_e0:" + nw.substring(nw.length() - 1, nw.length());
            if (nw.length() > 1) {
                feats += " nw_e1:" + nw.substring(nw.length() - 2, nw.length() - 1);
            } else {
                feats += " nw_e1:_";
            }
            if (nw.length() > 2) {
                feats += " nw_e2:" + nw.substring(nw.length() - 3, nw.length() - 2);
            } else {
                feats += " nw_e2:_";
            }
            if (nw.length() > 3) {
                feats += " nw_e3:" + nw.substring(nw.length() - 3, nw.length() - 3);
            } else {
                feats += " nw_e3:_";
            }
            String[] pp = feats.split(" ");
            List<String> f = new ArrayList<>();
            for (int k = 0; k < pp.length; k++) {
                f.add(pp[k]);
            }

            String tag = id3.classify(f).split(" ")[0];
            ppt = pt;
            pt = tag;
            tokens[i].tag = tag;
        }
    }
    ID3 id3;

    @Override
    public void loadModel(String folder) {
        try {
            id3 = ID3.createFromFile(folder + "/tag.id3");
        } catch (IOException ex) {
            Logger.getLogger(BasicLTS.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
