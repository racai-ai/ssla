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
import java.util.regex.Pattern;
import ro.racai.tts.WS.NLP.language.BaseProcessor;
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_LEMMATIZER;
import ro.racai.tts.WS.NLP.language.Token;
import ro.racai.tts.WS.NLP.machinelearning.ID3;

/**
 * Implements standard lemmatization functionality
 *
 * @author tibi
 */
public class BasicLemmatizerOld implements BaseProcessor {

    @Override
    public int getProcessorType() {
        return PROCESSOR_LEMMATIZER;
    }

    @Override
    public void processTokens(Token[] tokens) {
        List<String> feats = new ArrayList<>();
        List<String> labels = new ArrayList<>();
        //List<String> syllables = new ArrayList<>();

        for (int i = 0; i < tokens.length; i++) {
            Token t = tokens[i];
            if (Pattern.matches("\\p{Punct}", t.word)) {
                t.lemma = t.word;
                continue;
            }
            //t.phonemes = new String[t.word.length()];
            labels.clear();
            String plab = "_";
            for (int k = 0; k < t.word.length(); k++) {
                String pppl = "_";
                String ppl = "_";
                String pl = "_";
                String nl = "_";
                String nnl = "_";
                String nnnl = "_";
                //
                if (k > 0) {
                    pl = t.word.substring(k - 1, k);
                }
                if (k > 1) {
                    ppl = t.word.substring(k - 2, k - 1);
                }
                if (k > 2) {
                    pppl = t.word.substring(k - 3, k - 2);
                }
                String cl = t.word.substring(k, k + 1);
                if (k < t.word.length() - 1) {
                    nl = t.word.substring(k + 1, k + 2);
                }
                if (k < t.word.length() - 2) {
                    nnl = t.word.substring(k + 2, k + 3);
                }
                if (k < t.word.length() - 3) {
                    nnnl = t.word.substring(k + 3, k + 4);
                }
                feats.clear();
                feats.add("pppl:" + pppl.toLowerCase());
                feats.add("ppl:" + ppl.toLowerCase());
                feats.add("pl:" + pl.toLowerCase());
                feats.add("cl:" + cl.toLowerCase());
                feats.add("nl:" + nl.toLowerCase());
                feats.add("nnl:" + nnl.toLowerCase());
                feats.add("nnnl:" + nnnl.toLowerCase());
                feats.add("plab:" + plab);
                feats.add("pos:" + tokens[i].tag);
                String lab = id3.classify(feats).split(" ")[0];
                plab = lab;
                //t.phonemes[k] = phon;
                labels.add(lab);
            }
            //syllables.clear();
            StringBuilder sb = new StringBuilder();
            for (int k = 0; k < t.word.length(); k++) {
                String lab = labels.get(k);
                if (lab.equals("*")) {
                    sb.append(t.word.substring(k, k + 1));
                } else if (lab.equals("__nil__")) {
                    //nu adaugam nimic
                } else {
                    sb.append(lab);
                }
            }
            t.lemma = sb.toString();
        }
    }

    ID3 id3;

    @Override
    public void loadModel(String folder) {
        try {
            id3 = ID3.createFromFile(folder + "/lemma.id3");
        } catch (IOException ex) {
            Logger.getLogger(BasicLTS.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
