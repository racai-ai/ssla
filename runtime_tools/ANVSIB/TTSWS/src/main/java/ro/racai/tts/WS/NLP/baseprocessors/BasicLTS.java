/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.baseprocessors;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import ro.racai.tts.WS.NLP.language.BaseProcessor;
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_LTS;
import ro.racai.tts.WS.NLP.language.Token;
import ro.racai.tts.WS.NLP.machinelearning.ID3;

/**
 * Implements standard LTS functionality. Romanian accuracy on OOV words is
 * 94.1%
 *
 * @author tibi
 */
public class BasicLTS implements BaseProcessor {

    @Override
    public int getProcessorType() {
        return PROCESSOR_LTS;
    }

    @Override
    public void processTokens(Token[] tokens) {
        List<String> feats = new ArrayList<>();
        for (int i = 0; i < tokens.length; i++) {
            Token t = tokens[i];
            if (Pattern.matches("\\p{Punct}", t.word)) {
                if (t.word.equals(",")) {
                    t.phonemes = new String[]{"sp"};
                } else {
                    t.phonemes = new String[]{"pau"};
                }
                continue;
            }
            t.phonemes = new String[t.word.length()];
            String plab = "_";
            for (int k = 0; k < t.phonemes.length; k++) {
                int ndxP1 = k - 1;
                while (ndxP1 >= 0 && t.word.charAt(ndxP1) == '-') {
                    ndxP1--;
                }

                int ndxP2 = ndxP1 - 1;
                while (ndxP2 >= 0 && t.word.charAt(ndxP2) == '-') {
                    ndxP2--;
                }

                int ndxN1 = k + 1;
                while (ndxN1 < t.word.length() && t.word.charAt(ndxN1) == '-') {
                    ndxN1++;
                }

                int ndxN2 = ndxN1 + 1;
                while (ndxN2 < t.word.length() && t.word.charAt(ndxN2) == '-') {
                    ndxN2++;
                }

                String ppl = "_";
                String pl = "_";
                String nl = "_";
                String nnl = "_";

                if (ndxP1 >= 0) {
                    pl = t.word.substring(ndxP1, ndxP1 + 1);
                }
                if (ndxP2 >= 0) {
                    ppl = t.word.substring(ndxP2, ndxP2 + 1);
                }

                String cl = t.word.substring(k, k + 1);

                if (ndxN1 < t.word.length()) {
                    nl = t.word.substring(ndxN1, ndxN1 + 1);
                }
                if (ndxN2 < t.word.length()) {
                    nnl = t.word.substring(ndxN2, ndxN2 + 1);
                }

                feats.clear();
                if (!cl.equals("-")) {
                    feats.add("ppl:" + ppl.toLowerCase());
                    feats.add("pl:" + pl.toLowerCase());
                }
                feats.add("cl:" + cl.toLowerCase());
                if (!cl.equals("-")) {
                    feats.add("nl:" + nl.toLowerCase());
                    feats.add("nnl:" + nnl.toLowerCase());
                    feats.add("plab:" + plab);
                }
                String phon = id3.classify(feats).split(" ")[0];
                plab = phon;
                if (cl.equals("-")) {
                    t.phonemes[k] = "-";
                } else {
                    t.phonemes[k] = phon;
                }
            }
        }
    }

    ID3 id3;

    @Override
    public void loadModel(String folder) {
        try {
            id3 = ID3.createFromFile(folder + "/lts.id3");
            id3.saveGraph(new PrintStream(new File("lts.graph")));
        } catch (IOException ex) {
            Logger.getLogger(BasicLTS.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
