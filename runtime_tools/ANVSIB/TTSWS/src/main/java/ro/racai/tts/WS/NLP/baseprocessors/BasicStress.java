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
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_STRESS;
import ro.racai.tts.WS.NLP.language.Token;
import ro.racai.tts.WS.NLP.machinelearning.ID3;

/**
 * Provided basic functionality for lexical stress prediction on OOV words.
 * 97.44% accuracy on romanian OOV and 77.63% (could be a higher value - corpora
 * is not validated) accuracy on English OOV
 *
 * @author tibi
 */
public class BasicStress implements BaseProcessor {

    @Override
    public int getProcessorType() {
        return PROCESSOR_STRESS;
    }

    @Override
    public void processTokens(Token[] tokens) {
        List<String> feats = new ArrayList<>();
        List<String> labels = new ArrayList<>();
        //List<String> syllables = new ArrayList<>();

        for (int i = 0; i < tokens.length; i++) {
            Token t = tokens[i];
            if (Pattern.matches("\\p{Punct}", t.word)) {
                t.stress = new int[t.syllables.length];
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
                feats.add("pos:" + tokens[i].tag.substring(0, 1));
                String lab = id3.classify(feats).split(" ")[0];
                plab = lab;
                //t.phonemes[k] = phon;
                labels.add(lab);
            }

            t.stress = new int[t.syllables.length];
            int lp = 0;
            if (t.syllables.length == 1) {
                continue;
            }
            boolean has_stress = false;
            for (int k = 0; k < t.syllables.length; k++) {
                for (int l = 0; l < t.syllables[k].length(); l++) {
                    if (labels.get(lp).equals("PS")) {
                        t.stress[k] = 1;
                        has_stress = true;
                    }
                    lp++;
                }
            }
            //@All: daca avem un cuvant cu mai multe silabe, teoretic pentru romana si franceza puteam sa punem accent pe ultima silaba, insa pentru engleza nu tine
            //nu stiu ce sa fac.

        }

    }

    ID3 id3;

    @Override
    public void loadModel(String folder) {
        try {
            id3 = ID3.createFromFile(folder + "/stress.id3");
        } catch (IOException ex) {
            Logger.getLogger(BasicStress.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
