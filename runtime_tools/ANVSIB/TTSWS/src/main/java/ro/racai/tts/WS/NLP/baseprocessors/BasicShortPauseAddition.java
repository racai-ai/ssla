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
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_SHORTPAUSEADDITION;
import ro.racai.tts.WS.NLP.language.Token;
import ro.racai.tts.WS.NLP.machinelearning.ID3;

/**
 *
 * @author adrian
 */
public class BasicShortPauseAddition implements BaseProcessor{
     @Override
    public int getProcessorType() {
        return PROCESSOR_SHORTPAUSEADDITION;
    }

    @Override
    public void processTokens(Token[] tokens) {
        List<String> feats = new ArrayList<>();
        
//        //List<String> labels = new ArrayList<>();
//        //List<String> syllables = new ArrayList<>();
//
//        for (int i = 0; i < tokens.length; i++) {
//            Token t = tokens[i];
//            
//            // *** Nu sunt sigur ce face asta. Banuiesc ca e finalul si pun atributul pentru punct
//            if (Pattern.matches("\\p{Punct}", t.word)) {
//                t.shortPauseAdder = false;
//                continue;
//            }
//            
//            labels.clear();
//            
//            for (int k = 0; k < t.word.length(); k++) {
//                // de extras feature-ile finale
//                
//                feats.clear();
//                //feats.add("feat:" + value);
//                String lab = id3.classify(feats).split(" ")[0];
//                plab = lab;
//                //t.phonemes[k] = phon;
//                labels.add(lab);
//            }
//
//            t.stress = new int[t.syllables.length];
//            int lp = 0;
//            if (t.syllables.length == 1) {
//                continue;
//            }
//            boolean has_stress = false;
//            for (int k = 0; k < t.syllables.length; k++) {
//                for (int l = 0; l < t.syllables[k].length(); l++) {
//                    if (labels.get(lp).equals("PS")) {
//                        t.stress[k] = 1;
//                        has_stress = true;
//                    }
//                    lp++;
//                }
//            }
//            //@All: daca avem un cuvant cu mai multe silabe, teoretic pentru romana si franceza puteam sa punem accent pe ultima silaba, insa pentru engleza nu tine
//            //nu stiu ce sa fac.
//
//        }

    }

    ID3 id3;

    @Override
    public void loadModel(String folder) {
        try {
            id3 = ID3.createFromFile(folder + "/sp.id3");
        } catch (IOException ex) {
            Logger.getLogger(BasicLTS.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
