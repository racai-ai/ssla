/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.baseprocessors;

import ro.racai.tts.WS.NLP.language.BaseProcessor;
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_NORMALIZER;
import ro.racai.tts.WS.NLP.language.Token;

/**
 *
 * @author tibi
 */
public class BasicNormalizer implements BaseProcessor {

    @Override
    public int getProcessorType() {
        return PROCESSOR_NORMALIZER;
    }

    @Override
    public void processTokens(Token[] tokens) {
        
    }

    @Override
    public void loadModel(String folder) {
        
    }

}
