/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.baseprocessors;

import ro.racai.tts.WS.NLP.language.BaseProcessor;
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_PARSER;
import ro.racai.tts.WS.NLP.language.Token;

/**
 * Implements standard parsing functionality
 *
 * @author tibi
 */
public class BasicParser implements BaseProcessor {

    @Override
    public int getProcessorType() {
        return PROCESSOR_PARSER;
    }

    @Override
    public void processTokens(Token[] tokens) {
        for (int i = 0; i < tokens.length; i++) {
            Token t = tokens[i];
            t.dependencies = "";
        }
    }

    @Override
    public void loadModel(String folder) {
        //throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
