/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.formats;

import ro.racai.tts.WS.NLP.language.Token;



/**
 * This interface implements output formating functions. In practice it is used
 * to shape the output format from the MLPLA platform in various applications
 * such as basic text processing or speech synthesis (see HTS file format for
 * future details)
 *
 * @author tibi
 */
public interface BaseFeatureOutput {

    /**
     * Outputs a formatted feature-based representation of the processed tokens
     *
     * @param tokens Pre-processed tokens
     */
    public String print(Token[] tokens);

    /**
     * Creates a list of formatted feature-based representation of the processed
     * tokens
     *
     * @param tokens Pre-processed tokens
     * @return lines of feature sets at the desired granularity (ie. word or
     * phoneme) - depending on the applications
     */
    public String[] getFeatures(Token[] tokens);
}
