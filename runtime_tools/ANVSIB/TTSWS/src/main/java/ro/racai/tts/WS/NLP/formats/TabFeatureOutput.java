/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.formats;

import ro.racai.tts.WS.NLP.language.Token;


/**
 *
 * @author tibi
 */
public class TabFeatureOutput implements BaseFeatureOutput {

    /**
     *
     * @param tokens
     */
    @Override
    public String print(Token[] tokens) {
        StringBuilder sb=new StringBuilder();
        for (int i = 0; i < tokens.length; i++) {
            sb.append(tokens[i].word);
            if (tokens[i].tag != null) {
                sb.append("\t" + tokens[i].tag);
            }
            if (tokens[i].lemma != null) {
                sb.append("\t" + tokens[i].lemma);
            }
            if (tokens[i].ner != null) {
                sb.append("\t" + tokens[i].ner);
            }
            if (tokens[i].syllables != null) {
                sb.append("\t");
                for (int k = 0; k < tokens[i].syllables.length; k++) {
                    if (tokens[i].stress != null) {
                        if (tokens[i].stress[k] == 1) {
                            sb.append("'");
                        }
                    }
                    sb.append(tokens[i].syllables[k]);
                    if (k < tokens[i].syllables.length - 1) {
                        sb.append("-");
                    }
                }
            }
            if (tokens[i].phonemes != null) {
                sb.append("\t");
                for (int k = 0; k < tokens[i].phonemes.length; k++) {
                    sb.append(tokens[i].phonemes[k]);
                    if (k < tokens[i].phonemes.length - 1) {
                        sb.append(" ");
                    }
                }
            }
            if (tokens[i].sentimentWord != null) {
                sb.append("\t" + tokens[i].sentimentWord);
            }
            if (tokens[i].sentimentGroup != null) {
                sb.append("\t" + tokens[i].sentimentGroup);
            }
            if (tokens[i].chunk != null) {
                sb.append("\t" + tokens[i].chunk);
            }
            if (tokens[i].dependencies != null) {
                sb.append("\t" + tokens[i].dependencies);
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    /**
     *
     * @param tokens
     * @return
     */
    @Override
    public String[] getFeatures(Token[] tokens) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
