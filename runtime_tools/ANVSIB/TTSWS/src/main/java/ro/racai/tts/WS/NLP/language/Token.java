/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.language;

/**
 *
 * @author tibi
 */
public class Token {

    /**
     *
     */
    public String word;

    /**
     *
     */
    public String tag;
    /**
     * Contains a stacked list of chunks the word is part of (can be converted
     * into a tree)
     */
    public String chunk;

    /**
     * Contains a list of dependencies with associated types
     */
    public String dependencies;

    /**
     *
     */
    public String lemma;

    /**
     *
     */
    public String[] phonemes;
    public String[] phonemes_stress;

    /**
     *
     */
    public String[] syllables;
    public int[] stress;
    public boolean shortPauseAdder;
    
    public String ner="";
    
    public String sentimentGroup="";
    public String sentimentWord="";
    
    public String paraphrase="";
}
