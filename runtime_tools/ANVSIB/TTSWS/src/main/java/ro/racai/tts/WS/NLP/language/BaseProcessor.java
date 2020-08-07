/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.language;

/**
 * Each instance of a base processor must implement either a POS tagger, LTS
 * converter, Syllabifier or Lemmatizer. A single instance can perform multiple
 task (though this is not recommended). For example, a single BaseProcessor
 can fill in the information regarding both lemma and syllables, provided of
 course that the getProcessorType() method returns PROCESSOR_LEMMATIZER |
 PROCESSOR_SYLLABIFIER
 *
 * @author tibi
 */
public interface BaseProcessor {

    /**
     *
     */
    public static int PROCESSOR_TAGGER = 1;

    /**
     *
     */
    public static int PROCESSOR_LTS = 2;

    /**
     *
     */
    public static int PROCESSOR_SYLLABIFIER = 4;

    /**
     *
     */
    public static int PROCESSOR_LEMMATIZER = 8;

    public static int PROCESSOR_CHUNKER = 16;

    public static int PROCESSOR_PARSER = 32;

    public static int PROCESSOR_STRESS = 64;

    public static int PROCESSOR_NORMALIZER = 128;

    public static int PROCESSOR_SHORTPAUSEADDITION = 256;
    
    public static int PROCESSOR_TAGGER_2 = 512;

    /**
     * This method is used to control the completeness of the processing
     * pipeline.
     *
     * @return returns one or a combination of the BASE_PROCESSOR types,
     * according to its designated usage
     */
    public int getProcessorType();

    /**
     * This method receives a list of Tokens that contain sentence words and
     * their associated attributes computed up to the moment. For every token in
 the list, the processor must compute and fill in his own target attribute
 (i.e. a PROCESSOR_LEMMATIZER must fill in the lemma of the word)
     *
     * @param tokens - the list of tokens that will be operated upon
     */
    public void processTokens(Token[] tokens);

    /**
     * This method is invoked by the LanguagePipe. It is used to load the model
     * data and perform any additional operations on the model.
     *
     * @param folder - absolute or relative folder from which to load the
     * model(s). It does not contain the filename itself
     */
    public void loadModel(String folder);
}
