/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.language;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;
import ro.racai.tts.WS.NLP.formats.BaseFeatureOutput;
import ro.racai.tts.WS.NLP.machinelearning.ID3;
import ro.racai.tts.WS.NLP.utils.IO;
import ro.racai.tts.WS.NLP.utils.Log;
//ro
//--train models/ro/training_data/lts.train models/ro/training_data/lts.test models/ro/lts.id3
//--train models/ro/training_data/syl.train models/ro/training_data/syl.test models/ro/syl.id3
//--train models/ro/training_data/lemma.train models/ro/training_data/lemma.test models/ro/lemma.id3
//--train models/ro/training_data/tag.train models/ro/training_data/tag.test models/ro/tag.id3
//--train models/ro/training_data/stress.train models/ro/training_data/stress.test models/ro/stress.id3
//en
//--train models/en/training_data/lts.train models/en/training_data/lts.test models/en/lts.id3
//--train models/en/training_data/syl.train models/en/training_data/syl.test models/en/syl.id3
//--train models/en/training_data/tag.train models/en/training_data/tag.test models/en/tag.id3
//--train models/en/training_data/tag.train models/en/training_data/tag.test models/en/tag.id3
//--train models/en/training_data/stress.train models/en/training_data/stress.test mode

/**
 * This class is used to instantiate a custom language processing pipeline. It
 * reads data from a conf file that must contain references to baseprocessors
 * and an output formatter
 *
 * @author tibi
 */
public class LanguagePipe {

    private static void displayHelp() {
        System.out.println("MLPA v0.9 beta");
        System.out.println("authors: Adrian Zafiu and Tiberiu Boros");
        System.out.println("USAGE:");
        System.out.println("\t--test\t run a synthetic test using a simple NLP pipeline");
        System.out.println("\t--process <conf_pipeline_file> <models_folder> <input_file>\t run MLPLA on the input file");
        System.out.println("\t--word2lts <conf_pipeline_file> <models_folder> <input_file> <output_file>\t run MLPLA on the input file");
        System.out.println("\t--process_train <conf_pipeline_file> <models_folder> <input_file> <phs file>\t run MLPLA on the input file and force feature output to compensate for actual phonetic sequence");
        System.out.println("\t--train <train_file> <test_file> <output_file>\t build an ID3 model from the training file, display accuracy on the test file and store the model in the output file");
        System.out.println("\t--draw-tree <train_file> <test_file> <output_file>\t build an ID3 model and render it from the training file, display accuracy on the test file and store the model in the output file");
    }

    private static void buildID3(String train, String test, String output) throws IOException {
        ID3 id3 = ID3.trainFromFile(train, 0, true);
        testID3(id3, test);
        id3.saveModel(new PrintStream(output));
        id3.saveGraph(new PrintStream(output + ".dot"));

    }

    private static void drawID3(String train, String test, String output) throws IOException {
        ID3 id3 = ID3.trainFromFile(train, 0, true);
        testID3(id3, test);
        id3.saveGraph(new PrintStream(output));
    }

    private static void testID3(ID3 id3, String test) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(test));
        String line = null;
        int perr = 0;
        int werr = 0;
        int pt = 0;
        int wt = 0;
        boolean ok = true;
        while ((line = br.readLine()) != null) {
            if (line.trim().isEmpty()) {
                wt++;
                if (!ok) {
                    werr++;
                }
                ok = true;
            } else {
                pt++;
                String[] parts = line.split(" ");
                List<String> feats = new ArrayList<>();
                for (int i = 1; i < parts.length; i++) {
                    feats.add(parts[i]);
                }
                String output = id3.classify(feats);
                if (!output.split(" ")[0].equals(parts[0])) {
                    perr++;
                    ok = false;
                }
            }
        }
        br.close();
        if (wt == 0) {
            wt++;
        }
        System.out.println("Processed " + wt + " examples\n\twith " + werr + " errors (" + (1.0 - (float) werr / wt) + " accuracy) \n\twith " + pt + " individual states\n\t\tand " + perr + " errors (" + (1.0 - (float) perr / pt) + " accuracy)");
    }

    private static void process(String conf_pipeline_file, String models, String input) throws IOException, FileNotFoundException, ClassNotFoundException, InstantiationException, IllegalAccessException {
        LanguagePipe lp = new LanguagePipe(conf_pipeline_file, models);
        BufferedReader br = IO.openFile(input);
        String line = null;

        while ((line = br.readLine()) != null) {
            if (!line.trim().startsWith("#")) {
                /*String[] parts = line.split(" ");
             Token[] t = new Token[parts.length];
             for (int i = 0; i < parts.length; i++) {
             t[i] = new Token();
             t[i].word = parts[i];
             }*/
                List<List<Token>> tokList = lp.input.parseInput(line);
                for (int i = 0; i < tokList.size(); i++) {
                    Token[] temp = tokList.get(i).toArray(new Token[0]);
                    lp.process(temp);
                    lp.print(temp);
                }
            }
        }
        br.close();
    }

    private static void process_word2lts(String conf_pipeline_file, String models, String input, String output) throws IOException, FileNotFoundException, ClassNotFoundException, InstantiationException, IllegalAccessException {
        System.out.println("Conf " + conf_pipeline_file);
        System.out.println("Model " + models);
        System.out.println("Input " + input);
        System.out.println("Output " + output);

        LanguagePipe lp = new LanguagePipe(conf_pipeline_file, models);
        BufferedReader br = IO.openFile(input);
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));

        String line = null;

        while ((line = br.readLine()) != null) {
            if (!line.trim().startsWith("#")) {
                /*String[] parts = line.split(" ");
             Token[] t = new Token[parts.length];
             for (int i = 0; i < parts.length; i++) {
             t[i] = new Token();
             t[i].word = parts[i];
             }*/
                List<List<Token>> tokList = lp.input.parseInput(line);
                for (int i = 0; i < tokList.size(); i++) {
                    Token[] temp = tokList.get(i).toArray(new Token[0]);
                    lp.process(temp);
                    for (int j = 0; j < temp.length; j++) {
                        bw.write(temp[j].word);
                        bw.write("\t");
                        boolean first = true;
                        for (int s = 0; s < temp[j].phonemes.length; s++) {
                            if (first) {
                                first = false;
                            } else {
                                bw.write(" ");
                            }
                            bw.write(temp[j].phonemes[s]);
                        }
                        bw.write("\n");
                    }
                    //lp.print(temp);
                }
            }
        }
        br.close();
        bw.close();
    }

    private static void process_train(String conf_pipeline_file, String models, String input, String phs) throws IOException, FileNotFoundException, ClassNotFoundException, InstantiationException, IllegalAccessException {
        LanguagePipe lp = new LanguagePipe(conf_pipeline_file, models);
        BufferedReader br = new BufferedReader(new FileReader(input));
        String line = null;
        while ((line = br.readLine()) != null) {
            List<List<Token>> tokList = lp.input.parseInput(line);
            for (int i = 0; i < tokList.size(); i++) {
                Token[] temp = tokList.get(i).toArray(new Token[0]);
                lp.process(temp);
                temp = resolve_conflicts(temp, phs);
                if (temp != null) {
                    lp.print(temp);
                }
            }
        }
        br.close();
    }

    private static Token[] resolve_conflicts(Token[] t, String phs) throws FileNotFoundException, IOException {
        List<Token> tmp = new ArrayList<>();
        //citim phs
        List<Token> phs_tokens = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(phs));
        String line = null;
        String word = null;
        List<String> phonemes = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            String[] parts = line.split(" ");
            if (parts.length == 4) {
                if (word != null && !word.equals("sil")) {
                    Token tok = new Token();
                    tok.word = ASCII2UTF8(word);
                    tok.phonemes = phonemes.toArray(new String[0]);
                    phs_tokens.add(tok);
                }
                phonemes.clear();
                word = parts[3];
            }
            phonemes.add(parts[2]);
        }
        br.close();
        //facem armonizarea
        //trebuie ca intotdeauna sa avem mai multi tokeni in text decat in PHS
        if (t.length < phs_tokens.size()) {
            Log.e("Something may be wrong with the data. Unable to harmonize output. The number of tokens inside the PHS file exceedes that of the original sentence.");
        }

        int pos = 0;
        for (int i = 0; i < t.length; i++) {
            if (pos >= phs_tokens.size()) {
                Log.e("File " + phs + " has inapropiate content");
                return null;
            }
            if (ASCII2UTF8(t[i].word).toLowerCase().equals(phs_tokens.get(pos).word)) {
                tmp.add(t[i]);
                Token pt = phs_tokens.get(pos);
                if (pos < phs_tokens.size() - 1) {
                    pos++;
                }
                if (pt.phonemes[pt.phonemes.length - 1].equals("sp")) {
                    if (i < t.length - 1 && (t[i + 1].phonemes.length != 1 || t[i + 1].phonemes[0].equals("pau"))) {
                        Token vt = new Token();
                        vt.word = "<VIRTUAL_TOKEN>";
                        vt.phonemes = new String[]{"sp"};
                        vt.tag = "<VITRUAL_TOKEN>";
                        vt.syllables = new String[]{"P"};
                        vt.lemma = "sp";
                        vt.dependencies = "";
                        vt.stress = new int[1];
                        vt.chunk = "";
                        tmp.add(vt);
                    }
                }

            } else if ((t[i].phonemes.length != 1 || !(t[i].phonemes[0].equals("sp") || t[i].phonemes[0].equals("pau"))) && pos < phs_tokens.size() - 1) {
                pos++;
                i--;
            } else {
                tmp.add(t[i]);
            }
        }
        return tmp.toArray(new Token[0]);
    }

    private static String ASCII2UTF8(String word) throws UnsupportedEncodingException {
        StringBuilder sb = new StringBuilder();
        String s = word.replaceAll("[^\\x20-\\x7e]", "");
        for (int i = 0; i < s.length(); i++) {
            if (s.substring(i, i + 1).equals("\\")) {
                i += 3;
            } else {
                sb.append(s.substring(i, i + 1));
            }
        }
        return sb.toString().trim();
    }

    BaseProcessor[] bpList = null;
    public BaseFeatureOutput output = null;
    BaseInputProcessor input = null;

    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     * @throws java.lang.ClassNotFoundException
     * @throws java.lang.InstantiationException
     * @throws java.lang.IllegalAccessException
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, ClassNotFoundException, InstantiationException, IllegalAccessException {
//        if (true){
//            NeuralTagger nt=new NeuralTagger();
//            //nt.train("models/ro/training_data/1984.train", "models/ro/training_data/1984.test", word_embeddings_file);
//            nt.computeEmbeddings("models/ro/training_data/raw_text.txt", 300, 2, "models/ro/word.embeddings");
//            return;
//        }
//        //pregatim LTS si STRESS pentru CRF
//        clearFormat("models/en/training_data/lts.train", "models/en/training_data/lts.crf.train");
//        clearFormat("models/en/training_data/lts.test", "models/en/training_data/lts.crf.test");
//        clearFormat("models/en/training_data/stress.train", "models/en/training_data/stress.crf.train");
//        clearFormat("models/en/training_data/stress.test", "models/en/training_data/stress.crf.test");
//        clearFormat("models/en/training_data/syl.train", "models/en/training_data/syl.crf.train");
//        clearFormat("models/en/training_data/syl.test", "models/en/training_data/syl.crf.test");
//        if (true) {
//            return;
//        }

        /*        String s = "- Căpitanul spuse: 55*4=10!";
        for(String ss : Utils.simpleWhiteSpaceSplit(s)) 
            IO.outln(ss);
        if(true) return;
        
         
        boolean german = true;
        boolean train = false;
        if (german) {
            if (train) {
                // aici facem un test default
                //TaggerTrainer.CreateTrainingFileFromNegraCorpus( "res/de/tag/negra-corpus.tt", "res/de/tag/out.tag.de");
                //SyllabifierTrainer.CreateTrainingFileFrom("res/de/syl/german.syl", "res/de/syl/german.syl.feats");
                //DecompositionTrainer.CreateTrainingFileFrom("res/de/tok/sdaz.decompound", "res/de/tok/sdaz.decompound.feats");
            } else {
                String sent = "Am plecat la piata.";
                LanguagePipe lp = new LanguagePipe("etc/languagepipe.de.conf", "models/de");
                List<List<Token>> tmp = lp.input.parseInput(sent);
                for (int i = 0; i < tmp.size(); i++) {
                    Token[] temp = tmp.get(i).toArray(new Token[0]);
                    lp.process(temp);
                    lp.print(temp);
                }

            }
            return;
        }
         */
        if (args.length == 0) {
            displayHelp();
        } else if (args[0].equals("--test")) {

//            StringBuilder bal1000 = new StringBuilder();
//            String line;
//            BufferedReader br = new BufferedReader(new FileReader("test/bal1000.txt"));
//            while ((line = br.readLine()) != null) {
//                if (bal1000.length() > 0) {
//                    bal1000.append("\n");
//                }
//                bal1000.append(line);
//            }
//            br.close();
            LanguagePipe lp = new LanguagePipe("etc/languagepipe.hts.en.conf", "models/en");
            //String sent = "Georgel și Ion merg la cinema ! Cine ești tu ? Mă gândesc să ... acest test este pentru limba română . chintesența leoaicei tinere este cutezanța . exemplu și experiență . Ionel este laureat Nobel pentru realizări exorbitante .";
            String sent = "This is a simple test. \"We'll have to go out,\" she said with a shout, \"before all the butchers are closed.\"";
            List<List<Token>> tmp = lp.input.parseInput(sent);
            for (int i = 0; i < tmp.size(); i++) {
                Token[] temp = tmp.get(i).toArray(new Token[0]);
                lp.process(temp);
                lp.print(temp);
            }

        } else if (args.length == 4 && args[0].equals("--train")) {
            buildID3(args[1], args[2], args[3]);
        } else if (args.length == 4 && args[0].equals("--draw-tree")) {
            drawID3(args[1], args[2], args[3]);
        } else if (args.length == 4 && args[0].equals("--process")) {
            // --process etc/languagepipe.de.conf models/de input.de.txt
            process(args[1], args[2], args[3]);
        } else if (args.length == 5 && args[0].equals("--process_train")) {
            process_train(args[1], args[2], args[3], args[4]);
        } else if (args.length == 5 && args[0].equals("--word2lts")) {
            // --word2lts etc/languagepipe.de.conf models/de input output
            String conf = args[1];
            String model = args[2];
            String input = args[3];
            String output = args[4];
            process_word2lts(conf, model, input, output);

        } else {
            displayHelp();
        }
    }

    private static void clearFormat(String source, String destination) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(source));
        BufferedWriter bw = new BufferedWriter(new FileWriter(destination));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) {
                bw.write("\n");
            } else {
                String[] parts = line.split(" ");
                for (int i = 1; i < parts.length; i++) {
                    if (parts[i].startsWith("cl:")) {
                        bw.write(parts[i].substring(3) + " ");
                    }
                }
                bw.write(parts[0] + "\n");
            }
        }
        bw.close();
        br.close();
    }

    /**
     *
     * @param etclanguagepipeconf
     * @param models
     * @throws FileNotFoundException
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws InstantiationException
     * @throws IllegalAccessException
     */
    public LanguagePipe(String etclanguagepipeconf, String models) throws FileNotFoundException, IOException, ClassNotFoundException, InstantiationException, IllegalAccessException {
        //TODO: de implementat citirea cum trebuie (acum e doar de test)
        Log.i("Parsing configuration file [" + etclanguagepipeconf + "] with models in [" + models + "]");
        BufferedReader br = IO.openFile(etclanguagepipeconf);
        String line = null;
        List<String> processors = new ArrayList<>();

        while ((line = br.readLine()) != null) {
            if (line.startsWith("[") || line.startsWith("#")) {
                continue;
            }
            processors.add(line);
        }
        //procesatoare
        Log.i("Creating processing pipeline");
        int all_processors = 0;
        Log.i("Data input is handled by " + processors.get(0));
        Class<?> clazz = Class.forName(processors.get(0));
        input = (BaseInputProcessor) clazz.newInstance();
        Log.i("Pipeline length is " + (processors.size() - 2));
        bpList = new BaseProcessor[processors.size() - 2];
        for (int i = 0; i < bpList.length; i++) {
            Log.i("Step " + (i + 1) + " is performed using " + processors.get(i + 1));
            clazz = Class.forName(processors.get(i + 1));
            BaseProcessor bp = (BaseProcessor) clazz.newInstance();
            bp.loadModel(models);
            bpList[i] = bp;
            if ((all_processors & bp.getProcessorType()) != 0) {
                Log.w(processors.get(i) + " duplicates functionality");
            }
            all_processors = all_processors | bp.getProcessorType();
        }
        //output
        Log.i("Data output is handled by " + processors.get(processors.size() - 1));
        clazz = Class.forName(processors.get(processors.size() - 1));
        output = (BaseFeatureOutput) clazz.newInstance();

        //validare
        if ((BaseProcessor.PROCESSOR_CHUNKER & all_processors) == 0) {
            Log.w("Pipeline is missing a chunker");
        }
        if ((BaseProcessor.PROCESSOR_LEMMATIZER & all_processors) == 0) {
            Log.w("Pipeline is missing a lemmatizer");
        }
        if ((BaseProcessor.PROCESSOR_LTS & all_processors) == 0) {
            Log.w("Pipeline is missing a LTS converter");
        }
        if ((BaseProcessor.PROCESSOR_PARSER & all_processors) == 0) {
            Log.w("Pipeline is missing a parser");
        }
        if ((BaseProcessor.PROCESSOR_SYLLABIFIER & all_processors) == 0) {
            Log.w("Pipeline is missing a syllabifier");
        }
        if ((BaseProcessor.PROCESSOR_TAGGER & all_processors) == 0) {
            Log.w("Pipeline is missing a tagger");
        }
        if ((BaseProcessor.PROCESSOR_NORMALIZER & all_processors) == 0) {
            Log.w("Pipeline is missing a word normalizer");
        }
    }

    public Token[][] process(String text) {
        List<List<Token>> tmp = input.parseInput(text);
        Token[][] o = new Token[tmp.size()][];
        for (int i = 0; i < tmp.size(); i++) {
            Token[] t = tmp.get(i).toArray(new Token[0]);
            process(t);
            o[i] = t;
        }
        return o;
    }

    /**
     *
     * @param tokens
     */
    public void process(Token[] tokens) {
        for (int i = 0; i < bpList.length; i++) {
            bpList[i].processTokens(tokens);
        }
    }

    /**
     *
     * @param tokens
     */
    public void print(Token[] tokens) {
        output.print(tokens);
    }

}
