/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.baseprocessors;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import ro.racai.tts.WS.NLP.language.BaseProcessor;
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_LEMMATIZER;
import ro.racai.tts.WS.NLP.language.Token;
import ro.racai.tts.WS.NLP.machinelearning.ID3;

/**
 * Implements standard lemmatization functionality
 *
 * @author tibi
 */
public class BasicLemmatizer implements BaseProcessor {

    @Override
    public int getProcessorType() {
        return PROCESSOR_LEMMATIZER;
    }

    @Override
    public void processTokens(Token[] tokens) {
        List<String> feats = new ArrayList<>();
        //List<String> labels = new ArrayList<>();
        //List<String> syllables = new ArrayList<>();

        for (int i = 0; i < tokens.length; i++) {
            Token t = tokens[i];
            if (Pattern.matches("\\p{Punct}", t.word)) {
                t.lemma = t.word;
                continue;
            }

            int n = t.word.length();
            feats.clear();
            for (int k = 1; k < 8; k++) {
                if (n - k >= 0) {
                    feats.add("c" + k + ":" + t.word.charAt(n - k));
                } else {
                    feats.add("c" + k + ":_");
                }
            }
            if (t.tag != null) {
                feats.add("m:" + t.tag);
            }
            String[] outs = id3.classify(feats).split(" ")[0].split("\\|");
            int sz = Integer.parseInt(outs[0]);
            String rpl = "";
            if (outs.length > 1) {
                rpl = outs[1];
            }
            if (n-sz>0){
                t.lemma = t.word.substring(0, n - sz) + rpl;
            }else{
                t.lemma=t.word.toLowerCase();
            }
        }
    }

    ID3 id3;

    @Override
    public void loadModel(String folder) {
        try {
            id3 = ID3.createFromFile(folder + "/lemma2.id3");
        } catch (IOException ex) {
            Logger.getLogger(BasicLTS.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static class Model {

        public static Map<String, Map<String, String>> readAll(String fileName, String msd2ctagFile) throws FileNotFoundException, IOException {
            BufferedReader br;
            Map<String, String> msd2ctag = new HashMap<String, String>();
            br = new BufferedReader(new FileReader(msd2ctagFile));
            for (String line = br.readLine(); line != null; line = br.readLine()) {
                String[] toks = line.split("[ \t]");
                if (toks.length == 2) {
                    if (toks[0].length() > 0 && toks[1].length() > 0) {
                        msd2ctag.put(toks[0], toks[1]);
                    }
                }
            }
            br.close();
            br = new BufferedReader(new FileReader(fileName));
            Map<String, Map<String, String>> packs = new HashMap<String, Map<String, String>>();
            for (String line = br.readLine(); line != null; line = br.readLine()) {
                if (!line.startsWith("#")) {
                    String[] toks = line.split("[ \t]");
                    if (toks.length == 3) {
                        String word = toks[0];
                        String lemma = toks[1];
                        if (toks[1].equals("=")) {
                            lemma = toks[0];
                        }
                        String msd = toks[2];
                        if (msd2ctag.containsKey(msd)) {
                            msd = msd2ctag.get(msd);
                        }
                        if (!packs.containsKey(lemma)) {
                            packs.put(lemma, new HashMap<String, String>());
                        }
                        if (lemma.length() == 0 || word.length() == 0 || msd.length() == 0) {
                            System.out.println(String.format("W:%s\tL:%s\tC:%s", word, lemma, msd));
                        }
                        packs.get(lemma).put(word, msd);
                    }
                }
            }
            br.close();
            return packs;
        }

        public static void tenfold(Map<String, Map<String, String>> packs, Map<String, Map<String, String>> train, Map<String, Map<String, String>> test) {

            List<Integer> stat = new ArrayList<Integer>();

            int i = 0;
            for (String lemma : packs.keySet()) {
                Map<String, String> pack = new HashMap<String, String>();
                for (String word : packs.get(lemma).keySet()) {
                    String ctag = packs.get(lemma).get(word);
                    if (getPrefixSize(word, lemma, ctag) < 8 && word.indexOf("_") == -1) {
                        pack.put(word, ctag);
                    }
                }
                if (pack.size() > 0) {

                    if (i % 10 == 0) {
                        test.put(lemma, pack);
                    } else {
                        train.put(lemma, pack);
                    }
                    i++;

                    for (String word : pack.keySet()) {
                        int n = getPrefixSize(word, lemma, pack.get(word));
                        while (n >= stat.size()) {
                            stat.add(0);
                        }
                        stat.set(n, stat.get(n) + 1);
                    }
                }
            }

            for (int s = 0; s < stat.size(); s++) {
                System.out.println(s + " " + stat.get(s));
            }
        }

        public static void save(Map<String, Map<String, String>> packs, String fileName, boolean full) throws IOException {
            try (BufferedWriter bw = new BufferedWriter(new FileWriter(fileName))) {
                for (String lemma : packs.keySet()) {
                    for (String word : packs.get(lemma).keySet()) {
                        bw.write(extractFeatures(word, lemma, packs.get(lemma).get(word), full) + "\n");
                    }
                }
            }
        }

        public static int getPrefixSize(String word, String lemma, String ctag) {
            int sz = 0;
            while (sz < word.length() && sz < lemma.length() && word.charAt(sz) == lemma.charAt(sz)) {
                sz++;
            }
            return word.length() - sz;
        }

        public static String extractFeatures(String word, String lemma, String msd, boolean full) {
            final int max = 7;
            List<String> feats = new ArrayList<String>();
            int s = 0;
            while (s < word.length() && s < lemma.length() && word.charAt(s) == lemma.charAt(s)) {
                s++;
            }
            StringBuilder out = new StringBuilder((word.length() - s) + "|" + lemma.substring(s));
            if (word.length() - s > max) {
                s = word.length() - max;
            }
            for (int j = word.length() - max; j < word.length(); j++) {
                out.append(" c" + (word.length() - j) + ":" + (j < 0 ? "_" : word.charAt(j)));
            }
            char c = msd.charAt(0);
            out.append(" m:" + msd);
//        if (c == 'N' ) {
//            out.append(" s:" + msd.charAt(2));
//        }
            if (full) {
                out.append(" w:" + word);
                out.append(" wl:" + lemma);
            }
            return out.toString();
        }

        public static void create() throws NumberFormatException, IOException, FileNotFoundException {
            String lAllIn = "models/ro/training_data/tbl.wordform.ro.v81";
            String msd2tagFile = "models/ro/training_data/MSDtoTAG_ro.map";
            String lTrainFile = "models/ro/training_data/lemma2.train";
            String lTestFile = "models/ro/training_data/lemma2.test";
            String lID3File = "models/ro/lemma2.id3";

            Map<String, Map<String, String>> lallin = BasicLemmatizer.Model.readAll(lAllIn, msd2tagFile);
            Map<String, Map<String, String>> ltrainin = new HashMap<String, Map<String, String>>();
            Map<String, Map<String, String>> ltestin = new HashMap<String, Map<String, String>>();

            BasicLemmatizer.Model.tenfold(lallin, ltrainin, ltestin);
            BasicLemmatizer.Model.save(ltrainin, lTrainFile, false);
            BasicLemmatizer.Model.save(ltestin, lTestFile, true);

            //Utils.buildID3(lTrainFile, lTestFile, lID3File);
            ID3.shrink(lID3File);
            ID3 id3 = ID3.createFromFile(lID3File);

            BufferedReader br = new BufferedReader(new FileReader(lTestFile));
            String line = null;
            int perr = 0;
            int werr = 0;
            int pt = 0;
            int wt = 0;
            boolean ok = true;
            PrintWriter pw = new PrintWriter("tmp.txt");
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    wt++;
                    if (!ok) {
                        werr++;
                    }
                    ok = true;
                } else {
                    ok = true;
                    pt++;
                    String[] parts = line.split(" ");
                    List<String> feats = new ArrayList<>();
                    String w = "";
                    String wl = "";
                    for (int i = 1; i < parts.length; i++) {
                        feats.add(parts[i]);
                        if (parts[i].startsWith("w:")) {
                            w = parts[i].substring(2);
                        }
                        if (parts[i].startsWith("wl:")) {
                            wl = parts[i].substring(3);
                        }
                    }
                    String out = id3.classify(feats).split(" ")[0];
                    String[] output = out.split("\\|");
                    int sz = Integer.parseInt(output[0]);
                    String rpl = "";
                    if (output.length > 1) {
                        rpl = output[1];
                    }
                    String owl = "???";
                    if (w.length() < sz) {
                        perr++;
                        ok = false;
                        pw.write(w + " " + wl + " => ??? (" + sz + "_" + rpl + ")\n");
                    } else {
                        owl = w.substring(0, w.length() - sz) + rpl;
                        if (!owl.equals(wl)) {
                            perr++;
                            ok = false;
                            pw.write(w + "\torig: (" + parts[0] + ") => " + wl + " \tcod: (" + out + ") => " + owl + "\n");
                        }
                    }
                    if (!parts[0].equals(out) && ok == true) {
                    }
                }
            }
            pw.close();
            br.close();
            if (wt == 0) {
                wt++;
            }
            System.out.println("Processed " + wt + " examples\n\twith " + werr + " errors (" + (1.0 - (float) werr / wt) + " accuracy) \n\twith " + pt + " individual states\n\t\tand " + perr + " errors (" + (1.0 - (float) perr / pt) + " accuracy)");
        }

    }

}
