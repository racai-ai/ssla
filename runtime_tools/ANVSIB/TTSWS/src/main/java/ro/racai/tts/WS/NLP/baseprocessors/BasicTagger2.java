/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.baseprocessors;


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import ro.racai.tts.WS.NLP.language.BaseProcessor;
import static ro.racai.tts.WS.NLP.language.BaseProcessor.PROCESSOR_TAGGER;
import ro.racai.tts.WS.NLP.language.Token;
import ro.racai.tts.WS.NLP.machinelearning.ID3;

/**
 * Implements standard POS tagging functionality. Romanian accuracy is 94.19%
 *
 * @author tibi
 */
public class BasicTagger2 implements BaseProcessor {

    @Override
    public int getProcessorType() {
        return PROCESSOR_TAGGER;
    }

    @Override
    public void processTokens(Token[] tokens) {
        List<String> sentWords = new ArrayList<String>();
        List<String> id3Results = new ArrayList<String>();
        for (int i = 0; i < tokens.length; i++) {
            sentWords.add(tokens[i].word);
            id3Results.add(id3.classify(WordFeatures.create(tokens[i].word)));
        }
        Viterbi vb = new Viterbi(sentWords, id3Results, tranProb, tagCounts, wordTagsCounts, allTagsCount);
        String[] res = vb.run();
        for (int i = 0; i < tokens.length; i++) {
            tokens[i].tag = res[i];
        }
    }
    
    MapProbs tranProb;
    Map<String, Integer> tagCounts;
    int allTagsCount;
    Map<String, Map<String, Integer>> wordTagsCounts;
    ID3 id3;

    @Override
    public void loadModel(String folder) {
        try {
            // tranProb
            tranProb = MapProbs.createFromFile(folder + "/trans.prob");

            // tag counts & allTagCounts
            tagCounts = new HashMap<String, Integer>();
            allTagsCount = 0;
            BufferedReader br = new BufferedReader(new FileReader(folder+"/tagCounts.txt"));
            String line;
            while ((line = br.readLine()) != null) {
                String[] toks = line.split("[ \t]");
                if (toks.length == 2) {
                    int count = Integer.parseInt(toks[1]);
                    tagCounts.put(toks[0], count);
                    allTagsCount += count;
                }
            }

            //wordTagsCounts
            wordTagsCounts = new HashMap<String, Map<String, Integer>>();
            br = new BufferedReader(new FileReader(folder+"/wordTagsCounts.txt"));
            while ((line = br.readLine()) != null) {
                String[] toks = line.split("[ \t]");
                if (toks.length > 1) {
                    Map<String, Integer> tags = new HashMap<String, Integer>();
                    for (int i = 1; i < toks.length; i += 2) {
                        int count = Integer.parseInt(toks[i + 1]);
                        tags.put(toks[i], count);
                        allTagsCount += count;
                    }
                    wordTagsCounts.put(toks[0], tags);
                }
            }

            //id3
            id3 = ID3.createFromFile(folder + "/tag.7.id3");
        } catch (IOException ex) {
            Logger.getLogger(BasicLTS.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    // Probabilitati de tranzitie
    private static class MapProbs extends HashMap<String, Map<String, Double>> {

        public void add(String global, String part) {
            if (!containsKey(global)) {
                put(global, new HashMap<String, Double>());
            }
            Map<String, Double> parts = get(global);
            double value = 0;
            if (parts.containsKey(part)) {
                value = parts.get(part);
            }
            parts.put(part, value + 1);
        }

        public void normalize() {
            for (Map<String, Double> parts : this.values()) {
                double sum = 0;
                for (double value : parts.values()) {
                    sum += value;
                }
                for (String key : parts.keySet()) {
                    parts.put(key, parts.get(key) / sum);
                }
            }
        }

        public void saveToFile(String fileName) throws FileNotFoundException {
            PrintWriter pw = new PrintWriter(fileName);
            for (String glob : keySet()) {
                pw.print(glob);
                Map<String, Double> parts = get(glob);
                for (String part : parts.keySet()) {
                    pw.print(" " + part + " " + parts.get(part));
                }
                pw.println();
            }
            pw.close();
        }

        public static MapProbs createFromFile(String fileName) throws FileNotFoundException, IOException {
            MapProbs res = new MapProbs();
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line;
            while ((line = br.readLine()) != null) {
                String[] tokens = line.split("[ \t]");
                if (tokens.length > 0) {
                    String glob = tokens[0];
                    Map<String, Double> parts = new HashMap<String, Double>();
                    for (int i = 1; i < tokens.length; i += 2) {
                        parts.put(tokens[i], Double.parseDouble(tokens[i + 1]));
                    }
                    res.put(glob, parts);
                }
            }
            return res;
        }

        public static MapProbs computeTransitions(List<List<String[]>> posData) {
            MapProbs res = new MapProbs();
            for (List<String[]> data : posData) {
                res.add("_", data.get(0)[1]);
                for (int i = 0; i < data.size() - 1; i++) {
                    res.add(data.get(i)[1], data.get(i + 1)[1]);
                }
                res.add(data.get(data.size() - 1)[1], "_");
            }
            res.normalize();
            return res;
        }

        public static MapProbs computeTags(List<List<String[]>> posData) {
            MapProbs res = new MapProbs();
            for (List<String[]> data : posData) {
                for (String[] pair : data) {
                    res.add(pair[0], pair[1]);
                }
            }
            res.normalize();
            return res;
        }
    }

    private static class Viterbi {

        List<String> words;
        List<Map<String, Double>> tags;
        Map<String, Map<String, Double>> tranProbs;
        Map<String, Integer> tagsCounts;
        Map<String, Map<String, Integer>> wordTagsCount;
        double total;

        public Viterbi(List<String> words, List<String> id3Results, Map<String, Map<String, Double>> tranProbs, Map<String, Integer> tagsCounts, Map<String, Map<String, Integer>> wordTagsCount, int allTagCounts) {
            this.words = words;
            tags = new ArrayList<Map<String, Double>>();
            this.tranProbs = tranProbs;
            this.tagsCounts = tagsCounts;
            this.wordTagsCount = wordTagsCount;
            total = allTagCounts;
            for (String result : id3Results) {
                String[] toks = result.split("[ \t]");
                Map<String, Double> res = new HashMap<String, Double>();
                for (int i = 0; i < toks.length; i += 2) {
                    res.put(toks[i], Double.parseDouble(toks[i + 1]));
                }
                tags.add(res);
            }
        }

        public String[] run() {
            String[] res = new String[words.size()];

            List<List<Token>> tokens = new ArrayList<List<Token>>();

            for (int i = 0; i < res.length; i++) {
                List<Token> tmp = new ArrayList<>();
                if (!wordTagsCount.containsKey(words.get(i))) {
                    for (String tag : tags.get(i).keySet()) {
                        Token t = new Token();
                        t.tag = tag;
                        if (tagsCounts.containsKey(tag)) {
                            t.em = 1D / (tagsCounts.get(tag) + 1D);
                        } else {
                            t.em = 1D / (total + 1D);
                        }
                        if (t.em == 0) {
                            t.em = 0;
                        }
                        t.score = Math.log(t.em);
                        t.maxFrom = -1;
                        tmp.add(t);
                    }
                } else {
                    for (String tag : wordTagsCount.get(words.get(i)).keySet()) {
                        Token t = new Token();
                        t.tag = tag;
                        if (wordTagsCount.get(words.get(i)).containsKey(tag)) {
                            t.em = (float) wordTagsCount.get(words.get(i)).get(tag) / (tagsCounts.get(tag) + 1D);
                        } else {
                            t.em = 1D / (tagsCounts.get(tag) + 1D);
                        }
                        if (t.em == 0) {
                            t.em = 0;
                        }
                        t.maxFrom = -1;
                        t.score = Math.log(t.em);
                        tmp.add(t);
                    }
                }
                tokens.add(tmp);
            }

            for (int k = 1; k < res.length; k++) {
                List<Token> prev = tokens.get(k - 1);
                List<Token> cur = tokens.get(k);
                for (int i = 0; i < cur.size(); i++) {
                    double max = Double.NEGATIVE_INFINITY;
                    for (int j = 0; j < prev.size(); j++) {
                        String tagp = prev.get(j).tag;
                        String tagn = cur.get(i).tag;
                        double probBi = 1 / (total + 1);
                        if (tranProbs.containsKey(tagp) && tranProbs.get(tagp).containsKey(tagn)) {
                            probBi = tranProbs.get(tagp).get(tagn);
                        }
                        double probUni = 1 / (total + 1);
                        if (tagsCounts.containsKey(tagn)) {
                            probUni = (tagsCounts.get(tagn) + 1) / (total + 1);
                        }
                        double prob = Math.log(probBi);// Math.log(0.5 * probBi + 0.5 * probUni);
                        double score = prev.get(j).score + prob + Math.log(cur.get(i).em);
                        if (score == Double.NEGATIVE_INFINITY) {
                            score = 0;
                        }
                        if (score > max) {
                            max = score;
                            cur.get(i).maxFrom = j;
                            cur.get(i).score = max;
                        }
                    }
                }
            }

            //cautam maxim
            List<Token> last = tokens.get(tokens.size() - 1);
            int maxFrom = 0;
            for (int i = 1; i < last.size(); i++) {
                if (last.get(i).score > last.get(maxFrom).score) {
                    maxFrom = i;
                }
            }
            int pos = res.length - 1;
            while (pos >= 0) {
                res[pos] = last.get(maxFrom).tag;
                maxFrom = last.get(maxFrom).maxFrom;
                pos--;
                if (pos >= 0) {
                    last = tokens.get(pos);
                }
            }

            return res;
        }

        static class Token {

            public String tag;
            public double em;
            public double score;
            public int maxFrom;
        }
    }

    private static class WordFeatures {

        static int sz = 7;

        public static List<String> create(String word) {
            List<String> feats = new ArrayList<String>();

            for (int i = 0; i < sz; i++) {
                if (word.length() >= (i + 1)) {
                    feats.add("e" + i + ":" + word.charAt(word.length() - (i + 1)));
                } else {
                    feats.add("e" + i + ":_");
                }
            }
            return feats;
        }

        private static String compact(List<String> feats) {
            StringBuilder sb = new StringBuilder();
            for (String feat : feats) {
                if (sb.length() != 0) {
                    sb.append(" ");
                }
                sb.append(feat);
            }
            return sb.toString();

        }

        public static void saveToFile(List<String> feats, String fileName) throws FileNotFoundException {
            PrintWriter pw = new PrintWriter(fileName);
            for (String feat : feats) {
                pw.println(feat);
            }
            pw.close();
        }
    }
}
