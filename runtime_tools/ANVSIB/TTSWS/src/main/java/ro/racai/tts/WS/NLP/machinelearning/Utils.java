/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.machinelearning;


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author adrian
 */
public class Utils {

    static class FScore {

        int truePositive = 0, falsePositive = 0;
        int positives = 0;

        public void incPositives() {
            positives++;
        }

        public void incTruePositive() {
            truePositive++;
        }

        ;
        public void incFalsePositive() {
            falsePositive++;
        }

        ;
        
        public double getValue() {
            double p;
            if (truePositive + falsePositive == 0) {
                p = 1;
            } else {
                p = (double) truePositive / (truePositive + falsePositive);
            }

            double r = (double) truePositive / positives;
            if (p + r == 0) {
                return 0;
            }
            return 2 * p * r / (p + r);
        }
    }

    public static void buildID3(String train, String test, String output) throws IOException {
        ID3 id3 = ID3.trainFromFile(train, 0, true);
        testID3(id3, test);
        id3.saveModel(new PrintStream(output));
    }

    public static void buildID3(String train, String test, String output, String target) throws IOException {
        ID3 id3 = ID3.trainFromFile(train, 0, true);
        testID3(id3, test, target);
        id3.saveModel(new PrintStream(output));
    }

    public static void testID3(ID3 id3, String fileName) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        int perr = 0;
        int werr = 0;
        int pt = 0;
        int wt = 0;
        boolean ok = true;
        //PrintWriter pw = new PrintWriter("tmp.txt");
        while ((line = br.readLine()) != null) {
            if (line.trim().isEmpty()) {
                wt++;
                if (!ok) {
                    werr++;
                }
                ok = true;
            } else {
                pt++;
                String[] parts = line.split("[ \t]");
                List<String> feats = new ArrayList<>();
                for (int i = 1; i < parts.length; i++) {
                    feats.add(parts[i]);
                }
                String output = id3.classify(feats);
                if (!output.split(" ")[0].equals(parts[0])) {
                    perr++;
                    ok = false;
                }
                //pw.write(output.split(" ")[0] + "\n");
            }
        }
        //pw.close();
        br.close();
        if (wt == 0) {
            wt++;
        }
        System.out.println("Processed " + wt + " examples\n\twith " + werr + " errors (" + (1.0 - (float) werr / wt) + " accuracy) \n\twith " + pt + " individual states\n\t\tand " + perr + " errors (" + (1.0 - (float) perr / pt) + " accuracy)");

//        System.out.println("Processed " + fileName);
//        System.out.println("\t" + pt + " examples with " + perr + " errors (" + (1.0 - (float) perr / pt) + " accuracy)");
//        System.out.println("\t" + wt + " individual states with " + werr + " errors (" + (1.0 - (float) werr / wt) + " accuracy)");
    }

    public static void testID3(ID3 id3, String testFile, String target) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(testFile));
        String line = null;
        int perr = 0;
        int werr = 0;
        int pt = 0;
        int wt = 0;
        boolean ok = true;
        FScore fscore = new FScore();
        //PrintWriter pw = new PrintWriter("tmp.txt");
        while ((line = br.readLine()) != null) {
            if (line.trim().isEmpty()) {
                wt++;
                if (!ok) {
                    werr++;
                }
                ok = true;
            } else {
                pt++;
                String[] parts = line.split("[ \t]");
                List<String> feats = new ArrayList<>();
                for (int i = 1; i < parts.length; i++) {
                    feats.add(parts[i]);
                }
                String output = id3.classify(feats).split(" ")[0];
                if (!output.equals(parts[0])) {
                    perr++;
                    ok = false;
                }
                if (parts[0].equals(target)) {
                    fscore.incPositives();
                    if (output.equals(target)) {
                        fscore.incTruePositive();
                    } else {
                        fscore.incFalsePositive();
                    }
                }
                //pw.write(output + "\n");
            }
        }
        //pw.close();
        br.close();
        if (wt == 0) {
            wt++;
        }
//        System.out.println("Processed " + wt + " examples");
//        System.out.println("\twith " + werr + " errors (" + (1.0 - (float) werr / wt) + " accuracy)");
//        System.out.println("\twith " + pt + " individual states\n\t\tand " + perr + " errors (" + (1.0 - (float) perr / pt) + " accuracy)");
//        System.out.println("\tfscore = " + fscore.getValue());

        System.out.println("Processed " + testFile);
        System.out.println("\t" + pt + " examples with " + perr + " errors (" + (1.0 - (float) perr / pt) + " accuracy)");
        System.out.println("\t" + wt + " individual states with " + werr + " errors (" + (1.0 - (float) werr / wt) + " accuracy)");
    }

    public static boolean firstCaps(String s) {
        if (s.length() < 2) {
            return allCaps(s);
        }
        return s.toLowerCase().substring(1).equals(s.substring(1)) && Character.isUpperCase(s.charAt(0));
    }

    public static boolean allCaps(String s) {
        return s.toUpperCase().equals(s);
    }

    public static char firstChar(String s) {
        if (s.length() == 0) {
            return '_';
        }
        return s.charAt(0);
    }

    public static char lastChar(String s) {
        if (s.length() == 0) {
            return '_';
        }
        return s.charAt(s.length() - 1);
    }

    public static boolean hasLastPunctuation(String s) {
        if (s.length() == 0) {
            return false;
        }
        return isPunctuation(s.charAt(s.length() - 1));
    }

    public static boolean hasOnlyPunctuation(String s) {
        for (int i = 0; i < s.length(); i++) {
            if (!isPunctuation(s.charAt(i))) {
                return false;
            }
        }
        return true;
    }

    public static boolean hasPunctuation(String s) {
        for (int i = 0; i < s.length(); i++) {
            if (isPunctuation(s.charAt(i))) {
                return true;
            }
        }
        return false;
    }

    //private static List<String> punctuations = Arrays.asList(".", "'", "\"", "\\", "/", ",", "-", ":", ";", "`", "!", "?", "...", "(", ")", "[", "]", "{", "}", "<", ">", "@", "_","„","”");
    public static boolean isPunctuation(char s) {
        return !Character.isLetterOrDigit(s);
    }

    public static boolean hasManyCaps(String s) {
        int count = 0;
        for (int i = 0; i < s.length() && count < 2; i++) {
            if (("" + s.charAt(i)).toUpperCase().equals("" + s.charAt(i))) {
                count++;
            }
        }
        return count >= 2;
    }

    public static boolean isNumeric(String str) {
        String s = str.replaceAll("-", "").replaceAll(".", "").replaceAll(",", "");
        if (s.length() == 0) {
            return false;
        }
        for (char c : s.toCharArray()) {
            if (!Character.isDigit(c)) {
                return false;
            }
        }
        return true;
    }

    public static String[] simpleWhiteSpaceSplit(String sentence) {
        ArrayList<String> tokens = new ArrayList<String>();
        boolean isPunct;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < sentence.length(); i++) {
            isPunct = isPunctuation(sentence.charAt(i));

            if (isPunct) {
                // scrie ce era in buffer
                if (!sb.toString().equals("")) {
                    tokens.add(sb.toString());
                    //IO.outln(sb.toString());
                }
                // creaza un buffer nou pentru punctuatia curenta, daca nu e whitespace
                if (!Character.isWhitespace(sentence.charAt(i))) {
                    sb = new StringBuilder();
                    sb.append(sentence.charAt(i));
                    tokens.add(sb.toString());
                    //IO.outln(sb.toString());                
                }

                // creaza un buffer nou pentru urmatorul caracter
                sb = new StringBuilder();
            } else {
                sb.append(sentence.charAt(i));
            }
        }
        if (!sb.toString().equals("")) {
            tokens.add(sb.toString());
        }
        //IO.outln(sb.toString());;

        return tokens.toArray(new String[0]);
    }
}
