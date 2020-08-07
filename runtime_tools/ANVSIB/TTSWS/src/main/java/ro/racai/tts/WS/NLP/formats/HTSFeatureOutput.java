/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.formats;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import ro.racai.tts.WS.NLP.language.Token;

/**
 *
 * @author tibi
 */
public class HTSFeatureOutput implements BaseFeatureOutput {

    @Override
    public String print(Token[] tokens) {
        StringBuilder sb = new StringBuilder();
        String[] features = getFeatures(tokens);
        for (int i = 0; i < features.length; i++) {
            //System.out.println(features[i]);
            sb.append(features[i] + "\n");
        }
        return sb.toString();
    }

    @Override
    public String[] getFeatures(Token[] tokens) {

        List<String> feats = new ArrayList<>();
        //filtram tokenii
        List<String> allPhonemes = new ArrayList<>();
        List<Token> filtered = new ArrayList<>();
        boolean paraphrase = false;
        List<Integer> p2w = new ArrayList<>();
        List<Integer> p2s = new ArrayList<>();
        List<Integer> s2sc = new ArrayList<>();
        for (int i = 0; i < tokens.length; i++) {
            Token t = tokens[i];
            if (t.word.equals("\"") || t.word.equals("``") || t.word.equals("''")) {
                paraphrase = !paraphrase;
            }
            t.paraphrase = String.valueOf(paraphrase);
            if ((!t.word.equals("``") && !t.word.equals("\"") && !t.equals("'")) && (t.word.equals(".") || t.word.equals(",") || t.word.equals("!") || t.word.equals("?") || !Pattern.matches("\\p{Punct}", t.word))) {
                filtered.add(t);
                int k = 0;
                for (int l = 0; l < t.syllables.length; l++) {
                    int pis = 0;
                    for (int ll = 0; ll < t.syllables[l].length(); ll++) {

                        if (!t.phonemes[k].equals("-")) {
                            String[] parts = t.phonemes[k].split("\\.");
                            for (int zz = 0; zz < parts.length; zz++) {
                                allPhonemes.add(parts[zz]);
                                p2w.add(filtered.size() - 1);
                                p2s.add(l);
                                pis++;
                            }
                        }
                        k++;
                    }
                    for (int zz = 0; zz < pis; zz++) {
                        s2sc.add(pis);
                    }
                }
            }
        }
        //daca lipseste punctuatia la final adaugam noi un punct
        if (allPhonemes.size() == 0) {
            return feats.toArray(new String[0]);
        }
        if (!tokens[tokens.length - 1].equals(".") && !tokens[tokens.length - 1].equals("!") && !tokens[tokens.length - 1].equals("?")) {
            Token t = new Token();
            t.word = ".";
            t.phonemes = new String[1];
            t.phonemes[0] = "pau";
            filtered.add(t);
        }
        tokens = filtered.toArray(new Token[0]);
        int numSyllablesInSentence = 0;
        int numAccentedSyllablesInSentence = 0;
        int numWordsInSentence = 0;
        for (int i = 0; i < tokens.length; i++) {
            Token t = tokens[i];
            if (!Pattern.matches("\\p{Punct}", t.word)) {
                numWordsInSentence++;
                numSyllablesInSentence += t.syllables.length;
                for (int k = 0; k < t.syllables.length; k++) {
                    if (t.stress[k] == 1) {
                        numAccentedSyllablesInSentence++;
                    }
                }
            }
        }
        String sentType = "AFIRMATIVE";
        for (int i = tokens.length - 1; i >= 0; i--) {
            if (tokens[i].word.equals("?")) {
                sentType = "INTEROGATIVE";
            } else if (tokens[i].word.equals("!")) {
                sentType = "EXCLAMATIVE";
            }
        }
        //incepe distractia

        //avem un caz special pentru pauza de la inceput
        feats.add("#^#-pau+" + allPhonemes.get(0) + "=" + allPhonemes.get(1) + "@");
        String pp = "pau";
        String ppp = "#";
        int pindex = 0;
        int pindex_sil = 0;
        int tindex = 0;
        int sindex = 0;
        int npis = 0;
        while (pindex < allPhonemes.size()) {
            Token t = tokens[tindex];

            String cp = allPhonemes.get(pindex);
            String np = "#";
            String nnp = "#";
            if (pindex + 1 < allPhonemes.size()) {
                np = allPhonemes.get(pindex + 1);
            }
            if (pindex + 2 < allPhonemes.size()) {
                nnp = allPhonemes.get(pindex + 2);
            }
            StringBuilder context = new StringBuilder();
            //context fonetic
            context.append(ppp + "^" + pp + "-" + cp + "+" + np + "=" + nnp + "@");
            if (!cp.equals("pau") && !cp.equals("sp")) {
                //numarul de silabe in cuvant
                context.append("/NSYL:" + t.syllables.length);
                //pozitia silabei in cuvand
                String siw = "single";
                if (t.syllables.length > 1) {
                    if (sindex == 0) {
                        siw = "start";
                    } else if (sindex == t.syllables.length - 1) {
                        siw = "end";
                    } else {
                        siw = "mid";
                    }
                }
                context.append("/SIW:" + siw);
                //accent lexical
                if (t.stress.length <= sindex) {
                    System.err.println("PULA MEA!!!!");
                }
                if (t.stress[sindex] == 1) {
                    context.append("/STRESS:1");
                } else {
                    context.append("/STRESS:0");
                }
                context.append("/SIL:" + t.syllables[sindex]);
                //pozitia fonemei in silaba
                context.append("/PIS:" + pindex_sil);
                //numarul de foneme din silaba
                context.append("/NPIS:" + npis);
                //POS
                context.append("/POS:" + t.tag);
                //NER
                context.append("/NER:" + t.ner);
                //sentiment
                context.append("/SENTIMENT:" + t.sentimentGroup.replace(" ", "_"));
                //paraphrase
                context.append("/PARAPHRASE:" + t.paraphrase);
                //punctuatie si dinstantele fata de ea
                PhoneticDistance pd = findPunctuation(tindex, sindex, -1, tokens);
                context.append("/PREV_PUNCT:" + pd.type);
                context.append("/PREV_PUNCT_DIST_SYLS:" + pd.syllables);
                context.append("/PREV_PUNCT_DIST_WORDS:" + pd.words);
                pd = findPunctuation(tindex, sindex, +1, tokens);
                context.append("/NEXT_PUNCT:" + pd.type);
                context.append("/NEXT_PUNCT_DIST_SYLS:" + pd.syllables);
                context.append("/NEXT_PUNCT_DIST_WORDS:" + pd.words);
                pd = findAccentedSyllable(tindex, sindex, -1, tokens);
                context.append("/PREV_ACC_SYL_DIST_SYLS:" + pd.syllables);
                context.append("/PREV_ACC_SYL_DIST_WORDS:" + pd.words);
                pd = findAccentedSyllable(tindex, sindex, +1, tokens);
                context.append("/NEXT_ACC_SYL_DIST_SYLS:" + pd.syllables);
                context.append("/NEXT_ACC_SYL_DIST_WORDS:" + pd.words);
                //parsing

            }
            context.append("/TSYL:" + numSyllablesInSentence);
            context.append("/SENTTYPE:" + sentType);
            context.append("/TWORDS:" + tokens.length);
            feats.add(context.toString().trim());
            ppp = pp;
            pp = cp;
            pindex++;

            //de recalculat tindex sindex npis pindex_sil
            if (pindex < allPhonemes.size()) {
                pindex_sil++;
                int ltindex = tindex;
                tindex = p2w.get(pindex);
                int lsindex = sindex;
                sindex = p2s.get(pindex);
                if (lsindex != sindex || ltindex != tindex) {
                    pindex_sil = 0;
                }
                npis = s2sc.get(pindex);
            }
        }

        return feats.toArray(
                new String[0]);
    }

    private PhoneticDistance findPunctuation(int tindex, int sindex, int dir, Token[] tokens) {
        PhoneticDistance pd = new PhoneticDistance();
        pd.type = "sent_start";
        if (dir < 0) {
            pd.syllables = sindex;
        } else {
            pd.syllables = tokens[tindex].syllables.length - sindex - 1;
        }
        tindex += dir;
        while (tindex >= 0 && tindex < tokens.length) {
            if (Pattern.matches("\\p{Punct}", tokens[tindex].word)) {
                //tipul
                switch (tokens[tindex].word) {
                    case "\"":
                        pd.type = "paraph";
                        break;
                    case "``":
                        pd.type = "paraph";
                        break;
                    case "''":
                        pd.type = "paraph";
                        break;
                    case ".":
                        pd.type = "period";
                        break;
                    case "!":
                        pd.type = "excl";
                        break;
                    case "?":
                        pd.type = "quest";
                        break;
                    case ",":
                        pd.type = "comma";
                        break;
                    case "-":
                        pd.type = "dash";
                        break;
                    default:
                        pd.type = "other";
                        break;

                }

                break;
            }
            pd.words++;
            pd.syllables += tokens[tindex].syllables.length;
            tindex += dir;
        }
        return pd;
    }

    private PhoneticDistance findAccentedSyllable(int tindex, int sindex, int dir, Token[] tokens) {
        PhoneticDistance pd = new PhoneticDistance();
        boolean found = false;
        while (tindex >= 0 && tindex < tokens.length) {
            if (Pattern.matches("\\p{Punct}", tokens[tindex].word)) {
                break;
            }
            if (tokens[tindex].stress[sindex] == 1) {
                found = true;
                break;
            }
            sindex += dir;
            pd.syllables++;
            if (sindex < 0) {
                tindex--;
                if (tindex < 0) {
                    break;
                }
                pd.words++;
                if (!Pattern.matches("\\p{Punct}", tokens[tindex].word)) {
                    sindex = tokens[tindex].syllables.length - 1;
                }
            } else if (sindex >= tokens[tindex].syllables.length) {
                tindex++;
                pd.words++;
                sindex = 0;
            }
        }
        if (!found) {
            pd.words = -1;
            pd.syllables = -1;
        }
        return pd;
    }

    class PhoneticDistance {

        int words;
        int syllables;
        String type;
    }

    class Context {

        String phoneme;
        int windex;
        int sindex;
        int sstress;
        int numPhones = 0;
    }

}
