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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import ro.racai.tts.WS.NLP.language.BaseProcessor;
import ro.racai.tts.WS.NLP.language.BaseProcessor;
import ro.racai.tts.WS.NLP.language.Token;

/**
 * Implements standard chunking functionality
 *
 * @author tibi
 */
public class BasicChunker implements BaseProcessor {

    private static enum Card {

        MAX_ONE, ONLY_ONE, AT_LEAST_ONE, ANY;
    }

    private static class Rule {

        String name;
        String item = null;
        boolean isTerminal = false;
        List<Rule> orRules = null; // or rules
        List<Rule> andRules = null; // or rules
        Card card = Card.ONLY_ONE;

        public int apply(Token[] tokens, int pos, Map<String, Rule> rules) {
            if (pos>=tokens.length)
                return 0;
            if (item != null) {
                if (isTerminal) {
                    if (item.equals(tokens[pos].tag)) {
                        return 1;
                    } else {
                        return 0;
                    }

                } else {
                    Rule rule = rules.get(item);
                    return rule.apply(tokens, pos, rules);
                }
            }
            //avem doua variante or si and
            int pp = pos;
            if (andRules != null) {
                for (int i = 0; i < andRules.size(); i++) {
                    Rule r = andRules.get(i);
                    int count = r.apply(tokens, pp, rules);
                    //ANY si MAX_ONE
                    if (count == 0 && r.card != Card.ANY && r.card != Card.MAX_ONE) {
                        return 0;
                    }
                    pp += count;
                }
                return pp - pos;
            } else {

                //se pare ca daca avem or nu avem si and
                int count = 0;
                boolean ok = false;
                for (int i = 0; i < orRules.size(); i++) {
                    Rule r = orRules.get(i);
                    //Adrian, astea nu ar trebui sa aiba cardinalitate, nu?
                    count = r.apply(tokens, pp, rules);
                    if (count != 0) {
                        break;
                    }
                }
                pp += count;
                return pp - pos;
            }

            //return 0;
        }

        public int apply2(Token[] tokens, int pos, Map<String, Rule> rules) {

            if (pos == -1) {
                return -1;
            }

            if (tokens.length <= pos) {
                return -1;
            }

            if (name != null) {
                System.out.println(pos + " " + tokens[pos].word + " " + name);
            }

            int result = 0;
            if (item != null) {
                if (isTerminal) {
                    if (item.equals(tokens[pos].tag)) {
                        tokens[pos].chunk = item;
                        System.out.println("OK: " + pos + " " + tokens[pos].word + " " + item);
                        result = pos;
                    } else {
                        result = -1;
                    }
                    return result;
                } else {
                    Rule rule = rules.get(item);
                    result = rule.apply(tokens, pos, rules);
                }
            } else {

                if (orRules != null) {
                    result = -1;
                    for (int i = 0; i < orRules.size() && result == -1; i++) {
                        result = orRules.get(i).apply(tokens, pos, rules);
                    }

                } else if (andRules != null) {
                    result = pos;

                    for (int i = 0; i < andRules.size() && result != -1; i++) {
                        int oldResult = result;
                        result = andRules.get(i).apply(tokens, result, rules);
                        boolean last = (i + 1 == andRules.size());
                        if (result == -1 && andRules.get(i).card == Card.ANY) {
                            result = oldResult;
                        } else if (result != -1 && andRules.get(i).card == Card.ANY) {
                            i--;
                            if (!last) {
                                result++;
                            }
                        } else if (result != -1 && !last) {
                            result++;
                        }
                    }

                    if (result == pos) {
                        result = -1;
                    }
                }
            }

            if (result == -1) {
                //curatare nepotrivire
                for (int i = pos; i < tokens.length; i++) {
                    tokens[i].chunk = null;
                }
            }
            return result;
        }

        public static Rule createOrRule(List<String> tokens) {
            if (!tokens.isEmpty()) {

                int deep;
                if (tokens.get(0).equals("(") && tokens.get(tokens.size() - 1).startsWith(")")) {
                    // eliminare paranteze din capete
                    deep = -1;
                    boolean canDelete = true;
                    for (int i = 1; i < tokens.size() - 1 && canDelete; i++) {
                        if (tokens.get(i).equals("(")) {
                            deep--;
                        }
                        if (tokens.get(i).startsWith(")")) {
                            deep++;
                        }
                        if (deep == 0) {
                            canDelete = false;
                        }
                    }
                    if (canDelete) {
                        tokens.remove(0);
                        String card = tokens.get(tokens.size() - 1).substring(1);
                        tokens.remove(tokens.size() - 1);
                        Rule rule = createOrRule(tokens);
                        if (card.length() == 0) {
                            return rule;
                        }
                        Rule overRule = new Rule();
                        overRule.orRules = new ArrayList<Rule>();
                        overRule.orRules.add(rule);
                        if (card.equals("+")) {
                            overRule.card = Card.AT_LEAST_ONE;
                        }
                        if (card.equals("?")) {
                            overRule.card = Card.MAX_ONE;
                        }
                        if (card.equals("*")) {
                            overRule.card = Card.ANY;
                        }
                        return overRule;
                    }
                }

                // verificare or
                deep = 0;
                List<Integer> orPos = new ArrayList<Integer>();
                for (int i = 0; i < tokens.size(); i++) {
                    if (tokens.get(i).startsWith("(")) {
                        deep--;
                    }
                    if (tokens.get(i).startsWith(")")) {
                        deep++;
                    }
                    if (deep == 0 && tokens.get(i).equals("|")) {
                        orPos.add(i);
                    }
                }
                if (orPos.size() == 0) {
                    return createAndRule(tokens);
                }

                Rule rule = new Rule();
                rule.orRules = new ArrayList<Rule>();
                rule.card = Card.ONLY_ONE;

                int firstPos, lastPos;
                for (int i = 0; i <= orPos.size(); i++) {
                    if (i == 0) {
                        firstPos = 0;
                        lastPos = orPos.get(i) - 1;
                    } else if (i == orPos.size()) {
                        firstPos = orPos.get(i - 1) + 1;
                        lastPos = tokens.size() - 1;
                    } else {
                        firstPos = orPos.get(i - 1) + 1;
                        lastPos = orPos.get(i) - 1;
                    }
                    List<String> subRuleTokens = new ArrayList<String>();
                    for (int t = firstPos; t <= lastPos; t++) {
                        subRuleTokens.add(tokens.get(t));
                    }
                    rule.orRules.add(createOrRule(subRuleTokens));
                }
                return rule;
            }
            return null;
        }

        public static Rule createAndRule(List<String> tokens) {
            // Sigur nu contine paranteze pe capete
            if (tokens.size() == 1) {
                Rule rule = new Rule();

                if (tokens.get(0).endsWith("*")) {
                    rule.item = tokens.get(0).substring(0, tokens.get(0).length() - 1);
                    rule.card = Card.ANY;
                } else if (tokens.get(0).endsWith("+")) {
                    rule.item = tokens.get(0).substring(0, tokens.get(0).length() - 1);
                    rule.card = Card.AT_LEAST_ONE;
                } else if (tokens.get(0).endsWith("?")) {
                    rule.item = tokens.get(0).substring(0, tokens.get(0).length() - 1);
                    rule.card = Card.MAX_ONE;
                } else {
                    rule.item = tokens.get(0);
                    rule.card = Card.ONLY_ONE;
                }

                if (rule.item.startsWith("'<") && rule.item.endsWith(">'")) {
                    rule.item = rule.item.substring(2, rule.item.length() - 2);
                    rule.isTerminal = true;
                } else {
                    rule.isTerminal = false;
                }
//                if(rule.item.startsWith("TP")){
//                    rule=rule;
//                }

                return rule;
            }

            Rule rule = new Rule();
            rule.andRules = new ArrayList<Rule>();
            while (tokens.size() != 0) {
                if (!tokens.get(0).equals("(")) {
                    rule.andRules.add(createAndRule(Arrays.asList(tokens.get(0))));
                    tokens.remove(0);
                } else {
                    int deep = -1;
                    int pos = 1;
                    for (int i = 1; i < tokens.size() && deep != 0; i++) {
                        if (tokens.get(i).startsWith("(")) {
                            deep--;
                        }
                        if (tokens.get(i).startsWith(")")) {
                            deep++;
                        }
                        if (deep == 0) {
                            pos = i;
                        }
                    }
                    List<String> subTokensList = new ArrayList<String>();
                    for (int i = 0; i <= pos; i++) {
                        subTokensList.add(tokens.get(0));
                        tokens.remove(0);
                    }
                    rule.andRules.add(createOrRule(subTokensList));
                }
            }
            return rule;
        }
    }

    Set<String> nonterms;
    Set<String> meta;
    Set<String> start;
    Set<String> pair;
    Set<String> grules;
    Map<String, Rule> rules;// Pair(neterminal, rule description)

    @Override
    public int getProcessorType() {
        return PROCESSOR_CHUNKER;
    }

    @Override
    public void processTokens(Token[] tokens) {
        //chunking aici
        for (int i = 0; i < tokens.length; i++) {
            tokens[i].chunk = "";
        }
        Rule rule = rules.get("Np");
        int pos = 0;
        int np = 0;
        while (pos < tokens.length - 1) {
            int count = rule.apply(tokens, pos, rules);
            if (count > 1) {
                np++;
                for (int i = pos; i < pos + count; i++) {
                    tokens[i].chunk = "Np" + np;
                }
            }
            pos += count + 1;
        }

        rule = rules.get("Pp");
        pos = 0;
        np = 0;
        while (pos < tokens.length - 1) {
            int count = rule.apply(tokens, pos, rules);
            if (count > 1) {
                np++;
                for (int i = pos; i < pos + count; i++) {
                    tokens[i].chunk = tokens[i].chunk + " Pp" + np;
                }
            }
            pos += count + 1;
        }
        rule = rules.get("Ap");
        pos = 0;
        np = 0;
        while (pos < tokens.length - 1) {
            int count = rule.apply(tokens, pos, rules);
            if (count > 1) {
                np++;
                for (int i = pos; i < pos + count; i++) {
                    tokens[i].chunk = tokens[i].chunk + " Ap" + np;
                }
            }
            pos += count + 1;
        }
        rule = rules.get("Vp");
        pos = 0;
        np = 0;
        while (pos < tokens.length - 1) {
            int count = rule.apply(tokens, pos, rules);
            if (count > 1) {
                np++;
                for (int i = pos; i < pos + count; i++) {
                    tokens[i].chunk = tokens[i].chunk + " Vp" + np;
                }
            }
            pos += count + 1;
        }
        for (int i = 0; i < tokens.length; i++) {
            tokens[i].chunk = tokens[i].chunk.trim().replace(" ", ",");
        }

    }

    @Override
    public void loadModel(String folder) {
        //pe asta trebuie sa-l incarci
        String rules_file_name = folder + "/chunk.rgx";

        try {
            nonterms = null;
            meta = null;
            start = null;
            pair = null;
            rules = null;

            BufferedReader br = new BufferedReader(new FileReader(rules_file_name));
            List<String> gram = new ArrayList<String>();
            String line;
            boolean isRule = false;
            Set<String> set = null;
            while ((line = br.readLine()) != null) {
                if (line.length() != 0 && !line.trim().startsWith("#")) { // nu e comentariu
                    String[] tokens = line.split("[ \t]");
                    if (!isRule) {
                        // Incarcare multimi
                        for (int i = 0; i < tokens.length; i++) {
                            String token = tokens[i].trim();
                            if (token.length() != 0) {
                                if (token.equals("NONTERM:")) {
                                    set = nonterms = new HashSet<String>();
                                } else if (token.equals("META:")) {
                                    set = meta = new HashSet<String>();
                                } else if (token.equals("STARTSYM:")) {
                                    set = start = new HashSet<String>();
                                } else if (token.equals("PAIR:")) {
                                    set = pair = new HashSet<String>();
                                } else if (token.equals("RULES:")) {
                                    set = null;
                                    isRule = true;
                                    rules = new HashMap<String, Rule>();
                                } else {
                                    set.add(token);
                                }
                            }
                        }
                    } else { // incarcare regula
                        if (tokens.length < 2 || !tokens[1].trim().equals("->")) {
                            System.err.println("ERR:\tWRONG RULE: " + line);
                        } else {
                            String name = tokens[0].trim();
                            List<String> ruleTokens = new ArrayList<String>();
                            for (int i = 2; i < tokens.length; i++) {
                                ruleTokens.add(tokens[i].trim());
                            }
                            Rule rule = Rule.createOrRule(ruleTokens);
                            rule.name = tokens[0].trim();
                            rules.put(rule.name, rule);
                        }
                    }
                }
            }
            br.close();

        } catch (FileNotFoundException ex) {
            Logger.getLogger(BasicChunker.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(BasicChunker.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
