/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.language;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import ro.racai.tts.WS.NLP.utils.Log;

/**
 *
 * @author tibi
 */
public class BasicTokenizer implements BaseInputProcessor {

    abstract class Filter {
        public abstract List<String> match(String string);
    }
    
    List<Pattern> patterns;
    List<Filter> filters;
    Filter eosFilter = new Filter() {

        @Override
        public List<String> match(String string) {
            // patterns
            Matcher m = Pattern.compile(".*[\\.!?]$").matcher(string);
            if (m.matches()) {
                boolean ready = false;
                int i = 0;
                while (i < string.length() && !ready) {
                    char ch = string.charAt(i);
                    if (ch == '.' || ch == '!' || ch == '?') {
                        ready = true;
                    } else {
                        i++;
                    }
                }

                List<String> result = new ArrayList<String>();
                if (i == 0) {
                    result.add(string);
                } else {
                    result.add(string.substring(0, i));
                    result.add(string.substring(i));
                }
                return result;
            }
            return null;
        }
    };

    public BasicTokenizer() {


//        List<String> regexs = Arrays.asList(
//                "^/\\(?([0-9]{3})\\)?([ .-]?)([0-9]{3})\\2([0-9]{4})/$", // phone
//                "^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", // hexa
//                "^[_A-Za-z0-9-]+(\\\\.[_A-Za-z0-9-]+)*@[A-Za-z0-9]+(\\\\.[A-Za-z0-9]+)*(\\\\.[A-Za-z]{2,})$", // email                
//                "^([01]?\\\\d\\\\d?|2[0-4]\\\\d|25[0-5])\\\\.([01]?\\\\d\\\\d?|2[0-4]\\\\d|25[0-5])\\\\([01]?\\\\d\\\\d?|2[0-4]\\\\d|25[0-5])\\\\.([01]?\\\\d\\\\d?|2[0-4]\\\\d|25[0-5])$", // ip
//                "^(1[012]|[1-9]):[0-5][0-9](\\\\s)?(?i)(am|pm)$", // time 12                
//                "^([01]?[0-9]|2[0-3]):[0-5][0-9]$", // time 24
//                "^(0?[1-9]|[12][0-9]|3[01])/(0?[1-9]|1[012])/((19|20)\\\\d\\\\d)$", // date                
//                "^<(\"[^\"]*\"|'[^']*'|[^'\">])*>$", // html tag                
//                "^\\s*(?i)href\\s*=\\s*(\\\"([^\"]*\\\")|'[^']*'|([^'\">\\s]+));$" // link                
//        );        

        filters = new ArrayList<Filter>();

        filters.add(new Filter() {

            @Override
            public List<String> match(String string) {
                Matcher m = Pattern.compile("^etc\\.?$").matcher(string);
                if (m.matches()) {
                    return Arrays.asList(string);
                }
                return null;
            }
        });
        
        filters.add(new Filter() {

            @Override
            public List<String> match(String string) {
                Matcher m = Pattern.compile("^[a-zăâîșțA-ZĂÂÎȘȚ]*$").matcher(string);
                if (m.matches()) {
                    return Arrays.asList(string);
                }
                return null;
            }
        });
        
        filters.add(new Filter() {
            // ultimul caracter este , : - sau –
            @Override
            public List<String> match(String string) {
                char ch=string.charAt(string.length()-1);
                if(ch==',' || ch==':' || ch=='-' || ch=='–' || ch==';'){
                    if(string.length()==1)
                        return Arrays.asList(""+ch);                    
                    else
                        return Arrays.asList(string.substring(0,string.length()-1),""+ch);                    
                }
                return null;
            }
        });
        
        filters.add(new Filter() {
            // contine -
            @Override
            public List<String> match(String string) {
                
                if(string.contains("-")){
                    int pos=string.lastIndexOf("-");
                    if(pos==0){
                        return Arrays.asList(string);                    
                    } else {
                        return Arrays.asList(string.substring(0,pos),string.substring(pos));                    
                    }                                                            
                }
                return null;
            }
        });
        
        filters.add(new Filter() {
            // ultimul caracter e ,
            @Override
            public List<String> match(String string) {
                
                if(string.contains(",")){
                    int pos=string.lastIndexOf(",");
                    if(pos==0){
                        return Arrays.asList(",", string.substring(1));                    
                    } else {
                        return Arrays.asList(string.substring(0,pos+1),string.substring(pos+1));                    
                    }                                                            
                }
                return null;
            }
        });
        
        filters.add(new Filter() {
            // ultimul caracter e ,
            @Override
            public List<String> match(String string) {
                
                if(string.endsWith(".")){
                    return Arrays.asList(string);                                                                                                    
                }
                return null;
            }
        });

        filters.add(new Filter() {

            @Override
            public List<String> match(String string) {
                return null;
            }
        });
        
        

    }

    @Override

    public List<List<Token>> parseInput(String text) {
        List<List<Token>> tmp = new ArrayList<>();
        List<Token> sentence = null;

        //TODO: Adrian, baga tu aici codul pentru tokenizare simpla
        // eliminare taburi si spatii dublate. Se inlocuieste totul cu text
        for (String sen : text.split("\n")) {
            sentence = null;
            List<String> toks = new ArrayList<String>();
            for (String tok : sen.split("[ \t]")) {
                toks.add(tok);
            }
            for (int t = 0; t < toks.size(); t++) {
                String tok = toks.get(t);
                if (tok.length() > 0) {

                    if (sentence == null) { // Propozitie noua
                        sentence = new ArrayList<Token>();
                        tmp.add(sentence);
                    }

                    // verificare pentru ultimul cuvant din propozitie
                    List<String> eos;
                    if (t == toks.size() - 1) {
                        eos = eosFilter.match(tok);
                    } else {
                        eos = null;
                    }

                    if (eos != null) { // daca este ultimul cuvant

                        if (eos.size() == 1) {
                            Token token = new Token();
                            token.word = tok;
                            sentence.add(token);
                        } else { // inseram (cuvant, sfarsit) in locul lui tok si continuam de la cuvant
                            toks.add(t, eos.get(0));
                            toks.set(t + 1, eos.get(1));
                            t--;
                        }

                    } else {

                        // test filtre
                        List<String> result = null;
                        for (int i = 0; i < filters.size() && result == null; i++) {
                            result = filters.get(i).match(tok);
                        }

                        if (result != null) {
                            if (result.size() == 1) {
                                Token token = new Token();
                                token.word = tok;
                                sentence.add(token);
                            } else { // inseram (cuvant, sfarsit) in locul lui tok si continuam de la cuvant
                                toks.add(t, result.get(0));
                                toks.set(t + 1, result.get(1));
                                t--;
                            }
                        } else {  // daca nu a fost gasit
                            Token token = new Token();
                            token.word = tok;
                            token.tag="!";
                            sentence.add(token);                        
                        }
                    }
                }
            }
        }

//        String[] parts = text.split(" ");
//        for (int i = 0; i < parts.length; i++) {
//            Token t = new Token();
//            t.word = parts[i];
//            sentence.add(t);
//        }
//        tmp.add(sentence);
        return tmp;
    }

    @Override
    public void loadModel(String path) {
        Log.i("No input models required for tokenization");
    }

    
    
    public static String[] basicSplit (String sentence) {
        ArrayList<String> toks = new ArrayList<String>();
        Pattern p = Pattern.compile("\\w+", Pattern.UNICODE_CHARACTER_CLASS);
        Matcher m = p.matcher(sentence);
        while (m.find()) {
            
            
            System.out.println(m.group()+" "+m.start()+" "+m.end());
        } 
        
        
        return null;
    }
    
    public static String basicSplitJoin (String[] words) {
        if(words.length==0) return "";
        StringBuilder ret = new StringBuilder();
        for(String s:words){
            ret.append(s);
            ret.append(" ");
        }        
        return ret.toString().substring(0,ret.toString().length()-1);
    }
}
