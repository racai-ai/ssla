/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.TTS.machinelearning;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import ro.racai.tts.WS.TTS.utils.Log;


/**
 *
 * @author tibi
 */
public class DecissionTree {

    Map<String, String[]> rules = new TreeMap<>();
    List<String> leafs = new ArrayList<>();

    List<ArrayList<Node>> StateNodes = new ArrayList<>();
    float[][][] pdfs = null;
    boolean isMSD = false;
    int ssize = 0;

    public int getNumStates() {
        return StateNodes.size();
    }

    public DecissionTree(String tree, String pdf, String s_label, int num_windows) throws FileNotFoundException, IOException {
        Log.i("\tloading '" + tree + "'");
        BufferedReader br = new BufferedReader(new FileReader(tree));
        String line = null;
        ArrayList<Node> currentList = null;
        while ((line = br.readLine()) != null) {
            if (line.trim().isEmpty()) {
                continue;
            }

            String[] parts = line.replaceAll(" +", " ").trim().split(" ");
            if (parts.length == 0) {
                continue;
            }
            if (parts[0].equals("QS")) {
                String label = parts[1];
                String[] rules = parts[3].trim().replace(",", " ").replace("\"", "").split(" ");
                this.rules.put(label, rules);
                //Feature f = new Feature();
                //f.label = label;
                //f.rules = rules;
                //features.add(f);
            } else if (parts[0].startsWith("{*}")) {
                currentList = new ArrayList<>();
                StateNodes.add(currentList);
            } else if (!parts[0].equals("}") && !parts[0].equals("{")) {
                //caz special (nu avem decat o singura clasa peste tot)
                if (parts.length == 1) {
                    Node n = new Node();
                    n.feature = "";
                    n.yes = 0;
                    n.no = 0;
                    currentList.add(n);
                    leafs.add("single");
                } else {
                    String feature = parts[1];
                    String left = parts[2].replace("\"", "");
                    String right = parts[3].replace("\"", "");
                    int left_i;
                    int right_i;
                    if (!left.startsWith("-")) {
                        //avem o frunza
                        int index = -1;
                        if (!leafs.contains(left)) {
                            index = leafs.size();
                            leafs.add(left);
                        } else {
                            index = leafs.indexOf(left);
                        }
                        left_i = index;
                    } else {
                        left_i = Integer.parseInt(left);
                    }
                    if (!right.startsWith("-")) {
                        //avem o frunza
                        int index = -1;
                        if (!leafs.contains(right)) {
                            index = leafs.size();
                            leafs.add(right);
                        } else {
                            index = leafs.indexOf(right);
                        }
                        right_i = index;
                    } else {
                        right_i = Integer.parseInt(right);
                    }

                    Node n = new Node();
                    n.no = left_i;
                    n.yes = right_i;
                    n.feature = feature;
                    currentList.add(n);
                }

            }

        }
        br.close();
        //incarcam PDF-urile
        if (pdf == null) {
            return;
        }
        Log.i("\tloading PDF '" + pdf + "'");
        RandomAccessFile f = new RandomAccessFile(pdf, "r");
        byte[] bytes = new byte[(int) f.length()];
        f.read(bytes);

        int is_msd = ByteBuffer.wrap(bytes, 0, 4).order(ByteOrder.BIG_ENDIAN).getInt(); //asta e validat
        if (is_msd == 1) {
            isMSD = true;
        } else {
            isMSD = false;
        }
        ssize = ByteBuffer.wrap(bytes, 4, 4).order(ByteOrder.BIG_ENDIAN).getInt();//cica reprezinta stream_size
        int vector_length = ByteBuffer.wrap(bytes, 8, 4).order(ByteOrder.BIG_ENDIAN).getInt(); //si asta e validat
        ///validat tot ce e jos pana la ///////
        int[] nPdfs = new int[StateNodes.size()];
        int total_pdfs = 0;
        for (int i = 0; i < nPdfs.length; i++) {
            nPdfs[i] = ByteBuffer.wrap(bytes, 12 + i * 4, 4).order(ByteOrder.BIG_ENDIAN).getInt();
            total_pdfs += nPdfs[i];
        }
        if (total_pdfs != leafs.size()) {
            Log.e("The PDF and tree info are not in sync. The number of distinct leafs is " + leafs.size() + " and the number of PDFs is " + total_pdfs);
            return;
        }
        int fp = 12 + nPdfs.length * 4;
        pdfs = new float[nPdfs.length][][];
        /////////
        if (isMSD) {
            for (int i = 0; i < nPdfs.length; i++) {
                pdfs[i] = new float[nPdfs[i]][];
                for (int k = 0; k < nPdfs[i]; k++) {
                    float[] cpdf = new float[2 * vector_length + 1];
                    for (int l = 0; l < ssize; l++) {
                        for (int m = 0; m < vector_length / ssize; m++) {
                            int temp = ByteBuffer.wrap(bytes, fp, 4).order(ByteOrder.BIG_ENDIAN).getInt();
                            fp += 4;
                            cpdf[l * vector_length / ssize + m] = Float.intBitsToFloat(temp);
                            temp = ByteBuffer.wrap(bytes, fp, 4).order(ByteOrder.BIG_ENDIAN).getInt();
                            fp += 4;
                            cpdf[l * vector_length / ssize + m + vector_length] = Float.intBitsToFloat(temp);
                        }
                        int temp = ByteBuffer.wrap(bytes, fp, 4).order(ByteOrder.BIG_ENDIAN).getInt();
                        fp += 4;
                        if (l == 0) {
                            if (Float.intBitsToFloat(temp) < 0.0 || Float.intBitsToFloat(temp) > 1.0) {
                                Log.e("HTS_Model_load_pdf: MSD weight should be within 0.0 to 1.0.\n");
                            }
                            cpdf[2 * vector_length] = Float.intBitsToFloat(temp);
                        }
                        fp += 4;

                    }
                    pdfs[i][k] = cpdf;
                }
            }
        } else {
            for (int i = 0; i < nPdfs.length; i++) {
                pdfs[i] = new float[nPdfs[i]][];

                for (int k = 0; k < nPdfs[i]; k++) {
                    float[] cpdf = new float[2 * vector_length];
                    for (int m = 0; m < vector_length; m++) {
                        int temp = ByteBuffer.wrap(bytes, fp, 4).order(ByteOrder.BIG_ENDIAN).getInt();
                        fp += 4;
                        cpdf[m] = Float.intBitsToFloat(temp);
                        temp = ByteBuffer.wrap(bytes, fp, 4).order(ByteOrder.BIG_ENDIAN).getInt();
                        fp += 4;
                        cpdf[m + vector_length] = Float.intBitsToFloat(temp);
                    }
                    pdfs[i][k] = cpdf;
                }
            }
        }

        f.close();
    }

    public float[] getPdf(String label) {
        if (label.equals("single")) {
            return pdfs[0][0];
        } else {
            //incercam sa parsam
            String[] parts = label.split("_");
            if (StateNodes.size() > 1) {
                int i = Integer.parseInt(parts[1].substring(1)) - 2;
                int j = Integer.parseInt(parts[2]) - 1;
                return pdfs[i][j];
            } else {
                int j = Integer.parseInt(parts[parts.length - 1]) - 1;
                return pdfs[0][j];
            }
        }
    }

    public String parse(String features, int state) {
        List<Node> nodes = StateNodes.get(state);
        Node current = nodes.get(0);
        if (current.yes == current.no) {
            return leafs.get(0);
        }
        //altfel parcugem nodurile
        while (true) {
            String[] rule = rules.get(current.feature);
            //verificam daca se aplica ceva
            int next = current.no;
            for (int i = 0; i < rule.length; i++) {
                if (applies(rule[i], features)) {
                    next = current.yes;
                    break;
                }
            }
            if (next < 0) {
                current = nodes.get(-next);
            } else {
                return leafs.get(next);
            }
        }
    }

    private boolean applies(String rule, String features) {
        byte[] brule = rule.getBytes();
        byte[] bf = features.getBytes();
        for (int i = 0; i < bf.length - brule.length; i++) {
            boolean ok = true;
            for (int k = 0; k < brule.length; k++) {
                if (brule[k] != 42 && bf[i + k] != brule[k]) {//(42='*')
                    ok = false;
                    break;
                }
            }
            if (ok) {
                return true;
            }
        }
        return false;

    }

    class Node {

        public String feature;
        public int no;//presupun ca stanga (cu - inseamna nod cu + inseama frunza)
        public int yes;//presupun ca dreapta
    }
}
