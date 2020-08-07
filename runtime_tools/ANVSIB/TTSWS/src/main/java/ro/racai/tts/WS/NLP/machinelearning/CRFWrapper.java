/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.machinelearning;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.UUID;

/**
 *
 * @author tibi
 */
public class CRFWrapper {

    public static CRFResults train(String template_file, String train_file, String model_file) throws IOException {
        CRFResults tmp = new CRFResults();
        Process p = Runtime.getRuntime().exec("external_tools/crf/crf_learn " + template_file + " " + train_file + " " + model_file);
        InputStream is = p.getInputStream();
        byte[] b = new byte[255];
        while (true) {
            int rez = is.read(b);
            if (rez == -1) {
                break;
            }
            String s = new String(b, 0, rez);
            String[] lines = s.split("\n");
            for (int i = 0; i < lines.length; i++) {
                System.err.println("\t::" + lines[i]);
            }
        }
        return tmp;
    }

    public static CRFResults test(String model_file, String test_file) throws IOException {
        CRFResults tmp = new CRFResults();
        byte[] b = new byte[255];

        Process p = Runtime.getRuntime().exec("external_tools/crf/crf_test -v1 -m " + model_file + " " + test_file);
        InputStream is = p.getInputStream();
        while (true) {
            int rez = is.read(b);
            if (rez == -1) {
                break;
            }
            String s = new String(b, 0, rez);
//            String[] lines = s.split("\n");
//            for (int i = 0; i < lines.length; i++) {
//                //System.err.println("\t::" + lines[i]);
//                bw.write(lines[i]+"\n");
//            }
        }
        return tmp;
    }

    public static String[] labelSequence(String model_file, String[] sequence) throws IOException {
        //generam un fisier temporar
        String test_file = "tmp" + File.separator + UUID.randomUUID().toString();
        File f = new File("tmp");
        f.mkdirs();

        BufferedWriter bw = new BufferedWriter(new FileWriter(test_file));
        for (int i = 0; i < sequence.length; i++) {
            bw.write(sequence[i] + "\n");
        }
        bw.close();
        String[] output = new String[sequence.length];
        byte[] b = new byte[255];

        Process p = Runtime.getRuntime().exec("external_tools/crf/crf_test -v1 -m " + model_file + " " + test_file);
        InputStream is = p.getInputStream();
        StringBuilder sb = new StringBuilder();
        while (true) {
            int rez = is.read(b);
            if (rez == -1) {
                break;
            }
            String s = new String(b, 0, rez);
            sb.append(s);
        }
        //stergem fisierul
        f = new File(test_file);
        f.delete();

        String[] lines = sb.toString().trim().split("\n");
        int pos = 0;
        for (int i = 0; i < lines.length; i++) {
            if (!lines[i].startsWith("#") && !lines[i].isEmpty()) {
                String[] pp = lines[i].split("\t");
                String label = pp[pp.length - 1];
                output[pos++] = label.substring(0, label.lastIndexOf("/"));
            } else if (lines[i].isEmpty()) {
                pos++;
            }
        }
//        System.out.println("=======================================");
//        System.out.println(sb.toString().trim());
//        System.out.println("=======================================");
        return output;
    }

    public static class CRFResults {

        public float serr;
        public float terr;
    }
}
