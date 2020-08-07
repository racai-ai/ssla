/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ssla;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.List;
import ssla.machinelearning.RawFrame;
import ssla.machinelearning.SpeechFrame;
import ssla.machinelearning.SpeechModel;
import ssla.signalprocessing.Bap2Ap;
import ssla.signalprocessing.MLSAVocoder;
import ssla.signalprocessing.STRAIGHTVocoder;
import ssla.signalprocessing.WavWrite;
import ssla.utils.Log;

/**
 *
 * @author tibi
 */
public class SSLA {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, Exception {
//        importData("test/import_test", "test/import_test/models");
//        if (true) {
//            return;
//        }
        /*        Bap2Ap bap2ap = new Bap2Ap();
         double[][] bap = new double[2][2];
         bap[0][0] = 1;
         bap[0][1] = 2;
         bap[0][0] = 3;
         bap[0][1] = 4;
         double[][] ap = bap2ap.bap2ap(bap);
         if (true) {
         return;
         }*/
        //convertim din f0 in pitch
        if (args.length == 0) {
            displayHelp();
        } else if (args.length == 1) {
            //test("models/en", "test/blizzard_test_001.lab");
            synthesize("models/en", "test/blizzard_test_003.lab", "test.wav", 48000, true);
        } else if (args.length == 5 && args[0].equals("--synthesize")) {
            synthesize(args[1], args[2], args[3], 48000, Boolean.valueOf(args[4]));
        } else if (args.length == 3 && args[0].equals("--import")) {
            importData(args[1], args[2]);
        }
        // TODO code application logic here
    }

    private static void displayHelp() {
        System.out.println("Speech Synthesis for Lightwight Applications v0.9 beta.");
        System.out.println("Author(s) (as they come): Tiberiu 'Balls' Boros");
        System.out.println("USAGE:");
        System.out.println("\t--sythesize <model> <lab> <wave> <useGV>");
        System.out.println("\t--import <hts data folder> <output_folder>");
    }

    private static void test(String modelsro, String testtestlab) throws IOException {
        SpeechModel sm = new SpeechModel(modelsro);

        BufferedReader br = new BufferedReader(new FileReader(testtestlab));
        String line = null;
        List<String> parts = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            parts.add(line);
        }

        RawFrame[] tmp = sm.parse(parts.toArray(new String[0]), 1f, false);
        for (int i = 0; i < tmp.length; i++) {
            System.out.println(tmp[i].duration + "\t" + tmp[i].f0 + "\t\t" + tmp[i].lab);
        }

        //scriem in fisier sa vedem de kkt am facut
        DataOutputStream f0 = new DataOutputStream(new FileOutputStream("f0"));
        DataOutputStream mgc = new DataOutputStream(new FileOutputStream("mgc"));
        //smooth
        List<Float> fv = new ArrayList<>();
        List<float[]> mv = new ArrayList<>();

        for (int i = 0; i < tmp.length; i++) {
            for (int k = 0; k < tmp[i].duration; k++) {
                float f = 0;
                if (tmp[i].f0[0] != 0) {
                    f = 48000f / tmp[i].f0[0];
                }
                fv.add(f);
                float[] m = new float[50];
                //writeLE(f0, (float) (f + (Math.random() * 0.1) * f));

                for (int l = 0; l < 50; l++) {
                    f = tmp[i].mgc[l];
                    m[l] = f;
                    //writeLE(mgc, (float) (f + (Math.random() * 0.1) * f));
                }
                mv.add(m);
            }
        }
        //smooth
        for (int k = 0; k < 3; k++) {
            for (int i = 1; i < fv.size() - 1; i++) {
                if (fv.get(i - 1) != 0 && fv.get(i + 1) != 0) {
                    fv.set(i, (fv.get(i - 1) + fv.get(i + 1)) / 2);
                }
                //for (int l = 0; l < 50; l++) {
                //    mv.get(i)[l] = (mv.get(i - 1)[l] + mv.get(i + 1)[l]) / 2;
                //}
            }
        }
        for (int i = 0; i < fv.size(); i++) {
            writeLE(f0, fv.get(i));
            for (int l = 0; l < 50; l++) {
                float f = mv.get(i)[l];
                writeLE(mgc, (float) (f + (Math.random() * 0.05) * f));
            }
        }
        f0.close();
        mgc.close();
        br.close();
    }

    private static void writeLE(DataOutputStream out, float f) throws IOException {
        int value = Float.floatToRawIntBits(f);
        out.writeByte(value & 0xFF);
        out.writeByte((value >> 8) & 0xFF);
        out.writeByte((value >> 16) & 0xFF);
        out.writeByte((value >> 24) & 0xFF);
    }

    private static void synthesize(String models, String lab, String wav, int sr, boolean useGV) throws IOException, Exception {
        SpeechModel sm = new SpeechModel(models);
        BufferedReader br = new BufferedReader(new FileReader(lab));
        String line = null;
        List<String> parts = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            parts.add(line);
        }

        RawFrame[] tmp = sm.parse(parts.toArray(new String[0]), 0.5f, useGV);
//        //TODO: nu stiu daca este ok ceea ce fac acum: incerc sa gasesc anomalii in spectru ca sa elimin artefactele sonore
//        boolean[] anomaly = new boolean[tmp.length];

//        for (int l = 0; l < 30; l++) {
//            for (int i = 1; i < tmp.length - 1; i++) {
//                for (int k = 0; k < tmp[i].mgc.length; k++) {
//                    //if (i % 4 != 0 && i % 5 != 0) {
//                    tmp[i].mgc[k] = (tmp[i - 1].mgc[k] + tmp[i].mgc[k] + tmp[i + 1].mgc[k]) / 3;
//                    //}
//                }
//            }
//        }
        SpeechFrame[] speechFrames = sm.getActualFrames(tmp, useGV);
        float[] f0 = new float[speechFrames.length];
        float[][] mgc = new float[speechFrames.length][];
        float[][] ap = new float[speechFrames.length][];
        for (int i = 0; i < speechFrames.length; i++) {
            f0[i] = speechFrames[i].f0;
            mgc[i] = speechFrames[i].mgc;
            ap[i] = speechFrames[i].ap;
        }

        Log.i("Running vocoder to generate actual speech");
        //MLSAVocoder v = new MLSAVocoder();
        //double[] audio = v.runMLSAVocoder(mgc, f0, mgc[0].length, sr / 200, sr);
        STRAIGHTVocoder v = new STRAIGHTVocoder(models);
        String[] phonemes = new String[tmp.length / 5];
        int[] durations = new int[phonemes.length];
        for (int i = 0; i < tmp.length; i++) {
            phonemes[i / 5] = tmp[i].lab.substring(tmp[i].lab.indexOf("-") + 1);
            phonemes[i / 5] = phonemes[i / 5].substring(0, phonemes[i / 5].indexOf("+"));
            durations[i / 5] += tmp[i].duration;
        }

        double[] audio = v.runSTRAIGHTVocoder(mgc, ap, f0, sr, phonemes, durations);
        DataOutputStream dos = new DataOutputStream(new FileOutputStream(wav));
        Log.i("Creating output WAVE");
        WavWrite.Save2Wave(audio, sr, dos);
        dos.close();
    }

    private static List<String> listFiles(String source, String extension) {
        List<String> tmp = new ArrayList<>();
        File folder = new File(source);
        File[] listOfFiles = folder.listFiles();

        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile() && listOfFiles[i].getName().endsWith(extension)) {
                tmp.add(listOfFiles[i].getAbsolutePath());
            }
        }
        return tmp;
    }

    private static void importData(String source, String destination) throws FileNotFoundException, IOException {
        List<String> labs = listFiles(source + "/labels/mono", ".lab");
        DataOutputStream dosSpectrum = new DataOutputStream(new FileOutputStream(destination + "/hybrid.sp"));
        DataOutputStream dosF0 = new DataOutputStream(new FileOutputStream(destination + "/hybrid.f0"));
        DataOutputStream dosIF0 = new DataOutputStream(new FileOutputStream(destination + "/hybrid.if0"));
        DataOutputStream dosAp = new DataOutputStream(new FileOutputStream(destination + "/hybrid.ap"));
        BufferedWriter dosPhs = new BufferedWriter(new FileWriter(destination + "/hybrid.phs"));
        int current_frame = 0;
        for (int i = 0; i < labs.size(); i++) {
            System.out.print("Processing file " + (i + 1) + "/" + labs.size() + " '" + labs.get(i) + "'...");
            System.out.flush();
            current_frame += importFile(current_frame, source, labs.get(i), dosSpectrum, dosF0, dosIF0, dosAp, dosPhs);
            System.out.println("done");
        }
        dosSpectrum.close();
        dosF0.close();
        dosAp.close();
        dosPhs.close();
    }

    private static int importFile(int frame_offset, String source, String lab, DataOutputStream dosSpectrum, DataOutputStream dosF0, DataOutputStream dosIF0, DataOutputStream dosAp, BufferedWriter dosPhs) throws FileNotFoundException, IOException {
        String baseName = lab.substring(lab.lastIndexOf("/") + 1);
        baseName = baseName.replace(".lab", "");
        File f;
        f = new File(source + "/sp/" + baseName + ".sp");
        if (!f.exists()) {
            return 0;
        }
        f = new File(source + "/ap/" + baseName + ".ap");
        if (!f.exists()) {
            return 0;
        }
        f = new File(source + "/f0/" + baseName + ".f0");
        if (!f.exists()) {
            return 0;
        }
        //numaram cate frame-uri sunt. Exista un BUG in STRAIGHT care face ca numarul de frame-uri pentru SP sa fie diferit de AP si F0
        BufferedReader br;
        String line;
        int cnt;
        int min_frame = Integer.MAX_VALUE;
        br = new BufferedReader(new FileReader(source + "/sp/" + baseName + ".sp"));
        cnt = 0;
        while ((line = br.readLine()) != null) {
            if (!line.trim().isEmpty()) {
                cnt++;
            }
        }
        if (min_frame > cnt) {
            min_frame = cnt;
        }
        br = new BufferedReader(new FileReader(source + "/ap/" + baseName + ".ap"));
        cnt = 0;
        while ((line = br.readLine()) != null) {
            if (!line.trim().isEmpty()) {
                cnt++;
            }
        }
        if (min_frame > cnt) {
            min_frame = cnt;
        }
        br = new BufferedReader(new FileReader(source + "/f0/" + baseName + ".f0"));
        cnt = 0;
        while ((line = br.readLine()) != null) {
            while (!line.replace("  ", " ").equals(line)) {
                line = line.replace("  ", " ");
            }
            line = line.trim();
            cnt = line.split(" ").length;
        }
        if (min_frame > cnt) {
            min_frame = cnt;
        }
        //citim spectrul
        System.out.print("sp");
        System.out.flush();
        br = new BufferedReader(new FileReader(source + "/sp/" + baseName + ".sp"));
        cnt = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) {
                continue;
            }
            cnt++;
            if (cnt > min_frame) {
                cnt--;
                break;
            }
            while (!line.replace("  ", " ").equals(line)) {
                line = line.replace("  ", " ");
            }
            line = line.trim();
            String[] parts = line.split(" ");
            byte[] data = new byte[parts.length * 4];
            for (int i = 0; i < parts.length; i++) {
                float fff = Float.parseFloat(parts[i]);
                byte[] tmp = ByteBuffer.allocate(4).putFloat(fff).array();

                System.arraycopy(tmp, 0, data, i * 4, 4);
            }
            dosSpectrum.write(data);
        }
        System.out.print("(" + cnt + " frames)");
        br.close();
        br = new BufferedReader(new FileReader(source + "/ap/" + baseName + ".ap"));
        System.out.print(" ap");
        System.out.flush();
        cnt = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) {
                continue;
            }
            cnt++;
            if (cnt > min_frame) {
                cnt--;
                break;
            }
            while (!line.replace("  ", " ").equals(line)) {
                line = line.replace("  ", " ");
            }
            line = line.trim();
            String[] parts = line.split(" ");
            byte[] data = new byte[parts.length * 4];
            for (int i = 0; i < parts.length; i++) {
                float fff = Float.parseFloat(parts[i]);
                byte[] tmp = ByteBuffer.allocate(4).putFloat(fff).array();

                System.arraycopy(tmp, 0, data, i * 4, 4);
            }
            dosAp.write(data);

        }
        br.close();
        System.out.print("(" + cnt + " frames)");

        br = new BufferedReader(new FileReader(source + "/f0/" + baseName + ".f0"));
        System.out.print(" F0");
        System.out.flush();

        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) {
                continue;
            }
            while (!line.replace("  ", " ").equals(line)) {
                line = line.replace("  ", " ");
            }
            line = line.trim();
            String[] parts = line.split(" ");
            byte[] data = new byte[parts.length * 4];
            //cnt = parts.length;

            float[] f0 = new float[min_frame];
            float[] if0 = new float[min_frame];

            for (int i = 0; i < min_frame; i++) {
                float fff = Float.parseFloat(parts[i]);
                f0[i] = fff;
                byte[] tmp = ByteBuffer.allocate(4).putFloat(fff).array();

                System.arraycopy(tmp, 0, data, i * 4, 4);
            }

            float last_f0 = 0;
            float next_f0 = 0;
            int dist_last = 0;
            int dist_next = 0;
//            double[] nf0 = new double[hr.f0.length];
            for (int i = 0; i < f0.length; i++) {
                if (f0[i] == 0) {
                    last_f0 = 0;
                    next_f0 = 0;
                    for (int k = i - 1; k >= 0; k--) {
                        if (f0[k] != 0) {
                            last_f0 = f0[k];
                            break;
                        } else {
                            dist_last++;
                        }
                    }

                    for (int k = i + 1; k < f0.length; k++) {
                        if (f0[k] != 0) {
                            next_f0 = f0[k];
                            break;
                        } else {
                            dist_next++;
                        }
                    }
                    if (last_f0 == 0) {
                        if0[i] = next_f0;
                    } else if (next_f0 == 0) {
                        if0[i] = last_f0;
                    } else {
//                if (last_f0!=0 && next)
                        float total_dist = dist_last + dist_next;
                        float p0 = 1.0f - (float) (dist_last) / total_dist;
                        float p1 = 1.0f - (float) (dist_next) / total_dist;
                        if0[i] = p0 * last_f0 + p1 * next_f0;
                    }
                } else {
                    if0[i] = f0[i];
                }
            }
            cnt = min_frame;
            dosF0.write(data);
            for (int i = 0; i < min_frame; i++) {
                byte[] tmp = ByteBuffer.allocate(4).putFloat(if0[i]).array();

                System.arraycopy(tmp, 0, data, i * 4, 4);
            }
            dosIF0.write(data);
        }
        br.close();
        System.out.print("(" + cnt + " frames) ");

        br = new BufferedReader(new FileReader(source + "/labels/mono/" + baseName + ".lab"));
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) {
                continue;
            }
            String[] parts = line.replace("   ", " ").split(" ");
            int start = Integer.parseInt(parts[0]) / 10000;
            int stop = Integer.parseInt(parts[1]) / 10000;
            String phoneme = parts[2];
            start /= 5;
            stop /= 5;
            dosPhs.write(baseName + "\t" + phoneme + "\t" + (start + frame_offset) + "\t" + (stop + frame_offset) + "\n");
        }
        br.close();
        return min_frame;
    }

}
