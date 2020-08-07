/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package importcorpus;

import importcorpus.io.WavFile;
import importcorpus.io.WavFileException;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author boros
 */
public class ImportCorpus {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, WavFileException {
        // TODO code application logic here
        if (args.length == 0) {
            displayHelp();
        } else if (args.length == 4 && args[0].equals("--import")) {
            importCorpus(args[1], args[2], args[3]);
        }
    }

    private static void displayHelp() {
        System.out.println("SSLA Neural version corpus import tool");
        System.out.println("Usage: ");
        System.out.println("\t--import <train_base> <dev_base> <model_output_base>");
    }

    private static void importCorpus(String train_base, String dev_base, String model_output) throws IOException, WavFileException {
        System.out.println("\nImporting training set...");
        List<String> trainFiles = getCompleteFileList(train_base);
        BufferedWriter bw = new BufferedWriter(new FileWriter(model_output + "/train"));
        for (int i = 0; i < trainFiles.size(); i++) {
            String path = trainFiles.get(i);
            System.out.print("\rImporting file " + (i + 1) + "/" + trainFiles.size() + ": " + path);
            importExample(path, bw);
        }
        bw.close();

        System.out.println("\n\nImporting development set...");
        List<String> devFiles = getCompleteFileList(dev_base);
        bw = new BufferedWriter(new FileWriter(model_output + "/dev"));

        for (int i = 0; i < devFiles.size(); i++) {
            String path = devFiles.get(i);
            System.out.print("\rImporting file " + (i + 1) + "/" + devFiles.size() + ": " + path);
            importExample(path, bw);
        }

        bw.close();
        System.out.println("\n\nDone");
    }

    private static List<String> getCompleteFileList(String train_base) {
        List<String> tmp = new ArrayList<>();

        File folder = new File(train_base);
        File[] listOfFiles = folder.listFiles();
        for (File listOfFile : listOfFiles) {
            if (listOfFile.isFile()) {
                String path = listOfFile.getPath();
                if (path.endsWith(".wav")) {
                    File lab = new File(path.replace(".wav", ".lab"));
                    //File phs = new File(path.replace(".wav", ".phs"));
                    if (lab.exists()) {
                        tmp.add(path.replace(".wav", ""));
                    }
                }
            }
        }
        return tmp;
    }

    static int ulawEncode(int sample) {
        double x = (double) sample / 32768;
        double mean = 255;
        double y = Math.signum(x) * Math.log(1.0 + mean * Math.abs(x)) / Math.log(1.0 + mean);
        int discrete = (int) (y * 127);//value is between -127 and +127
        return discrete;
    }

    static int ulawDecode(int sample) {
        double y = (double) sample / 127;
        double mean = 255;
        double x = Math.signum(y) * (1.0 / mean) * (Math.pow(1 + mean, Math.abs(y)) - 1);
        return (int) (x * 32768);
    }

    private static void importExample(String path, BufferedWriter bw) throws IOException, WavFileException {
        File wav = new File(path.concat(".wav"));
        File phs = new File(path.concat(".phs"));
        File lab = new File(path.concat(".lab"));
        WavFile wavFile = WavFile.openWavFile(wav);
        int[] buffer = new int[(int) wavFile.getNumFrames()];
        wavFile.readFrames(buffer, buffer.length);
        wavFile.close();
        int last_sample = 0;
        boolean okToWrite = false;
        for (int x = 0; x < buffer.length; x++) {
            int enc = ulawEncode(buffer[x]);
            if (Math.abs(enc)>0)
                okToWrite=true;
            if (okToWrite) {
                bw.write(last_sample + " " + enc + "\n");
            }
            last_sample = enc;
        }
        bw.write("\n");
    }

    private List<PHONEME_INFO> readLab(String lab) {
        List<PHONEME_INFO> tmp = new ArrayList<>();
        return tmp;
    }

    private static class PHONEME_INFO {

        int start;
        int stop;
        String label;
    }

}
