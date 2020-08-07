/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ssla.signalprocessing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import static java.util.Arrays.sort;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import ssla.straight.straightSynth;
import static java.util.Arrays.sort;

/**
 *
 * @author tibi
 */
public class STRAIGHTVocoder {

    String models = null;

    public STRAIGHTVocoder(String hybrid_data) {
        models = hybrid_data;
    }

    public double[] runSTRAIGHTVocoder(float[][] mgcf, float[][] apf, float[] f0Rawf, int sr, String[] phonemes, int[] durations) throws Exception {
        double[][] ap = new double[apf.length][apf[0].length];
        double[] f0Raw = new double[f0Rawf.length];
        double[][] mgc = new double[mgcf.length][mgcf[0].length];
        for (int i = 0; i < ap.length; i++) {
            f0Raw[i] = f0Rawf[i];
            for (int j = 0; j < ap[0].length; j++) {
                ap[i][j] = apf[i][j];
            }
            for (int j = 0; j < mgcf[0].length; j++) {
                mgc[i][j] = mgcf[i][j];
            }
        }
        double[][] sp = mgc2sp(mgc);
        double shiftm = 5;
        double fs = sr;
        double pconv = 1;
        double fconv = 1;
        double sconv = 1;
        double gdbw = 70;
        double delfrac = 0.2000f;
        double delsp = 0.5000f;
        double cornf = 4000f;
        double delfracind = 0;
        double imap = 1;
        double lowestF0 = 50;

        //LOL: codul de mai jos va purta numele Phane Transpose
        double[][] tsp = new double[sp[0].length][sp.length];
        for (int i = 0; i < sp.length; i++) {
            for (int j = 0; j < sp[0].length; j++) {
                tsp[j][i] = sp[i][j];
            }
        }

        double[][] tap = new double[ap[0].length][ap.length];
        for (int i = 0; i < ap.length; i++) {
            for (int j = 0; j < ap[i].length; j++) {
                tap[j][i] = ap[i][j];
            }
        }

        HybridResult hr = hybridify(tsp, tap, f0Raw, phonemes, durations);

        drawPPM(tsp, "artif.ppm");

        double[] signal = straightSynth.straightSynth(hr.tsp, hr.f0, shiftm, fs, pconv, fconv, sconv, gdbw, delfrac, delsp, cornf, delfracind, hr.tap, imap, lowestF0);
        //signal = remodelUsingNaturalUnits(signal, hr);
        return signal;
    }

    private double computeTransition(ViterbiToken src, double[][] spDst, double[] f0Dst, double[][] tsp, double[] tf0) throws IOException {
//        if (true) {
//            return 0;
//        }
//        //double[][] spSrc = readArrays(models + "/hybrid.sp", src.stopFrame - 1, src.stopFrame);
        double[][] spSrc;
        double[] f0Src;
        if (src.startFrame > 0) {
            spSrc = readArrays(models + "/hybrid.sp", src.stopFrame - 1, src.stopFrame);
            f0Src = readArray(models + "/hybrid.f0", src.stopFrame - 4, src.stopFrame);
        } else {
            f0Src = new double[5];
            spSrc = new double[2][1025];
            for (int i = 0; i < 1025; i++) {
                spSrc[0][i] = tsp[i][-src.startFrame];
                spSrc[1][i] = tsp[i][-src.stopFrame];
            }
            for (int i = 0; i < f0Src.length; i++) {
                f0Src[i] = tf0[-src.stopFrame - f0Src.length + i];
            }
        }
        //double[][] spDst = readArrays(models + "/hybrid.sp", dst.startFrame, dst.startFrame + 1);
        double cost = 0;//dtw(f0Src, f0Dst);
        //facem doar pentru banda de 6khz;
        for (int i = 0; i < 1025; i++) {
            cost += Math.abs(spSrc[spSrc.length - 1][i] - spDst[0][i]);
        }
        return cost;
    }

    public double dtw(double[] a, double[] b) {
        for (int i = 0; i < a.length; i++) {
            if (a[i] != 0) {
                a[i] = Math.log(a[i]);
            }
        }
        for (int i = 0; i < b.length; i++) {
            if (b[i] != 0) {
                b[i] = Math.log(b[i]);
            }
        }
        double[][] x = new double[a.length + 1][b.length + 1];
        for (int i = 1; i < x.length; i++) {
            x[i][0] = Math.pow(a[i - 1], 2) + x[i - 1][0];
        }
        for (int i = 1; i < x[0].length; i++) {
            x[0][i] = Math.pow(b[i - 1], 2) + x[0][i - 1];
        }

        for (int i = 1; i < x.length; i++) {
            for (int j = 1; j < x[0].length; j++) {
                double min = x[i - 1][j - 1];
                if (x[i - 1][j] < min) {
                    min = x[i - 1][j];
                }
                if (x[i][j - 1] < min) {
                    min = x[i][j - 1];
                }
                double cost = Math.pow(a[i - 1] - b[j - 1], 2);
                x[i][j] = min + cost;
            }
        }
        return x[x.length - 1][x[0].length - 1];
    }

    public double dtw(double[][] a, double[][] b) {
        double[][] x = new double[a.length + 1][b.length + 1];
        double cc = 0;
        for (int i = 1; i < x.length; i++) {
            cc = 0;
            for (int zz = 0; zz < b[0].length / 2; zz++) {
                cc += Math.abs(a[i - 1][zz]);
            }
            x[i][0] = cc + x[i - 1][0];
        }
        for (int i = 1; i < x[0].length; i++) {
            cc = 0;
            for (int zz = 0; zz < b[0].length / 2; zz++) {
                cc += Math.abs(b[i - 1][zz]);
            }
            x[0][i] = cc + x[0][i - 1];
        }

        for (int i = 1; i < x.length; i++) {
            for (int j = 1; j < x[0].length; j++) {
                double min = x[i - 1][j - 1];
                if (x[i - 1][j] < min) {
                    min = x[i - 1][j];
                }
                if (x[i][j - 1] < min) {
                    min = x[i][j - 1];
                }
                double cost = 0;
                for (int zz = 0; zz < b[0].length / 2; zz++) {
                    cost += Math.abs(a[i - 1][zz] - b[j - 1][zz]);
                }
                x[i][j] = min + cost;
            }
        }
        return x[x.length - 1][x[0].length - 1];
    }

    private double computeEmissionScore(int startFrame, int stopFrame, double[][] tsp, double[] tf0, int pos, int length) throws IOException {
        double[][] sp = readArrays(models + "/hybrid.sp", startFrame, stopFrame);
        double[] f0 = readArray(models + "/hybrid.f0", startFrame, stopFrame);
        double[] tff0 = new double[length];
        for (int i = 0; i < length; i++) {
            tff0[i] = tf0[i + pos];
        }
        double f0Cost = 0;//dtw(f0, tff0);
        //doar banda de 6000 hz
        double[][] tssp = new double[length][1025];
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < 1025; j++) {
                tssp[i][j] = tsp[j][i + pos];
            }
        }
        double spCost = dtw(sp, tssp);
        double score = (f0Cost + spCost) / length;

        return score;
    }

    private Map<String, List<PhoneInfo>> readPhones(String file) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        Map<String, List<PhoneInfo>> map = new TreeMap<>();
        String line;
        String lastLabelFile = "__";
        int ofs = 0;
        while ((line = br.readLine()) != null) {
            String[] parts = line.split("\t");
            if (!lastLabelFile.equals(parts[0])) {
                lastLabelFile = parts[0];
                ofs = Integer.parseInt(parts[2]);
            }
            String cp = parts[1];
            String phoneme = cp;

            int strFrame = Integer.parseInt(parts[2]);
            int stoFrame = Integer.parseInt(parts[3]);
            int startFrame = strFrame;
            int stopFrame = stoFrame;

            if (!phoneme.startsWith("#")) {
                if (!map.containsKey(phoneme)) {
                    map.put(phoneme, new ArrayList<PhoneInfo>());
                }
                PhoneInfo di = new PhoneInfo();
                di.wav = parts[0];
                di.startFrame = startFrame;
                di.stopFrame = stopFrame;
                di.wavStartFrame = startFrame - ofs;
                di.lab = phoneme;
                List<PhoneInfo> value = map.get(phoneme);
                value.add(di);
            }
        }
        br.close();
        return map;
    }
    double EMISSION_THRESHOLD = 3;
    double DURATION_THRESHOLD = 3;
    double LAMBDA_EMISSION = 1;
    double LAMBDA_TRANSITION = 1;

    private ViterbiToken[] getCompatibleTokens(ViterbiToken vt, double[][] tsp, double[] f0Raw, Map<String, List<PhoneInfo>> realPhones) throws IOException {
        if (true) {
            System.out.println("Hybrid synthesis disabled for " + vt.lab);
            ViterbiToken[] tmp = new ViterbiToken[1];
            tmp[0] = vt;
            vt.signalStartFrame = -vt.startFrame;
            vt.signalStopFrame = -vt.stopFrame;
            vt.totalCost = Double.POSITIVE_INFINITY;

            return tmp;

        }
        String lab = vt.lab;
        List<PhoneInfo> ru = realPhones.get(lab);
        if (ru == null || lab.equals("pau") || lab.equals("sp")) {
            System.out.println("No candidates found for " + lab + " using synthetic unit instead");
            ViterbiToken[] tmp = new ViterbiToken[1];
            tmp[0] = vt;
            vt.signalStartFrame = -vt.startFrame;
            vt.signalStopFrame = -vt.stopFrame;
            vt.totalCost = Double.POSITIVE_INFINITY;

            return tmp;
        } else {
            int cnt = 0;
            double[] ems_scores = new double[ru.size()];
            for (int i = 0; i < ru.size(); i++) {
                PhoneInfo xx = ru.get(i);
                int ddur = Math.abs(vt.stopFrame - vt.startFrame);
                int rdur = xx.stopFrame - xx.startFrame;
                if ((rdur >= ddur - DURATION_THRESHOLD) && (rdur <= ddur + DURATION_THRESHOLD)) {
                    ems_scores[i] = computeEmissionScore(xx.startFrame, xx.stopFrame, tsp, f0Raw, -vt.startFrame, -(vt.stopFrame - vt.startFrame));
                    if (ems_scores[i] < EMISSION_THRESHOLD) {
                        cnt++;
                    }
                }
            }
            if (cnt == 0) {
                System.out.println("No candidates found for " + lab + " using synthetic unit instead");
                ViterbiToken[] tmp = new ViterbiToken[1];
                tmp[0] = vt;
                vt.signalStartFrame = -vt.startFrame;
                vt.signalStopFrame = -vt.stopFrame;
                vt.totalCost = Double.POSITIVE_INFINITY;
                return tmp;
            } else {
                System.out.println("Found " + cnt + " natural units for " + lab);
                ViterbiToken[] tmp = new ViterbiToken[cnt];
                int pos = 0;
                double tot_ems_cost = 0;
                for (int i = 0; i < ru.size(); i++) {
                    PhoneInfo xx = ru.get(i);
                    int ddur = Math.abs(vt.stopFrame - vt.startFrame);
                    int rdur = xx.stopFrame - xx.startFrame;
                    if ((rdur >= ddur - DURATION_THRESHOLD) && (rdur <= ddur + DURATION_THRESHOLD)) {
                        if (ems_scores[i] < EMISSION_THRESHOLD) {
                            ViterbiToken nvt = new ViterbiToken();
                            nvt.bestFrom = -1;
                            nvt.startFrame = xx.startFrame;
                            nvt.stopFrame = xx.stopFrame;
                            nvt.totalCost = Double.POSITIVE_INFINITY;
                            nvt.lab = lab;
                            nvt.emissionCost = ems_scores[i];
                            tot_ems_cost += ems_scores[i];
                            nvt.signalStartFrame = -vt.startFrame;
                            nvt.signalStopFrame = -vt.stopFrame;
                            nvt.wav = xx.wav;
                            nvt.wavStartFrame = xx.wavStartFrame;
                            tmp[pos++] = nvt;
                        }
                    }
                }
//                tmp[tmp.length - 1] = vt;
//                vt.signalStartFrame = -vt.startFrame;
//                vt.signalStopFrame = -vt.stopFrame;
//                vt.totalCost = Double.POSITIVE_INFINITY;
//                vt.emissionCost = 64;
                return tmp;
            }
        }

    }

    private void drawPPM(double[][] tsp, String file) throws IOException {
        double min = tsp[0][0];
        double max = tsp[0][0];
        for (int i = 0; i < tsp.length; i++) {
            for (int j = 0; j < tsp[0].length; j++) {
                if (min > tsp[i][j]) {
                    min = tsp[i][j];
                }
                if (max < tsp[i][j]) {
                    max = tsp[i][j];
                }
            }
        }
        max -= 0.5 * max;
        min += 0.5 * min;
        BufferedWriter bw = new BufferedWriter(new FileWriter(file));
        bw.write("P3\n");
        bw.write(tsp.length + " " + tsp[0].length + "\n");
        bw.write("255\n");
        for (int j = 0; j < tsp[0].length; j++) {
            for (int i = 0; i < tsp.length; i++) {
                int col = (int) (((tsp[i][j] - min) / (max - min)) * 255);
                if (col < 0) {
                    col = 0;
                }
                if (col > 255) {
                    col = 255;
                }

                bw.write(col + " " + col + " " + col + " ");
            }
            bw.write("\n");
        }
        bw.close();
    }

    private double[] remodelUsingNaturalUnits(double[] signal, ViterbiToken[] hr) throws FileNotFoundException, IOException {
        List<Double> sign = new ArrayList<Double>();
        for (int i = 0; i < hr.length; i++) {
            if (hr[i].startFrame < 0 || hr[i].wav == null) {
                for (int k = hr[i].signalStartFrame * 240; k < hr[i].signalStopFrame * 240; k++) {
                    sign.add(signal[k]);
                }
            } else {
                long skip = hr[i].wavStartFrame * 2 * 240 + 46;
                System.out.println("replacing with natural from " + hr[i].signalStartFrame * 240 + " to " + hr[i].signalStopFrame * 240 + " lab=" + hr[i].wav + " ws=" + hr[i].wavStartFrame + " skip=" + skip);
                DataInputStream dis = new DataInputStream(new FileInputStream("normalized_corpus/" + hr[i].wav + ".wav"));
                dis.skip(skip);
                byte[] tmp = new byte[(hr[i].signalStopFrame - hr[i].signalStartFrame + 1) * 240 * 2];
                dis.read(tmp);
                dis.close();
                for (int z = 0; z < tmp.length / 2; z++) {
                    byte b1 = tmp[z * 2];
                    byte b2 = tmp[z * 2 + 1];
                    int val = (b2 << 8) + b1;
                    double d = (double) val / 32767;
                    sign.add(d * 0.5);
                }
            }
        }
        double[] tmp = new double[sign.size()];
        for (int i = 0; i < tmp.length; i++) {
            tmp[i] = sign.get(i);
        }
        return tmp;
    }
    double ARITFICIAL_WEIGHT = 0.1;

    private HybridResult combineUnits(ViterbiToken artif, ViterbiToken natural, double[][] tsp, double[][] tap, double[] f0Raw) throws IOException {
        HybridResult hr = new HybridResult();
        int aStart = -artif.startFrame;
        int aStop = -artif.stopFrame;
        double[] f0Artif = new double[aStop - aStart + 1];
        double[][] spArtif = new double[aStop - aStart + 1][1025];
        double[][] apArtif = new double[aStop - aStart + 1][1025];
        for (int i = aStart; i <= aStop; i++) {
            for (int j = 0; j < 1025; j++) {
                spArtif[i - aStart][j] = tsp[j][i];
                apArtif[i - aStart][j] = tap[j][i];
            }
            f0Artif[i - aStart] = f0Raw[i];
        }
        double[] f0Natural = readArray(models + "/hybrid.f0", natural.startFrame, natural.stopFrame);
        double[][] spNatural = readArrays(models + "/hybrid.sp", natural.startFrame, natural.stopFrame);
        double[][] apNatural = readArrays(models + "/hybrid.ap", natural.startFrame, natural.stopFrame);
        //int[][] f0Pairs = align(f0Natural, f0Artif);
        int[][] pairs = align(spArtif, spNatural);

        for (int i = 0; i < pairs.length; i++) {
            int iArtif = pairs[i][0];
            int iNatural = pairs[i][1];
            for (int j = 0; j < 1025; j++) {
                spNatural[iNatural][j] = ARITFICIAL_WEIGHT * spArtif[iArtif][j] + (1.0 - ARITFICIAL_WEIGHT) * spNatural[iNatural][j];
                apNatural[iNatural][j] = ARITFICIAL_WEIGHT * apArtif[iArtif][j] + (1.0 - ARITFICIAL_WEIGHT) * apNatural[iNatural][j];
            }
            if (f0Artif[iArtif] != 0 && f0Natural[iNatural] != 0) {
                f0Natural[iNatural] = ARITFICIAL_WEIGHT * f0Artif[iArtif] + (1.0 - ARITFICIAL_WEIGHT) * f0Natural[iNatural];
            }
        }
//        for (int i = 0; i < f0Artif.length; i++) {
//            f0Artif[i] += 250;
//        }
        hr.f0 = f0Natural;
        hr.tap = apNatural;
        hr.tsp = spNatural;
        return hr;
    }

    private int[][] align(double[][] a, double[][] b) {
        double[][] x = new double[a.length + 1][b.length + 1];
        double cc = 0;
        for (int i = 1; i < x.length; i++) {
            cc = 0;
            for (int zz = 0; zz < b[0].length; zz++) {
                cc += Math.abs(a[i - 1][zz]);
            }
            x[i][0] = cc + x[i - 1][0];
        }
        for (int i = 1; i < x[0].length; i++) {
            cc = 0;
            for (int zz = 0; zz < b[0].length; zz++) {
                cc += Math.abs(b[i - 1][zz]);
            }
            x[0][i] = cc + x[0][i - 1];
        }

        for (int i = 1; i < x.length; i++) {
            for (int j = 1; j < x[0].length; j++) {
                double min = x[i - 1][j - 1];
                if (x[i - 1][j] < min) {
                    min = x[i - 1][j];
                }
                if (x[i][j - 1] < min) {
                    min = x[i][j - 1];
                }
                double cost = 0;
                for (int zz = 0; zz < b[0].length; zz++) {
                    cost += Math.abs(a[i - 1][zz] - b[j - 1][zz]);
                }
                x[i][j] = min + cost;
            }
        }
        int i = x.length - 1;
        int j = x[0].length - 1;
        List<int[]> pairs = new ArrayList<>();
        while (i != 1 || j != 1) {
            int[] pair = new int[2];
            pair[0] = i - 1;
            pair[1] = j - 1;
            pairs.add(0, pair);
            if (i == 1) {
                j--;
            } else if (j == 1) {
                i--;
            } else if (x[i - 1][j - 1] <= x[i - 1][j] && x[i - 1][j - 1] <= x[i][j - 1]) {
                i--;
                j--;
            } else if (x[i - 1][j] <= x[i - 1][j - 1] && x[i - 1][j] <= x[i][j - 1]) {
                i--;
            } else {
                j--;
            }
        }

        return pairs.toArray(new int[0][]);
    }

    private int[][] align(double[] a, double[] b) {
//        for (int i = 0; i < a.length; i++) {
//            if (a[i] != 0) {
//                a[i] = Math.log(a[i]);
//            }
//        }
//        for (int i = 0; i < b.length; i++) {
//            if (b[i] != 0) {
//                b[i] = Math.log(b[i]);
//            }
//        }
        double[][] x = new double[a.length + 1][b.length + 1];
        for (int i = 1; i < x.length; i++) {
            x[i][0] = Math.pow(a[i - 1], 2) + x[i - 1][0];
        }
        for (int i = 1; i < x[0].length; i++) {
            x[0][i] = Math.pow(b[i - 1], 2) + x[0][i - 1];
        }

        for (int i = 1; i < x.length; i++) {
            for (int j = 1; j < x[0].length; j++) {
                double min = x[i - 1][j - 1];
                if (x[i - 1][j] < min) {
                    min = x[i - 1][j];
                }
                if (x[i][j - 1] < min) {
                    min = x[i][j - 1];
                }
                double cost = Math.pow(a[i - 1] - b[j - 1], 2);
                x[i][j] = min + cost;
            }
        }

        int i = x.length - 1;
        int j = x[0].length - 1;
        List<int[]> pairs = new ArrayList<>();
        while (i != 1 || j != 1) {
            int[] pair = new int[2];
            pair[0] = i - 1;
            pair[1] = j - 1;
            pairs.add(0, pair);
            if (i == 1) {
                j--;
            } else if (j == 1) {
                i--;
            } else if (x[i - 1][j - 1] <= x[i - 1][j] && x[i - 1][j - 1] <= x[i][j - 1]) {
                i--;
                j--;
            } else if (x[i - 1][j] <= x[i - 1][j - 1] && x[i - 1][j] <= x[i][j - 1]) {
                i--;
            } else {
                j--;
            }
        }

        return pairs.toArray(new int[0][]);
    }

    private class HybridResult {

        double[][] tsp;
        double[][] tap;
        double[] f0;
    }

    int l = 1024;

    private double[][] mgc2sp(double[][] mgc) {
        float a = 0.55f;
        float g = 0f;
        double[][] sp = new double[mgc.length][l + 1];
        for (int i = 0; i < mgc.length; i++) {
            double[] x = new double[l * 2];
            double[] y = new double[l * 2];
            double[] c = new double[l / 2 + 1];//aici merge si cu l direct
            //MLSAVocoder.ignorm(mgc[i], mgc[i], mgc[i].length - 1, g);
            //incercam niste postfiltering pe taran
            //MLSAVocoder.postfilter_mgc(mgc[i], mgc[i].length - 1, 1.4, a);
            MLSAVocoder.mgc2mgc(mgc[i], mgc[i].length - 1, a, g, c, c.length - 1, 0, 0);
            c2sp(c, c.length - 1, x, y, l / 2);

            for (int k = 0; k < l / 2 + 1; k++) {
                sp[i][k] = Math.exp(x[k]);
                sp[i][l - k] = sp[i][k];
            }
            //sp[i][0] = 0;
        }
        return sp;
    }
    FFT fft = FFT.get(l * 2);

    void c2sp(double[] c, int m, double[] x, double[] y, int l) {
        //movem(c, x, sizeof( * c), m1);
        //fillz(x + m1, sizeof( * x), l - m1);
        System.arraycopy(c, 0, x, 0, Math.min(x.length, c.length));
        fft.fft(x, y);

        //fftr(x, y, l);
    }

    private HybridResult hybridify(double[][] tsp, double[][] tap, double[] f0Raw, String[] phonemes, int[] durations) {
        HybridResult hr = new HybridResult();
        try {
            if (new File(models + "/hybrid.sp").exists()) {
                BufferedWriter bw = new BufferedWriter(new FileWriter("frames.trace"));
//                //megatest
//                double[][] xsp = readArrays(models + "/hybrid.sp", 3307, 5620);
//                double[][] xap = readArrays(models + "/hybrid.ap", 3307, 5620);
//                double[] xf0 = readArray(models + "/hybrid.f0", 3307, 5620);
//
//                hr.tsp = new double[1025][xsp.length];
//                hr.tap = new double[1025][xsp.length];
//                hr.f0 = new double[xsp.length];
//                for (int i = 0; i < xsp.length; i++) {
//                    for (int j = 0; j < 1025; j++) {
//                        hr.tsp[j][i] = xsp[i][j];
//                        hr.tap[j][i] = xap[i][j];
//                    }
//                    hr.f0[i] = xf0[i];
//                }
//                if (true) {
//                    return hr;
//                }

                System.out.println("Found hybrid models...");
                Map<String, List<PhoneInfo>> realPhones = readPhones(models + "/hybrid.phs");
                ViterbiToken[][] vtList;
                ViterbiToken[] vtArtif;
                vtList = new ViterbiToken[phonemes.length][];
                vtArtif = new ViterbiToken[phonemes.length];

                int lastStart = -1;
                for (int i = 0; i < phonemes.length; i++) {
                    String cp = phonemes[i];
                    String lab = cp;
                    ViterbiToken vt = new ViterbiToken();
                    vt.emissionCost = 0;
                    if (i == 0) {
                        vt.totalCost = vt.emissionCost;
                    } else {
                        vt.totalCost = Double.POSITIVE_INFINITY;
                    }
                    vt.lab = lab;
                    vt.bestFrom = -1;
                    vt.startFrame = lastStart;
                    vt.stopFrame = lastStart - durations[i];
                    if (i == 0) {
                        vt.stopFrame += 1;
                    }
                    lastStart = vt.stopFrame;
                    ViterbiToken[] tokens = getCompatibleTokens(vt, tsp, f0Raw, realPhones);
                    vtArtif[i] = vt;
                    vtList[i] = new ViterbiToken[tokens.length];
                    for (int j = 0; j < tokens.length; j++) {
                        vtList[i][j] = tokens[j];
                        if (i == 0) {
                            tokens[j].totalCost = tokens[j].emissionCost;
                        }
                    }
                }

                for (int i = 1; i < vtList.length; i++) {
                    System.out.println(i + "/" + (vtList.length - 1));
                    for (int k = 0; k < vtList[i].length; k++) {

                        //double[][] spDst = readArrays(models + "/hybrid.sp", vtList[i][k].startFrame, vtList[i][k].startFrame + 1);
                        double[][] spDst;
                        double[] f0Dst;
                        if (vtList[i][k].startFrame > 0) {
                            spDst = readArrays(models + "/hybrid.sp", vtList[i][k].startFrame, vtList[i][k].startFrame + 1);
                            f0Dst = readArray(models + "/hybrid.f0", vtList[i][k].startFrame, vtList[i][k].startFrame + 4);
                        } else {
                            spDst = new double[2][1025];
                            f0Dst = new double[5];
                            for (int zz = 0; zz < 1025; zz++) {
                                spDst[0][zz] = tsp[zz][-vtList[i][k].startFrame];
                                spDst[1][zz] = tsp[zz][-vtList[i][k].stopFrame - 1];
                            }
                            for (int zz = 0; zz < f0Dst.length; zz++) {
                                f0Dst[zz] = f0Raw[-vtList[i][k].startFrame + zz];
                            }
                        }

                        for (int l = 0; l < Math.min(10000, vtList[i - 1].length); l++) {
                            if (vtList[i - 1][l].totalCost < vtList[i][k].totalCost) {
                                double trans = computeTransition(vtList[i - 1][l], spDst, f0Dst, tsp, f0Raw);
                                double cost = trans * LAMBDA_TRANSITION + vtList[i][k].emissionCost * LAMBDA_EMISSION + vtList[i - 1][l].totalCost;
                                if (cost < vtList[i][k].totalCost) {
                                    vtList[i][k].totalCost = cost;
                                    vtList[i][k].bestFrom = l;
                                }
                            }
                        }
                    }
//                    sort(vtList[i], new Comparator<ViterbiToken>() {
//                        @Override
//                        public int compare(ViterbiToken o1, ViterbiToken o2) {
//                            Double d1 = o1.totalCost;
//                            return d1.compareTo(o2.totalCost);
//                        }
//                    });
                }
                //merge invers
                double bestScore = Double.POSITIVE_INFINITY;
                int bestIndex = 0;
                for (int i = 0; i < vtList[vtList.length - 1].length; i++) {
                    if (vtList[vtList.length - 1][i].totalCost < bestScore) {
                        bestScore = vtList[vtList.length - 1][i].totalCost;
                        bestIndex = i;
                    }
                }
                System.out.println("Best decoding path score: " + bestScore);
                int pos = vtList.length - 1;
                List<ViterbiToken> tmp = new ArrayList<>();
                int total_frames = 0;
                while (pos >= 0) {
                    System.out.println(pos + " " + bestIndex);
                    total_frames += Math.abs(vtList[pos][bestIndex].stopFrame - vtList[pos][bestIndex].startFrame) + 1;
                    tmp.add(0, vtList[pos][bestIndex]);
                    bestIndex = vtList[pos][bestIndex].bestFrom;
                    pos--;
                }
                hr.tap = new double[1025][total_frames];
                hr.tsp = new double[1025][total_frames];
                hr.f0 = new double[total_frames];
                pos = 0;
                List<Integer> smoothPoints = new ArrayList<>();
                for (int i = 0; i < tmp.size(); i++) {
                    double[][] sp;
                    double[][] ap;
                    double[] f0;
                    ViterbiToken vt = tmp.get(i);
                    boolean synth = false;
                    if (vt.startFrame >= 0 && !vt.lab.equals("pau")) {
//                        sp = readArrays(models + "/hybrid.sp", vt.startFrame, vt.stopFrame);
//                        ap = readArrays(models + "/hybrid.ap", vt.startFrame, vt.stopFrame);
//                        f0 = readArray(models + "/hybrid.f0", vt.startFrame, vt.stopFrame);
                        HybridResult combined = combineUnits(vtArtif[i], vt, tsp, tap, f0Raw);
                        sp = combined.tsp;
                        ap = combined.tap;
                        f0 = combined.f0;
                        synth = false;
                        for (int k = 0; k < f0.length; k++) {
                            bw.write(k + "\t" + vt.lab + "\t" + f0[k] + "\t" + (vt.wavStartFrame + k) * 480 + "\n");
                        }
                    } else {
                        synth = true;
                        int start = -vtArtif[i].startFrame;
                        int stop = -vtArtif[i].stopFrame;
                        sp = new double[stop - start][1025];
                        ap = new double[stop - start][1025];
                        f0 = new double[stop - start];
                        for (int zz = 0; zz < sp.length; zz++) {
                            for (int ll = 0; ll < 1025; ll++) {
                                sp[zz][ll] = tsp[ll][start + zz];
                                ap[zz][ll] = tap[ll][start + zz];
                            }
                            f0[zz] = f0Raw[start + zz];
                        }
                        for (int k = 0; k < f0.length; k++) {
                            bw.write(k + "\t" + vt.lab + "\t" + f0[k] + "\t" + (vt.wavStartFrame + k) * 480 + "\t");
                            for (int l = 0; l < 1025; l++) {
                                bw.write(tsp[l][start + k] + ((l == 1024) ? "" : " "));
                            }
                            bw.write("\n");
                        }
                    }
                    for (int k = 0; k < sp.length; k++) {
                        for (int l = 0; l < 1025; l++) {
                            hr.tap[l][pos + k] = ap[k][l];
                            hr.tsp[l][pos + k] = sp[k][l];
                        }
                        if (pos + k > 0) {
                            hr.f0[pos + k] = f0[k];
                        }
                    }
                    //smootify
                    if (synth) {
                        if (!smoothPoints.contains(pos)) {
                            smoothPoints.add(pos);
                        }
                    }
                    pos += sp.length;
                    if (synth) {
                        smoothPoints.add(pos);
                    }
                }
                for (int z = 0; z < smoothPoints.size(); z++) {
                    pos = smoothPoints.get(z);
                    if (pos > 0) {
                        for (int k = 0; k < 1; k++) {
                            for (int frame = pos - 1; frame <= pos + 1; frame++) {
                                if (frame > 1 && frame < hr.tsp.length) {
                                    for (int ll = 0; ll < 1025; ll++) {
                                        hr.tap[ll][frame] = (hr.tap[ll][frame - 1] + hr.tap[ll][frame + 1]) / 2;
                                        hr.tsp[ll][frame] = (hr.tsp[ll][frame - 1] + hr.tsp[ll][frame + 1]) / 2;
                                    }
                                    if (hr.f0[frame - 1] != 0 && hr.f0[frame + 1] != 0) {
                                        hr.f0[frame] = (hr.f0[frame - 1] + hr.f0[frame + 1]) / 2;
                                    }
                                }
                            }
                        }
                    }
                }
                bw.close();
            } else {
                hr.f0 = f0Raw;
                hr.tap = tap;
                hr.tsp = tsp;
                return hr;
            }
        } catch (Exception e) {
            e.printStackTrace();
            //System.out.println();
            hr.f0 = f0Raw;
            hr.tap = tap;
            hr.tsp = tsp;
            return hr;
            //return null;
        }
        //interpolam f0 pentru toate frame-urile
        double last_f0 = 0;
        double next_f0 = 0;
        int dist_last = 0;
        int dist_next = 0;
        double[] nf0 = new double[hr.f0.length];
        for (int i = 0; i < hr.f0.length; i++) {
            if (hr.f0[i] == 0) {
                last_f0 = 0;
                next_f0 = 0;
                for (int k = i - 1; k >= 0; k--) {
                    if (hr.f0[k] != 0) {
                        last_f0 = hr.f0[k];
                        break;
                    } else {
                        dist_last++;
                    }
                }

                for (int k = i + 1; k < hr.f0.length; k++) {
                    if (hr.f0[k] != 0) {
                        next_f0 = hr.f0[k];
                        break;
                    } else {
                        dist_next++;
                    }
                }
                if (last_f0 == 0) {
                    hr.f0[i] = next_f0;
                } else if (next_f0 == 0) {
                    hr.f0[i] = last_f0;
                } else {
//                if (last_f0!=0 && next)
                    double total_dist = dist_last + dist_next;
                    double p0 = 1.0 - (double) (dist_last) / total_dist;
                    double p1 = 1.0 - (double) (dist_next) / total_dist;
                    nf0[i] = p0 * last_f0 + p1 * next_f0;
                }
            } else {
                nf0[i] = hr.f0[i];
            }
        }
        hr.f0 = nf0;
        return hr;
    }

    private class ViterbiToken {

        String wav;
        String lab;
        int startFrame;
        int stopFrame;
        double emissionCost;
        double totalCost;
        int bestFrom;
        int signalStartFrame = 0;
        int signalStopFrame = 0;
        private long wavStartFrame;
    }

    private double[][] readArrays(String file, int start, int stop) throws FileNotFoundException, IOException {
        double[][] tmp = new double[stop - start + 1][1025];
        FileInputStream fis = new FileInputStream(file);
        DataInputStream dis = new DataInputStream(fis);
        int ofs = start * 1025 * 4;
        dis.skipBytes(ofs);
        int len = (tmp.length) * 4 * 1025;
        byte[] data = new byte[len];
        if (dis.read(data) != data.length) {
            System.out.println("errrrrr");
        }
        ByteBuffer bb = ByteBuffer.wrap(data);
        for (int i = 0; i < tmp.length; i++) {
            for (int j = 0; j < 1025; j++) {
                float f = bb.getFloat();
                tmp[i][j] = f;
            }
        }
        dis.close();
        fis.close();

        return tmp;
    }

    private double[] readArray(String file, int start, int stop) throws FileNotFoundException, IOException {
        double[] tmp = new double[stop - start + 1];
        FileInputStream fis = new FileInputStream(file);
        DataInputStream dis = new DataInputStream(fis);
        int ofs = start * 4;
        dis.skipBytes(ofs);
        int len = (tmp.length) * 4;
        byte[] data = new byte[len];
        dis.read(data);
        ByteBuffer bb = ByteBuffer.wrap(data);
        for (int i = 0; i < tmp.length; i++) {
            float f = bb.getFloat();
            tmp[i] = f;
        }
        dis.close();
        fis.close();
        return tmp;
    }

//    private Map<String, List<int[]>> readPhones(String file) throws IOException {
//        BufferedReader br = new BufferedReader(new FileReader(file));
//        Map<String, List<int[]>> map = new TreeMap<>();
//        String line;
//        while ((line = br.readLine()) != null) {
//            String[] parts = line.split("\t");
//            String phoneme = parts[1];
//            int startFrame = Integer.parseInt(parts[2]);
//            int stopFrame = Integer.parseInt(parts[3]);
//            if (!map.containsKey(phoneme)) {
//                map.put(phoneme, new ArrayList<int[]>());
//            }
//            List<int[]> value = map.get(phoneme);
//            int[] x = new int[2];
//            x[0] = startFrame;
//            x[1] = stopFrame;
//            value.add(x);
//        }
//        br.close();
//        return map;
//    }
    private int bestMatch(double[][] tsp, double[][] tap, double[] f0Raw, int pos, int duration, Map<String, List<int[]>> realPhones, String phoneme, String models) throws IOException {
        //convert to array spectrum
        if (!realPhones.containsKey(phoneme)) {
            return 0;
        }
        List<int[]> indices = realPhones.get(phoneme);
        double[][] spect = new double[duration][1025];
        for (int i = 0; i < duration; i++) {
            for (int j = 0; j < 1025; j++) {
                spect[i][j] = tsp[j][i + pos];
            }
        }
        int bestMatch = -1;
        double bestCost = Double.POSITIVE_INFINITY;
        for (int i = 0; i < indices.size(); i++) {
            if ((indices.get(i)[1] - indices.get(i)[0]) != spect.length) {
                continue;
            }
            double[][] realSpect = readArrays(models + "/hybrid.sp", indices.get(i)[0], indices.get(i)[1]);
            double cost = 0;
            for (int k = 0; k < spect.length; k++) {
                for (int l = 0; l < 1025; l++) {
                    cost += Math.abs(realSpect[k][l] - spect[k][l]);
                }
            }
            if (cost < bestCost) {
                bestCost = cost;
                bestMatch = i;
            }
        }

        if (bestMatch == -1) {
            return 0;
        }
        System.out.print("..." + (bestCost / duration) + "...");
        if ((bestCost / duration) > 30) {
            return 0;
        }
        int startF = indices.get(bestMatch)[0];
        int stopF = indices.get(bestMatch)[0] + duration;
        double[][] sp = readArrays(models + "/hybrid.sp", startF, stopF);
        double[][] ap = readArrays(models + "/hybrid.ap", startF, stopF);
        double[] f0 = readArray(models + "/hybrid.f0", startF, stopF);
        for (int i = 0; i < duration; i++) {
            for (int j = 0; j < 1025; j++) {
                tsp[j][i + pos] = sp[i][j];
                tap[j][i + pos] = ap[i][j];
            }
            f0Raw[i + pos] = (f0Raw[i + pos] + f0[i]) / 2;
        }
        return duration;
    }

    class PhoneInfo {

        String lab;
        String wav;
        int startFrame;
        int stopFrame;
        int wavStartFrame;
    }
}
