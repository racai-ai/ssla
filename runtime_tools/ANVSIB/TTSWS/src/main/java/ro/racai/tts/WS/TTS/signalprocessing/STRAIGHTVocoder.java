/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.TTS.signalprocessing;

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
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import ro.racai.tts.WS.TTS.straight.straightSynth;

/**
 *
 * @author tibi
 */
public class STRAIGHTVocoder {

    public STRAIGHTVocoder() {
    }

    public double[] runSTRAIGHTVocoder(float[][] mgcf, float[][] apf, float[] f0Rawf, int sr) throws Exception {
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

        double[] signal = straightSynth.straightSynth(tsp, f0Raw, shiftm, fs, pconv, fconv, sconv, gdbw, delfrac, delsp, cornf, delfracind, tap, imap, lowestF0);
        //signal = remodelUsingNaturalUnits(signal, hr);
        return signal;
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
}
