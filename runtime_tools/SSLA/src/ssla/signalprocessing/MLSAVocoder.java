/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ssla.signalprocessing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author tibi
 */
public class MLSAVocoder {

    public static final int IPERIOD = 1; /* interpolation period */

    public static final int SEED = 1;
    public static final int PADEORDER = 5; /* pade order for MLSA filter */

    public static final int IRLENG = 96; /* length of impulse response */

    public static final double ZERO = 1.0e-10; /* ~(0) */

    public static final double LZERO = (-1.0e+10); /* ~log(0) */
    /* ppade is a copy of pade in mlsadf() function : ppade = &( pade[pd*(pd+1)/2] ); */

    static final double[] pade = new double[]{ /* used in mlsadf */
        1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0.4999273, 0.1067005, 0.01170221, 0.0005656279, 1, 0.4999391, 0.1107098, 0.01369984,
        0.0009564853, 0.00003041721};
    static final int ppade = PADEORDER * (PADEORDER + 1) / 2;
    private Random rand;
    private int stage; /* Gamma=-1/stage : if stage=0 then Gamma=0 */

    private double gamma; /* Gamma */

    private boolean use_log_gain; /* log gain flag (for LSP) */

    private int fprd; /* frame shift */

    private double p1; /* used in excitation generation */

    private double pc; /* used in excitation generation */

    private double C[]; /* used in the MLSA/MGLSA filter */

    private double CC[]; /* used in the MLSA/MGLSA filter */

    private double CINC[]; /* used in the MLSA/MGLSA filter */

    private double D1[]; /* used in the MLSA/MGLSA filter */

    private double rate;
    int pt2; /* used in mlsadf2 */

    private final int pt3[] = new int[PADEORDER + 1]; /* used in mlsadf2 */

    private boolean mixedExcitation = false;
    private boolean fourierMagnitudes = false;

    /**
     * The initialisation of VocoderSetup should be done when there is already
     * information about the number of feature vectors to be processed, size of
     * the mcep vector file, etc.
     */
    private void initVocoder(int mcep_order, int mcep_vsize, int p, int sr) {
        stage = 0;
        gamma = 1;
        use_log_gain = false;
        fprd = p;
        rate = sr;
        rand = new Random(SEED);
        C = new double[mcep_order];
        CC = new double[mcep_order];
        CINC = new double[mcep_order];
        if (stage == 0) {
            int vector_size = (mcep_vsize * (3 + PADEORDER) + 5 * PADEORDER + 6) - (3 * (mcep_order));
            D1 = new double[vector_size];
            pt2 = (2 * (PADEORDER + 1)) + (PADEORDER * (mcep_order + 1));
            for (int i = PADEORDER; i >= 1; i--) {
                pt3[i] = (2 * (PADEORDER + 1)) + ((i - 1) * (mcep_order + 1));
            }
        } else { /* for LSP */

            int vector_size = ((mcep_vsize + 1) * (stage + 3)) - (3 * (mcep_order));
            D1 = new double[vector_size];
        }
        /* excitation initialisation */
        p1 = -1;
        pc = 0.0;
    } /* method initVocoder */


    public double[] runMLSAVocoder(float[][] mgc, float[] pitch, int mgcOrder, int p, int sr) throws Exception {
        //initVocoder(mgcOrder, mgcOrder - 1, p, sr);
        ArrayList<Double> audio = new ArrayList<>();
        double inc, x;
        double xp = 0.0, xn = 0.0, fxp, fxn, mix; /* samples for pulse and for noise and the filtered ones */

        int k, m, mcepframe, lf0frame;
        double alpha = 0.55f;
        //double beta = 1.0f - alpha;//probabil
        double[] magPulse = null; /* pulse generated from Fourier magnitudes */

        int magSample, magPulseSize;

        m = mgcOrder;
        initVocoder(m, mgcOrder - 1, p, sr);
        mixedExcitation = false;// de verificat
        fourierMagnitudes = true; //de verificat

        /* _______________________Synthesize speech waveforms_____________________ */
        /* generate Nperiod samples per mcepframe */
        int s = 0; /* number of samples */

        int s_double = 0;
        int audio_size = pitch.length * p; /* audio size in samples, calculated as num frames * frame period */

        double[] audio_double = new double[audio_size]; /* initialise buffer for audio */

        magSample = 1;
        magPulseSize = 0;
        for (mcepframe = 0; mcepframe < mgc.length; mcepframe++) { /* for each mcep frame */

            lf0frame = mcepframe;
            double mc[] = new double[m]; /* feature vector for a particular frame */

            for (int i = 0; i < m; i++) {
                mc[i] = mgc[mcepframe][i];
            }
            double f0 = pitch[lf0frame];

            double str = 0.0;
            if (f0 != 0.0) {
                f0 = rate / f0;
            }
            /* p1 is initialised in -1, so this will be done just for the first frame */
            if (p1 < 0) {
                p1 = f0;
                pc = p1;
                /* for LSP */
                if (stage != 0) {
                    C[0] = (use_log_gain) ? LZERO : ZERO;
                    double PI_m = Math.PI / m;
                    for (int i = 0; i < m; i++) {
                        C[i] = i * PI_m;
                    }
                    /* LSP -> MGC */
                    lsp2mgc(C, C, (m - 1), alpha);
                    mc2b(C, C, (m - 1), alpha);
                    gnorm(C, C, (m - 1), gamma);
                    for (int i = 1; i < m; i++) {
                        C[i] *= gamma;
                    }
                }
            }
            if (stage == 0) {
                mc2b(mc, CC, (m - 1), alpha);
                for (int i = 0; i < m; i++) {
                    CINC[i] = (CC[i] - C[i]) * IPERIOD / fprd;
                }
            } else {
                lsp2mgc(mc, CC, (m - 1), alpha);
                mc2b(CC, CC, (m - 1), alpha);
                gnorm(CC, CC, (m - 1), gamma);
                for (int i = 1; i < m; i++) {
                    CC[i] *= gamma;
                }
                for (int i = 0; i < m; i++) {
                    CINC[i] = (CC[i] - C[i]) * IPERIOD / fprd;
                }
            }
            /* p=f0 in c code!!! */
            if (p1 != 0.0 && f0 != 0.0) {
                inc = (f0 - p1) * (double) IPERIOD / (double) fprd;
            } else {
                inc = 0.0;
                pc = f0;
                p1 = 0.0;
            }

            for (int j = fprd - 1, i = (IPERIOD + 1) / 2; j >= 0; j--) {
                if (p1 == 0.0) {
                    x = uniformRand(); /* returns 1.0 or -1.0 uniformly distributed */

                    if (mixedExcitation) {
                        xn = x;
                        xp = 0.0;
                    }
                } else {
                    if ((pc += 1.0) >= p1) {

                        x = Math.sqrt(p1);
                        pc = pc - p1;
                    } else {
                        if (fourierMagnitudes) {
                            if (magSample >= magPulseSize) {
                                x = 0.0;
                            } else {
                                x = magPulse[magSample];
                            }
                            magSample++;
                        } else {
                            x = 0.0;
                        }
                    }
                    if (mixedExcitation) {
                        xp = x;
                        xn = uniformRand();
                    }
                }

                if (stage == 0) {
                    if (x != 0.0) {
                        x *= Math.exp(C[0]);
                    }
                    x = mlsadf(x, C, m, alpha, D1, pt2, pt3);
                } else {
                    x *= C[0];
                    x = mglsadf(x, C, (m - 1), alpha, stage, D1);
                }
// System.out.format("%f ", x);
                audio_double[s_double] = x;
                //if (audioProducer != null) {
                //    audioProducer.putOneDataPoint(x);
                //}
                s_double++;
                if ((--i) == 0) {
                    p1 += inc;
                    for (k = 0; k < m; k++) {
                        C[k] += CINC[k];
                    }
                    i = IPERIOD;
                }
            } /* for each sample in a period fprd */

            p1 = f0;
            System.arraycopy(CC, 0, C, 0, m);
        }
        return (audio_double);
    }

    /**
     * Compute the audio size, in samples, that this vocoder is going to produce
     * for the given data.
     *     
* @param mcepPst
     * @param htsData
     * @return
     */
    private void printVector(String val, int m, double vec[]) {
        int i;
        System.out.println(val);
        for (i = 0; i < m; i++) {
            System.out.println("v[" + i + "]=" + vec[i]);
        }
    }

    /**
     * mlsafir: sub functions for MLSA filter
     */
    private static double mlsafir(double x, double b[], int m, double a, double d[], int _pt3) {
        d[_pt3 + 0] = x;
        d[_pt3 + 1] = (1 - a * a) * d[_pt3 + 0] + (a * d[_pt3 + 1]);
        for (int i = 2; i <= m; i++) {
            d[_pt3 + i] += a * (d[_pt3 + i + 1] - d[_pt3 + i - 1]);
        }
        double y = 0.0;
        for (int i = 2; i <= m; i++) {
            y += d[_pt3 + i] * b[i];
        }
        for (int i = m + 1; i > 1; i--) {
            d[_pt3 + i] = d[_pt3 + i - 1];
        }
        return y;
    }

    /**
     * mlsdaf1: sub functions for MLSA filter
     */
    private static double mlsadf1(double x, double b[], int m, double a, double d[]) {
// pt1 --> pt = &d1[pd+1]
        double out = 0.0;
        for (int i = PADEORDER; i > 0; i--) {
            d[i] = (1 - a * a) * d[PADEORDER + i] + a * d[i];
            d[PADEORDER + 1 + i] = d[i] * b[1];
            double v = d[PADEORDER + 1 + i] * pade[ppade + i];
            x += ((1 & i) == 1) ? v : -v;
            /*
             * if(i == 1 || i == 3 || i == 5) x += v; else x += -v;
             */
            out += v;
        }
        d[PADEORDER + 1] = x;
        out += x;
        return out;
    }

    /**
     * mlsdaf2: sub functions for MLSA filter
     */
    private static double mlsadf2(double x, double b[], int m, double a, double d[], int pt2, int pt3[]) {
        double out = 0.0;
// pt2 --> pt = &d1[pd * (m+2)]
// pt3 --> pt = &d1[ 2*(pd+1) ]
        for (int i = PADEORDER; i > 0; i--) {
            int pt2_plus_i = pt2 + i;
            d[pt2_plus_i] = mlsafir(d[pt2_plus_i - 1], b, m, a, d, pt3[i]);
            double v = d[pt2_plus_i] * pade[ppade + i];
            x += ((1 & i) == 1) ? v : -v;
            /*
             * if(i == 1 || i == 3 || i == 5) x += v; else x += -v;
             */
            out += v;
        }
        d[pt2 /* +0 */] = x;
        out += x;
        return out;
    }

    /**
     * mlsadf: HTS Mel Log Spectrum Approximation filter
     */
    public static double mlsadf(double x, double b[], int m, double a, double d[], int pt2, int pt3[]) {
        x = mlsadf1(x, b, m, a, d);
        x = mlsadf2(x, b, m - 1, a, d, pt2, pt3);
        return x;
    }

    /**
     * uniform_rand: generate uniformly distributed random numbers 1 or -1
     */
    public double uniformRand() {
        return (rand.nextBoolean()) ? 1.0 : -1.0;
    }

    /**
     * mc2b: transform mel-cepstrum to MLSA digital filter coefficients
     */
    public static void mc2b(double mc[], double b[], int m, double a) {
        b[m] = mc[m];
        for (m--; m >= 0; m--) {
            b[m] = mc[m] - a * b[m + 1];
        }
    }

    /**
     * b2mc: transform MLSA digital filter coefficients to mel-cepstrum
     */
    public static void b2mc(double b[], double mc[], int m, double a) {
        double d = mc[m] = b[m];
        for (int i = m--; i >= 0; i--) {
            double o = b[i] + (a * d);
            d = b[i];
            mc[i] = o;
        }
    }

    /**
     * freqt: frequency transformation
     */
    public static void freqt(double c1[], int m1, double c2[], int m2, double a) {
        double b = 1 - a * a;
        double freqt_buff[] = new double[(m2 + m2 + 2)]; /* used in freqt */

        int g = m2 + 1; /* offset of freqt_buff */

        for (int i = -m1; i <= 0; i++) {
            if (0 <= m2) {
                freqt_buff[g + 0] = c1[-i] + a * (freqt_buff[0] = freqt_buff[g + 0]);
            }
            if (1 <= m2) {
                freqt_buff[g + 1] = b * freqt_buff[0] + a * (freqt_buff[1] = freqt_buff[g + 1]);
            }
            for (int j = 2; j <= m2; j++) {
                freqt_buff[g + j] = freqt_buff[j - 1] + a * ((freqt_buff[j] = freqt_buff[g + j]) - freqt_buff[g + j - 1]);
            }
        }
        /* move memory */
        System.arraycopy(freqt_buff, g, c2, 0, m2);
    }

    /**
     * c2ir: The minimum phase impulse response is evaluated from the minimum
     * phase cepstrum
     */
    public static void c2ir(double c[], int nc, double hh[], int leng) {
        hh[0] = Math.exp(c[0]);
        for (int n = 1; n < leng; n++) {
            double d = 0;
            int upl = (n >= nc) ? nc - 1 : n;
            for (int k = 1; k <= upl; k++) {
                d += k * c[k] * hh[n - k];
            }
            hh[n] = d / n;
        }
    }

    /**
     * b2en: functions for postfiltering
     */
    public static double b2en(double b[], int m, double a) {
        double cep[], ir[];
        int arrayLength = (m + 1) + 2 * IRLENG;
        double[] spectrum2en_buff = new double[arrayLength];
        cep = new double[arrayLength]; /* CHECK! these sizes!!! */

        ir = new double[arrayLength];
        b2mc(b, spectrum2en_buff, m, a);
        /* freqt(vs->mc, m, vs->cep, vs->irleng - 1, -a); */
        freqt(spectrum2en_buff, m, cep, IRLENG - 1, -a);
        /* HTS_c2ir(vs->cep, vs->irleng, vs->ir, vs->irleng); */
        c2ir(cep, IRLENG, ir, IRLENG);
        double en = 0.0;
        for (int i = 0; i < IRLENG; i++) {
            en += ir[i] * ir[i];
        }
        return en;
    }

    /**
     * ignorm: inverse gain normalization
     */
    public static void ignorm(double c1[], double c2[], int m, double ng) {
        if (ng != 0.0) {
            double k = Math.pow(c1[0], ng);
            for (int i = m; i >= 1; i--) {
                c2[i] = k * c1[i];
            }
            c2[0] = (k - 1.0) / ng;
        } else {
            /* movem */
            System.arraycopy(c1, 1, c2, 1, m - 1);
            c2[0] = Math.log(c1[0]);
        }
    }

    /**
     * ignorm: gain normalization
     */
    public static void gnorm(double c1[], double c2[], int m, double g) {
        if (g != 0.0) {
            double k = 1.0 + g * c1[0];
            for (; m >= 1; m--) {
                c2[m] = c1[m] / k;
            }
            c2[0] = Math.pow(k, 1.0 / g);
        } else {
            /* movem */
            System.arraycopy(c1, 1, c2, 1, m - 1);
            c2[0] = Math.exp(c1[0]);
        }
    }

    /**
     * lsp2lpc: transform LSP to LPC. lsp[1..m] --> a=lpc[0..m] a[0]=1.0
     */
    public static void lsp2lpc(double lsp[], double a[], int m) {
        int i, k, mh1, mh2, flag_odd;
        double xx, xf, xff;
        int p, q; /* offsets of lsp2lpc_buff */

        int a0, a1, a2, b0, b1, b2; /* offsets of lsp2lpc_buff */

        flag_odd = 0;
        if (m % 2 == 0) {
            mh1 = mh2 = m / 2;
        } else {
            mh1 = (m + 1) / 2;
            mh2 = (m - 1) / 2;
            flag_odd = 1;
        }
        double[] lsp2lpc_buff = new double[(5 * m + 6)];
        int lsp2lpc_size = m;
        /* offsets of lsp2lpcbuff */
        p = m;
        q = p + mh1;
        a0 = q + mh2;
        a1 = a0 + (mh1 + 1);
        a2 = a1 + (mh1 + 1);
        b0 = a2 + (mh1 + 1);
        b1 = b0 + (mh2 + 1);
        b2 = b1 + (mh2 + 1);
        /* move lsp -> lsp2lpc_buff */
        System.arraycopy(lsp, 1, lsp2lpc_buff, 0, m);
        for (i = 0; i < mh1 + 1; i++) {
            lsp2lpc_buff[a0 + i] = 0.0;
        }
        for (i = 0; i < mh1 + 1; i++) {
            lsp2lpc_buff[a1 + i] = 0.0;
        }
        for (i = 0; i < mh1 + 1; i++) {
            lsp2lpc_buff[a2 + i] = 0.0;
        }
        for (i = 0; i < mh2 + 1; i++) {
            lsp2lpc_buff[b0 + i] = 0.0;
        }
        for (i = 0; i < mh2 + 1; i++) {
            lsp2lpc_buff[b1 + i] = 0.0;
        }
        for (i = 0; i < mh2 + 1; i++) {
            lsp2lpc_buff[b2 + i] = 0.0;
        }
        /* lsp filter parameters */
        for (i = k = 0; i < mh1; i++, k += 2) {
            lsp2lpc_buff[p + i] = -2.0 * Math.cos(lsp2lpc_buff[k]);
        }
        for (i = k = 0; i < mh2; i++, k += 2) {
            lsp2lpc_buff[q + i] = -2.0 * Math.cos(lsp2lpc_buff[k + 1]);
        }
        /* impulse response of analysis filter */
        xx = 1.0;
        xf = xff = 0.0;
        for (k = 0; k <= m; k++) {
            if (flag_odd == 1) {
                lsp2lpc_buff[a0 + 0] = xx;
                lsp2lpc_buff[b0 + 0] = xx - xff;
                xff = xf;
                xf = xx;
            } else {
                lsp2lpc_buff[a0 + 0] = xx + xf;
                lsp2lpc_buff[b0 + 0] = xx - xf;
                xf = xx;
            }
            for (i = 0; i < mh1; i++) {
                lsp2lpc_buff[a0 + i + 1] = lsp2lpc_buff[a0 + i] + lsp2lpc_buff[p + i] * lsp2lpc_buff[a1 + i]
                        + lsp2lpc_buff[a2 + i];
                lsp2lpc_buff[a2 + i] = lsp2lpc_buff[a1 + i];
                lsp2lpc_buff[a1 + i] = lsp2lpc_buff[a0 + i];
            }
            for (i = 0; i < mh2; i++) {
                lsp2lpc_buff[b0 + i + 1] = lsp2lpc_buff[b0 + i] + lsp2lpc_buff[q + i] * lsp2lpc_buff[b1 + i]
                        + lsp2lpc_buff[b2 + i];
                lsp2lpc_buff[b2 + i] = lsp2lpc_buff[b1 + i];
                lsp2lpc_buff[b1 + i] = lsp2lpc_buff[b0 + i];
            }
            if (k != 0) {
                a[k - 1] = -0.5 * (lsp2lpc_buff[a0 + mh1] + lsp2lpc_buff[b0 + mh2]);
            }
            xx = 0.0;
        }
        for (i = m - 1; i >= 0; i--) {
            a[i + 1] = -a[i];
        }
        a[0] = 1.0;
    }

    /**
     * gc2gc: generalized cepstral transformation
     */
    public static void gc2gc(double c1[], int m1, double g1, double c2[], int m2, double g2) {
        double[] gc2gc_buff = Arrays.copyOf(c1, m1 + 1);
        c2[0] = gc2gc_buff[0];
        for (int i = 1; i <= m2; i++) {
            double ss1 = 0.0;
            double ss2 = 0.0;
            int min = m1 < i ? m1 : i - 1;
            for (int k = 1; k <= min; k++) {
                int mk = i - k;
                double cc = gc2gc_buff[k] * c2[mk];
                ss2 += k * cc;
                ss1 += mk * cc;
            }
            if (i <= m1) {
                c2[i] = gc2gc_buff[i] + (g2 * ss2 - g1 * ss1) / i;
            } else {
                c2[i] = (g2 * ss2 - g1 * ss1) / i;
            }
        }
    }

    /**
     * mgc2mgc: frequency and generalized cepstral transformation
     */
    public static void mgc2mgc(double c1[], int m1, double a1, double g1, double c2[], int m2, double a2, double g2) {
        if (a1 == a2) {
            gnorm(c1, c1, m1, g1);
            gc2gc(c1, m1, g1, c2, m2, g2);
            ignorm(c2, c2, m2, g2);
        } else {
            double a = (a2 - a1) / (1 - a1 * a2);
            freqt(c1, m1, c2, m2, a);
            gnorm(c2, c2, m2, g1);
            gc2gc(c2, m2, g1, c2, m2, g2);
            ignorm(c2, c2, m2, g2);
        }
    }

    /**
     * lsp2mgc: transform LSP to MGC. lsp=C[0..m] mgc=C[0..m]
     */
    public void lsp2mgc(double lsp[], double mgc[], int m, double alpha) {
        /* lsp2lpc */
        lsp2lpc(lsp, mgc, m); /* lsp starts in 1! lsp[1..m] --> mgc[0..m] */

        if (use_log_gain) {
            mgc[0] = Math.exp(lsp[0]);
        } else {
            mgc[0] = lsp[0];
        }
        /* mgc2mgc */
        ignorm(mgc, mgc, m, gamma);
        for (int i = m; i >= 1; i--) {
            mgc[i] *= -stage;
        }
        mgc2mgc(mgc, m, alpha, gamma, mgc, m, alpha, gamma); /* input and output is in mgc=C */

    }

    /**
     * mglsadff: sub functions for MGLSA filter
     */
    public static double mglsadf(double x, double b[], int m, double a, int n, double d[]) {
        for (int i = 0; i < n; i++) {
            x = mglsadff(x, b, m, a, d, (i * (m + 1)));
        }
        return x;
    }

    /**
     * mglsadf: sub functions for MGLSA filter
     */
    private static double mglsadff(double x, double b[], int m, double a, double d[], int d_offset) {
        double y = d[d_offset + 0] * b[1];
        for (int i = 1; i < m; i++) {
            d[d_offset + i] += a * (d[d_offset + i + 1] - d[d_offset + i - 1]);
            y += d[d_offset + i] * b[i + 1];
        }
        x -= y;
        for (int i = m; i > 0; i--) {
            d[d_offset + i] = d[d_offset + i - 1];
        }
        d[d_offset + 0] = a * d[d_offset + 0] + (1 - a * a) * x;
        return x;
    }

    /**
     * posfilter: postfilter for mel-cepstrum. It uses alpha and beta defined in
     * HMMData
     */
    public static void postfilter_mgc(double mgc[], int m, double alpha, double beta) {
        if (beta > 0.0 && m > 1) {
            double[] postfilter_buff = new double[m + 1];
            mc2b(mgc, postfilter_buff, m, alpha);
            double e1 = b2en(postfilter_buff, m, alpha);
            postfilter_buff[1] -= beta * alpha * mgc[2];
            for (int k = 2; k < m; k++) {
                postfilter_buff[k] *= (1.0 + beta);
            }
            double e2 = b2en(postfilter_buff, m, alpha);
            postfilter_buff[0] += Math.log(e1 / e2) / 2;
            b2mc(postfilter_buff, mgc, m, alpha);
        }
    }

//    public static double[] genPulseFromFourierMag(HTSPStream mag, int n, double f0) {
//        return genPulseFromFourierMag(mag.getParVec(n), f0);
//    }
    /**
     * Generate one pitch period from Fourier magnitudes
     */
    public static double[] genPulseFromFourierMag(double[] mag, double f0) {
        int numHarm = mag.length;
        int currentF0 = (int) Math.round(f0);
        int T;
        if (currentF0 < 512) {
            T = 512;
        } else {
            T = 1024;
        }
        int T2 = 2 * T;
        /* since is FFT2 no aperiodicFlag or jitter of 25% is applied */
        /* get the pulse */
        double[] pulse = new double[T];
        double[] real = new double[T2];
        double[] imag = new double[T2];
        /* copy Fourier magnitudes (Wai C. Chu "Speech Coding algorithms foundation and evolution of standardized coders" pg. 460) */
        real[0] = real[T] = 0.0; /* DC component set to zero */

        for (int i = 1; i <= numHarm; i++) {
            real[i] = real[T - i] = real[T + i] = real[T2 - i] = mag[i - 1]; /* Symetric extension */

            imag[i] = imag[T - i] = imag[T + i] = imag[T2 - i] = 0.0;
        }
        for (int i = (numHarm + 1); i < (T - numHarm); i++) { /* Default components set to 1.0 */

            real[i] = real[T - i] = real[T + i] = real[T2 - i] = 1.0;
            imag[i] = imag[T - i] = imag[T + i] = imag[T2 - i] = 0.0;
        }
        /* Calculate inverse Fourier transform */
        FFT fft = FFT.get(real.length);
        fft.ifft(real, imag);
        /* circular shift and normalise multiplying by sqrt(F0) */
        double sqrt_f0 = Math.sqrt(currentF0);
        for (int i = 0; i < T; i++) {
            pulse[i] = real[(i - numHarm) % T] * sqrt_f0;
        }
        return pulse;
    }

    private void circularShift(double y[], int T, int n) {
        double x[] = new double[T];
        for (int i = 0; i < T; i++) {
            x[i] = y[modShift(i - n, T)];
        }
        for (int i = 0; i < T; i++) {
            y[i] = x[i];
        }
    }

    private int modShift(int n, int N) {
        if (n < 0) {
            while (n < 0) {
                n = n + N;
            }
        } else {
            while (n >= N) {
                n = n - N;
            }
        }
        return n;
    }

    
}
