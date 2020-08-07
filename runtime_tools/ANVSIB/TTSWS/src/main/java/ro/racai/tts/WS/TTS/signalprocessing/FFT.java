/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.TTS.signalprocessing;

import java.util.HashMap;
import java.util.Map;

public class FFT {

    int n, m;

    // Lookup tables.  Only need to recompute when size of FFT changes.
    double[] cos;
    double[] sin;

    double[] window;

    static Map<Integer, FFT> ffts=new HashMap<Integer, FFT>();
    
    public static FFT get(int n){
        if(!ffts.containsKey(n)){
            ffts.put(n,new FFT(n));
        }
        return ffts.get(n);
    }
    FFT(int n) {
        System.out.println("FFT "+n);
        this.n = n;
        this.m = (int) (Math.log(n) / Math.log(2));
        
        // Make sure n is a power of 2
        if (n != (1 << m)) {
            throw new RuntimeException("FFT length must be power of 2");
        }

        cos = new double[n / 2];
        sin = new double[n / 2];

        for (int i = 0; i < n / 2; i++) {
            cos[i] = Math.cos(-2.0 * Math.PI * i / n);
            sin[i] = Math.sin(-2.0 * Math.PI * i / n);
        }
    }

    public void ifft(double[] x, double[] y) {
        //conjugam numarul complex
        for (int i = 0; i < y.length; i++) {
            y[i] = -y[i];
        }
        fft(x, y);
        for (int i = 0; i < x.length; i++) {
            x[i] = x[i] / x.length;
            y[i] = -y[i] / y.length;
        }
    }

    public double[] abs(double[] real, double[] imag) {
        double[] ret = new double[real.length];
        for (int i = 0; i < real.length; i++) {
            ret[i] = Math.sqrt(real[i] * real[i] + imag[i] * imag[i]);
        }
        return ret;
    }

    public void fft(double[] x, double[] y) {
        int i, j, k, n1, n2, a;
        double c, s, e, t1, t2;

        // Bit-reverse
        j = 0;
        n2 = n / 2;
        for (i = 1; i < n - 1; i++) {
            n1 = n2;
            while (j >= n1) {
                j = j - n1;
                n1 = n1 / 2;
            }
            j = j + n1;

            if (i < j) {
                t1 = x[i];
                x[i] = x[j];
                x[j] = t1;
                t1 = y[i];
                y[i] = y[j];
                y[j] = t1;
            }
        }

        // FFT
        n1 = 0;
        n2 = 1;

        for (i = 0; i < m; i++) {
            n1 = n2;
            n2 = n2 + n2;
            a = 0;

            for (j = 0; j < n1; j++) {
                c = cos[a];
                s = sin[a];
                a += 1 << (m - i - 1);

                for (k = j; k < n; k = k + n2) {
                    t1 = c * x[k + n1] - s * y[k + n1];
                    t2 = s * x[k + n1] + c * y[k + n1];
                    x[k + n1] = x[k] - t1;
                    y[k + n1] = y[k] - t2;
                    x[k] = x[k] + t1;
                    y[k] = y[k] + t2;
                }
            }
        }
    }

//    int fftr(double[] x, double[] y) {
//        int i, j;
//        int xp, yp, xq;
//        int yq;
//        int mv2, n, tblsize;
//        double xt, yt;
//        int sinp, cosp;
//        double arg;
//        int m = x.length;
//
//        mv2 = m / 2;
//
//        /* separate even and odd  */
//        xq = xp = 0;
//        yp = 0;
//        for (i = mv2; --i >= 0;) {
//            x[xp++] = x[xq++];
//            y[yp++] = x[xq++];
//        }
//        double[] xtemp = new double[mv2];
//        double[] ytemp = new double[mv2];
//        System.arraycopy(x, 0, xtemp, 0, mv2);
//        System.arraycopy(y, 0, ytemp, 0, mv2);
//        fft(xtemp, ytemp);
//
//        n = maxfftsize / m;
//        sinp = _sintbl;
//        cosp = _sintbl + maxfftsize / 4;
//
//        xp = x;
//        yp = y;
//        xq = xp + m;
//        yq = yp + m;
//         * (xp + mv2) =  * xp -  * yp;
//         * xp =  * xp +  * yp;
//         * (yp + mv2) =  * yp = 0;
//
//        for (i = mv2, j = mv2 - 2; --i; j -= 2) {
//            ++xp;
//            ++yp;
//            sinp += n;
//            cosp += n;
//            yt =  * yp +  * (yp + j);
//            xt =  * xp -  * (xp + j);
//             * (--xq) = ( * xp +  * (xp + j) +  * cosp * yt -  * sinp * xt) * 0.5;
//             * (--yq) = ( * (yp + j) -  * yp +  * sinp * yt +  * cosp * xt) * 0.5;
//        }
//
//        xp = x + 1;
//        yp = y + 1;
//        xq = x + m;
//        yq = y + m;
//
//        for (i = mv2; --i;) {
//             * xp++ =  * (--xq);
//             * yp++ = -( * (--yq));
//        }
//
//        return (0);
//    }

}
