/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.TTS.straight;



import java.io.BufferedWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import org.apache.commons.math3.complex.Complex;
import ro.racai.tts.WS.TTS.utils.IO;

/**
 *
 * @author Echo
 */
public class M {

    static public double[][] transpose_lineArray(double[] a) {
        double[][] ret = new double[a.length][1];
        for (int i = 0; i < a.length; i++) {
            ret[i][0] = a[i];
        }
        return ret;
    }

    static public double[] transpose_columnArray(double[][] a) {
        double[] ret = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            ret[i] = a[i][0];
        }
        return ret;
    }

    static public double[][] transpose_matrix(double[][] a) {
        double[][] ret = new double[a[0].length][a.length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                ret[j][i] = a[i][j];
            }
        }
        return ret;
    }

    static public double[][] multiply_matrix(double[][] a, double[][] b) {
        double[][] c = new double[a.length][b[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b[0].length; j++) {
                for (int k = 0; k < a[0].length; k++) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return c;
    }

    static public double[] scalarProduct_lineArray(double[] a, double scalar) {
        double[] ret = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            ret[i] = a[i] * scalar;
        }
        return ret;
    }

    static public double[][] scalarProduct_matrix(double[][] a, double scalar) {
        double[][] ret = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                ret[i][j] = a[i][j] * scalar;
            }
        }
        return ret;
    }

    public static double round(double value) {
        return round(value, 0);
    }

    public static double round(double value, int places) {
        if (places < 0) {
            throw new IllegalArgumentException();
        }
        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

    static public double[] ComputeLP(int FilterOrder) {
        double[] NumCoeffs;
        int m;
        int i;

        NumCoeffs = new double[FilterOrder + 1];

        NumCoeffs[0] = 1;
        NumCoeffs[1] = FilterOrder;
        m = FilterOrder / 2;
        for (i = 2; i <= m; ++i) {
            NumCoeffs[i] = (double) (FilterOrder - i + 1) * NumCoeffs[i - 1] / i;
            NumCoeffs[FilterOrder - i] = NumCoeffs[i];
        }
        NumCoeffs[FilterOrder - 1] = FilterOrder;
        NumCoeffs[FilterOrder] = 1;

        return NumCoeffs;
    }

    static public double[] ComputeHP(int FilterOrder) {
        double[] NumCoeffs;
        int i;

        NumCoeffs = ComputeLP(FilterOrder);

        for (i = 0; i <= FilterOrder; ++i) {
            if (i % 2 != 0) {
                NumCoeffs[i] = -NumCoeffs[i];
            }
        }

        return NumCoeffs;
    }

    public static double[] TrinomialMultiply(int FilterOrder, double[] b,
            double[] c) {
        int i, j;
        double[] RetVal;

        RetVal = new double[4 * FilterOrder];

        RetVal[2] = c[0];
        RetVal[3] = c[1];
        RetVal[0] = b[0];
        RetVal[1] = b[1];

        for (i = 1; i < FilterOrder; ++i) {
            RetVal[2 * (2 * i + 1)] += c[2 * i] * RetVal[2 * (2 * i - 1)]
                    - c[2 * i + 1] * RetVal[2 * (2 * i - 1) + 1];
            RetVal[2 * (2 * i + 1) + 1] += c[2 * i]
                    * RetVal[2 * (2 * i - 1) + 1] + c[2 * i + 1]
                    * RetVal[2 * (2 * i - 1)];

            for (j = 2 * i; j > 1; --j) {
                RetVal[2 * j] += b[2 * i] * RetVal[2 * (j - 1)] - b[2 * i + 1]
                        * RetVal[2 * (j - 1) + 1] + c[2 * i]
                        * RetVal[2 * (j - 2)] - c[2 * i + 1]
                        * RetVal[2 * (j - 2) + 1];
                RetVal[2 * j + 1] += b[2 * i] * RetVal[2 * (j - 1) + 1]
                        + b[2 * i + 1] * RetVal[2 * (j - 1)] + c[2 * i]
                        * RetVal[2 * (j - 2) + 1] + c[2 * i + 1]
                        * RetVal[2 * (j - 2)];
            }

            RetVal[2] += b[2 * i] * RetVal[0] - b[2 * i + 1] * RetVal[1]
                    + c[2 * i];
            RetVal[3] += b[2 * i] * RetVal[1] + b[2 * i + 1] * RetVal[0]
                    + c[2 * i + 1];
            RetVal[0] += b[2 * i];
            RetVal[1] += b[2 * i + 1];
        }

        return RetVal;
    }

    public static double[] ComputeDenCoeffs(int FilterOrder, double Lcutoff,
            double Ucutoff) {
        int k; // loop variables
        double theta; // PI * (Ucutoff - Lcutoff) / 2.0
        double cp; // cosine of phi
        double st; // sine of theta
        double ct; // cosine of theta
        double s2t; // sine of 2*theta
        double c2t; // cosine 0f 2*theta
        double[] RCoeffs; // z^-2 coefficients
        double[] TCoeffs; // z^-1 coefficients
        double[] DenomCoeffs; // dk coefficients
        double PoleAngle; // pole angle
        double SinPoleAngle; // sine of pole angle
        double CosPoleAngle; // cosine of pole angle
        double a; // workspace variables

        cp = (double) Math.cos(Math.PI * (Ucutoff + Lcutoff) / 2.0f);
        theta = (double) Math.PI * (double) (Ucutoff - Lcutoff) / 2.0f;
        st = (double) Math.sin(theta);
        ct = (double) Math.cos(theta);
        s2t = 2.0f * st * ct; // sine of 2*theta
        c2t = 2.0f * ct * ct - 1.0f; // cosine of 2*theta

        RCoeffs = new double[2 * FilterOrder];
        TCoeffs = new double[2 * FilterOrder];

        for (k = 0; k < FilterOrder; ++k) {
            PoleAngle = (double) Math.PI * (double) (2 * k + 1)
                    / (double) (2 * FilterOrder);
            SinPoleAngle = (double) Math.sin(PoleAngle);
            CosPoleAngle = (double) Math.cos(PoleAngle);
            a = 1.0f + s2t * SinPoleAngle;
            RCoeffs[2 * k] = c2t / a;
            RCoeffs[2 * k + 1] = s2t * CosPoleAngle / a;
            TCoeffs[2 * k] = -2.0f * cp * (ct + st * SinPoleAngle) / a;
            TCoeffs[2 * k + 1] = -2.0f * cp * st * CosPoleAngle / a;
        }

        DenomCoeffs = TrinomialMultiply(FilterOrder, TCoeffs, RCoeffs);

        DenomCoeffs[1] = DenomCoeffs[0];
        DenomCoeffs[0] = 1.0f;
        for (k = 3; k <= 2 * FilterOrder; ++k) {
            DenomCoeffs[k] = DenomCoeffs[2 * k - 2];
        }

        return DenomCoeffs;
    }

    static public double[] ComputeNumCoeffs(int FilterOrder, double Lcutoff,
            double Ucutoff, double[] DenC) {
        double[] TCoeffs;
        double[] NumCoeffs;
        double[] Numbers = new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        int i;

        NumCoeffs = new double[2 * FilterOrder + 1];

        Complex[] NormalizedKernel = new Complex[2 * FilterOrder + 1];
        // double[] NormalizedKernel = new double[2*FilterOrder+1];

        TCoeffs = ComputeHP(FilterOrder);

        for (i = 0; i < FilterOrder; ++i) {
            NumCoeffs[2 * i] = TCoeffs[i];
            NumCoeffs[2 * i + 1] = 0.0f;
        }
        NumCoeffs[2 * FilterOrder] = TCoeffs[FilterOrder];
        double[] cp = new double[2];
        double Bw, Wn;
        cp[0] = 2 * 2.0f * (double) Math.tan(Math.PI * Lcutoff / 2.0);
        cp[1] = 2 * 2.0f * (double) Math.tan(Math.PI * Ucutoff / 2.0);

        Bw = cp[1] - cp[0];
        // center frequency
        Wn = (double) Math.sqrt(cp[0] * cp[1]);
        Wn = 2 * (double) Math.atan2(Wn, 4);
        double kern;
        Complex result = new Complex(-1, 0);
        // double result = -1;

        for (int k = 0; k < 11; k++) {
            NormalizedKernel[k] = ((Complex) (result.sqrt().multiply(Wn
                    * Numbers[k]))).exp();
        }
        double b = 0;
        double den = 0;
        for (int d = 0; d < 11; d++) {
            b += (NormalizedKernel[d].multiply(NumCoeffs[d])).getReal();
            den += (NormalizedKernel[d].multiply(DenC[d])).getReal();
        }
        for (int c = 0; c < 11; c++) {
            NumCoeffs[c] = (double) (NumCoeffs[c] * den) / b;
        }

        return NumCoeffs;
    }

    public static double[] filtera(double[] b, double[] a, double[] x) {
        double[] y = new double[x.length];

        for (int n = 0; n < y.length; n++) {
            if (n - 1 < 0) {
                y[n] = b[0] * x[n];
            } else {
                y[n] = b[0] * x[n] + b[1] * x[n - 1] - a[1] * y[n - 1];
            }

        }
        return y;
    }

    static public double[] filter2(double[] b, double[] a, double[] x) {
        int i, j;
        int ord = 1;
        double[] y = new double[x.length];
        y[0] = b[0] * x[0];
        for (i = 1; i < ord + 1; i++) {
            y[i] = 0.0f;
            for (j = 0; j < i + 1; j++) {
                y[i] = y[i] + b[j] * x[i - j];
            }
            for (j = 0; j < i; j++) {
                y[i] = y[i] - a[j + 1] * y[i - j - 1];
            }
        }
        /* end of initial part */
        for (i = ord + 1; i < x.length; i++) {
            y[i] = 0.0f;
            for (j = 0; j < ord + 1; j++) {
                y[i] = y[i] + b[j] * x[i - j];
            }
            for (j = 0; j < ord; j++) {
                y[i] = y[i] - a[j + 1] * y[i - j - 1];
            }
        }
        return y;
    }

    static public double[] filter(double[] b, double[] a, double[] input) {
        ArrayList<Double> inputVector = new ArrayList<Double>();
        ArrayList<Double> outputVector = new ArrayList<Double>();
        for (int i = 0; i < input.length; i++) {
            inputVector.add(input[i]);
            filter_in(b, a, inputVector, outputVector);
        }

        double[] ret = new double[outputVector.size()];

        for (int i = 0; i < outputVector.size(); i++) {
            ret[i] = outputVector.get(i);
        }

        return ret;
    }

    static public void filter_in(double[] b, double[] a,
            ArrayList<Double> inputVector, ArrayList<Double> outputVector) {
        double rOutputY = 0.0f;
        int j = 0;
        for (int i = 0; i < inputVector.size(); i++) {
            if (j < b.length) {
                rOutputY += b[j] * inputVector.get(inputVector.size() - i - 1);
            }
            j++;
        }
        j = 1;
        for (int i = 0; i < outputVector.size(); i++) {
            if (j < a.length) {
                rOutputY -= a[j]
                        * outputVector.get(outputVector.size() - i - 1);
            }
            j++;
        }
        outputVector.add(rOutputY);
    }

    static public double std(double[] x) {
        double mean = 0;
        double sum = 0;
        for (int i = 0; i < x.length; i++) {
            mean += x[i];
        }
        mean /= x.length;
        for (int i = 0; i < x.length; i++) {
            sum += (x[i] - mean) * (x[i] - mean);
        }
        return (double) Math.sqrt(1.0f / (x.length - 1) * sum);
    }

    static public double[] randn_1row_ncolumns(int cols) {
        double[] ret = new double[cols];
        Random r = new Random();
        for (int i = 0; i < cols; i++) {
            ret[i] = (double) r.nextGaussian();
        }
        return ret;
    }

    static public double[] optimumsmoothing(double eta, double pc) {
        /*
         * function ovc=optimumsmoothing(eta,pc) % ovc=optimumsmoothing(eta,pc)
         * % Calculate the optimum smoothing function % ovc : coefficients for
         * 2nd order cardinal B-spline % eta : temporal stretch factor % pc :
         * power exponent for nonlinearity
         */

        // fx=-8:0.05:8;
        // cb=max(0,1-abs(fx));
        double[] fx = new double[321];// asta e
        double[] cb = new double[fx.length];
        int cnt = 0;
        for (double i = -8; i <= 8.0005; i = i + 0.05f) {
            fx[cnt] = i;
            cb[cnt] = Math.max(0.0f, 1 - Math.abs(i));
            cnt++;
        }

        // gw=exp(-pi*(fx*eta).^2).^pc;
        double[] gw = new double[fx.length];
        for (int i = 0; i < fx.length; i++) {
            gw[i] = (double) Math.pow(
                    Math.exp(-Math.PI * (fx[i] * eta * fx[i] * eta)), pc);
        }

        // cmw=conv(cb,gw);
        double[] cmw = multiplyPolynomials(cb, gw);

        // bb=(1:length(cb));
        // bbc=bb+(length(cb)-1)/2;
        double[] bbc = new double[cb.length];
        for (int i = 0; i < cb.length; i++) {
            bbc[i] = (double) Math.ceil((i + 1.0f) + (cb.length - 1) / 2.0f);
        }

        // cmw=cmw(bbc)/max(cmw);
        double max_cmw = cmw[0];
        for (int i = 1; i < cmw.length; i++) {
            if (cmw[i] > max_cmw) {
                max_cmw = cmw[i];
            }
        }
        ArrayList<Double> cwtemp = new ArrayList<Double>(bbc.length);
        for (int i = 0; i < bbc.length; i++) {
            // System.out.println("bbc "+i+":"+bbc[i]);
            cwtemp.add(cmw[(int) bbc[i]] / max_cmw);
        }
        cmw = new double[cwtemp.size()];
        for (int i = 0; i < cmw.length; i++) {
            cmw[i] = cwtemp.get(i);
        }

        // ss=(abs(fx-round(fx))<0.025).*(1:length(cb));
        double[] ss = new double[cb.length];
        for (int i = 0; i < cb.length; i++) {
            ss[i] = (Math.abs(fx[i] - M.round(fx[i])) < 0.025) ? (i) : 0.0f;
            // aici e i nu i+1 ??!?
        }
        // ss=ss(ss>0);
        // cmws=cmw(ss);
        ArrayList<Double> cmws_array = new ArrayList<Double>();
        cmws_array.add(new Double(cmw[0]));
        for (int i = 0; i < ss.length; i++) {
            if (ss[i] > 0) {
                cmws_array.add(new Double(cmw[(int) ss[i] - 1]));
            }
        }
        double[] cmws = new double[cmws_array.size()];
        for (int i = 0; i < cmws_array.size(); i++) {
            cmws[i] = cmws_array.get(i).doubleValue();
        }

        // nn=length(cmws);
        int nn = cmws.length;

        // idv=1:nn;
        int[] idv = new int[nn];

        // hh=zeros(2*nn,nn);
        double hh[][] = new double[2 * nn][nn];

        // for ii=1:nn
        // hh((ii-1)+idv,ii)=cmws';
        // end;
        for (int ii = 0; ii < nn; ii++) {
            for (int line = 0; line < cmws.length; line++) {
                hh[line + ii][ii] = cmws[line];
            }
        }

        // bv=zeros(2*nn,1);
        // bv(nn+1)=1; % This is the original unit impulse.
        // %bv(nn)=0.04; % You can design the target function as you wish
        // %bv(nn+2)=0.04;
        double[][] bv = new double[2 * nn][1];
        bv[nn][0] = 1; // aici e nn nu nn+1 !

        // h=hh'*hh;
        double[][] h = new double[hh[0].length][hh[0].length];
        h = multiply_matrix(transpose_matrix(hh), hh);

        // ov=inv(h)*(hh'*bv); % This is the optimum coefficient vector.
        double[][] ov;
        ov = multiply_matrix(invert_matrix(h),
                multiply_matrix(transpose_matrix(hh), bv));

        // idc=(nn-1)/2+2;
        int idc = (nn - 1) / 2 + 2;

        // ovc=ov(idc+(0:3));
        idc--; // pt ca e zero index
        double[] ovc = new double[4];
        /*
         * for(int i = 0;i<ovc.length;i++) { int line = (int)((idc+i) /
         * ov[0].length); int col = (idc+i) % ov[0].length; ovc[i] =
         * ov[col][line]; }
         */
        ovc[0] = ov[idc][0];
        ovc[1] = ov[idc + 1][0];
        ovc[2] = ov[idc + 2][0];
        ovc[3] = ov[idc + 3][0];

        return ovc;
    }

    public static double[] multiplyPolynomials(double A[], double B[]) {
        double[] prod = new double[A.length + B.length - 1];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < B.length; j++) {
                prod[i + j] += A[i] * B[j];
            }
        }
        return prod;
    }

    public static double[][] invert_matrix(double a[][]) {
        int n = a.length;
        double x[][] = new double[n][n];
        double b[][] = new double[n][n];
        int index[] = new int[n];
        for (int i = 0; i < n; ++i) {
            b[i][i] = 1;
        }

        gaussian(a, index);

        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    b[index[j]][k] -= a[index[j]][i] * b[index[i]][k];
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            x[n - 1][i] = b[index[n - 1]][i] / a[index[n - 1]][n - 1];
            for (int j = n - 2; j >= 0; --j) {
                x[j][i] = b[index[j]][i];
                for (int k = j + 1; k < n; ++k) {
                    x[j][i] -= a[index[j]][k] * x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }

    public static void gaussian(double a[][], int index[]) {
        int n = index.length;
        double c[] = new double[n];

        for (int i = 0; i < n; ++i) {
            index[i] = i;
        }

        for (int i = 0; i < n; ++i) {
            double c1 = 0;
            for (int j = 0; j < n; ++j) {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) {
                    c1 = c0;
                }
            }
            c[i] = c1;
        }

        int k = 0;
        for (int j = 0; j < n - 1; ++j) {
            double pi1 = 0;
            for (int i = j; i < n; ++i) {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1) {
                    pi1 = pi0;
                    k = i;
                }
            }

            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i = j + 1; i < n; ++i) {
                double pj = a[index[i]][j] / a[index[j]][j];
                a[index[i]][j] = pj;
                for (int l = j + 1; l < n; ++l) {
                    a[index[i]][l] -= pj * a[index[j]][l];
                }
            }
        }
    }

    static public double[] fft(double[] x) {
        int n = x.length; // assume n is a power of 2
        int nu = (int) (Math.log(n) / Math.log(2));
        int n2 = n / 2;
        int nu1 = nu - 1;
        double[] xre = new double[n];
        double[] xim = new double[n];
        double[] mag = new double[n2];
        double tr, ti, p, arg, c, s;
        for (int i = 0; i < n; i++) {
            xre[i] = x[i];
            xim[i] = 0.0f;
        }
        int k = 0;

        for (int l = 1; l <= nu; l++) {
            while (k < n) {
                for (int i = 1; i <= n2; i++) {
                    p = bitrev(k >> nu1, nu);
                    arg = 2 * (double) Math.PI * p / n;
                    c = (double) Math.cos(arg);
                    s = (double) Math.sin(arg);
                    tr = xre[k + n2] * c + xim[k + n2] * s;
                    ti = xim[k + n2] * c - xre[k + n2] * s;
                    xre[k + n2] = xre[k] - tr;
                    xim[k + n2] = xim[k] - ti;
                    xre[k] += tr;
                    xim[k] += ti;
                    k++;
                }
                k += n2;
            }
            k = 0;
            nu1--;
            n2 = n2 / 2;
        }
        k = 0;
        int r;
        while (k < n) {
            r = bitrev(k, nu);
            if (r > k) {
                tr = xre[k];
                ti = xim[k];
                xre[k] = xre[r];
                xim[k] = xim[r];
                xre[r] = tr;
                xim[r] = ti;
            }
            k++;
        }

        double[] ret = new double[xre.length * 2];
        for (int i = 0; i < xre.length; i++) {
            ret[i] = xre[i];
            ret[i + xre.length] = xim[i];
            // System.out.println(": "+xre[i]+" "+xim[i]+"i");
        }
        return ret;
    }

    static public double[] abs_fft(double[] fft) {
        int half = fft.length / 2;
        double[] ret = new double[half];
        for (int i = 0; i < half; i++) {
            ret[i] = (double) Math.sqrt(fft[i] * fft[i] + fft[i + half]
                    * fft[i + half]);
        }
        return ret;
    }

    static public int bitrev(int j, int nu) {

        int j2;
        int j1 = j;
        int k = 0;
        for (int i = 1; i <= nu; i++) {
            j2 = j1 / 2;
            k = 2 * k + j1 - 2 * j2;
            j1 = j2;
        }
        return k;
    }

    static public double mean(double[] x) {
        double mean = 0;
        for (int i = 0; i < x.length; i++) {
            mean += x[i];
        }
        return mean / x.length;
    }

    public static final double[] interpLinear(double[] x, double[] y, double[] xi)
            throws Exception {

        if (x.length != y.length) {
            throw new Exception("lungime X diferita de Y");
        }

        double[] dx = new double[x.length - 1];
        double[] dy = new double[x.length - 1];
        double[] slope = new double[x.length - 1];
        double[] intercept = new double[x.length - 1];

        for (int i = 0; i < x.length - 1; i++) {
            dx[i] = x[i + 1] - x[i];
            if (dx[i] == 0) {
                throw new Exception("X contine un duplicat");
            }
            if (dx[i] < 0) {
                throw new Exception("X trebuie sa fie sortat");
            }
            dy[i] = y[i + 1] - y[i];
            slope[i] = dy[i] / dx[i];
            intercept[i] = y[i] - x[i] * slope[i];
        }

        double[] yi = new double[xi.length];
        for (int i = 0; i < xi.length; i++) {
            if ((xi[i] > x[x.length - 1]) || (xi[i] < x[0])) {
                yi[i] = Double.NaN;
            } else {
                int loc = Arrays.binarySearch(x, xi[i]);
                if (loc < -1) {
                    loc = -loc - 2;
                    yi[i] = slope[loc] * xi[i] + intercept[loc];
                } else {
                    yi[i] = y[loc];
                }
            }
        }

        return yi;
    }

    public static double[] hanning(int N) {
        double[] data = new double[N];
        for (int i = 0; i < N; i++) {
            double n = i;
            //hamming
            //double p1 = (double)(0.54f - 0.46f * Math.cos(2.0 * Math.PI * n / (N - 1)));
            //hanning
            double p1 = (double) (0.5 * (1.0 - Math.cos(2.0 * Math.PI * (n + 1) / (N + 1))));
            //double p1 = (double)(2595 * Math.Log(i * 700 + 1));
            data[i] = p1;
        }
        return data;
    }

    public static double[] fftfilt(double[] H, double[] signal) {
        double[] result = new double[signal.length];

        for (int i = 0; i < signal.length; i++) {
            double val = 0;
            for (int k = 0; k < H.length; k++) {
                if (i - k > 0) {
                    if (i + k - H.length >= 0) {
                        val += (double) signal[i - k] * H[k];
                    }
                }
            }
            result[i] = (double) val;

        }
        return result;
    }

    public static double[][] smax(double[][] x, double a, double b) {
        double[][] y = new double[x.length][x[0].length];
        double y0 = 1.0f / (1.0f + (double) Math.exp(-a * (0.0 - b)));
        double y1 = 1.0f / (1.0f + (double) Math.exp(-a * (1.0 - b)));
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < x[0].length; j++) {
                y[i][j] = (double) (1.0 / (1.0 + Math.exp(-a * (x[i][j] - b))) - y0) / (y1 - y0);
            }
        }
        return y;
    }

    //%	Phase rotator for fractional pitch
    public static double[] fractpitch2(int fftl) {
        int amp = 15;
        double[] t = new double[fftl];

        //t=([1:fftl]-fftl/2-1)/fftl*2;
        for (int i = 0; i < t.length; i++) {
            t[i] = (double) (((double) i - (double) fftl / 2.0f - 0.0) / (double) fftl * 2.0f);
        }

        //phs=t+(1-exp(amp*t))./(1+exp(amp*t)) -(1+(1-exp(amp))/(1+exp(amp)))*t;\
        double[] phs = new double[t.length];
        double exp_amp, exp_amp_t;

        for (int i = 0; i < t.length; i++) {
            exp_amp = (double) Math.exp(amp);
            exp_amp_t = (double) Math.exp(amp * t[i]);
            phs[i] = t[i] + (1.0f - exp_amp_t) / (1.0f + exp_amp_t) - (1.0f + (1.0f - exp_amp) / (1 + exp_amp)) * t[i];
            phs[i] *= Math.PI;
        }
        //phs(1)=0;
        phs[0] = 0;

        //phs=phs*pi;
        return phs;
    }

    public static double[] fftshift_array(double[] x) {
        int half = x.length / 2;
        double[] y = new double[x.length];
        System.arraycopy(x, half, y, 0, half);
        System.arraycopy(x, 0, y, half, half);
        //for (int i = 0; i < half; i++) {
        //    y[i] = x[half + i];
        //    y[half + i] = x[i];
        //}
        return y;
    }

    public static double[] cumsum(double[] x) {
        double[] y = new double[x.length];
        double sum;
        for (int i = 0; i < x.length; i++) {
            sum = 0;
            for (int j = 0; j <= i; j++) {
                sum += x[j];
            }
            y[i] = sum;
        }
        return y;
    }

    public static void multiplyComplex(double[] a, double[] b, double[] c, double[] d, double[] rezx, double[] rezy) {
        for (int i = 0; i < a.length; i++) {
            rezx[i] = a[i] * c[i] - b[i] * d[i];
            rezy[i] = a[i] * d[i] + b[i] * c[i];
        }
    }

    public static void exponentComplex(double[] b, double[] c, double[] rezx, double[] rezy) {
        for (int i = 0; i < b.length; i++) {
            double e = Math.pow(Math.E, b[i]);
            rezx[i] = e * Math.cos(c[i]);
            rezy[i] = e * Math.sin(c[i]);
        }
    }

    public static void log(String s) {
        if (true) {
            return;
        }
        try {
            BufferedWriter bw = IO.openFileForWriting("D:\\Dropbox\\STRAIGHTsrc\\log.txt", true);
            bw.write(s + "\n\r");
            bw.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static void compare(double[] x, String nume) throws Exception {
        if (true) {
            return;
        }
        ArrayList<String> d = new ArrayList<String>();
        String[] data = null;
        try {
            data = IO.readAllLines("D:\\Dropbox\\STRAIGHTsrc\\temp-" + nume);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            //System.err.println("\tCompare skipping: "+nume);
            return;
            //e.printStackTrace();
        }
        int i = 0;
        System.out.println("Variabila " + nume + " are " + data.length + " randuri in Matlab si " + x.length + " la mine.");
        if (data.length != x.length) //throw new Exception ("Variabila " + nume + " are "+data.length+" randuri in Matlab si "+x.length+" la mine.");
        {
            return;
        }
        //System.err.println("Variabila " + nume + " are "+data.length+" randuri in Matlab si "+x.length+" la mine.");
        int dif = 0;
        for (String n : data) {
            //d.add(double.valueOf(n));
            String wrt = i + "\t" + n + "\t" + String.format("%f", x[i]);
            double diff = Math.abs(Double.parseDouble(n) - x[i]);
            if (diff > 0.00001) {
                wrt += "\t DIFF!";
                if (nume == "sy_temp" && dif == 0) {
                    System.err.println("DIFF sytemp");
                    dif = 1;
                }
            }
            d.add(wrt);
            i++;
        }
        IO.writeAllLines(d.toArray(new String[0]), "D:\\Dropbox\\STRAIGHTsrc\\temp-" + nume + ".txt");
    }

    public static double zpowerchk(double[] x, double fs, double segms) {
        /*%	Calculate average power of voiced portion
         %	pow=powerchk(x,fs,segms)
         %		x	: signal
         %		fs	: sampling frequency (Hz)
         %		segms	: segment length (ms)
         */
        double pow;
        /*
         x1=x(:);
         iv=(1:length(x1))';
         x1(isnan(x1))=iv(isnan(x1))*0+0.0000000001;
         x2=x1.*x1;
         %n=100; % 23/Sept./1999
         n=round(segms/1000*fs); % 23/Sept./1999
         nw=ceil(length(x)/n);
         if rem(length(x),n)>0
         %  x2=[x2;zeros(n*nw-length(x),1)];
         x2=[x2;0.000001*randn(n*nw-length(x),1).^2]; % 23/Sept./1999
         end;
         x2(x2==0)=x2(x2==0)+0.000001;

         pw=sum(reshape(x2,n,nw))/n;

         pow=10*log10(mean(pw(pw>(mean(pw)/30))));
         */
        pow = 0;
        return pow;
    }
}
