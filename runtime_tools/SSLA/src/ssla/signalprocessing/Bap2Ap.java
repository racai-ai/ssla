/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ssla.signalprocessing;

/**
 *
 * @author tibi
 */
public class Bap2Ap {

    public double[][] bap2ap(double[][] bap) {
        bufp = 0;
        double[] a = new double[]{1};
        double[] b = new double[]{-1, 1};
        double[][] ap = new double[bap.length][1025];
        double[] stream1 = new double[bap.length * 2];
        int spos = 0;
        getAp(bap, 0, 2, ap, 0);
        getAp(bap, 1, 3, ap, 2);
        getAp(bap, 2, 4, ap, 5);
        getAp(bap, 3, 4, ap, 9);
        getAp(bap, 4, 4, ap, 13);
        getAp(bap, 5, 5, ap, 17);
        getAp(bap, 6, 5, ap, 22);
        getAp(bap, 7, 6, ap, 27);
        getAp(bap, 8, 6, ap, 33);
        getAp(bap, 9, 7, ap, 39);
        getAp(bap, 10, 8, ap, 46);
        getAp(bap, 11, 9, ap, 54);
        getAp(bap, 12, 10, ap, 63);
        getAp(bap, 13, 12, ap, 73);
        getAp(bap, 14, 14, ap, 85);
        getAp(bap, 15, 16, ap, 99);
        getAp(bap, 16, 19, ap, 115);
        getAp(bap, 17, 24, ap, 134);
        getAp(bap, 18, 29, ap, 158);
        getAp(bap, 19, 37, ap, 187);
        getAp(bap, 20, 49, ap, 224);
        getAp(bap, 21, 68, ap, 273);
        getAp(bap, 22, 99, ap, 341);
        getAp(bap, 23, 160, ap, 440);
        getAp(bap, 24, 300, ap, 600);
        getAp(bap, 25, 125, ap, 900);
//        for (int i = 0; i < ap.length; i++) {
//            for (int k = 1; k < ap[i].length - 1; k++) {
//                ap[i][k] = (ap[i][k - 1] + ap[i][k + 1]) / 2;
//            }
//        }
        return ap;
    }
    int bufp = 0;

    double dfs(double x, double[] a, double[] b, double[] buf) {
        double y = 0.0;
        int i, p;
        int max;

        int m = a.length;
        int n = b.length;

        if (m < n) {
            max = n;
        } else {
            max = m;
        }

        x = x * a[0];
        for (i = 1; i < m; i++) {
            p = bufp + i;
            if (p >= max) {
                p -= max;
            }
            x -= buf[p] * a[i];
        }
        buf[bufp] = x;
        for (i = 0; i < n; i++) {
            p = bufp + i;
            if (p >= max) {
                p -= max;
            }
            y += buf[p] * b[i];
        }

        if ((--bufp) < 0) {
            bufp += max;
        }

        return (y);
    }
    double[] buf1 = new double[3];
    double[] buf2 = new double[3];

    private void getAp(double[][] bap, int bpos, int p, double[][] ap, int aps) {
        int bufp1 = 0;
        int bufp2 = 0;
        int apl = p;
        double[] a = new double[]{1};
        double[] b = new double[]{1, -1};

        for (int i = 0; i < bap.length; i++) {
            //dfs
            bufp = bufp1;
            double x = dfs(bap[i][bpos], a, b, buf1);
            bufp1 = bufp;
            //interpolate
            //aici teoretic bufp=0
            bufp = bufp2;
            double y = dfs(x, b, a, buf2);
            ap[i][aps] = y;

            for (int k = 0; k < apl - 1; k++) {
                y = dfs(0, b, a, buf2);
                ap[i][aps + k + 1] = y;
            }
            bufp2 = bufp;
        }
    }
}
