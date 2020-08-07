package ro.racai.tts.WS.TTS.straight;

import ro.racai.tts.WS.TTS.signalprocessing.FFT;

public class straightSynth {

    /*
     * % sy : synthsized speech % n2sgram : amplitude spectrogram % f0raw :
     * pitch pattern (Hz) % f0var : expected F0 variation with fricative
     * modification % f0varL : expected F0 variation % shiftm : frame shift (ms)
     * for spectrogram % fs : sampling freqnency (Hz) % pcnv : pitch stretch
     * factor % fconv : freqnency stretch factor % sconv : speaking duratin
     * stretch factor (overridden if || imap || >1 ) % gdbw : finest resolution
     * in group delay (Hz) % delfrac : ratio of standard deviation of group
     * delay in terms of F0 % delsp : standard deviation of group delay (ms) %
     * cornf : lower corner frequency for phase randomization (Hz) % delfracind
     * : selector of fixed and proportional group delay % ap : aperiodicity
     * measure % imap : arbirtary mapping from new time (sample) to old time
     * (frame) % imgi : display indicator, 1: display on (default), 0: off %
     * lowestF0 : lower limit of the resynthesized fundamental frequency (Hz)
     */
    /*
     VALORI DEFAULT, tot ce trebuie este n2sgram, f0raw si ap
     double shiftm = 1;
     double fs = 48000;
     double pconv =1;
     double fconv =1;
     double sconv =1;
     double gdbw = 70;
     double delfrac = 0.2000f;
     double delsp = 0.5000f;
     double cornf = 4000f;
     double delfracind = 0;
     double imap = 1;
     double lowestF0 = 50;
     */
    public static boolean debug = false;

    public static double[] straightSynth(double[][] n2sgram,
            double[] f0raw, double shiftm, double fs, double pcnv,
            double fconv, double sconv, double gdbw, double delfrac,
            double delsp, double cornf, double delfracind, double[][] ap,
            double imap_scalar, double lowestF0) throws Exception {

        //de proba:P
        double[][] newap = new double[513][f0raw.length];
        double[][] newsp = new double[513][f0raw.length];
        for (int i = 0; i < newap.length; i++) {
            for (int j = 0; j < newap[0].length; j++) {
                newap[i][j] = ap[i * 2][j];
                newsp[i][j] = n2sgram[i * 2][j];
            }
        }
        n2sgram = newsp;
        ap = newap;
        long start_time = System.currentTimeMillis();
        
        int ii = 0;

        int nii = n2sgram.length;
        int njj = n2sgram[0].length < f0raw.length ? n2sgram[0].length : f0raw.length;

        double[] f0l = new double[njj];
        System.arraycopy(f0raw, 0, f0l, 0, njj);

        if (debug) {
            for (int i = 0; i < f0l.length; i++) {
                if (f0l[i] > 0 && f0l[i] * pcnv < lowestF0) {
                    throw new Exception("Minimum synthesized F0 exceeded the lower limit " + lowestF0 + " Hz");
                }
            }
        }

        double fftLengthForLowestF0 = (double) Math.pow(2.0, Math.ceil(Math.log(2.0 * Math.round(fs / lowestF0))));

        int fftl = 2 * (nii - 1);
        double[] fft_real = new double[fftl];        
        double[] fft_imag = new double[fftl];
        

        if (debug) {
            if (fftl < fftLengthForLowestF0) {
                throw new Exception("The FFT length was inconsistent and replaced. XXX de refacut interpolare aici daca este nevoie si da eroarea asta, suna-l pe stefan.");
            }
        }

        double fftl2 = fftl / 2.0;
        double[] lft = new double[fftl];
        double[] hanning_fftl = M.hanning(fftl);

        double[][] aprms = new double[ap.length][ap[0].length];
        for (int i = 0; i < ap.length; i++) {
            for (int j = 0; j < ap[0].length; j++) {
                aprms[i][j] = (double) Math.pow(10, ap[i][j] / 20.0);
            }
        }

        double f_temp;
        double[][] aprm = new double[ap.length][ap[0].length];
        for (int i = 0; i < aprms.length; i++) {
            for (int j = 0; j < aprms[0].length; j++) {
                f_temp = aprms[i][j] * 1.6 - 0.015;
                if (f_temp < 0.001) {
                    f_temp = 0.001;
                }
                if (f_temp > 1.000) {
                    f_temp = 1;
                }
                aprm[i][j] = f_temp;
            }
        }

        double[] idcv = new double[fftl / 2 + 1];// e de la mine +1
        for (int i = 0; i < idcv.length; i++) {
            idcv[i] = Math.min((double) i / fconv + 1.0, fftl / 2.0 + 1.0);
        } // posibil sa trebuiasca cu minimul fara +1.0f !! ca depaseste indexul

        double[][] sy = new double[(int) ((njj * shiftm / 1000 * fs) * sconv + 3 * fftl + 1)][1];
        double[] imap = new double[sy.length];
        for (int i = 0; i < imap.length; i++) {
            imap[i] = Math.min(f0l.length, ((double) i / fs * 1000.0 / shiftm / sconv + 1.0));
        }

        // imap=[imap ones(1,round(fs*0.2))*length(f0l)]; % safe guard
        double[] imap_temp = new double[(int) (imap.length + (int) fs * 0.2)];
        for (int i = 0; i < imap_temp.length; i++) {
            imap_temp[i] = f0l.length;
        }
        System.arraycopy(imap, 0, imap_temp, 0, imap.length);
        //imap = new double[imap_temp.length];
        imap = imap_temp.clone();

        // ix=min(find(imap>=length(f0l)));
        int ix = -1;
        for (int i = 0; i < imap.length; i++) {
            if (imap[i] >= f0l.length) {
                ix = i;
                break;
            }
        }

        // rmap=interp1(imap(1:ix),[1:ix],[1:length(f0l)]);
        imap_temp = new double[ix];
        System.arraycopy(imap, 0, imap_temp, 0, imap_temp.length);
        double[] temp_array = new double[ix];
        for (int i = 0; i < temp_array.length; i++) {
            temp_array[i] = i + 1;
        }
        double[] temp_array2 = new double[f0l.length];
        for (int i = 0; i < temp_array2.length; i++) {
            temp_array2[i] = i + 1;
        }
        double[] rmap;//=new double[f0l.length];
        rmap = M.interpLinear(imap_temp, temp_array, temp_array2);
        // hack tampit ca ultima valoare da NaN si nu stiu de ce:
        rmap[rmap.length - 1] = ix + 1;

        double[] f = new double[fftl / 2];
        for (int i = 0; i < f.length; i++) {
            f[i] = (double) i / fftl * fs;
        }

        double[] lf = new double[f.length];
        for (int i = 0; i < lf.length; i++) {
            lf[i] = (double) Math.log10(f[i] + 0.1f);
        }

        double[] phs = M.fractpitch2(fftl);

        double[] a = new double[(int) (fftl2 * 2.0)];
        for (int i = 0; i < fftl2; i++) {
            a[i] = (double) i / fftl2;
        }
        for (int i = (int) fftl2 + 1; i < a.length; i++) {
            a[i] = (double) (-(fftl2 - (i - fftl2)) / fftl2);
        }

        fftl2 = fftl / 2;

        int nsyn = sy.length;

        double idx = 1;

        int[] bb = new int[fftl];
        for (int i = 0; i < bb.length; i++) {
            bb[i] = i;
        }

        double[] rbb2 = new double[fftl / 2 - 1];
        for (int i = 0; i < rbb2.length; i++) {
            rbb2[i] = fftl / 2.0 - i;
        }

        double[] fxa = new double[(int) fftl2];
        for (int i = 0; i < fxa.length; i++) {
            fxa[i] = i / fftl * fs;
        }

        double lowcutf = 0;
        double counter_mean = 0;
        for (int i = 0; i < f0l.length; i++) {
            if (f0l[i] > 0) {
                lowcutf += f0l[i];
                counter_mean++;
            }
        }
        lowcutf = lowcutf / counter_mean * 0.7 * pcnv;

        //double lowcutfav = lowcutf;

        //double[] wlcutav = new double[fxa.length];
        //for (int i = 0; i < fxa.length; i++) {
        //    wlcutav[i] = 1.0 / (1 + Math.exp(-14 * (fxa[i] - lowcutfav) / lowcutfav));
        //}
        //M.compare(wlcutav, "wlcutav");
        //double[] waprm = new double[fxa.length];
        //for (int i = 0; i < fxa.length; i++) {
        //    waprm[i] = 1.0 / (1 + Math.exp(-14 * (fxa[i] - 1000) / 200));
        //}
        //M.compare(waprm, "waprm");
        double[] t = new double[fftl];
        for (int i = 0; i < t.length; i++) {
            t[i] = (i - fftl / 2.0) / fftl * 2.0;
        }

        double[] adjd = new double[t.length];
        for (int i = 0; i < adjd.length; i++) {
            adjd[i] = 1.0 / (1.0 + Math.exp(-20 * t[i]));
        }

        double[] gw = new double[t.length];
        double sum_gw = 0;
        for (int i = 0; i < gw.length; i++) {
            gw[i] = Math.exp(-0.25 * Math.PI
                    * Math.pow(fs * (t[i] / 2.0) / gdbw, 2));
            sum_gw += gw[i];
        }

        for (int i = 0; i < gw.length; i++) {
            gw[i] /= sum_gw;
        }

        double[] fgw;//= new double[gw.length];
        fgw = M.fftshift_array(gw);
        FFT fft_object = FFT.get(gw.length);
        //double[] fft_real = new double[gw.length];
        //System.arraycopy(fgw, 0, fgw, 0, fgw.length);
        //double[] fft_imag = new double[gw.length];
        fft_object.fft(fgw, fft_imag);
        //System.arraycopy(fft_real, 0, fgw, 0, fgw.length);

        // df=fs/fftl*2*pi; % normalization constant for integration and
        // differentiation
        double df = fs / fftl * 2.0 * (double) Math.PI;

        // fw=(1:fftl2+1)/fftl*fs; % frequency axis
        double[] fw = new double[(int) (fftl2 + 1)];
        for (int i = 0; i < fw.length; i++) {
            fw[i] = (i + 1.0) / fftl * fs;
        }

        // trbw=300; % width of transition area
        double trbw = 300;

        // rho=1.0./(1+exp(-(fw-cornf)/trbw)); % rondom group delay weighting
        // function
        double[] rho = new double[fw.length];
        for (int i = 0; i < rho.length; i++) {
            rho[i] = 1.0 / (1.0 + (double) Math.exp(-(fw[i] - cornf) / trbw));
        }

        double[] nz = new double[(int) (fftl2 + 1)];
        nz = M.randn_1row_ncolumns(nz.length);
        for (int i = 0; i < nz.length; i++) {
            nz[i] *= rho[i];
        }

        double[] nz_temp = new double[nz.length + rbb2.length];
        System.arraycopy(nz, 0, nz_temp, 0, nz.length);
        for (int i = 0; i < rbb2.length; i++) {
            nz_temp[nz.length + i] = nz[(int) rbb2[i] - 1];// pizda masii.
        }

        nz = new double[nz_temp.length];
        System.arraycopy(nz_temp, 0, nz, 0, nz_temp.length);
        fft_object = FFT.get(nz_temp.length);
        fft_real = new double[nz_temp.length];
        System.arraycopy(nz_temp, 0, fft_real, 0, nz_temp.length);
        fft_imag = new double[nz_temp.length];
        fft_object.fft(fft_real, fft_imag);
        for (int i = 0; i < fgw.length; i++) {
            fft_real[i] *= fgw[i];
            fft_imag[i] *= fgw[i];
        }
        fft_object.ifft(fft_real, fft_imag);
        System.arraycopy(fft_real, 0, nz, 0, nz.length);

        double temp = (double) Math.sqrt(fftl * gdbw / fs);
        for (int i = 0; i < nz.length; i++) {
            nz[i] *= temp * delsp * df / 1000.0;
        }

        nz_temp = new double[(int) (fftl2 + 1 + rbb2.length)];
        System.arraycopy(nz, 0, nz_temp, 0, (int) fftl2 + 1);
        for (int i = 0; i < rbb2.length; i++) {
            nz_temp[(int) (fftl2 + 1 + i)] = nz[(int) rbb2[i] - 1];
        }

        double[] mz;//= new double[nz_temp.length];
        mz = M.cumsum(nz_temp);
        for (int i = 0; i < mz.length; i++) {
            mz[i] -= nz[0];
        }

        double[] mmz = new double[mz.length];
        for (int i = 0; i < mz.length; i++) {
            mmz[i] = (double) -(mz[i] - adjd[i]
                    * (((mz[fftl - 1] + mz[1]) % (2.0 * (double) Math.PI)) - 2.0 * Math.PI));

        }

        // pz=exp(-i*mmz)'; %.*[wlcut wlcut(rbb2)]';
        // nefolosit!!
        // %---------------------------------------------------------
        // [snn,smm]=size(n2sgram);
        int snn = n2sgram.length;
        int smm = n2sgram[0].length;

        // fqx=(0:snn-1)/snn*fs/2;
        double[] fqx = new double[snn];
        for (int i = 0; i < fqx.length; i++) {
            fqx[i] = (double) ((double) i / (double) snn * fs / 2.0);
        }

        // chigh=1.0./(1+exp(-(fqx-600)/100))';
        double[] chigh = new double[fqx.length];
        for (int i = 0; i < chigh.length; i++) {
            chigh[i] = (double) (1.0 / (1.0 + Math
                    .exp(-(fqx[i] - 600.0) / 100.0)));
        }

        for (int i = 0; i < fftl; i++) {
            lft[i] = (double) (1.0 / (1.0 + Math
                    .exp(-((1.0 - hanning_fftl[i]) - 0.5) * 60.0)));
        }

        double[] ww = new double[fftl];
        for (int i = 0; i < fftl; i++) {
            ww[i] = (double) (1.0 / (1.0 + Math
                    .exp(-(hanning_fftl[i] - 0.3) * 23.0)));
        }

        int iin = 1;

        int icntr = 0;

        double dmx = Double.MIN_VALUE;
        for (int i = 0; i < snn; i++) {
            for (int j = 0; j < smm; j++) {
                if (n2sgram[i][j] > dmx) {
                    dmx = n2sgram[i][j];
                }
            }
        }
        double dmx_1000000 = dmx/1000000.0;
        double dmx_100000 = dmx/100000.0;
        
        double[] gh = null;
        
        
        double[] ff = new double[idcv.length + rbb2.length];
        double[] ccp = new double[ff.length];
        double[] ccp2 = new double[fftl];
        double[] ffx = new double[ccp2.length];
        double[] ffx_imag = new double[ccp2.length];
        double[] i_real = new double[phs.length];
        double[] i_imag = new double[phs.length];
        double[] frtz = new double[phs.length];
        double[] frtz_imag = new double[phs.length];
        nz_temp = new double[(int) (fftl2 + 1 + rbb2.length)];            
        mmz = new double[mz.length];            
        double[] pz = new double[mmz.length];
        double[] pz_imag = new double[mmz.length];
        i_real = new double[mmz.length];
        i_imag = new double[mmz.length];
        double[] ep = new double[ffx.length];
        double[] epf = new double[ep.length];
        double[] epf_imag = new double[ep.length];
        
        //double[] tx_imag = new double[frtz.length];
        double[] temp_real = new double[frtz.length];
        double[] temp_imag = new double[frtz.length];
        double[] temp2_real = new double[frtz.length];
        double[] temp2_imag = new double[frtz.length];

        


        while ((idx < nsyn - fftl - 10) && (iin < f0l.length)) {
            //M.log("------------------");
            icntr++;
            int iix = (int) Math.round(imap[(int) (Math.round(idx) - 1)]);
            ii = Math.min(Math.min(Math.max(1, iix), njj), f0l.length);
            //M.log(" ii : " + ii + " iix " + iix + " idx " + idx);
            double f0 = f0l[Math.round(ii) - 1];
            //M.log(ii + "\t f0 este " + f0);
            if (f0 == 0) {
                f0 = 200;
            } else {
                f0 = Math.max(lowestF0 / pcnv, f0l[Math.round(ii) - 1]);
            }
            //M.log(ii + "\t f0 devine " + f0);

            f0 = f0 * pcnv;

            double tnf0 = fs / f0;

            double tidx = idx + tnf0;

            int tiix = (int) Math.round(imap[(int) Math.round(tidx) - 1]);

            int tii = Math.min(Math.min(Math.max(1, tiix), njj), f0l.length);

            double tf0 = f0l[tii - 1];
            //M.log("tidx=" + tidx + " tiix=" + tiix + " tii=" + tii + " tf0=" + tf0 + " f0=" + f0);

            if (tf0 > 0 && f0l[ii - 1] > 0) {
                //M.log("in if 2. f0l = " + f0l[ii - 1]);
                if (f0l[(int) Math.round((ii + tii) / 2.0) - 1] > 0.0) {
                    f0 = Math.max(lowestF0 / pcnv,
                            f0l[(int) Math.round((ii + tii) / 2.0) - 1]);
                    //M.log("  caz 1. f0 = " + f0 + " (tii = " + tii + ") index = " + ((int) Math.round((ii + tii) / 2.0)));
                } else {
                    f0 = f0l[ii - 1];
                    //M.log("  caz 2. f0 = " + f0);
                }
                f0 = f0 * pcnv;
            }
            //M.log(ii + "\t f0 NOU " + f0);

            
            for (int i = 0; i < idcv.length; i++) {
                ff[i] = n2sgram[(int) idcv[i] - 1][ii - 1];
            }
            for (int i = idcv.length; i < idcv.length + rbb2.length; i++) {
                ff[i] = n2sgram[(int) idcv[(int) rbb2[i - idcv.length] - 1] - 1][ii - 1];
            }
            //M.compare(ff, "ff" + ii);

            
            for (int i = 0; i < ccp.length; i++) {
                ccp[i] = Math.log(ff[i] + dmx_1000000);
            }
            //double[] delme_ccp = ccp.clone();
            //M.compare(delme_ccp, "delme_ccp" + ii);

            fft_object = FFT.get(ccp.length);
            //fft_real = new double[ccp.length];
            System.arraycopy(ccp, 0, fft_real, 0, ccp.length);
            fft_imag = new double[ccp.length];
            fft_object.fft(ccp, fft_imag);
            //System.arraycopy(fft_real, 0, ccp, 0, ccp.length);
            //M.compare(ccp, "ccp" + ii);

            ccp2[0] = ccp[0];
            for (int i = 1; i < fftl / 2; i++) {
                ccp2[i] = ccp[i] * 2.0;
            }

            //M.compare(ccp2, "ccp2" + ii);
            for (int i = 0; i < ffx.length; i++) {
                ffx[i] = ccp2[i] * lft[i] / fftl;
            }
            fft_object = FFT.get(ffx.length);
            //fft_real = new double[ffx.length];
            //System.arraycopy(ffx, 0, fft_real, 0, ffx.length);
            //fft_imag = new double[ffx.length];
            ffx_imag = new double[ffx.length];
            fft_object.fft(ffx, ffx_imag);
            //System.arraycopy(fft_real, 0, ffx, 0, ffx.length);
            //System.arraycopy(fft_imag, 0, ffx_imag, 0, ffx.length);

            //M.compare(ffx, "ffx" + ii);
            //M.compare(idcv, "idcv" + ii);
            //M.compare(ffx_imag, "ffx_imag" + ii);
            double nidx = idx;

            double nf0 = fs / f0;

            double frt = idx - nidx;

            
            for (int i = 0; i < phs.length; i++) {
                i_real[i] = 0;
                i_imag[i] = phs[i] * frt;
            }
            M.exponentComplex(i_real, i_imag, frtz, frtz_imag);

            //M.compare(frtz, "frtz" + ii);
            nz = new double[(int) (fftl2 + 1)];
            nz = M.randn_1row_ncolumns(nz.length);
            for (int i = 0; i < nz.length; i++) {
                nz[i] *= rho[i];
            }

            //nz_temp = new double[nz.length + rbb2.length];
            System.arraycopy(nz, 0, nz_temp, 0, nz.length);
            for (int i = 0; i < rbb2.length; i++) {
                nz_temp[nz.length + i] = nz[(int) rbb2[i] - 1];// pizda masii 2.
            }
            System.arraycopy(nz, 0, nz_temp, 0, nz.length);
            
            fft_object = FFT.get(nz_temp.length);
            fft_real = new double[nz_temp.length];
            System.arraycopy(nz_temp, 0, fft_real, 0, nz_temp.length);
            fft_imag = new double[nz_temp.length];
            fft_object.fft(fft_real, fft_imag);
            for (int i = 0; i < fgw.length; i++) {
                fft_real[i] *= fgw[i];
                fft_imag[i] *= fgw[i];
            }
            fft_object.ifft(fft_real, fft_imag);
            System.arraycopy(fft_real, 0, nz, 0, nz.length);

            temp = (double) Math.sqrt(fftl * gdbw / fs);
            for (int i = 0; i < nz.length; i++) {
                nz[i] *= temp * delsp * df / 1000.0;
            }

            //M.compare(nz, "nz" + ii);
            // mz=cumsum([nz(1:fftl2+1),nz(rbb2)])-nz(1);
            System.arraycopy(nz, 0, nz_temp, 0, (int) fftl2 + 1);
            for (int i = 0; i < rbb2.length; i++) {
                nz_temp[(int) (fftl2 + 1 + i)] = nz[(int) rbb2[i] - 1];
            }

            //mz = new double[nz_temp.length];
            mz = M.cumsum(nz_temp);
            for (int i = 0; i < mz.length; i++) {
                mz[i] -= nz[0];
            }

            //M.compare(mz, "mz" + ii);
           for (int i = 0; i < mz.length; i++) {
                mmz[i] = (double) -(mz[i] - adjd[i]
                        * (((mz[fftl - 1] + mz[1]) % (2.0 * (double) Math.PI)) - 2.0 * Math.PI));

            }

            //M.compare(mmz, "mmz" + ii);
            for (int i = 0; i < pz_imag.length; i++) {
                i_real[i] = 0;
                i_imag[i] = mmz[i];
            }
            M.exponentComplex(i_real, i_imag, pz, pz_imag);

            //M.compare(pz, "pz" + ii);
            //M.compare(pz_imag, "pz_imag" + ii);
            double[] wnz = new double[aprm.length];
            for (int i = 0; i < idcv.length; i++) {
                wnz[i] = aprm[(int) idcv[i] - 1][ii - 1];
            }

            //M.compare(wnz, "wnz" + ii);
            double[] wpr = new double[wnz.length];
            for (int i = 0; i < wnz.length; i++) {
                wpr[i] = Math.sqrt(Math.max(0, 1.0 - wnz[i] * wnz[i]));
            }
            //M.compare(wpr, "wpr" + ii);

            double zt0 = nf0 / fs;

            double ztc = 0.01;

            double[] ztp = new double[(int) nf0];
            for (int i = 0; i < ztp.length; i++) {
                ztp[i] = i / fs;
            }

            double[] nev = new double[ztp.length];
            for (int i = 0; i < nev.length; i++) {
                nev[i] = Math.sqrt(2.0 * zt0 / ztc
                        / (1.0 - Math.exp(-2.0 * zt0 / ztc)))
                        * Math.exp(-ztp[i] / ztc);
            }

            //M.compare(nev, "nev" + ii);
            double[] rx = new double[(int) Math.round(nf0)];
            rx = M.randn_1row_ncolumns(rx.length);

            //M.compare(rx, "rx" + (int) nf0);
            double[] wfv = new double[nev.length];
            double[] wfv_imag = new double[fftl];

            double mean_rx = 0;
            for (int i = 0; i < rx.length; i++) {
                mean_rx += rx[i];
            }
            mean_rx /= rx.length;

            for (int i = 0; i < wfv.length; i++) {
                wfv[i] = (rx[i] - mean_rx) * nev[i];
            }

            fft_object = FFT.get(fftl);
            fft_real = new double[fftl];
//            if (fft_real.length < wfv.length) {
//                System.out.print(".");
//            }
            //TIBI: aici am schimbat eu (trebuie scos daca dispare reducerea rezolutiei pe spectru de la inceputul metdei 
            System.arraycopy(wfv, 0, fft_real, 0, wfv.length < fft_real.length ? wfv.length : fft_real.length);
            fft_imag = new double[fftl];
            fft_object.fft(fft_real, fft_imag);
            wfv = new double[fftl];
            System.arraycopy(fft_real, 0, wfv, 0, fft_real.length);
            System.arraycopy(fft_imag, 0, wfv_imag, 0, fft_imag.length);

            //M.compare(wfv, "wfv" + ii);
            
            int nf0n = (int) nf0;
            //TIBI: idem cu ce am scris mai sus
            if (nf0n > fftl) {
                nf0n = fftl;
            }

            //M.compare(gh, "gh" + ii);
            if (gh == null || gh.length != nf0n * 2) {
                gh = M.hanning(nf0n * 2);
            }
            for (int i = 0; i < nf0n; i++) {
                ep[i] = gh[nf0n - 1 - i];
            }

            //M.compare(ep, "ep" + ii);
            int temp_int = 0;
            for (int i = ep.length - 1; i > ep.length - nf0n; i--) {
                ep[i] = ep[++temp_int];
            }

            //M.compare(ep, "ep" + ii);
            double ep_sum = 0;
            for (int i = 0; i < ep.length; i++) {
                ep_sum += ep[i];
            }
            for (int i = 0; i < ep.length; i++) {
                ep[i] = -ep[i] / ep_sum;
            }

            ep[0]++;

            //M.compare(ep, "ep" + ii);
            fft_object = FFT.get(ep.length);
            fft_real = new double[ep.length];
            System.arraycopy(ep, 0, fft_real, 0, ep.length);
            fft_imag = new double[ep.length];
            fft_object.fft(fft_real, fft_imag);
            System.arraycopy(fft_real, 0, epf, 0, fft_real.length);
            System.arraycopy(fft_imag, 0, epf_imag, 0, fft_imag.length);

            //M.compare(epf, "epf" + ii);
            double[] wpr_temp = new double[wpr.length + rbb2.length];
            System.arraycopy(wpr, 0, wpr_temp, 0, wpr.length);
            for (int i = 0; i < rbb2.length; i++) {
                wpr_temp[wpr.length + i] = wpr[(int) rbb2[i] - 1];// pizda masii
                // #3.
            }
            //M.compare(wpr_temp, "wpr_temp");

            double[] exp_ffx_real = new double[frtz.length];
            double[] exp_ffx_imag = new double[frtz.length];
            for (int i = 0; i < exp_ffx_real.length; i++) {
                exp_ffx_real[i] = Math.exp(ffx[i]) * Math.cos(ffx_imag[i]);
                exp_ffx_imag[i] = Math.exp(ffx[i]) * Math.sin(ffx_imag[i]);
            }

            
            M.multiplyComplex(epf, epf_imag, exp_ffx_real, exp_ffx_imag,
                    temp2_real, temp2_imag);

            M.multiplyComplex(temp2_real, temp2_imag, pz, pz_imag,
                    temp_real, temp_imag);
            
            double[] tx = new double[wpr_temp.length];
            for (int i = 0; i < tx.length; i++) {
                temp_real[i] *= frtz[i] * wpr_temp[i];
                temp_imag[i] *= frtz[i] * wpr_temp[i];
            }
            //double[] btar_real; //= new double[temp_real.length];
            //double[] btar_imag;// = new double[temp_real.length];
            //btar_real = temp_real.clone();
            //btar_imag = temp_imag.clone();

            //M.compare(btar_real, "btar_real" + ii);
            //M.compare(btar_imag, "btar_imag" + ii);
            fft_object = FFT.get(tx.length);
            //fft_real = new double[tx.length];
            //System.arraycopy(temp_real, 0, fft_real, 0, tx.length);
            //fft_imag = new double[tx.length];
            //System.arraycopy(temp_imag, 0, fft_imag, 0, tx.length);
            fft_object.ifft(temp_real, temp_imag);
            //temp_real = fft_real.clone();
            //temp_imag = fft_imag.clone();
            //M.compare(temp_real, "temp_real" + ii);
            //M.compare(temp_imag, "temp_imag" + ii);

            tx = M.fftshift_array(temp_real);

            for (int i = 0; i < tx.length; i++) {
                tx[i] *= ww[i];
            }

            //M.compare(tx, "tx-" + ii);
            double[] wnz_temp = new double[wnz.length + rbb2.length];
            System.arraycopy(wnz, 0, wnz_temp, 0, wnz.length);
            for (int i = 0; i < rbb2.length; i++) {
                wnz_temp[wnz.length + i] = wnz[(int) rbb2[i] - 1];// pizda masii
            }

            M.multiplyComplex(wfv, wfv_imag, exp_ffx_real, exp_ffx_imag,
                    temp_real, temp_imag);
            double[] tx2 = new double[tx.length];
            double[] tx2_imag = new double[tx.length];
            for (int i = 0; i < tx2.length; i++) {
                tx2[i] = temp_real[i] * frtz[i] * wnz_temp[i];
                tx2_imag[i] = temp_imag[i] * frtz[i] * wnz_temp[i];
            }
            //M.compare(tx2, "tx2");
            //M.compare(tx2_imag, "tx2_imag");

            fft_object = FFT.get(tx2.length);
            //fft_real = new double[tx2.length];
            //System.arraycopy(tx2, 0, fft_real, 0, tx2.length);
            //fft_imag = new double[tx2.length];
            //System.arraycopy(tx2_imag, 0, fft_imag, 0, tx2.length);
            fft_object.ifft(tx2, tx2_imag);

            tx2 = M.fftshift_array(tx2);

            for (int i = 0; i < tx2.length; i++) {
                tx2[i] *= ww[i];
            }

            //M.compare(tx2, "tx2-" + ii);
            //M.compare(tx2, "tx2");
            double sqrt_nf0 = Math.sqrt(nf0);
            //double[] sy_temp = new double[bb.length];
            for (int i = 0; i < bb.length; i++) { // cu +1?? XXX de testat pe
                sy[(int) (bb[i] + nidx)][0] += (tx[i] * sqrt_nf0 + tx2[i]) * ((f0raw[ii] > 0) ? 1.0 : 0.0);
                //sy_temp[i] = sy[(int) (bb[i] + nidx)][0];
            }
            //M.compare(sy_temp, "sy_temp" + ii);
            //M.compare(sy_temp, "sy_temp");

            idx = (idx + nf0);
            //M.log("idx este idx +nf0 (" + nf0 + ") = " + idx);

            iin = (int) Math.min(
                    Math.max(1, Math.round(imap[(int) (Math.round(idx)) - 1])),
                    Math.min(njj, f0raw.length));
            //M.log("iin inainte de if este: " + iin);

            if (f0raw[ii - 1] == 0 && f0raw[iin - 1] > 0) {
                double idxo = idx;

                //M.log("in if. ii = " + ii + " iin = " + iin);
                String ss = "";
                for (int q = ii - 1; q < iin; q++) {
                    ss += " f0raw[" + q + "]=" + f0raw[q];
                }
                //M.log(ss);

                int ipos = -1;
                for (int i = ii - 1; i < iin; i++) {
                    if (f0raw[i] > 0) {
                        ipos = i + 1;
                        break;
                    }
                }
                //M.log("   ipos dupa find : " + ipos);

                if (ipos == -1) // idx=idxo;
                {
                    idx = idxo;
                } // else
                else {
                    idx = Math.max(idxo - nf0 + 1, rmap[ipos - 1]);
                    //M.log("\t idx devine max idxo " + idxo + " - nf0 " + nf0 + " +1 , rmap[ipos-1] " + rmap[ipos - 1]);
                }

            }
            //M.log("  noul idx este: " + idx);
        }

        ii = 1;

        idx = 1;

        double f0 = 1000;

        //double[] wlcutfric = new double[fxa.length];
        //for (int i = 0; i < fxa.length; i++) {
        //    wlcutfric[i] = 1.0 / (1.0 + Math.exp(-14 * (fxa[i] - lowcutfav) / lowcutfav));
        //}
        icntr = 0;

        double nf0 = Double.MIN_VALUE;
        double[] cccp = new double[ff.length];
        double[] cccp2 = new double[fftl];
        double[] fffx = new double[cccp2.length];
        double[] fffx_imag = new double[cccp2.length];
        double[] exp_fffx_real = new double[fffx.length];
        double[] exp_fffx_imag = new double[fffx.length];
        double[] ttx = new double[fffx.length];
        double[] tnx;
        
        while ((idx < nsyn - fftl) && (ii < f0l.length)) {

            //M.log("-----------------");
            icntr++;

            //iin = (int) Math.min(Math.max(1, Math.round(imap[(int) Math.round(idx)])), Math.min(njj, f0raw.length));
            double nidx = Math.round(idx);

            //M.log(">START ii " + ii + " idx " + idx + " iin " + iin + " nidx " + nidx);
            if (f0raw[ii - 1] == 0) {

                //double[] ff = new double[idcv.length + rbb2.length];
                for (int i = 0; i < idcv.length; i++) {
                    ff[i] = n2sgram[(int) idcv[i] - 1][ii - 1];
                }
                for (int i = idcv.length; i < idcv.length + rbb2.length; i++) {
                    ff[i] = n2sgram[(int) idcv[(int) rbb2[i - idcv.length] - 1] - 1][ii - 1];
                }
                //M.compare(ff, "ff" + ii);

                
                for (int i = 0; i < cccp.length; i++) {
                    cccp[i] = Math.log(ff[i] + dmx_100000);
                }
                //double[] delme_cccp = cccp.clone();
                //M.compare(delme_cccp, "delme_cccp" + ii);

                fft_object = FFT.get(cccp.length);
                fft_real = new double[cccp.length];
                System.arraycopy(cccp, 0, fft_real, 0, cccp.length);
                fft_imag = new double[cccp.length];
                fft_object.fft(fft_real, fft_imag);
                System.arraycopy(fft_real, 0, cccp, 0, cccp.length);
                //M.compare(cccp, "cccp" + ii);

                cccp2[0] = cccp[0];
                for (int i = 1; i < fftl / 2; i++) {
                    cccp2[i] = cccp[i] * 2.0;
                }
                //M.compare(cccp2, "cccp2" + ii);

                for (int i = 0; i < fffx.length; i++) {
                    fffx[i] = cccp2[i] * lft[i] / fftl;
                }
                fft_object = FFT.get(fffx.length);
                fft_real = new double[fffx.length];
                System.arraycopy(fffx, 0, fft_real, 0, fffx.length);
                fft_imag = new double[fffx.length];
                fft_object.fft(fft_real, fft_imag);
                System.arraycopy(fft_real, 0, fffx, 0, fffx.length);
                System.arraycopy(fft_imag, 0, fffx_imag, 0, fffx.length);

                //M.compare(fffx, "fffx" + ii);
                //M.compare(fffx_imag, "fffx_imag" + ii);
                nf0 = fs / f0;

                for (int i = 0; i < exp_fffx_real.length; i++) {
                    exp_fffx_real[i] = Math.exp(fffx[i]) * Math.cos(fffx_imag[i]);
                    exp_fffx_imag[i] = Math.exp(fffx[i]) * Math.sin(fffx_imag[i]);
                }

                //double[] ttx_imag = new double[fffx.length];

                fft_object = FFT.get(ttx.length);
                fft_real = new double[ttx.length];
                System.arraycopy(exp_fffx_real, 0, fft_real, 0, ttx.length);
                fft_imag = new double[ttx.length];
                System.arraycopy(exp_fffx_imag, 0, fft_imag, 0, ttx.length);
                fft_object.ifft(fft_real, fft_imag);
                ttx = M.fftshift_array(fft_real);
                //M.compare(ttx, "ttx" + ii);

                double[] rx = new double[(int) Math.round(nf0)];
                rx = M.randn_1row_ncolumns(rx.length);

                double mean_rx = M.mean(rx);
                for (int i = 0; i < rx.length; i++) {
                    rx[i] -= mean_rx;
                }
                tnx = M.fftfilt(rx, ttx);

                //double[] ssy_temp = new double[bb.length];
                for (int i = 0; i < bb.length; i++) {
                    sy[(int) (bb[i] + nidx)][0] += tnx[bb[i]];
                    //ssy_temp[i] = sy[(int) (bb[i] + nidx)][0];
                }
                //M.compare(ssy_temp, "ssy_temp" + ii);
            }

            idx = idx + nf0;
            //M.log(" nf0 = " + nf0);
            ii = (int) Math.round(imap[(int) Math.round(idx) - 1]);
            //M.log(">END ii " + ii + " idx " + idx + " iin " + iin + " nidx " + nidx);
        }

        //M.log("FINAL ix = " + ix);
        //double[][] sy2 = new double[ix][1];
        double[] fin = new double[ix + 1];
        double max_abs = Math.abs(fin[0]);
        
        for (int i = 0; i < ix; i++) {
            //sy2[i][0] = sy[fftl / 2 + i][0];
            //fin[i] = sy2[i][0];
        	fin[i] = sy[fftl / 2 + i][0];
            if (Math.abs(fin[i]) > max_abs) 
                max_abs = Math.abs(fin[i]);
        }
        //M.compare(fin, "fin");
        
        /*for (int i = 0; i < fin.length; i++) {
            if (Math.abs(fin[i]) > max_abs) {
                max_abs = Math.abs(fin[i]);
            }
        }
        */
        
        for (int i = 0; i < fin.length; i++) {
            fin[i] = fin[i] / max_abs;
        }

        long end_time = System.currentTimeMillis();
        long difference = end_time - start_time;
        System.out.println("DONE (" + difference + "ms)");

        return fin;
    }

}
