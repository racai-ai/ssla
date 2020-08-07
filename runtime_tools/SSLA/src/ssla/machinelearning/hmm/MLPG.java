/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ssla.machinelearning.hmm;

import java.io.BufferedWriter;
import java.io.FileWriter;
import ssla.machinelearning.RawFrame;
import ssla.machinelearning.SpeechFrame;
import ssla.signalprocessing.Bap2Ap;

/**
 *
 * Maximum Likelihood Parameter Generation
 *
 * @author tibi
 */
public class MLPG {

    /**
     * Genereaza parametri vocali pe baza criteriului Maximum Likelihood
     * Parameter Generation (MLPG)
     *
     * @param frames - setul de obeservatii, compus din media valorilor, a
     * vitezei si a accelaratiei precum si a varianta acestora
     * @return un set de cadre de parametri vocali care maximizeaza
     * probabilitatea observatiilor tinand cont de medie, viteza, acceleratie,
     * media lor si a variantei
     */
    public SpeechFrame[] mlpg(RawFrame[] frames, boolean useGV) throws Exception {
        //aflam numarul total de cadre de voce
        int total = 0;
        int totalVoiced = 0;
        for (int i = 0; i < frames.length; i++) {
            if (frames[i].f0[0] != 0) {
                totalVoiced += frames[i].duration;
            }
            total += frames[i].duration;
        }
        RawFrame[] copy = new RawFrame[total];
        //facem o copie a cadrelor de voce neprocesate, tinand cont de durata
        int cnt = 0;
        for (int i = 0; i < frames.length; i++) {
            for (int t = 0; t < frames[i].duration; t++) {
                copy[cnt++] = new RawFrame(frames[i]);
            }
        }

        //niste magie
        PStream mgc = new PStream(copy[0].mgc.length, total, 100);
        PStream f0 = new PStream(copy[0].f0.length, totalVoiced, 100);
        PStream bap = new PStream(copy[0].bap.length, total, 100);

        int voicedI = 0;
        for (int i = 0; i < copy.length; i++) {
            double[] tempM = new double[copy[i].mgc.length];
            for (int k = 0; k < tempM.length; k++) {
                tempM[k] = copy[i].mgc[k];
            }
            double[] tempV = new double[copy[i].mgcVar.length];
            for (int k = 0; k < tempV.length; k++) {
                tempV[k] = copy[i].mgcVar[k];
            }
            mgc.setGvSwitch(i, copy[i].gv);
            mgc.setMseq(i, tempM);
            mgc.setVseq(i, tempV);
            //TODO: de refacut asta

            tempM = new double[copy[i].gvMGC.length];
            for (int k = 0; k < tempM.length; k++) {
                tempM[k] = copy[i].gvMGC[k];
            }
            tempV = new double[copy[i].gvMGCVar.length];
            for (int k = 0; k < tempV.length; k++) {
                tempV[k] = copy[i].gvMGCVar[k];
            }

            mgc.setGvMeanVar(tempM, tempV);
            //bap
            tempM = new double[copy[i].bap.length];
            for (int k = 0; k < tempM.length; k++) {
                tempM[k] = copy[i].bap[k];
            }
            tempV = new double[copy[i].bapVar.length];
            for (int k = 0; k < tempV.length; k++) {
                tempV[k] = copy[i].bapVar[k];
            }

            bap.setMseq(i, tempM);
            bap.setVseq(i, tempV);
            bap.setGvSwitch(i, copy[i].gv);
            tempM = new double[copy[i].gvBap.length];
            for (int k = 0; k < tempM.length; k++) {
                tempM[k] = copy[i].gvBap[k];
            }
            tempV = new double[copy[i].gvBapVar.length];
            for (int k = 0; k < tempV.length; k++) {
                tempV[k] = copy[i].gvBapVar[k];
            }
            bap.setGvMeanVar(tempM, tempV);

            //end bap
            if (copy[i].f0[0] != 0) {
                tempM = new double[copy[i].f0.length];
                for (int k = 0; k < tempM.length; k++) {
                    tempM[k] = copy[i].f0[k];
                }
                tempV = new double[copy[i].f0Var.length];
                for (int k = 0; k < tempV.length; k++) {
                    tempV[k] = copy[i].f0Var[k];
                }
                f0.setVseq(voicedI, tempV);

                f0.setMseq(voicedI, tempM);

                tempM = new double[copy[i].gvF0.length];
                for (int k = 0; k < tempM.length; k++) {
                    tempM[k] = copy[i].gvF0[k];
                }
                tempV = new double[copy[i].gvF0Var.length];
                for (int k = 0; k < tempV.length; k++) {
                    tempV[k] = copy[i].gvF0Var[k];
                }

                f0.setGvSwitch(voicedI++, copy[i].gv);
                f0.setGvMeanVar(tempM, tempV);
            }

        }

        mgc.fixDynFeatOnBoundaries();
        f0.fixDynFeatOnBoundaries();
        bap.fixDynFeatOnBoundaries();
        mgc.mlpg(useGV);
        f0.mlpg(useGV);
        bap.mlpg(useGV);
        //BufferedWriter bw = new BufferedWriter(new FileWriter("speech_raw.txt"));
        //rezultat final
        SpeechFrame[] sfs = new SpeechFrame[copy.length];
        voicedI = 0;
        double[][] bapVect = new double[sfs.length][];
        for (int i = 0; i < bapVect.length; i++) {
            bapVect[i] = new double[copy[i].bap.length / 3];
            for (int k = 0; k < bapVect[i].length; k++) {
                bapVect[i][k] = bap.getPar(i, k);
            }
        }
        Bap2Ap bap2ap = new Bap2Ap();
        double[][] ap = bap2ap.bap2ap(bapVect);
        for (int i = 0; i < sfs.length; i++) {
            sfs[i] = new SpeechFrame(copy[i]);
            for (int k = 0; k < sfs[i].mgc.length; k++) {
                sfs[i].mgc[k] = (float) mgc.getPar(i, k);
                //bw.write(sfs[i].mgc[k] + " ");
            }
            sfs[i].ap = new float[ap[i].length];
            for (int k = 0; k < sfs[i].ap.length; k++) {
                sfs[i].ap[k] = (float) ap[i][k];
            }
            //bw.write("\n");
            if (sfs[i].f0 != 0) {
                float f = (float) Math.exp(f0.getPar(voicedI++, 0));
                sfs[i].f0 = f;
            }
        }
        //bw.close();
        return sfs;
    }
}
