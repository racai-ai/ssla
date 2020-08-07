/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.TTS.machinelearning;

import java.io.File;
import java.io.IOException;
import java.util.Random;
import ro.racai.tts.WS.TTS.hmm.MLPG;
import ro.racai.tts.WS.TTS.utils.Log;


/**
 *
 * @author tibi
 */
public class SpeechModel {

    DecissionTree treeBap;
    DecissionTree treeDur;
    DecissionTree treeMGC;
    DecissionTree treeLF0;
    DecissionTree treeGVBap;
    DecissionTree treeGVLF0;
    DecissionTree treeGVMGC;
    DecissionTree treeGVSwitch;
    int maxStates = 0;

    public SpeechModel(String path) throws IOException {
        Log.i("Loading decision trees and PDFs.");
        long start = System.currentTimeMillis();
        treeDur = new DecissionTree(path.concat(File.separator).concat("tree-dur.inf"), path.concat(File.separator).concat("dur.pdf"), "dur", 1);
        treeLF0 = new DecissionTree(path.concat(File.separator).concat("tree-lf0.inf"), path.concat(File.separator).concat("lf0.pdf"), "fl0", 3);
        treeMGC = new DecissionTree(path.concat(File.separator).concat("tree-mgc.inf"), path.concat(File.separator).concat("mgc.pdf"), "mgc", 3);
        treeBap = new DecissionTree(path.concat(File.separator).concat("tree-bap.inf"), path.concat(File.separator).concat("bap.pdf"), "bap", 3);
        treeGVBap = new DecissionTree(path.concat(File.separator).concat("tree-gv-bap.inf"), path.concat(File.separator).concat("gv-bap.pdf"), "gv-bap", 1);
        treeGVLF0 = new DecissionTree(path.concat(File.separator).concat("tree-gv-lf0.inf"), path.concat(File.separator).concat("gv-lf0.pdf"), "gv-lf0", 1);
        treeGVMGC = new DecissionTree(path.concat(File.separator).concat("tree-gv-mgc.inf"), path.concat(File.separator).concat("gv-mgc.pdf"), "gv-mgc", 1);
        treeGVSwitch = new DecissionTree(path.concat(File.separator).concat("gv-switch.inf"), null, "gv-switch", 1);
        long stop = System.currentTimeMillis();
        maxStates = Math.max(maxStates, treeDur.getNumStates());
        maxStates = Math.max(maxStates, treeLF0.getNumStates());
        maxStates = Math.max(maxStates, treeMGC.getNumStates());
        maxStates = Math.max(maxStates, treeBap.getNumStates());
        maxStates = Math.max(maxStates, treeGVBap.getNumStates());
        maxStates = Math.max(maxStates, treeGVLF0.getNumStates());
        maxStates = Math.max(maxStates, treeGVMGC.getNumStates());
        Log.i("Loading trees finished in " + (stop - start) + "ms");
    }

    public RawFrame[] parse(String[] features, float vuv, boolean useGV) {
        int cnt = 0;
        RawFrame[] frames = new RawFrame[features.length * 5];
        for (int i = 0; i < features.length; i++) {
            //pentru durata

            float[] tmp;
            for (int k = 0; k < maxStates; k++) {
                RawFrame f = new RawFrame();
                f.lab = features[i];
                String pdf = getPDF(treeDur, k, features[i]);
                tmp = treeDur.getPdf(pdf);
                f.duration = (int) ((tmp[k] + 0.5));
                pdf = getPDF(treeLF0, k, features[i]);
                f.f0Label = pdf;
                tmp = treeLF0.getPdf(pdf);
                /*f.f0 = (float) Math.exp(tmp[0]);
                 if (tmp[6] < vuv) {
                 f.f0 = 0;
                 }*/
                f.f0 = new float[(tmp.length - 1) / 2];
                f.f0Var = new float[f.f0.length];
                for (int z = 0; z < (tmp.length - 1) / 2; z++) {
                    f.f0[z] = tmp[z];
                    f.f0Var[z] = tmp[z + f.f0.length];
                }

                if (tmp[6] < vuv) {
                    f.f0[0] = 0;
                }
                pdf = getPDF(treeMGC, k, features[i]);
                tmp = treeMGC.getPdf(pdf);
                f.mgc = new float[tmp.length / 2];
                f.mgcVar = new float[tmp.length / 2];
                for (int l = 0; l < tmp.length / 2; l++) {
                    f.mgc[l] = tmp[l];
                    f.mgcVar[l] = tmp[f.mgcVar.length + l];
                }

                pdf = getPDF(treeBap, k, features[i]);
                tmp = treeBap.getPdf(pdf);
                f.bap = new float[tmp.length / 2];
                f.bapVar = new float[tmp.length / 2];
                for (int l = 0; l < tmp.length / 2; l++) {
                    f.bap[l] = tmp[l];
                    f.bapVar[l] = tmp[f.bapVar.length + l];
                }

                pdf = getPDF(treeGVBap, k, features[i]);
                tmp = treeGVBap.getPdf(pdf);

                f.gvBap = new float[tmp.length / 2];
                for (int l = 0; l < f.gvBap.length; l++) {
                    f.gvBap[l] = tmp[l];
                }
                f.gvBapVar = new float[tmp.length / 2];
                for (int l = 0; l < f.gvBapVar.length; l++) {
                    f.gvBapVar[l] = tmp[l + f.gvBap.length];
                }

                pdf = getPDF(treeGVLF0, k, features[i]);
                tmp = treeGVLF0.getPdf(pdf);

                f.gvF0 = new float[tmp.length / 2];
                for (int l = 0; l < f.gvF0.length; l++) {
                    f.gvF0[l] = tmp[l];
                }
                f.gvF0Var = new float[tmp.length / 2];
                for (int l = 0; l < f.gvF0Var.length; l++) {
                    f.gvF0Var[l] = tmp[l + f.gvF0Var.length];
                }

                pdf = getPDF(treeGVMGC, k, features[i]);
                tmp = treeGVMGC.getPdf(pdf);
                f.gvMGC = new float[tmp.length / 2];
                for (int l = 0; l < f.gvMGC.length; l++) {
                    f.gvMGC[l] = tmp[l];
                }
                f.gvMGCVar = new float[tmp.length / 2];
                for (int l = 0; l < f.gvMGCVar.length; l++) {
                    f.gvMGCVar[l] = tmp[l + f.gvMGCVar.length];
                }

                //si acum pentru gv-switch
                String use_GV = treeGVSwitch.parse(features[i], 0);
                if (use_GV.equals("gv-switch_2")) {
                    f.gv = true;
                } else {
                    f.gv = false;
                }

                frames[i * 5 + k] = f;

            }
        }
        return frames;
    }

    private String getPDF(DecissionTree tree, int state, String features) {
        if (state >= tree.getNumStates()) {
            state = tree.getNumStates() - 1;
        }
        return tree.parse(features, state);
    }

    /**
     * Genereaza parametrii pentru cadrele de voce intermediare, pa baza
     * vitezei, acceleratiei si a variantelor date
     *
     * @param raw - cadre de voce generate din arbori
     * @return - o lista de cadre reale de voce ce trebuie folosite in Vocoder
     */
    public SpeechFrame[] getActualFrames(RawFrame[] raw, boolean useGV) throws Exception {
        //TODO: @Adrian: aici codul :D:D:D
        MLPG mlpg = new MLPG();
        return mlpg.mlpg(raw, useGV);
    }
    boolean generate = false;
    Random rand = new Random(1234);
    double z0 = 0, z1 = 0;

    float generateGaussianNoise(float mu, float sigma) {
        double epsilon = Float.MIN_NORMAL;
        double two_pi = 2.0 * 3.14159265358979323846;

        generate = !generate;

        if (!generate) {
            return (float) (z1 * sigma + mu);
        }

        double u1, u2;

        do {
            u1 = rand.nextDouble();
            u2 = rand.nextDouble();
        } while (u1 <= epsilon);

        z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(two_pi * u2);
        z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(two_pi * u2);
        return (float) (z0 * sigma + mu);
    }
}
