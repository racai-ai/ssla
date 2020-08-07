/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.TTS.machinelearning;

/**
 *
 * @author tibi
 */
public class RawFrame {

    public boolean gv;

    public int duration; //se inmulteste cu 5ms pentru durata reala
    public float[] f0;
    public float[] f0Var;
    public float[] mgc;
    public float[] mgcVar;
    public float[] bap;
    public float[] bapVar;

    public float[] gvF0;
    public float[] gvF0Var;
    public float[] gvBap;
    public float[] gvBapVar;
    public float[] gvMGC;
    public float[] gvMGCVar;

    public String f0Label;
    public String MGCLabel;
    public String BAPLabel;
    public String lab;

    public RawFrame(RawFrame frame) {
        f0Label = frame.f0Label;
        MGCLabel = frame.MGCLabel;
        BAPLabel = frame.BAPLabel;
        lab = frame.lab;
        duration = frame.duration;
        f0 = new float[frame.f0.length];
        System.arraycopy(frame.f0, 0, f0, 0, f0.length);
        f0Var = new float[frame.f0Var.length];
        System.arraycopy(frame.f0Var, 0, f0Var, 0, f0Var.length);
        mgc = new float[frame.mgc.length];
        System.arraycopy(frame.mgc, 0, mgc, 0, mgc.length);
        mgcVar = new float[frame.mgcVar.length];
        System.arraycopy(frame.mgcVar, 0, mgcVar, 0, mgcVar.length);
        bap = new float[frame.bap.length];
        System.arraycopy(frame.bap, 0, bap, 0, bap.length);
        bapVar = new float[frame.bapVar.length];
        System.arraycopy(frame.bapVar, 0, bapVar, 0, bapVar.length);

        gvMGC = new float[frame.gvMGC.length];
        System.arraycopy(frame.gvMGC, 0, gvMGC, 0, gvMGC.length);
        gvMGCVar = new float[frame.gvMGCVar.length];
        System.arraycopy(frame.gvMGCVar, 0, gvMGCVar, 0, gvMGCVar.length);
        gvF0 = new float[frame.gvF0.length];
        System.arraycopy(frame.gvF0, 0, gvF0, 0, gvF0.length);
        gvF0Var = new float[frame.gvF0Var.length];
        System.arraycopy(frame.gvF0Var, 0, gvF0Var, 0, gvF0Var.length);
        gvBap = new float[frame.gvBap.length];
        System.arraycopy(frame.gvBap, 0, gvBap, 0, gvBap.length);
        gvBapVar = new float[frame.gvBapVar.length];
        System.arraycopy(frame.gvBapVar, 0, gvBapVar, 0, gvBapVar.length);
        this.gv = frame.gv;

    }

    public RawFrame() {
        //constructor gol
    }
}
