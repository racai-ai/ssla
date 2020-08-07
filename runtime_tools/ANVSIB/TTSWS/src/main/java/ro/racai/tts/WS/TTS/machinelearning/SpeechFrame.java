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
public class SpeechFrame {

    public float f0;
    public float[] mgc;
    public float[] ap;

    public SpeechFrame(RawFrame copy) {
        this.f0 = copy.f0[0];
        this.mgc = new float[copy.mgc.length / 3];
        System.arraycopy(copy.mgc, 0, mgc, 0, mgc.length);
        this.ap = new float[copy.bap.length / 3];
        //TODO: aici trebuie convertit din BAP in AP (este exemplu in antrenarea HTS despre cum se face) - ei ruleaza "merge" din SPTK
        System.arraycopy(copy.bap, 0, this.ap, 0, ap.length);
    }
}
