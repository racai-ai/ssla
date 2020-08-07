/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.TTS;

import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import ro.racai.tts.WS.TTS.machinelearning.RawFrame;
import ro.racai.tts.WS.TTS.machinelearning.SpeechFrame;
import ro.racai.tts.WS.TTS.machinelearning.SpeechModel;
import ro.racai.tts.WS.TTS.signalprocessing.MLSAVocoder;
import ro.racai.tts.WS.TTS.signalprocessing.STRAIGHTVocoder;
import ro.racai.tts.WS.TTS.signalprocessing.WavWrite;
import ro.racai.tts.WS.TTS.utils.Log;

/**
 *
 * @author tibi
 */
public class SSLARuntime {

    SpeechModel sm;
    MLSAVocoder mlsa = new MLSAVocoder();
    STRAIGHTVocoder straight;
    int sr = 48000;

    public SSLARuntime(String models) throws IOException {
        sm = new SpeechModel(models);
        straight = new STRAIGHTVocoder();
    }

    public byte[] synthesize(VocoderType vocoder, String[] labs) throws Exception {
        RawFrame[] tmp = sm.parse(labs, 0.5f, true);

        SpeechFrame[] speechFrames = sm.getActualFrames(tmp, true);
        float[] f0 = new float[speechFrames.length];
        float[][] mgc = new float[speechFrames.length][];
        float[][] ap = new float[speechFrames.length][];
        for (int i = 0; i < speechFrames.length; i++) {
            f0[i] = speechFrames[i].f0;
            mgc[i] = speechFrames[i].mgc;
            ap[i] = speechFrames[i].ap;
        }

        double[] audio = null;
        if (vocoder == VocoderType.MLSA) {
            audio = mlsa.runMLSAVocoder(mgc, f0, mgc[0].length, sr / 200, sr);
        }else{
            audio = straight.runSTRAIGHTVocoder(mgc, ap, f0, sr);
        }
        ByteArrayOutputStream bos=new ByteArrayOutputStream();
        DataOutputStream dos = new DataOutputStream(bos);
        Log.i("Creating output WAVE");
        WavWrite.Save2Wave(audio, sr, dos);
        dos.close();
        return bos.toByteArray();
    }

    public static enum VocoderType {
        STRAIGHT,
        MLSA
    }
}
