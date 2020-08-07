/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ssla.signalprocessing;

import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author tibi
 */
public class WavWrite {

    public static void Save2Wave(double[] audio, int sr, DataOutputStream dos) throws IOException {
        dos.write("RIFF".getBytes());
        int len = 44 - 8 + audio.length * 2;
        dos.writeByte(len & 0xFF);
        dos.writeByte((len >> 8) & 0xFF);
        dos.writeByte((len >> 16) & 0xFF);
        dos.writeByte((len >> 24) & 0xFF);
        dos.write("WAVE".getBytes());
        dos.write("fmt ".getBytes());
        dos.writeByte((byte) 16);
        dos.writeByte((byte) 0);
        dos.writeByte((byte) 0);
        dos.writeByte((byte) 0);
        dos.writeByte((byte) 1); //PCM
        dos.writeByte((byte) 0);
        dos.writeByte((byte) 1); //1 canal
        dos.writeByte((byte) 0);
        dos.writeByte(sr & 0xFF); //sample rate
        dos.writeByte((sr >> 8) & 0xFF);
        dos.writeByte((sr >> 16) & 0xFF);
        dos.writeByte((sr >> 24) & 0xFF);
        int x = sr * 16 * 1 / 8;
        dos.writeByte((byte) (x & 0xFF)); //bps
        dos.writeByte((byte) ((x >> 8) & 0xFF));
        dos.writeByte((byte) ((x >> 16) & 0xFF));
        dos.writeByte((byte) ((x >> 24) & 0xFF));
        dos.writeByte(2);
        dos.writeByte(0);
        dos.writeByte((byte) 16); //16 bit integer   
        dos.writeByte((byte) 0);
        dos.write("data".getBytes());
        len = audio.length * 2;
        dos.writeByte(len & 0xFF);
        dos.writeByte((len >> 8) & 0xFF);
        dos.writeByte((len >> 16) & 0xFF);
        dos.writeByte((len >> 24) & 0xFF);
        byte[] temp = new byte[audio.length * 2];
        for (int i = 0; i < audio.length; i++) {
            int samp = (int) (audio[i] * 32767);
//            if (samp > 32768) {
//                samp = 32768;
//            } else if (samp < -32760) {
//                samp = -32760;
//            }
            int index = i * 2;
            temp[index + 1] = (byte) ((samp >> 8) & 0xFF);
            temp[index] = (byte) (samp & 0xFF);
        }
        dos.write(temp);
    }
}
