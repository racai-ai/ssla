/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.TTS.utils;

public class Log {

    public static void w(String message) {
        System.err.println("WARN: " + message);
    }

    public static void i(String message) {
        System.err.println("INFO: " + message);
    }

    public static void e(String message) {
        System.err.println("ERROR: " + message);
    }

}
