/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.NLP.language;

import java.util.List;

/**
 *
 * @author tibi
 */
public interface BaseInputProcessor {
    public List<List<Token>> parseInput(String text);
    public void loadModel(String path);
}
