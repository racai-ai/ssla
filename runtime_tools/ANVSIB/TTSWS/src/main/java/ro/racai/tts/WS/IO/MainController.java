/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ro.racai.tts.WS.IO;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import javax.servlet.http.HttpServletResponse;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpRequest;
import org.springframework.http.MediaType;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.RestController;
import ro.racai.tts.WS.NLP.language.LanguagePipe;
import ro.racai.tts.WS.NLP.language.Token;
import ro.racai.tts.WS.TTS.SSLARuntime;

/**
 *
 * @author tibi
 */
@RestController
public class MainController {

    LanguagePipe mlpla = null;
    SSLARuntime ssla = null;

    @RequestMapping("/version")
    public String version(Model model) {
        return "RACAI TTS Server beta version";
    }

    @RequestMapping("/synth.php")
    public HttpEntity<byte[]> synthesize(@RequestParam(value = "lang", defaultValue = "ro") String lang, @RequestParam(value = "coder", defaultValue = "mlsa") String coder, @RequestParam(value = "voice", defaultValue = "anca") String voice, @RequestParam(value = "text", defaultValue = "Acesta este un test") String text, HttpServletResponse response) throws IOException, FileNotFoundException, ClassNotFoundException, InstantiationException, IllegalAccessException, Exception {
        long start=System.currentTimeMillis();
        System.out.println("Received text :'"+text+"'");
        System.out.println("Coder: "+coder);
        if (mlpla == null) {
            System.out.println("Models are not loaded. Attempting to read from disk");
            mlpla = new LanguagePipe("etc/languagepipe.conf", "models/ro");
            ssla = new SSLARuntime("models/ro/" + voice);
        }
        HttpHeaders header = new HttpHeaders();
        header.setContentType(new MediaType("audio", "wav"));
        header.setExpires(0);
        header.setCacheControl("must-revalidate");

        //NLP
        Token[][] tokens = mlpla.process(text);
        StringBuilder sb=new StringBuilder();
        for (int i = 0; i < tokens.length; i++) {
            
            sb.append(mlpla.output.print(tokens[i]));
        }
        byte []wave=ssla.synthesize(coder.equals("mlsa")?SSLARuntime.VocoderType.MLSA:SSLARuntime.VocoderType.STRAIGHT,sb.toString().split("\n"));
        long stop=System.currentTimeMillis();
        System.out.println("Done in "+(stop-start)+"ms");
        return new HttpEntity<>(wave, header);
    }

}
