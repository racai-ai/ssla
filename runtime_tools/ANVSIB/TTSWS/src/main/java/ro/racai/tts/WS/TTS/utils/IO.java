package ro.racai.tts.WS.TTS.utils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import ro.racai.tts.WS.TTS.utils.UnicodeBOMInputStream.BOM;


public class IO {
	
	private static boolean verbose = true;
	
	static public BufferedReader openFile (String filePath) throws IOException {		
		if(!((File)new File(filePath)).exists()) throw new IOException("File "+filePath+" does not exist!");

		FileInputStream f1 = new FileInputStream(filePath);
		UnicodeBOMInputStream ubis = new UnicodeBOMInputStream(f1);
		//skip BOM if present;
		if(ubis.getBOM()!=BOM.NONE) ubis.skipBOM();		
		InputStreamReader f2 = new InputStreamReader(ubis,"UTF-8");
		BufferedReader f3 = new BufferedReader(f2);		
		return f3;
	}
	
	static public BufferedReader openFile (InputStream stream) throws IOException {
		UnicodeBOMInputStream ubis = new UnicodeBOMInputStream(stream);
		//skip BOM if present;
		if(ubis.getBOM()!=BOM.NONE) ubis.skipBOM();		
		InputStreamReader f2 = new InputStreamReader(ubis,"UTF-8");
		BufferedReader f3 = new BufferedReader(f2);		
		return f3;
	}
	
	public static boolean isVerbose() {
		return verbose;
	}

	public static void setVerbose(boolean verbose) {
		IO.verbose = verbose;
	}

	static public DataInputStream readBinaryFile (String filePath) throws IOException {		
		if(!((File)new File(filePath)).exists()) throw new IOException("File "+filePath+" does not exist!");

		FileInputStream f1 = new FileInputStream(filePath);				
		BufferedInputStream f2 = new BufferedInputStream(f1);
		DataInputStream f3 = new DataInputStream(f2);		
		return f3;
	}

	static public DataOutputStream writeBinaryFile (String filePath) throws IOException {
		FileOutputStream f1 = new FileOutputStream(filePath);				
		BufferedOutputStream f2 = new BufferedOutputStream(f1);
		DataOutputStream f3 = new DataOutputStream(f2);		
		return f3;
	}

	static public String [] readAllLines(String filePath) throws IOException{
		BufferedReader br =openFile(filePath);
		List<String> lines=new ArrayList<String>();
		String line;
		while ((line=br.readLine())!=null){
			lines.add(line);
		}
		br.close();
		return lines.toArray(new String[0]);
	}
	
	static public void writeAllLines(String[] data, String filePath) throws IOException{
		BufferedWriter br = new BufferedWriter(new FileWriter(filePath,false));			
		for(String line:data) {
			br.write(line);
			if(!line.endsWith("\n")) br.write("\n");
		}	
		br.close();		
	}

	public static BufferedWriter openFileForWriting(String filePath, boolean append) throws IOException {
		return new BufferedWriter (new OutputStreamWriter(new FileOutputStream(filePath,append),"UTF-8"));		
	}

	
	public static void outln (String out) {
		if(!verbose) return;
		System.out.println(out);	
	}
	
	public static void out (String out) {
		if(!verbose) return;
		System.out.print(out);	
	}
	
	public static void outln () {
		if(!verbose) return;
		System.out.println();	
	}
	
	public static void outln (double val) {
		if(!verbose) return;
		System.out.println(val);	
	}
	

	public static void errln (String err) {
		if(!verbose) return;
		System.err.println(err);	
	}
	
	public static void err (String err) {
		if(!verbose) return;
		System.err.print(err);	
	}
	
	public static void errln () {
		if(!verbose) return;
		System.err.println();	
	}
	
}
