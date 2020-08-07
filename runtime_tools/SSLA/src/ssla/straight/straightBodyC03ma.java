/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ssla.straight;


import ssla.utils.IO;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import ssla.signalprocessing.FFT;

/**
 *
 * @author Echo
 */
public class straightBodyC03ma {
    // scos f0var si f0varL - meaningless
    public static double jstraightBodyC03ma (double[] x,double fs,double shiftm,int fftl,double[] f0raw, double eta, double pc) throws Exception {
    /*
      function [n2sgram,nsgram]=straightBodyC03ma(x,fs,shiftm,fftl,f0raw,f0var,f0varL,eta,pc,imgi);
%  [n2sgram,nsgram]=straightBodyC03ma(x,fs,shiftm,fftl,f0raw,f0var,f0varL,eta,pc,imgi)
%	n2sgram		: smoothed spectrogram
%	nsgram		: isometric spectrogram
%	x		: input waveform
%	fs		: sampling frequency (Hz)
%	shiftm		: frame shift (ms)
%	fftl		: length of FFT
%	f0raw		: Pitch information to gude analysis (TEMPO) assumed
%	f0var		: expected f0 variance including zerocross information
%	f0varL		: expected f0 variance 
%	eta		: 
%	pc		:
%	imgi	: display indicator 1: display on (default), 0: off

%	f0shiftm	: frame shift (ms) for F0 analysis
*/
        
        
//if nargin==9; imgi=1; end; % 10/Sept./2005
//f0l=f0raw(:);
    	double f0 = 0; // e definit mai jos in for si trebuie sa ramana dupa
        double[] f0l=f0raw.clone();
//framem=40;
        int framem=40;        
//fx=((1:fftl)-fftl/2-1)/fftl*fs;
        double[] fx = new double[fftl];
        for(int i = 0;i<fftl;i++) {
            fx[i] = ((i+1)-fftl/2-1)/fftl*fs;
        }
//framel=round(framem*fs/1000);        
        int framel=(int) M.round(framem*fs/1000);
                
/*
if fftl<framel
  disp('Warning! fftl is too small.');
  fftl=2^ceil(log(framel)/log(2) );
  disp(['New length:' num2str(fftl) ' is used.']);
end;
x=x(:)';
        */
       if(fftl<framel) {
           System.out.println("Warning! fftl is too small.");
           fftl = (int) Math.pow(2.0,Math.ceil(Math.log(framel)/Math.log(2)));
           System.out.println("New length:"+fftl);
       } 
       double[][] x_transposed = new double[x.length][1];
       x_transposed = M.transpose_lineArray(x);
       
        

/* smoothinid=1;       

if shiftm >1.5
  disp('Frame shift has to be small enough. (Less than 1 ms is recommended.)');
  disp('Temporal smoothing will be disabled');
  smoothinid=0;
end;
smoothinid=0; % modificatin for ICSLP'2002
*/
// pt ca oricum e 0:
       int smoothinid=0;

//shiftl=shiftm*fs/1000;
       double shiftl=shiftm*fs/1000;

//   High-pass filtering using 70Hz cut-off butterworth filter
//[b,a]=butter(6,70/fs*2,'high');
  
       double[] b=null;
       double[] a=null;
       
       // simulare butterworth trecesus 
		if (fs == 8000) {
			b = new double[] { 0.8992f, -5.3953f, 13.4883f, -17.9844f, 13.4883f,
					-5.3953f, 0.8992f };
			a = new double[] { 1.0000f, -5.7876f, 13.9604f, -17.9641f, 13.0061f,
					-5.0234f, 0.8086f };
		}
		if (fs == 16000) {
			b = new double[] { 0.9483f, -5.6897f, 14.2242f, -18.9656f, 14.2242f,
					-5.6897f, 0.9483f };
			a = new double[] { 1.0000f, -5.8938f, 14.4746f, -18.9602f, 13.9711f,
					-5.4909f, 0.8992f };
		}
		if (fs == 24000) {
			b = new double[] { 0.9652f, -5.7913f, 14.4782f, -19.3043f, 14.4782f,
					-5.7913f, 0.9652f };
			a = new double[] { 1.0000f, -5.9292f, 14.6485f, -19.3019f, 14.3068f,
					-5.6558f, 0.9316f };
		}
		if (fs == 48000) {
			b = new double[] { 0.9825f, -5.8947f, 14.7368f, -19.6491f, 14.7368f,
					-5.8947f, 0.9825f };
			a = new double[] { 1.0000f, -5.9646f, 14.8236f, -19.6485f, 14.6497f,
					-5.8255f, 0.9652f };
		}
		if (fs == 96000) {
			b = new double[] { 0.9912f, -5.9471f, 14.8678f, -19.8238f, 14.8678f,
					-5.9471f, 0.9912f };
			a = new double[] { 1.0000f, -5.9823f, 14.9116f, -19.8236f, 14.8239f,
					-5.9121f, 0.9825f };
		}
//xh=filter(b,a,x);
// XXX de pus filtrul, acum e hardcodat;
		
	double[] xh=M.filter2(b,a,x);
	
	// YYY
    ArrayList<Float> d = new ArrayList<Float>();
    String[] data = IO.readAllLines("D:\\Dropbox\\STRAIGHTsrc\\xh");
    for(String n : data) {
    	d.add(Float.valueOf(n));
    }
    xh = new double[d.size()];
    int iiq = 0;
    for (Float f : d) {
        xh[iiq++] = (f != null ? f : Float.NaN); // Or whatever default you want.
    }
	
//rmsp=std(xh);
	double rmsp=M.std(xh);
	
//%xh=x;  % mod for ICSLP 2002

//[b,a]=butter(6,300/fs*2,'high');  % 08/Sept./1999
// XXX de pus filtrul si aici
	
	if (fs == 8000) {
		b = new double[] {  0.6335f, -3.8009f,  9.5022f, -12.6697f,  9.5022f, -3.8009f,  0.6335f };
		a = new double[] {  1.0000f, -5.0900f, 10.8550f, -12.4084f,  8.0152f, -2.7730f,  0.4013f };
	}
	if (fs == 16000) {
		b = new double[] {  0.7963f, -4.7779f, 11.9448f, -15.9263f, 11.9448f, -4.7779f,  0.7963f };
		a = new double[] { 1.0000f, -5.5449f, 12.8268f, -15.8439f, 11.0213f, -4.0933f,  0.6341f };
	}
	if (fs == 24000) {
		b = new double[] { 0.8592f, -5.1551f, 12.8877f, -17.1836f, 12.8877f, -5.1551f,  0.8592f };
		a = new double[] {  1.0000f, -5.6966f, 13.5285f, -17.1441f, 12.2271f, -4.6531f,  0.7382f };
	}
	if (fs == 48000) {
		b = new double[] { 0.9269f, -5.5616f, 13.9041f, -18.5387f, 13.9041f, -5.5616f,  0.9269f };
		a = new double[] {  1.0000f, -5.8483f, 14.2528f, -18.5281f, 13.5499f, -5.2856f,  0.8592f };
	}
	if (fs == 96000) {
		b = new double[] { 0.9628f, -5.7767f, 14.4417f, -19.2556f, 14.4417f, -5.7767f,  0.9628f };
		a = new double[] { 1.0000f, -5.9241f, 14.6236f, -19.2528f, 14.2584f, -5.6320f,  0.9269f };
	}
	
	
//xh2=filter(b,a,x);
	
	double[] xh2 = M.filter(b,a,x);
    d.clear(); d = new ArrayList<Float>();
     data = IO.readAllLines("D:\\Dropbox\\STRAIGHTsrc\\xh2");
    for(String n : data) {
    	d.add(Float.valueOf(n));
    }
    xh2 = new double[d.size()];
    iiq = 0;
    for (Float f : d) {
    	xh2[iiq++] = (f != null ? f : Float.NaN); // Or whatever default you want.
    }
/// YYY de scos, nu cade cum trebuie filter pe valori

	
//%	High-pass filter using 3000Hz cut-off butterworth filter
//[b,a]=butter(6,3000/fs*2,'high');    
    if (fs == 8000) {
		b = new double[] {  0.0011f, -0.0063f,  0.0158f, -0.0210f,  0.0158f, -0.0063f,  0.0011f };
		a = new double[] {  1.0000f,  2.9785f,  4.1361f,  3.2598f,  1.5173f,  0.3911f,  0.0434f };
	}
	if (fs == 16000) {
		b = new double[] {  0.0852f, -0.5114f,  1.2785f, -1.7047f,  1.2785f, -0.5114f,  0.0852f };
		a = new double[] {  1.0000f, -1.4851f,  1.6036f, -0.9241f,  0.3592f, -0.0756f,  0.0073f };
	}
	if (fs == 24000) {
		b = new double[] {  0.2082f, -1.2493f,  3.1233f, -4.1644f,  3.1233f, -1.2493f,  0.2082f };
		a = new double[] {  1.0000f, -2.9785f,  4.1361f, -3.2598f,  1.5173f, -0.3911f,  0.0434f };
	}
	if (fs == 48000) {
		b = new double[] {  0.4654f, -2.7923f,  6.9808f, -9.3077f,  6.9808f, -2.7923f,  0.4654f };
		a = new double[] {  1.0000f, -4.4846f,  8.5290f, -8.7791f,  5.1476f, -1.6277f,  0.2166f };
	}
	if (fs == 96000) {
		b = new double[] {  0.6838f, -4.1028f, 10.2570f, -13.6760f, 10.2570f, -4.1028f,  0.6838f };
		a = new double[] {  1.0000f, -5.2416f, 11.4903f, -13.4798f,  8.9237f, -3.1601f,  0.4676f };
	}
    
//xhh=filter(b,a,x);
	
	double[] xhh = M.filter(b,a,x);
	d.clear(); d = new ArrayList<Float>();
    data = IO.readAllLines("D:\\Dropbox\\STRAIGHTsrc\\xhh");
   for(String n : data) {
   	d.add(Float.valueOf(n));
   }
   xhh = new double[d.size()];
   iiq = 0;
   for (Float f : d) {
   	xhh[iiq++] = (f != null ? f : Float.NaN); // Or whatever default you want.
   }
	/// YYY
	
//tx=[randn(1,framel/2)*rmsp/4000,xh,randn(1,framel)*rmsp/4000];
//txs=tx;
	
	double[] tx = new double[(int) (framel/2+xh.length+framel)];
	double[] txs;
	
	tx = M.randn_1row_ncolumns(tx.length);
	for(int i = 0;i<tx.length;i++) 
		tx[i] *= rmsp/4000;
	for(int i = 0;i<xh.length;i++)
		tx[i+(int)framel/2] = xh[i]; 

	// scos YYY
	d.clear(); d = new ArrayList<Float>();
    data = IO.readAllLines("D:\\Dropbox\\STRAIGHTsrc\\tx");
   for(String n : data) {
   	d.add(Float.valueOf(n));
   }
   tx = new double[tx.length];
   iiq = 0;
   for (Float f : d) {
	   tx[iiq++] = (f != null ? f : Float.NaN); // Or whatever default you want.
   }
   // pana aici YYYY sus
	
	
	txs = tx.clone();
	

//%datalength=length(tx);
//%nframe=floor((datalength-framel)/shiftl);
//nframe=min(length(f0l),round(length(x)/shiftl));
	
	int nframe = (int) Math.min(f0l.length,M.round(x.length/shiftl));
	
//nsgram=zeros(fftl/2+1,nframe);  % adaptive spectrogram
//n2sgram=zeros(fftl/2+1,nframe);
	double[][] nsgram = new double[fftl/2+1][nframe]; // java le face by default zero
	double[][] n2sgram = new double[fftl/2+1][nframe];// era +1 la fftl/2 +1
	
//t=([1:framel]-framel/2)/framel*2;
	
	double[] t= new double[framel];
	for(int i = 0;i<framel;i++) 
		t[i] = (i + 1.0 - framel/2.0)/framel*2.0;
	
//wx=(exp(-(4*t).^2/2));
	double[] wx = new double[t.length];
	double refw = 0;
	for(int i = 0;i<t.length; i++) {
		wx[i]=(double)Math.exp(-1.0*(double)Math.pow(4*t[i],2)/2.0);
		refw+=wx[i];
	}
	
	
//tt=([1:framel]-framel/2)/fs;
//refw=sum(wx);
	double[] tt= new double[framel];
	for(int i = 0;i<framel;i++) 
		tt[i] = (i + 1 - framel/2)/fs;
	


//bbase=1:fftl/2+1;
	double[] bbase = new double[fftl/2+1];
	for(int i = 0;i<bbase.length;i++)
		bbase[i] = i+1;
	// asta ar trebui scos pt ca e un counter chior

//zerobase=zeros(1,fftl);
	double[][] zerobase = new double[1][fftl];
	
//ist=1; % ii=1;
	int ist = 1;
	
//f0x=f0l*0;
	double[] f0x = new double[f0l.length]; //?? wtf

//% Optimum blending table for interference free spec.
//cfv=[1.03 0.83 0.67 0.54 0.43 0.343 0.2695];
//muv=[1    1.1  1.2  1.3  1.4  1.5   1.6];
	double[] cfv= new double[]{1.03f,0.83f,0.67f,0.54f,0.43f,0.343f,0.2695f};
	double[] muv= new double[]{1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f};
	
//bcf=spline(muv,cfv,eta);
	SplineInterpolator spline = new SplineInterpolator();
	PolynomialSplineFunction splineFunction;
	splineFunction = spline.interpolate(muv, cfv);
	double bcf = (double) splineFunction.value((double)eta);
	
//% Calculate the optimum smoothing function coefficients

//ovc=optimumsmoothing(eta,pc);
	double[] ovc = new double[4];
	ovc = M.optimumsmoothing(eta,pc);

//% Design windows for unvoiced
//fc= 300;  % default frequency resolution for unvoiced signal
//fc= 160;  % default frequency resolution for unvoiced signal (13/April/1999)
//c0=4*exp(-pi/2);
//t0=1/fc;
	
	double fc = 160.0;
	double c0=(double) (4*Math.exp(-Math.PI/2.0));
	double t0=1.0/fc;
	
//wxec=exp(-pi*(tt*fc).^2);
	double[] wxec = new double[tt.length];
	for(int i = 0;i<tt.length;i++) {
		wxec[i] = (double)Math.exp(-Math.PI*(tt[i]*fc*tt[i]*fc));		
	}

//cfc=sqrt(sum(wxec.^2));
	double cfc = 0;
	for(int i = 0; i<wxec.length; i++) {
		cfc += wxec[i]*wxec[i];
	}
	cfc = (double)Math.sqrt(cfc);
	
//wxec=wxec/cfc;
	for(int i = 0;i<wxec.length;i++) {
		wxec[i]/=cfc;
	}
	
//wxdc=bcf*wxec.*sin(pi*tt/t0);
	double[] wxdc = new double[wxec.length];
	for(int i = 0;i<wxdc.length;i++) {
		wxdc[i] = (double) (bcf*wxec[i]*Math.sin(Math.PI*tt[i]/t0));		
	}

//ttm=[0.00001 1:fftl/2 -fftl/2+1:-1]/fs;
	int fftl_pe_doi = (int) Math.ceil(fftl/2); // pentru ca java.
	double[] ttm = new double[1+2*fftl_pe_doi - 1]; // pentru ca automobile.
	ttm[0] = 0.00001f;
	for(int i = 1;i<=fftl_pe_doi;i++) 
		ttm[i] = i/fs;
	int fftl_cnter = fftl_pe_doi;
	for(int i = -fftl_pe_doi+1;i<=-1;i++) 
		ttm[++fftl_cnter] = i/fs;
	
//bfq=fftl:-1:2;
	double[] bfq = new double[fftl-2+1];
	for(int i = 0;i<bfq.length;i++)
		bfq[i] = fftl-i;
	
//ffq=2:fftl;
	double[] ffq = new double[fftl-2+1];
	for(int i = 0;i<ffq.length;i++)
		ffq[i] = 2.0+i;

//lft=1.0./(1+exp(-(abs((1:fftl)-fftl/2-1)-fftl/30)/2)); % safeguard 05/Oct./2005 by HK
	double[] lft = new double[fftl];
	for(int i = 0;i<lft.length;i++) {
		lft[i] = (double) (1.0/(1.0+(double)Math.exp( (Math.abs( (i+1)-(fftl/2)-1)-fftl/30)/-2.0)));
	}	
	

//if imgi==1; hpg=waitbar(0,'F0 adaptive time-frequency analysis.'); end;% 10/Aug./2005
//%while (ist+framel<datalength) & (ii<=min(length(f0l),nframe))
//for ii=1:nframe

	for(int ii = 0;ii<nframe-1;ii++) {
		
//%  tic;
//  if imgi==1 & rem(ii,10)==0 % 10/Aug./2005
//%    fprintf('o')
//%    if rem(ii,200)==0
//%       fprintf('\n')
//%    end;
//     waitbar(ii/nframe);
//  end;

//  f0=f0l(max(1,ii));
		f0 = f0l[0];
		for (int i = 1;i<=ii;i++)
			if(f0<f0l[i]) f0 = f0l[i];
		
//  if f0==0 
//     f0=200; %
//     f0=160; % 09/Sept./1999
//  end;
		if(f0 == 0)
			f0 = 160.0;
  
  //f0x(ii)=f0;
  //t0=1/f0;		
		f0x[ii] = f0;
		t0 = 1.0/f0;

//  wxe=exp(-pi*(tt/t0/eta).^2);
		double[] wxe = new double[tt.length];
		for(int i = 0;i<tt.length;i++) {
			wxe[i] = (double)Math.exp(-1.0 * (double)Math.PI * ((tt[i]/t0/eta) * (tt[i]/t0/eta)) );
		}
		
//  wxe=wxe/sqrt(sum(wxe.^2));
		double temp_wxe = 0;
		for(int i = 0;i<tt.length;i++) {
			temp_wxe += wxe[i]*wxe[i];
		}
		temp_wxe = (double)Math.sqrt(temp_wxe);
		for(int i = 0;i<tt.length;i++) {
			wxe[i] = wxe[i]/temp_wxe;
		}
		
//  wxd=bcf*wxe.*sin(pi*tt/t0);
		double[] wxd = new double[tt.length];
		for(int i = 0;i<tt.length;i++) {
			wxd[i] = bcf*wxe[i]*(double)Math.sin((double)Math.PI*tt[i]/t0);
		}

//  iix=round(ist:ist+framel-1);
		double[] iix = new double[framel-1];
		for(int i = 0;i<iix.length; i++) {
			iix[i] = i+ist; // wtf mai e nevoie de round???
		}
		
//  pw=sqrt(abs(fft( (tx(iix)-mean(tx(iix))).*wxe,fftl)).^2+ ...
//          abs(fft((tx(iix)-mean(tx(iix))).*wxd,fftl)).^2).^pc;
		double[] tx_iix = new double[iix.length];
		double[] tx_iix1 = new double[iix.length];
		double[] tx_iix2 = new double[iix.length];
		for(int i = 0;i<iix.length;i++) {
			tx_iix[i] = tx[(int)Math.round(iix[i]-1.0)]; // sa ma suga spiridusul pentru scaderea asta
			// aici ar trebui direct i, si scapat de vectorul iix, de studiat pt viitor pt ca matlab suge.
		}
		double mean_tx_iix = M.mean(tx_iix);
		for(int i = 0;i<iix.length;i++) {
			tx_iix1[i] = (tx_iix[i]-mean_tx_iix)*wxe[i];
			tx_iix2[i] = (tx_iix[i]-mean_tx_iix)*wxd[i];
		}
		//fft aici
		//double[] tx_fft1 = new double[tx_iix.length*2];
		//double[] tx_fft2 = new double[tx_iix.length*2];
		//tx_fft1 = M.fft(tx_iix1); // aici unde intra fftl???
		//tx_fft2 = M.fft(tx_iix2); // aici unde intra fftl???
		
		FFT fft_object = FFT.get(fftl);
		double[] fft_iix1_real = new double[fftl];
		System.arraycopy(tx_iix1, 0, fft_iix1_real, 0, tx_iix1.length);
		double[] fft_iix1_imag = new double[fftl];
		fft_object.fft(fft_iix1_real, fft_iix1_imag);
		
		double[] fft_iix2_real = new double[fftl];
		System.arraycopy(tx_iix2, 0, fft_iix2_real, 0, tx_iix2.length);
		double[] fft_iix2_imag = new double[fftl];
		fft_object.fft(fft_iix2_real, fft_iix2_imag);
		
		// abs aici
		tx_iix1 = fft_object.abs(fft_iix1_real, fft_iix1_imag);
		tx_iix2 = fft_object.abs(fft_iix2_real, fft_iix2_imag);
		
		
		double[] pw = new double[fftl];
		for(int i = 0;i<fftl;i++) {
			pw[i] = (double)Math.pow((double)Math.sqrt((tx_iix1[i]*tx_iix1[i]+tx_iix2[i]*tx_iix2[i])),pc);
		}
		
// nsgram(:,ii)=pw(bbase)';
		for(int i = 0;i<bbase.length;i++) {
			nsgram[i][ii] = pw[(int) bbase[i]-1];
		}
		
/*		
% mod for ICSLP2002  
%% Warning!! this part introduced additional low frequency noise!!
% It introduced quantization effects
%% (3/Feb./2003) found by Hideki Kawahara
%   Notes for the following 9 lines. Low to middle frequency aperiodic
%   noisy clusters were observed when a very high quality alto recording
%   was resynthesized using STRAIGHTv30kr205. The follwing lines were
%   originally designed to circumvent group delay fluctuations due to
%   unreliable estimate of low frequency power spectrum. (Situation was
%   made worse, because the lack of the DC component introduces an apparent
%   sharp spectral transition which results into very long group delay.
%   This huge group delay acts as an ampifier of group delay induced jitter
%   due to small spectral fluctuations. The basic idea was to add
%   relatively stabel DC component by mirroring the first component to the
%   DC component. Generally it worked reasonable, but due to the sloppy
%   imiplementation, it suffered from quantization in making mirror image.
%   The following codes use a linear interpolation, instead of the sloppy
%   indexing approximation. (3/Feb./2003 by H.K.)
%
%- -------------------------------------
%  f0p2=round((f0/fs*fftl)/2+1); % removed by H.K. on 3/Feb./2003
%  f0p=round((f0/fs*fftl)+1); % removed by H.K. on 3/Feb./2003
*/		
		
//  f0p2=floor((f0/fs*fftl)/2+1); % modified by H.K. on 3/Feb./2003		
		double f0p2 = (double)Math.floor((f0/fs*fftl)/2.0+1.0);
				
//  f0p=ceil((f0/fs*fftl)+1); % modified by H.K. on 3/Feb./2003
		double f0p = (double)Math.ceil((f0/fs*fftl)+1.0);
		
//  f0pr=f0/fs*fftl+1; % added by H.K. on 3/Feb./2003
		double f0pr = (double) (f0/fs*fftl+1.0);
		
//  tmppw=interp1(1:f0p,pw(1:f0p),f0pr-((1:f0p2)-1)); % added by H.K. on 3/Feb./2003
		double[] tmppw = new double[(int)f0p2];
		double[] temp_x, temp_y, temp_xd;
		temp_x = new double[(int)f0p];
		temp_y = new double[(int)f0p];
		temp_xd = new double[(int)f0p2];
		for(int i = 0;i<(int)f0p;i++) {
			temp_x[i] = i+1;
			temp_y[i] = pw[i]; // aici de verificat sa nu fie i+1			
		}
		for(int i = 0;i<(int)f0p2;i++) {
			temp_xd[i] = f0pr-i;
		}
		tmppw = M.interpLinear(temp_x, temp_y, temp_xd);
  
//  pw(1:f0p2)=tmppw; % modified by H.K. on 3/Feb./2003
		for(int i = 0;i<f0p2;i++) {
			pw[i] = tmppw[i];
		}
		
//  %%pw(1:f0p2)=pw(f0p:-1:f0p-f0p2+1); % removed by H.K. on 3/Feb./2003
//  pw(fftl:-1:fftl-f0p2+2)=pw(2:f0p2);
	int poz = 0;
	for(int i = fftl-1;i>=fftl-f0p2+1;i--) {
		pw[i] = pw[++poz];  
	}	
	
//%  local level equalization
//  ww2t=(sin(ttm/(t0/3)*pi)./(ttm/(t0/3)*pi)).^2;
	double[] ww2t = new double[ttm.length];
	for(int i = 0;i<ttm.length;i++) {
		ww2t[i] = (double)Math.pow(Math.sin(ttm[i]/(t0/3)*Math.PI)/(ttm[i]/(t0/3)*Math.PI), 2);
	}		
	
//  spw2=real(ifft(ww2t.*fft(pw).*lft));
	fft_object = FFT.get(pw.length);
	double[] fft_pw_real = new double[pw.length];
	System.arraycopy(pw, 0, fft_pw_real, 0, fft_pw_real.length);
	double[] fft_pw_imag = new double[pw.length];
	fft_object.fft(fft_pw_real, fft_pw_imag);	
	for(int i = 0;i<pw.length;i++) {
		fft_pw_real[i] = ww2t[i]*fft_pw_real[i]*lft[i];
		fft_pw_imag[i] = ww2t[i]*fft_pw_imag[i]*lft[i];
	}
	fft_object.ifft(fft_pw_real, fft_pw_imag);
	double[] spw2 = new double[fft_pw_real.length];
	System.arraycopy(fft_pw_real,0,spw2,0,fft_pw_real.length);
	
//  spw2(spw2==0)=spw2(spw2==0)+eps;  %%% safe guard added on 15/Jan./2003	
	for(int i = 0; i<spw2.length;i++)
		if(spw2[i]==0.0) spw2[i]=0.0000000000000001f;
	
//%	Optimum weighting

//  wwt=(sin(ttm/t0*pi)./(ttm/t0*pi)).^2.*(ovc(1)+ovc(2)*2*cos(ttm/t0*2*pi) ...
//     +ovc(3)*2*cos(ttm/(t0/2)*2*pi));
	
	double[] wwt = new double[ttm.length];
	for (int i = 0;i<wwt.length;i++) {
		wwt[i] = (double) (Math.sin(ttm[i]/t0*Math.PI) / (ttm[i]/t0*Math.PI));// * (double) (Math.sin(ttm[i]/t0*Math.PI)/(ttm[i]/t0*Math.PI));
		wwt[i] *= wwt[i];
		wwt[i] *= (double) (ovc[0]+ovc[1]*2*Math.cos(ttm[i]/t0*2*Math.PI)+ovc[2]*2*Math.cos(ttm[i]/(t0/2)*2*Math.PI));
	}
	
//  spw=real(ifft(wwt.*fft(pw./spw2)))/wwt(1);
	fft_object = FFT.get(pw.length);
	fft_pw_real = new double[pw.length];
	for(int i = 0;i<fft_pw_real.length;i++) {
		fft_pw_real[i] = pw[i]/spw2[i];
	}
	fft_pw_imag = new double[pw.length];
	fft_object.fft(fft_pw_real, fft_pw_imag);	
	for(int i = 0;i<pw.length;i++) {
		fft_pw_real[i] = wwt[i]*fft_pw_real[i];
		fft_pw_imag[i] = wwt[i]*fft_pw_imag[i];
	}
	fft_object.ifft(fft_pw_real, fft_pw_imag);
	double[] spw = new double[fft_pw_real.length];
	for(int i = 0;i<fft_pw_real.length;i++) {
		spw[i] = fft_pw_real[i]/wwt[0];
	}

//%   smooth half wave rectification
//  n2sgram(:,ii) = (spw2(bbase).*(0.25*(log(2*cosh(spw(bbase)*4/1.4))*1.4+spw(bbase)*4)/2))';	
	
	double[] spw_bbase = new double[bbase.length];
	System.arraycopy(spw, 0, spw_bbase, 0, spw_bbase.length);
	double[] spw2_bbase = new double[bbase.length];
	System.arraycopy(spw2, 0, spw2_bbase, 0, spw2_bbase.length);
	/*for(int i = 0;i<bbase.length;i++) {
		spw2_bbase[i] = spw2[(int)bbase[i]];
	}*/
	
	for(int i = 0;i<n2sgram.length;i++) {
		n2sgram[i][ii] = (double) (spw2_bbase[i]*(0.25*(Math.log(2.0*Math.cosh(spw_bbase[i]*4.0/1.4))*1.4+spw_bbase[i]*4.0)/2.0));
		System.out.println("["+i+"]["+ii+"] "+n2sgram[i][ii]);
	}	
	
//%  ii=ii+1;
//  ist=ist+shiftl;
	ist = ist+(int)shiftl;
	
//%  toc*1000
//end; 
	}
	
//if imgi==1; close(hpg); end; % added 06/Dec./2002% 10/Aug./2005
//f0x(ii:length(f0l))=f0+f0x(ii:length(f0l))*0;	
	for(int i = nframe;i<f0l.length;i++) {
		f0x[i] = f0;
	}
	
//if imgi==1; fprintf('\n'); end;% 10/Aug./2005	
//nsgram=nsgram.^(1/pc);
	for(int i = 0;i<nsgram.length;i++){
		for(int j = 0;j<nsgram[0].length;j++)
			nsgram[i][j] = (double)Math.pow(nsgram[i][j], 1.0/pc);
	}
	
//n2sgram=n2sgram.^(2/pc);	
	for(int i = 0;i<n2sgram.length;i++){
		for(int j = 0;j<n2sgram[0].length;j++) {
			//n2sgram[i][j] = (double)Math.pow(n2sgram[i][j], 2.0/pc);
			System.out.println("["+i+"]["+j+"] "+n2sgram[i][j]);
		}
		System.out.println("linie");
	}
	System.out.println("asd");
//[snn,smm]=size(n2sgram);
//[nii,njj]=size(n2sgram);
	int snn, nii; snn = nii = n2sgram.length;
	int smm, njj; smm = njj = n2sgram[0].length;

//lamb=0.25./(f0var+0.25);
//lambL=0.25./(f0varL+0.25);
// neutilizate deoarece f0var si f0varL sunt .. neutilizate.	

//fqx=(0:snn-1)/snn*fs/2;
	double[] fqx = new double[snn];
	for(int i = 0;i<snn;i++) {
		 fqx[i] = (double) ((double)(i/(double)snn)*fs/2.0);
	}
		
//chigh=1.0./(1+exp(-(fqx-600)/100))';
	double[] chigh = new double[snn];
	for(int i = 0;i<snn;i++) {
		chigh[i] = (double) (1.0/(1.0+Math.exp(-(fqx[i]-600.0)/100.0)));
	}
	
//clow=1.0-chigh;
	double[] clow = new double[snn];
	for(int i = 0;i<snn;i++) {
		clow[i] = (double) (1.0-chigh[i]);
	}
//tx=txs;
System.arraycopy(txs, 0, txs, 0, txs.length);


/* pentru ca e zero smoothinid 
if smoothinid
  nssgram=n2sgram;

  lowestf0=40;

  tunitw=ceil(1.1*(1000/shiftm)/lowestf0);
  tx=(-tunitw:tunitw)';
  cumfreq=cumsum(f0x)/(1000/shiftm);
  t0x=(1000/shiftm)./f0x;

  bb=(1:length(tx))';
  for jj=1:min(length(f0l),njj)
    txx=cumfreq(min(njj,max(1,jj+tx)))-cumfreq(jj);
    txt=(tx+jj>0).*(tx+jj<=njj);
    idx=(abs(txx)<=1.1).*bb.*txt;
    idx=idx(idx>0);
    txx=txx(idx);
    wt=max(0,1-abs(txx));
    wt=wt.*f0x(min(njj,max(1,jj+tx(idx)+1)));
    wt=wt.*(abs(f0x(jj)-f0x(jj+tx(idx)))/f0x(jj)<0.25);
    n2sgram(:,jj)=nssgram(:,min(njj,max(1,jj+tx(idx))))*wt/sum(wt);
  end;
end;
*/

/*
%-----------------------------------------------------
%	Dirty hack for controling time constant in
%	unvoiced part analysis
%-----------------------------------------------------
*/

//if imgi==1; hpg=waitbar(0,'spline-based F0 adaptive smooting'); end;% 10/Aug./2005

//ttlv=sum(sum(n2sgram));
	double ttlv = 0;
	for (int i = 0;i<n2sgram.length;i++)
		for (int j=0;j<n2sgram[0].length;j++)
			ttlv+=n2sgram[i][j];

//ncw=round(2*fs/1000);
	double ncw = Math.round(2.0*fs/1000.0);

//lbb=round(300/fs*fftl);  % 22/Sept./1999
	double lbb = Math.round(300.0/fs*fftl);

//%pwc=fftfilt(hanning(ncw*2+1),abs([xh, zeros(1,ncw*10)]).^2);
//h3=(conv(hanning(round(fs/1000)),exp(-1400/fs*(0:ncw*2))));  % 30/July/1999
	double[] h3;
	double[] temp_h = M.hanning((int) Math.round(fs/1000));
	double[] temp_exp = new double[(int) (ncw*2)+1];
	for(int i=0;i<temp_exp.length;i++) {
		temp_exp[i] = (double) Math.exp(-1400.0/fs*i);
	}
	h3 = M.multiplyPolynomials(temp_h, temp_exp);
	
//pwc=fftfilt(h3,abs([xh2, zeros(1,ncw*10)]).^2); % 30/July/1999,   % 08/Sept./1999
	double[] temp_abs = new double[(int) (ncw*10+xh2.length)];
	for(int i = 0;i<xh2.length;i++)
	temp_abs[i] = xh2[i]*xh2[i];
	
	double[] pwc = M.fftfilt(h3, temp_abs); 
	
//if imgi==1; waitbar(0.1); end; % 08/Dec./2002% 10/Aug./2005
//pwc=pwc(round(1:fs/(1000/shiftm):length(pwc)));
	double[] temp_pwc = new double[(int) (pwc.length/(fs/(1000/shiftm)))+1];	
	for(int i = 0;i<temp_pwc.length;i++) {
		temp_pwc[i] = i*fs/(1000/shiftm);
	}
	double[] temp2_pwc = new double[temp_pwc.length];
	for(int i = 0;i<temp_pwc.length;i++) {
		temp2_pwc[i] = pwc[(int) temp_pwc[i]];
	}
	pwc = new double[temp2_pwc.length];
	System.arraycopy(temp2_pwc, 0, pwc, 0, n2sgram[0].length);
	
//[nn,mm]=size(n2sgram);
	int nn = n2sgram.length;
	int mm = n2sgram[0].length;
	
//pwc=pwc(1:mm);
	// scris la array copy 3 linii mai sus
	
//pwc=pwc/sum(pwc)*sum(sum(n2sgram(lbb:nn,:)));
	double temp_sum = 0;	
	for (int i = (int) lbb;i<n2sgram.length;i++) {
		for (int j=0;j<n2sgram[0].length;j++){
			temp_sum+=n2sgram[i][j];
			System.out.println ("["+j+"] :"+n2sgram[i][j]);
		}
		System.out.println(">>> "+temp_sum);
	}
	double temp_pwc_sum = 0;
	for(int i = 0; i<pwc.length; i++) {
		temp_pwc_sum+=pwc[i];
	}		
	for(int i = 0; i<pwc.length; i++) {
		pwc[i] = pwc[i]/temp_pwc_sum*temp_sum;
	}
	
//if imgi==1; waitbar(0.2); end; % 08/Dec./2002% 10/Aug./2005

//%pwch=fftfilt(hanning(ncw*2+1),abs([xhh, zeros(1,ncw*10)]).^2);
//pwch=fftfilt(h3,abs([xhh, zeros(1,ncw*10)]).^2);% 30/July/1999
	temp_exp = new double[(int) (ncw*10)+xhh.length];
	for(int i=0;i<xhh.length;i++) {
		temp_exp[i] = xhh[i]*xhh[i];
	}
	//h3 = M.fftfilt(xhh, temp_exp);
	
	
	
	
//if imgi==1; waitbar(0.3); end; % 08/Dec./2002% 10/Aug./2005
//pwch=pwch(round(1:fs/(1000/shiftm):length(pwch)));
/*	double[] pwch = new double[];
[nn,mm]=size(n2sgram);
pwch=pwch(1:mm);
pwch=pwch/sum(pwch)*ttlv;

ipwm=7;	% impact detection window size
ipl=round(ipwm/shiftm);
ww=hanning(ipl*2+1);
ww=ww/sum(ww);
apwt=fftfilt(ww,[pwch(:)' zeros(1,length(ww)*2)]);
apwt=apwt((1:length(pwch))+ipl);
dpwt=fftfilt(ww,[diff(pwch(:)').^2 zeros(1,length(ww)*2)]);
dpwt=dpwt((1:length(pwch))+ipl);
mmaa=max(apwt);
apwt(apwt<=0)=apwt(apwt<=0)*0+mmaa;  % bug fix 03/Sept./1999
rr=(sqrt(dpwt)./apwt);
lmbd=(1.0./(1+exp(-(sqrt(rr)-0.75)*20)));

pwc=pwc.*lmbd+(1-lmbd).*sum(n2sgram);  %  time constant controller

%	Shaping amplitude envelope

for ii=1:mm
	if f0raw(ii)==0
		n2sgram(:,ii)=pwc(ii)*n2sgram(:,ii)/sum(n2sgram(:,ii));
	end;
	if imgi==1 & rem(ii,10)==0% 10/Aug./2005
		waitbar(0.4+0.5*ii/mm); % 08/Dec./2002
	end;
end;

n2sgram=abs(n2sgram+0.0000000001);
n2sgram=sqrt(n2sgram);
if imgi==1; waitbar(1); end; % 08/Dec./2002% 10/Aug./2005
if imgi==1; fprintf('\n'); end;
if imgi==1; close(hpg); end;   
        
        
    */        
    
	
	return 0;
    }
}
