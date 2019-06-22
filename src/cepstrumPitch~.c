/*

cepstrumPitch~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *cepstrumPitch_tilde_class;

typedef struct _cepstrumPitch_tilde
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
	t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_bool x_powerSpectrum;
    t_bool x_powerCepstrum;
    t_bool x_spectrumOffset;
    t_float x_loFreq;
    t_float x_hiFreq;
    t_float x_thresh;
    double x_lastDspTime;
    t_sample *x_signalBuffer;
    t_float *x_fftwIn;
    fftwf_complex *x_fftwOut;
	fftwf_plan x_fftwForwardPlan;
	fftwf_plan x_fftwBackwardPlan;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_float x_f;
    t_outlet *x_pitch;

} t_cepstrumPitch_tilde;


/* ------------------------ cepstrumPitch~ -------------------------------- */

static void cepstrumPitch_tilde_bang(t_cepstrumPitch_tilde *x)
{
    t_sampIdx i, j, binCount, window, windowHalf, bangSample;
    t_binIdx loFreqBin, hiFreqBin, maxValIdx;
    t_float *windowFuncPtr, nRecip, maxVal, pitch, loFreqBinFloat, hiFreqBinFloat, mean, std, sum;
	double currentTime;

    window = x->x_window;
    windowHalf = x->x_windowHalf;
	
	nRecip = 1.0/window;
	
	currentTime = clock_gettimesince(x->x_lastDspTime);
	bangSample = roundf((currentTime/1000.0)*x->x_sr);
	
	loFreqBinFloat = roundf(x->x_sr/x->x_loFreq);
	hiFreqBinFloat = roundf(x->x_sr/x->x_hiFreq);

	hiFreqBin = (hiFreqBinFloat<0)?0:hiFreqBinFloat;
	hiFreqBin = (hiFreqBin>windowHalf)?windowHalf:hiFreqBin;

	loFreqBin = (loFreqBinFloat<0)?0:loFreqBinFloat;
	loFreqBin = (loFreqBin>windowHalf)?windowHalf:loFreqBin;
		
	if(bangSample >= x->x_n)
        bangSample = x->x_n-1;
            
	// construct analysis window using bangSample as the end of the window
	for(i=0, j=bangSample; i<window; i++, j++)
		x->x_fftwIn[i] = x->x_signalBuffer[j];

	switch(x->x_windowFunction)
	{
		case rectangular:
			break;
		case blackman:
			windowFuncPtr = x->x_blackman;
			break;
		case cosine:
			windowFuncPtr = x->x_cosine;
			break;
		case hamming:
			windowFuncPtr = x->x_hamming;
			break;
		case hann:
			windowFuncPtr = x->x_hann;
			break;
		default:
			windowFuncPtr = x->x_blackman;
			break;
	};

	// if windowFunction == 0, skip the windowing (rectangular)
	if(x->x_windowFunction!=rectangular)
		for(i=0; i<window; i++, windowFuncPtr++)
			x->x_fftwIn[i] *= *windowFuncPtr;

	fftwf_execute(x->x_fftwForwardPlan);

	// put the result of power calc back in x_fftwIn
	tIDLib_power(windowHalf+1, x->x_fftwOut, x->x_fftwIn);

	if(!x->x_powerSpectrum)
		tIDLib_mag(windowHalf+1, x->x_fftwIn);

	// add 1.0 to power or magnitude spectrum before taking the log and then IFT. Avoid large negative values from log(negativeNum). MPM (McCleod Pitch Method)
	if(x->x_spectrumOffset)
		for(i=0; i<windowHalf+1; i++)
			x->x_fftwIn[i] += 1.0;
 
	tIDLib_log(windowHalf+1, x->x_fftwIn);

	// copy forward DFT magnitude result into real part of backward DFT complex input buffer, and zero out the imaginary part. fftwOut is only N/2+1 points long, while fftwIn is N points long
	for(i=0; i<windowHalf+1; i++)
	{
		x->x_fftwOut[i][0] = x->x_fftwIn[i];
		x->x_fftwOut[i][1] = 0.0;
	}

	fftwf_execute(x->x_fftwBackwardPlan);

	for(i=0; i<windowHalf+1; i++)
		x->x_fftwIn[i] *= nRecip;

	// optionally square the cepstrum results for power cepstrum
	if(x->x_powerCepstrum)
		for(i=0; i<windowHalf+1; i++)
			x->x_fftwIn[i] = x->x_fftwIn[i]*x->x_fftwIn[i];

	maxVal = 0;
	maxValIdx = 0;
	
	sum=mean=std=0.0;
	binCount=0;

	// traverse from hiFreq to loFreq because the high frequency cepstrum bin is lower than the low frequency cepstrum bin
	for(i=hiFreqBin; i<=loFreqBin; i++, binCount++)
	{
		// check that loFreqBin doesn't go above Nyquist bin
		if(i>=windowHalf)
			break;
		
		// accumulate a sum to get the mean below
		sum += x->x_fftwIn[i];

		if(x->x_fftwIn[i]>maxVal)
		{
			maxVal = x->x_fftwIn[i];
			maxValIdx = i;
		}	
	}

	mean = sum/binCount;
	sum = 0.0;
	binCount = 0;

	// center & square the data
	for(i=hiFreqBin; i<=loFreqBin; i++, binCount++)
	{
		x->x_fftwIn[i] -= mean;
		x->x_fftwIn[i] *= x->x_fftwIn[i];
		sum += x->x_fftwIn[i];
	}

	// get standard deviation
	std = sum/(binCount-1);
	std = sqrt(std);

	// see if maxVal is above the mean by more than x_thresh standard deviations
	if( fabs(maxVal-mean) > (x->x_thresh*std) )
		pitch = ftom(x->x_sr/((t_float)maxValIdx));
	else
		pitch = -1500.0;

 	outlet_float(x->x_pitch, pitch);
}


static void cepstrumPitch_tilde_print(t_cepstrumPitch_tilde *x)
{
	post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr/x->x_overlap));
	post("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
	post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
	post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
	post("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
	post("%s power cepstrum: %i", x->x_objSymbol->s_name, x->x_powerCepstrum);
	post("%s spectrum offset: %i", x->x_objSymbol->s_name, x->x_spectrumOffset);
	post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
	post("%s pitch range: %0.2f to %0.2f", x->x_objSymbol->s_name, ftom(x->x_loFreq), ftom(x->x_hiFreq));
	post("%s peak threshold: %0.4f standard deviations above mean", x->x_objSymbol->s_name, x->x_thresh);
}


static void cepstrumPitch_tilde_window(t_cepstrumPitch_tilde *x, t_floatarg w)
{
	t_sampIdx i, window, windowHalf;
	
	window = w;
	
	if(window<MINWINDOWSIZE)
	{
		window = WINDOWSIZEDEFAULT;
		post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
	}
	
	windowHalf = window*0.5;
	
	x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window+x->x_n) * sizeof(t_sample), (window+x->x_n) * sizeof(t_sample));
	x->x_fftwIn = (t_float *)t_resizebytes(x->x_fftwIn, x->x_window * sizeof(t_float), window * sizeof(t_float));

	x->x_blackman = (t_float *)t_resizebytes(x->x_blackman, x->x_window*sizeof(t_float), window*sizeof(t_float));
	x->x_cosine = (t_float *)t_resizebytes(x->x_cosine, x->x_window*sizeof(t_float), window*sizeof(t_float));
	x->x_hamming = (t_float *)t_resizebytes(x->x_hamming, x->x_window*sizeof(t_float), window*sizeof(t_float));
	x->x_hann = (t_float *)t_resizebytes(x->x_hann, x->x_window*sizeof(t_float), window*sizeof(t_float));

	x->x_window = window;
	x->x_windowHalf = windowHalf;

	// free the FFTW output buffer, and re-malloc according to new window
	fftwf_free(x->x_fftwOut);
	
	// destroy old plan, which depended on x->x_window
	fftwf_destroy_plan(x->x_fftwForwardPlan); 
	fftwf_destroy_plan(x->x_fftwBackwardPlan); 
	
	// allocate new fftwf_complex memory for the plan based on new window size
	x->x_fftwOut = (fftwf_complex *) fftwf_alloc_complex(windowHalf+1);

	// create a new DFT plan based on new window size
	x->x_fftwForwardPlan = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);
	x->x_fftwBackwardPlan = fftwf_plan_dft_c2r_1d(x->x_window, x->x_fftwOut, x->x_fftwIn, FFTWPLANNERFLAG);

	// we're supposed to initialize the input array after we create the plan
 	for(i=0; i<x->x_window; i++)
		x->x_fftwIn[i] = 0.0;
		
	// initialize signal buffer
	for(i=0; i<x->x_window+x->x_n; i++)
		x->x_signalBuffer[i] = 0.0;
		
	// re-init window functions
	tIDLib_blackmanWindow(x->x_blackman, x->x_window);
	tIDLib_cosineWindow(x->x_cosine, x->x_window);
	tIDLib_hammingWindow(x->x_hamming, x->x_window);
	tIDLib_hannWindow(x->x_hann, x->x_window);
 
	post("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void cepstrumPitch_tilde_overlap(t_cepstrumPitch_tilde *x, t_floatarg o)
{
	// this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr/x->x_overlap;
	x->x_overlap = (o<1)?1:o;

    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void cepstrumPitch_tilde_windowFunction(t_cepstrumPitch_tilde *x, t_floatarg f)
{
    f = (f<0)?0:f;
    f = (f>4)?4:f;
	x->x_windowFunction = f;

	switch(x->x_windowFunction)
	{
		case rectangular:
			post("%s window function: rectangular.", x->x_objSymbol->s_name);
			break;
		case blackman:
			post("%s window function: blackman.", x->x_objSymbol->s_name);
			break;
		case cosine:
			post("%s window function: cosine.", x->x_objSymbol->s_name);
			break;
		case hamming:
			post("%s window function: hamming.", x->x_objSymbol->s_name);
			break;
		case hann:
			post("%s window function: hann.", x->x_objSymbol->s_name);
			break;
		default:
			break;
	};
}


// magnitude spectrum == 0, power spectrum == 1
static void cepstrumPitch_tilde_powerSpectrum(t_cepstrumPitch_tilde *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
	x->x_powerSpectrum = spec;

	if(x->x_powerSpectrum)
		post("%s using power spectrum", x->x_objSymbol->s_name);
	else
		post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void cepstrumPitch_tilde_powerCepstrum(t_cepstrumPitch_tilde *x, t_floatarg power)
{
    power = (power<0)?0:power;
    power = (power>1)?1:power;
	x->x_powerCepstrum = power;

	if(x->x_powerCepstrum)
		post("%s using power cepstrum", x->x_objSymbol->s_name);
	else
		post("%s using magnitude cepstrum", x->x_objSymbol->s_name);
}


static void cepstrumPitch_tilde_spectrumOffset(t_cepstrumPitch_tilde *x, t_floatarg offset)
{
    offset = (offset<0)?0:offset;
    offset = (offset>1)?1:offset;
	x->x_spectrumOffset = offset;

	if(x->x_spectrumOffset)
		post("%s spectrum offset ON", x->x_objSymbol->s_name);
	else
		post("%s spectrum offset OFF", x->x_objSymbol->s_name);
}


static void cepstrumPitch_tilde_pitchRange(t_cepstrumPitch_tilde *x, t_floatarg low, t_floatarg hi)
{
    low = (low<0)?0:low;
    low = (low>20000)?20000:low;

    hi = (hi<0)?0:hi;
    hi = (hi>20000)?20000:hi;
    
    if(low>hi)
    {
    	t_float tmp;
    	tmp = hi;
    	hi = low;
    	low = tmp;
    }

	x->x_loFreq = mtof(low);
	x->x_hiFreq = mtof(hi);
	
	post("%s pitch range: %0.2f to %0.2f", x->x_objSymbol->s_name, ftom(x->x_loFreq), ftom(x->x_hiFreq));
}


static void cepstrumPitch_tilde_threshold(t_cepstrumPitch_tilde *x, t_floatarg thresh)
{
	x->x_thresh = (thresh<0)?0:thresh;

	post("%s peak threshold: %0.4f standard deviations above mean", x->x_objSymbol->s_name, x->x_thresh);
}


static void *cepstrumPitch_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_cepstrumPitch_tilde *x = (t_cepstrumPitch_tilde *)pd_new(cepstrumPitch_tilde_class);
	t_sampIdx i;
	
	x->x_pitch = outlet_new(&x->x_obj, &s_float);

	// store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
	x->x_objSymbol = s;
	
	switch(argc)
	{
		case 3:
			x->x_window = atom_getfloat(argv);
			if(x->x_window<MINWINDOWSIZE)
			{
				x->x_window = WINDOWSIZEDEFAULT;
				post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
			}
			x->x_loFreq = atom_getfloat(argv+1);
			x->x_hiFreq = atom_getfloat(argv+2);
			x->x_loFreq = mtof(x->x_loFreq);
			x->x_hiFreq = mtof(x->x_hiFreq);
			x->x_loFreq = (x->x_loFreq<0)?0:x->x_loFreq;
			x->x_hiFreq = (x->x_hiFreq<0)?0:x->x_hiFreq;
			x->x_loFreq = (x->x_loFreq>20000)?20000:x->x_loFreq;
			x->x_hiFreq = (x->x_hiFreq>20000)?20000:x->x_hiFreq;
			break;

		case 2:
			x->x_window = atom_getfloat(argv);
			if(x->x_window<MINWINDOWSIZE)
			{
				x->x_window = WINDOWSIZEDEFAULT;
				post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
			}
			x->x_loFreq = atom_getfloat(argv+1);
			x->x_hiFreq = x->x_loFreq+12; // if no upper limit specified, make it an octave higher
			x->x_loFreq = mtof(x->x_loFreq);
			x->x_hiFreq = mtof(x->x_hiFreq);
			x->x_loFreq = (x->x_loFreq<0)?0:x->x_loFreq;
			x->x_hiFreq = (x->x_hiFreq<0)?0:x->x_hiFreq;
			x->x_loFreq = (x->x_loFreq>20000)?20000:x->x_loFreq;
			x->x_hiFreq = (x->x_hiFreq>20000)?20000:x->x_hiFreq;
			break;

		case 1:
			x->x_window = atom_getfloat(argv);
			if(x->x_window<MINWINDOWSIZE)
			{
				x->x_window = WINDOWSIZEDEFAULT;
				post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
			}
			x->x_loFreq = 50;
			x->x_hiFreq = 500;
			break;

		case 0:
			x->x_window = WINDOWSIZEDEFAULT;
			x->x_loFreq = 50;
			x->x_hiFreq = 500;
			break;
						
		default:
			post("%s WARNING: Too many arguments supplied. Using default window size of %i.", x->x_objSymbol->s_name, WINDOWSIZEDEFAULT);
			x->x_window = WINDOWSIZEDEFAULT;
			x->x_loFreq = 50;
			x->x_hiFreq = 500;
			break;
	}	

    if(x->x_loFreq>x->x_hiFreq)
    {
    	t_float tmp;
    	tmp = x->x_hiFreq;
    	x->x_hiFreq = x->x_loFreq;
    	x->x_loFreq = tmp;
    }
    
	x->x_windowHalf = x->x_window*0.5;	
	x->x_sr = SAMPLERATEDEFAULT;
	x->x_n = BLOCKSIZEDEFAULT;
	x->x_overlap = 1;
	x->x_windowFunction = rectangular;
	x->x_powerSpectrum = true;
	x->x_powerCepstrum = false;
	x->x_spectrumOffset = true;
	x->x_thresh = 0.0;
	x->x_lastDspTime = clock_getlogicaltime();
	
	x->x_signalBuffer = (t_sample *)t_getbytes((x->x_window+x->x_n) * sizeof(t_sample));
	x->x_fftwIn = (t_float *)t_getbytes(x->x_window * sizeof(t_float));

 	for(i=0; i<(x->x_window+x->x_n); i++)
		x->x_signalBuffer[i] = 0.0;

  	x->x_blackman = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_cosine = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_hamming = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_hann = (t_float *)t_getbytes(x->x_window*sizeof(t_float));

 	// initialize signal windowing functions
	tIDLib_blackmanWindow(x->x_blackman, x->x_window);
	tIDLib_cosineWindow(x->x_cosine, x->x_window);
	tIDLib_hammingWindow(x->x_hamming, x->x_window);
	tIDLib_hannWindow(x->x_hann, x->x_window);

	// set up the FFTW output buffer.
	x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

	// Forward DFT plan
	x->x_fftwForwardPlan = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

	// Backward DFT plan
	x->x_fftwBackwardPlan = fftwf_plan_dft_c2r_1d(x->x_window, x->x_fftwOut, x->x_fftwIn, FFTWPLANNERFLAG);
	
	// we're supposed to initialize the input array after we create the plan
 	for(i=0; i<x->x_window; i++)
		x->x_fftwIn[i] = 0.0;
    
    return (x);
}


static t_int *cepstrumPitch_tilde_perform(t_int *w)
{
	t_uShortInt n;
	t_sampIdx i;

    t_cepstrumPitch_tilde *x = (t_cepstrumPitch_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];
 			
 	// shift signal buffer contents back.
	for(i=0; i<x->x_window; i++)
		x->x_signalBuffer[i] = x->x_signalBuffer[i+n];
	
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->x_signalBuffer[x->x_window+i] = in[i];
		
	x->x_lastDspTime = clock_getlogicaltime();

    return (w+4);
}


static void cepstrumPitch_tilde_dsp(t_cepstrumPitch_tilde *x, t_signal **sp)
{	
	dsp_add(
		cepstrumPitch_tilde_perform,
		3,
		x,
		sp[0]->s_vec,
		sp[0]->s_n
	); 

// compare sr to stored sr and update if different
	if(sp[0]->s_sr != (x->x_sr*x->x_overlap))
	{
		x->x_sr = sp[0]->s_sr/x->x_overlap;
		x->x_lastDspTime = clock_getlogicaltime();
	};

// compare n to stored n and update/resize buffer if different
	if(sp[0]->s_n != x->x_n)
	{
		t_sampIdx i;

		x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window+x->x_n) * sizeof(t_sample), (x->x_window+sp[0]->s_n) * sizeof(t_sample));

		x->x_n = sp[0]->s_n;
		x->x_lastDspTime = clock_getlogicaltime();

		// init signal buffer
		for(i=0; i<(x->x_window+x->x_n); i++)
			x->x_signalBuffer[i] = 0.0;
			
    	post("%s block size: %i", x->x_objSymbol->s_name, x->x_n);
	};
};

static void cepstrumPitch_tilde_free(t_cepstrumPitch_tilde *x)
{
	// free the input buffer memory
    t_freebytes(x->x_signalBuffer, (x->x_window+x->x_n)*sizeof(t_sample));

	// free FFTW stuff
    t_freebytes(x->x_fftwIn, (x->x_window)*sizeof(t_float));
	fftwf_free(x->x_fftwOut);
	fftwf_destroy_plan(x->x_fftwForwardPlan); 
	fftwf_destroy_plan(x->x_fftwBackwardPlan); 
	
	// free the window memory
	t_freebytes(x->x_blackman, x->x_window*sizeof(t_float));
	t_freebytes(x->x_cosine, x->x_window*sizeof(t_float));
	t_freebytes(x->x_hamming, x->x_window*sizeof(t_float));
	t_freebytes(x->x_hann, x->x_window*sizeof(t_float));
}

void cepstrumPitch_tilde_setup(void)
{
    cepstrumPitch_tilde_class = 
    class_new(
    	gensym("cepstrumPitch~"),
    	(t_newmethod)cepstrumPitch_tilde_new,
    	(t_method)cepstrumPitch_tilde_free,
        sizeof(t_cepstrumPitch_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

	class_addcreator(
		(t_newmethod)cepstrumPitch_tilde_new,
		gensym("timbreIDLib/cepstrumPitch~"),
		A_GIMME,
		0
	);

    CLASS_MAINSIGNALIN(cepstrumPitch_tilde_class, t_cepstrumPitch_tilde, x_f);

	class_addbang(cepstrumPitch_tilde_class, cepstrumPitch_tilde_bang);

	class_addmethod(
		cepstrumPitch_tilde_class, 
        (t_method)cepstrumPitch_tilde_print,
		gensym("print"),
		0
	);

	class_addmethod(
		cepstrumPitch_tilde_class, 
        (t_method)cepstrumPitch_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		cepstrumPitch_tilde_class, 
        (t_method)cepstrumPitch_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		cepstrumPitch_tilde_class,
        (t_method)cepstrumPitch_tilde_windowFunction,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		cepstrumPitch_tilde_class, 
        (t_method)cepstrumPitch_tilde_powerSpectrum,
		gensym("power_spectrum"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		cepstrumPitch_tilde_class,
        (t_method)cepstrumPitch_tilde_powerCepstrum,
		gensym("power_cepstrum"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		cepstrumPitch_tilde_class,
        (t_method)cepstrumPitch_tilde_spectrumOffset,
		gensym("spectrum_offset"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		cepstrumPitch_tilde_class,
        (t_method)cepstrumPitch_tilde_pitchRange,
		gensym("pitch_range"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		cepstrumPitch_tilde_class,
        (t_method)cepstrumPitch_tilde_threshold,
		gensym("threshold"),
		A_DEFFLOAT,
		0
	);
	
    class_addmethod(
    	cepstrumPitch_tilde_class,
    	(t_method)cepstrumPitch_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}
