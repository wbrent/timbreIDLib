/*

specHarmonicity~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* specHarmonicity_tilde_class;

typedef struct _specHarmonicity_tilde
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_bool x_powerSpectrum;
    t_bool x_inputFund;
    t_float x_fundFreq;
    t_float x_minFund;
    t_float x_maxFund;
    t_float x_threshPct;
    t_uShortInt x_maxPeaks;
    double x_lastDspTime;
    t_sample* x_signalBuffer;
    t_sample* x_fftwIn;
    fftwf_complex* x_fftwOut;
    fftwf_plan x_fftwPlan;
    t_float* x_blackman;
    t_float* x_cosine;
    t_float* x_hamming;
    t_float* x_hann;
    t_outlet* x_harm;
    t_outlet* x_inHarm;
    t_float x_f;
} t_specHarmonicity_tilde;


/* ------------------------ specHarmonicity~ -------------------------------- */

static void specHarmonicity_tilde_bang (t_specHarmonicity_tilde* x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float* windowFuncPtr, *flagsBuf, *peakFreqs, *peakAmps, minPeakVal, maxPeakVal, thresh, fund, harmSpacing, halfHarmSpacing, harm, inHarm, harmDividend, inHarmDividend, divisor;
    t_uShortInt numPeaks;
    double currentTime;

    window = x->x_window;
    windowHalf = x->x_windowHalf;

    currentTime = clock_gettimesince (x->x_lastDspTime);
    bangSample = roundf ((currentTime / 1000.0) * x->x_sr);

    if (bangSample >= x->x_n)
        bangSample = x->x_n - 1;

    // construct analysis window using bangSample as the end of the window
    for (i = 0, j = bangSample; i < window; i++, j++)
        x->x_fftwIn[i] = x->x_signalBuffer[j];

    windowFuncPtr = x->x_blackman;

    switch (x->x_windowFunction)
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
    if (x->x_windowFunction != rectangular)
        for (i = 0; i < window; i++, windowFuncPtr++)
            x->x_fftwIn[i] *= *windowFuncPtr;

    fftwf_execute (x->x_fftwPlan);

    // put the result of power calc back in x_fftwIn
    tIDLib_power (windowHalf + 1, x->x_fftwOut, x->x_fftwIn);

    if ( !x->x_powerSpectrum)
        tIDLib_mag (windowHalf + 1, x->x_fftwIn);

    // find all significant peaks
    flagsBuf = (t_float *)t_getbytes ((windowHalf + 1) * sizeof (t_float));
    minPeakVal = FLT_MAX;
    maxPeakVal = -FLT_MAX;
    numPeaks = 0;

    tIDLib_peaksValleys (windowHalf + 1, x->x_fftwIn, flagsBuf, &minPeakVal, &maxPeakVal);

    thresh = maxPeakVal * (x->x_threshPct/100.0);
    peakFreqs = (t_float *)t_getbytes (0);
    peakAmps = (t_float *)t_getbytes (0);

    for (i = 0; i <= windowHalf; i++)
    {
        // 0.5 in the flagsBuf means a half peak, which we'll ignore
        if (flagsBuf[i]>0.5)
        {
            t_float thisAmp;

            thisAmp = x->x_fftwIn[i];

            if (thisAmp>=thresh)
            {
                peakFreqs = (t_float *)t_resizebytes (peakFreqs, numPeaks * sizeof (t_float), (numPeaks+1) * sizeof (t_float));
                peakAmps = (t_float *)t_resizebytes (peakAmps, numPeaks * sizeof (t_float), (numPeaks+1) * sizeof (t_float));

                peakAmps[numPeaks] = thisAmp;
                peakFreqs[numPeaks] = tIDLib_bin2freq (i, window, x->x_sr);
                numPeaks++;

                if (numPeaks>=x->x_maxPeaks)
                    break;
            }
        }
    }

    t_freebytes (flagsBuf, (windowHalf + 1) * sizeof (t_float));

    harm = inHarm = harmDividend = inHarmDividend = divisor = 0.0;

    if (x->x_inputFund)
        fund = x->x_fundFreq;
    else
    {
        if (peakFreqs[0]==0.0)
            fund = peakFreqs[1];
        else
            fund = peakFreqs[0];
    }

    if (fund<x->x_minFund || fund>x->x_maxFund)
    {
        harmDividend = -numPeaks; // to make harm value  -1.0
        goto earlyExit;
    }

    harmSpacing = fund;
    halfHarmSpacing = harmSpacing * 0.5;

    for (i = 0; i < numPeaks; i++)
    {
        t_float thisAmp;

        thisAmp = peakAmps[i];

        if (thisAmp>0.0)
        {
            t_float thisFreq, deviation;
            t_uShortInt roundedHarm;

            thisFreq = peakFreqs[i];
            roundedHarm = roundf (thisFreq/fund);
            deviation = fabs (thisFreq - (roundedHarm*fund));
            inHarmDividend += deviation * thisAmp;
            divisor += thisAmp;

            harmDividend += (halfHarmSpacing-deviation)/halfHarmSpacing;
        }
    }

    earlyExit:

    t_freebytes (peakAmps, numPeaks * sizeof (t_float));
    t_freebytes (peakFreqs, numPeaks * sizeof (t_float));

    if (divisor <= 0.0 || fund <= 0.0)
        inHarm = -1.0;
    else
        inHarm = (2*inHarmDividend)/(divisor*fund);

    if (numPeaks <= 0)
        harm = -1.0;
    else
        harm = harmDividend/numPeaks;

    outlet_float (x->x_inHarm, inHarm);
    outlet_float (x->x_harm, harm);
}


static void specHarmonicity_tilde_print (t_specHarmonicity_tilde* x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_sampIdx)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post ("%s use input fundamental: %i", x->x_objSymbol->s_name, x->x_inputFund);
    post ("%s current fundamental frequency: %f", x->x_objSymbol->s_name, x->x_fundFreq);
    post ("%s minimum fundamental frequency: %f", x->x_objSymbol->s_name, x->x_minFund);
    post ("%s maximum fundamental frequency: %f", x->x_objSymbol->s_name, x->x_maxFund);
    post ("%s spectral peak threshold: %f", x->x_objSymbol->s_name, x->x_threshPct);
    post ("%s maximum spectral peaks to consider: %i", x->x_objSymbol->s_name, x->x_maxPeaks);
}


static void specHarmonicity_tilde_window (t_specHarmonicity_tilde* x, t_floatarg w)
{
    t_sampIdx i, window, windowHalf;

    window = w;

    if (window < TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
    }

    windowHalf = window * 0.5;

    x->x_signalBuffer = (t_sample *)t_resizebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample), (window + x->x_n) * sizeof (t_sample));
    x->x_fftwIn = (t_sample *)t_resizebytes (x->x_fftwIn, x->x_window * sizeof (t_sample), window * sizeof (t_sample));

    x->x_blackman = (t_float *)t_resizebytes (x->x_blackman, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_cosine = (t_float *)t_resizebytes (x->x_cosine, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_hamming = (t_float *)t_resizebytes (x->x_hamming, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_hann = (t_float *)t_resizebytes (x->x_hann, x->x_window * sizeof (t_float), window * sizeof (t_float));

    x->x_window = window;
    x->x_windowHalf = windowHalf;

    // free the FFTW output buffer, and re-malloc according to new window
    fftwf_free (x->x_fftwOut);

    // destroy old plan, which depended on x->x_window
    fftwf_destroy_plan (x->x_fftwPlan);

    // allocate new fftwf_complex memory for the plan based on new window size
    x->x_fftwOut = (fftwf_complex *) fftwf_alloc_complex (windowHalf + 1);

    // create a new DFT plan based on new window size
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_window; i++)
        x->x_fftwIn[i] = 0.0;

    // initialize signal buffer
    for (i = 0; i < x->x_window + x->x_n; i++)
        x->x_signalBuffer[i] = 0.0;

    // re-init window functions
    tIDLib_blackmanWindow (x->x_blackman, x->x_window);
    tIDLib_cosineWindow (x->x_cosine, x->x_window);
    tIDLib_hammingWindow (x->x_hamming, x->x_window);
    tIDLib_hannWindow (x->x_hann, x->x_window);

    post ("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void specHarmonicity_tilde_overlap (t_specHarmonicity_tilde* x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o<1.0)?1.0:o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void specHarmonicity_tilde_windowFunction (t_specHarmonicity_tilde* x, t_floatarg f)
{
    f = (f<0.0)?0.0:f;
    f = (f>4.0)?4.0:f;
    x->x_windowFunction = f;

    switch (x->x_windowFunction)
    {
        case rectangular:
            post ("%s window function: rectangular.", x->x_objSymbol->s_name);
            break;
        case blackman:
            post ("%s window function: blackman.", x->x_objSymbol->s_name);
            break;
        case cosine:
            post ("%s window function: cosine.", x->x_objSymbol->s_name);
            break;
        case hamming:
            post ("%s window function: hamming.", x->x_objSymbol->s_name);
            break;
        case hann:
            post ("%s window function: hann.", x->x_objSymbol->s_name);
            break;
        default:
            break;
    };
}


static void specHarmonicity_tilde_powerSpectrum (t_specHarmonicity_tilde* x, t_floatarg spec)
{
    spec = (spec<0.0)?0.0:spec;
    spec = (spec>1.0)?1.0:spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void specHarmonicity_tilde_maxPeaks (t_specHarmonicity_tilde* x, t_floatarg max)
{
    max = (max<1.0)?1.0:max;
    max = (max>(x->x_window/4.0))?(x->x_window/4.0):max;
    x->x_maxPeaks = max;

    post ("%s maximum spectral peaks to consider: %i.", x->x_objSymbol->s_name, x->x_maxPeaks);
}


static void specHarmonicity_tilde_inputFund (t_specHarmonicity_tilde* x, t_floatarg useFund)
{
    useFund = (useFund<0)?0:useFund;
    useFund = (useFund>1)?1:useFund;
    x->x_inputFund = useFund;

    if (x->x_inputFund)
        post ("%s using incoming fundamental.", x->x_objSymbol->s_name);
    else
        post ("%s using first spectral peak as fundamental.", x->x_objSymbol->s_name);
}


static void specHarmonicity_tilde_peakThresh (t_specHarmonicity_tilde* x, t_floatarg thresh)
{
    thresh = (thresh<0.0)?0.0:thresh;
    thresh = (thresh>100.0)?100.0:thresh;
    x->x_threshPct = thresh;

    post ("%s spectral peak thresh: %0.2f%% of maximum peak amplitude.", x->x_objSymbol->s_name, x->x_threshPct);
}


static void specHarmonicity_tilde_fundFreq (t_specHarmonicity_tilde* x, t_floatarg fund)
{
    if (fund <= 0.0)
        x->x_fundFreq = 0.0;
    else
        x->x_fundFreq = fund;
}


static void specHarmonicity_tilde_minFund (t_specHarmonicity_tilde* x, t_floatarg min)
{
    if (min < 0.0 || min > 20000.0)
        pd_error (x, "%s: minimum fundamental frequency must be between 0 and 20kHz.", x->x_objSymbol->s_name);
    else
        x->x_minFund = min;
}


static void specHarmonicity_tilde_maxFund (t_specHarmonicity_tilde* x, t_floatarg max)
{
    if (max < 0.0 || max > 20000.0)
        pd_error (x, "%s: maximum fundamental frequency must be between 0 and 20kHz.", x->x_objSymbol->s_name);
    else
        x->x_maxFund = max;
}


static void* specHarmonicity_tilde_new (t_symbol* s, int argc, t_atom* argv)
{
    t_specHarmonicity_tilde* x = (t_specHarmonicity_tilde *)pd_new (specHarmonicity_tilde_class);
    t_sampIdx i;

    inlet_new (&x->x_obj, &x->x_obj.ob_pd, gensym ("float"), gensym ("fund"));
    x->x_harm = outlet_new (&x->x_obj, &s_float);
    x->x_inHarm = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 1:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            break;

        default:
            post ("%s WARNING: Too many arguments supplied. Using default window size of %i.", x->x_objSymbol->s_name, TID_WINDOWSIZEDEFAULT);
            x->x_window = TID_WINDOWSIZEDEFAULT;
            break;
    }

    x->x_windowHalf = x->x_window * 0.5;
    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_windowFunction = blackman;
    x->x_powerSpectrum = false;
    x->x_inputFund = false;
    x->x_fundFreq = 100.0;
    x->x_minFund = 30.0;
    x->x_maxFund = 4000.0;
    x->x_threshPct = 5.0;
    x->x_maxPeaks = 24;
    x->x_lastDspTime = clock_getlogicaltime();

    x->x_signalBuffer = (t_sample *)t_getbytes ((x->x_window + x->x_n) * sizeof (t_sample));
    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));

     for (i = 0; i < x->x_window + x->x_n; i++)
        x->x_signalBuffer[i] = 0.0;

      x->x_blackman = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
      x->x_cosine = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
      x->x_hamming = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
      x->x_hann = (t_float *)t_getbytes (x->x_window * sizeof (t_float));

     // initialize signal windowing functions
    tIDLib_blackmanWindow (x->x_blackman, x->x_window);
    tIDLib_cosineWindow (x->x_cosine, x->x_window);
    tIDLib_hammingWindow (x->x_hamming, x->x_window);
    tIDLib_hannWindow (x->x_hann, x->x_window);

    // set up the FFTW output buffer.
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);

    // DFT plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_window; i++)
        x->x_fftwIn[i] = 0.0;

    return (x);
}


static t_int *specHarmonicity_tilde_perform (t_int *w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_specHarmonicity_tilde* x = (t_specHarmonicity_tilde *)(w[1]);

    t_sample* in = (t_float *)(w[2]);
    n = w[3];

     // shift signal buffer contents back.
    for (i = 0; i < x->x_window; i++)
        x->x_signalBuffer[i] = x->x_signalBuffer[i + n];

    // write new block to end of signal buffer.
    for (i = 0; i < n; i++)
        x->x_signalBuffer[x->x_window + i] = in[i];

    x->x_lastDspTime = clock_getlogicaltime();

    return (w + 4);
}


static void specHarmonicity_tilde_dsp (t_specHarmonicity_tilde* x, t_signal **sp)
{
    dsp_add (
        specHarmonicity_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if (sp[0]->s_sr != x->x_sr * x->x_overlap)
    {
        x->x_sr = sp[0]->s_sr / x->x_overlap;
        x->x_lastDspTime = clock_getlogicaltime();
    };

// compare n to stored n and update/resize buffer if different
    if (sp[0]->s_n != x->x_n)
    {
        t_sampIdx i;

        x->x_signalBuffer = (t_sample *)t_resizebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample), (x->x_window + sp[0]->s_n) * sizeof (t_sample));

        x->x_n = sp[0]->s_n;
        x->x_lastDspTime = clock_getlogicaltime();

        // init signal buffer
        for (i = 0; i < x->x_window + x->x_n; i++)
            x->x_signalBuffer[i] = 0.0;

        post ("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    };
};


static void specHarmonicity_tilde_free (t_specHarmonicity_tilde* x)
{
    // free the input buffer memory
    t_freebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample));

    // free FFTW stuff
    t_freebytes (x->x_fftwIn, x->x_window * sizeof (t_sample));
    fftwf_free (x->x_fftwOut);
    fftwf_destroy_plan (x->x_fftwPlan);

    // free the window memory
    t_freebytes (x->x_blackman, x->x_window * sizeof (t_float));
    t_freebytes (x->x_cosine, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hamming, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hann, x->x_window * sizeof (t_float));
}

void specHarmonicity_tilde_setup (void)
{
    specHarmonicity_tilde_class =
    class_new (
        gensym ("specHarmonicity~"),
        (t_newmethod)specHarmonicity_tilde_new,
        (t_method)specHarmonicity_tilde_free,
        sizeof (t_specHarmonicity_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)specHarmonicity_tilde_new,
        gensym ("timbreIDLib/specHarmonicity~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (specHarmonicity_tilde_class, t_specHarmonicity_tilde, x_f);

    class_addbang (specHarmonicity_tilde_class, specHarmonicity_tilde_bang);

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_maxPeaks,
        gensym ("max_peaks"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_peakThresh,
        gensym ("peak_thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_inputFund,
        gensym ("input_fund"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_fundFreq,
        gensym ("fund"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_minFund,
        gensym ("min_fund"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_maxFund,
        gensym ("max_fund"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specHarmonicity_tilde_class,
        (t_method)specHarmonicity_tilde_dsp,
        gensym ("dsp"),
        0
    );
}
