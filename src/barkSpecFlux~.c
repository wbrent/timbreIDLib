/*

barkSpecFlux~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *barkSpecFlux_tilde_class;

typedef struct _barkSpecFlux_tilde
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_bool x_normalize;
    t_bool x_powerSpectrum;
    t_bool x_logSpectrum;
    double x_lastDspTime;
    t_sample *x_signalBuffer;
    t_sample *x_fftwInForwardWindow;
    t_sample *x_fftwInBackWindow;
    fftwf_complex *x_fftwOutForwardWindow;
    fftwf_complex *x_fftwOutBackWindow;
    fftwf_plan x_fftwPlanForwardWindow;
    fftwf_plan x_fftwPlanBackWindow;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_fluxMode x_mode;
    t_bool x_squaredDiff;
    t_uInt x_separation;
    t_filterIdx x_sizeFilterFreqs;
    t_filterIdx x_numFilters;
    t_float x_barkSpacing;
    t_float *x_filterFreqs;
    t_filter *x_filterbank;
    t_bool x_specBandAvg;
    t_bool x_filterAvg;
    t_atom *x_listOut;
    t_outlet *x_fluxList;
    t_outlet *x_flux;
    t_float x_f;

} t_barkSpecFlux_tilde;


/* ------------------------ barkSpecFlux~ -------------------------------- */

static void barkSpecFlux_tilde_bang(t_barkSpecFlux_tilde *x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float flux, *windowFuncPtr;
    double currentTime;

    window = x->x_window;
    windowHalf = x->x_windowHalf;

    currentTime = clock_gettimesince(x->x_lastDspTime);
    bangSample = roundf((currentTime/1000.0)*x->x_sr);

    if(bangSample >= x->x_n)
        bangSample = x->x_n-1;

    // construct forward analysis window
    for(i=0, j=bangSample; i<window; i++, j++)
        x->x_fftwInForwardWindow[i] = x->x_signalBuffer[window + j];

    // construct back analysis window x->x_separation frames earlier
    for(i=0, j=bangSample; i<window; i++, j++)
        x->x_fftwInBackWindow[i] = x->x_signalBuffer[window - x->x_separation + j];

    windowFuncPtr = x->x_blackman;

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

    // if x_windowFunction == 0, skip the windowing (rectangular)
    if(x->x_windowFunction!=rectangular)
        for(i=0; i<window; i++, windowFuncPtr++)
            x->x_fftwInForwardWindow[i] *= *windowFuncPtr;

    fftwf_execute(x->x_fftwPlanForwardWindow);

    // put the result of power calc back in x_fftwIn
    tIDLib_power(windowHalf+1, x->x_fftwOutForwardWindow, x->x_fftwInForwardWindow);

    if(!x->x_powerSpectrum)
        tIDLib_mag(windowHalf+1, x->x_fftwInForwardWindow);

    if(x->x_specBandAvg)
        tIDLib_specFilterBands(windowHalf+1, x->x_numFilters, x->x_fftwInForwardWindow, x->x_filterbank, x->x_normalize);
    else
        tIDLib_filterbankMultiply(x->x_fftwInForwardWindow, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

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

    // if x_windowFunction == 0, skip the windowing (rectangular)
    if(x->x_windowFunction!=rectangular)
        for(i=0; i<window; i++, windowFuncPtr++)
            x->x_fftwInBackWindow[i] *= *windowFuncPtr;

    fftwf_execute(x->x_fftwPlanBackWindow);

    // put the result of power calc back in x_fftwIn
    tIDLib_power(windowHalf+1, x->x_fftwOutBackWindow, x->x_fftwInBackWindow);

    if(!x->x_powerSpectrum)
        tIDLib_mag(windowHalf+1, x->x_fftwInBackWindow);

    if(x->x_specBandAvg)
        tIDLib_specFilterBands(windowHalf+1, x->x_numFilters, x->x_fftwInBackWindow, x->x_filterbank, x->x_normalize);
    else
        tIDLib_filterbankMultiply(x->x_fftwInBackWindow, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

    flux=0.0;

    for(i=0; i<x->x_numFilters; i++)
    {
        t_float diff, val;

        if (x->x_logSpectrum)
        {
            t_float logForward, logBack;

            if (x->x_fftwInForwardWindow[i] == 0.0)
                logForward = 0.0;
            else
                logForward = log (x->x_fftwInForwardWindow[i]);

            if (x->x_fftwInBackWindow[i] == 0.0)
                logBack = 0.0;
            else
                logBack = log (x->x_fftwInBackWindow[i]);

            diff = logForward - logBack;
        }
        else
            diff = x->x_fftwInForwardWindow[i] - x->x_fftwInBackWindow[i];

        switch(x->x_mode)
        {
            case mGrowth:
                diff = (diff<0)?0:diff;
                break;
            case mDecay:
                diff = (diff>0)?0:diff;
                break;
            default:
                break;
        }

        if(x->x_squaredDiff)
            val = diff*diff;
        else
            val = fabs(diff);

        SETFLOAT(x->x_listOut+i, diff);
        flux += val;
    }

     outlet_list(x->x_fluxList, 0, x->x_numFilters, x->x_listOut);
    outlet_float(x->x_flux, flux);
}


static void barkSpecFlux_tilde_createFilterbank(t_barkSpecFlux_tilde *x, t_floatarg bs)
{
    t_filterIdx oldNumFilters;

    x->x_barkSpacing = bs;

    if(x->x_barkSpacing<TID_MINBARKSPACING || x->x_barkSpacing>TID_MAXBARKSPACING)
    {
        x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
        post("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, TID_MINBARKSPACING, TID_MAXBARKSPACING, TID_BARKSPACINGDEFAULT);
    }

    oldNumFilters = x->x_numFilters;

    x->x_sizeFilterFreqs = tIDLib_getBarkBoundFreqs(&x->x_filterFreqs, x->x_sizeFilterFreqs, x->x_barkSpacing, x->x_sr);

    x->x_numFilters = x->x_sizeFilterFreqs-2;

    tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, oldNumFilters, x->x_numFilters, x->x_window, x->x_sr);

    // resize listOut memory
    x->x_listOut = (t_atom *)t_resizebytes(x->x_listOut, oldNumFilters*sizeof(t_atom), x->x_numFilters*sizeof(t_atom));
}


static void barkSpecFlux_tilde_spec_band_avg(t_barkSpecFlux_tilde *x, t_floatarg avg)
{
    avg = (avg<0)?0:avg;
    avg = (avg>1)?1:avg;
    x->x_specBandAvg = avg;

    if(x->x_specBandAvg)
        post("%s: averaging energy in spectrum bands.", x->x_objSymbol->s_name);
    else
        post("%s: using triangular filterbank.", x->x_objSymbol->s_name);
}


static void barkSpecFlux_tilde_filter_avg(t_barkSpecFlux_tilde *x, t_floatarg avg)
{
    avg = (avg<0)?0:avg;
    avg = (avg>1)?1:avg;
    x->x_filterAvg = avg;

    if(x->x_filterAvg)
        post("%s: averaging energy in triangular filters.", x->x_objSymbol->s_name);
    else
        post("%s: summing energy in triangular filters.", x->x_objSymbol->s_name);
}


static void barkSpecFlux_tilde_print(t_barkSpecFlux_tilde *x)
{
    post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr/x->x_overlap));
    post("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
    post("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post("%s log spectrum: %i", x->x_objSymbol->s_name, x->x_logSpectrum);
    post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post("%s Bark spacing: %f", x->x_objSymbol->s_name, x->x_barkSpacing);
    post("%s number of filters: %i", x->x_objSymbol->s_name, x->x_numFilters);
    post("%s separation: %i", x->x_objSymbol->s_name, x->x_separation);
    post("%s squared difference: %i", x->x_objSymbol->s_name, x->x_squaredDiff);
    post("%s spectrum band averaging: %i", x->x_objSymbol->s_name, x->x_specBandAvg);
    post("%s triangular filter averaging: %i", x->x_objSymbol->s_name, x->x_filterAvg);

    switch(x->x_mode)
    {
        case mFlux:
            post("%s mode: flux", x->x_objSymbol->s_name);
            break;
        case mGrowth:
            post("%s mode: growth", x->x_objSymbol->s_name);
            break;
        case mDecay:
            post("%s mode: decay", x->x_objSymbol->s_name);
            break;
        default:
            post("%s mode: flux", x->x_objSymbol->s_name);
            break;
    }
}


static void barkSpecFlux_tilde_window(t_barkSpecFlux_tilde *x, t_floatarg w)
{
    t_sampIdx i, window, windowHalf;

    window = w;

    if(window<TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
    }

    if(x->x_separation > window)
    {
        post("%s change in window size caused frame separation to be less than %i samples apart. Setting frame separation to half of current window size instead.", x->x_objSymbol->s_name, window);
        x->x_separation = window*0.5;
    }

    windowHalf = window*0.5;

    x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window*2+x->x_n) * sizeof(t_sample), (window*2+x->x_n) * sizeof(t_sample));
    x->x_fftwInForwardWindow = (t_sample *)t_resizebytes(x->x_fftwInForwardWindow, x->x_window*sizeof(t_sample), window*sizeof(t_sample));
    x->x_fftwInBackWindow = (t_sample *)t_resizebytes(x->x_fftwInBackWindow, x->x_window*sizeof(t_sample), window*sizeof(t_sample));

    x->x_blackman = (t_float *)t_resizebytes(x->x_blackman, x->x_window*sizeof(t_float), window*sizeof(t_float));
    x->x_cosine = (t_float *)t_resizebytes(x->x_cosine, x->x_window*sizeof(t_float), window*sizeof(t_float));
    x->x_hamming = (t_float *)t_resizebytes(x->x_hamming, x->x_window*sizeof(t_float), window*sizeof(t_float));
    x->x_hann = (t_float *)t_resizebytes(x->x_hann, x->x_window*sizeof(t_float), window*sizeof(t_float));

    x->x_window = window;
    x->x_windowHalf = windowHalf;

    // free the FFTW output buffer, and re-malloc according to new window
    fftwf_free(x->x_fftwOutForwardWindow);
    fftwf_free(x->x_fftwOutBackWindow);

    // destroy old plan, which depended on x->x_window
    fftwf_destroy_plan(x->x_fftwPlanForwardWindow);
    fftwf_destroy_plan(x->x_fftwPlanBackWindow);

    // allocate new FFTW output buffer memory
    x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
    x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

    // create a new plan
    x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
    x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

     // we're supposed to initialize the input array after we create the plan
    for(i=0; i<x->x_window; i++)
     {
        x->x_fftwInForwardWindow[i] = 0.0;
        x->x_fftwInBackWindow[i] = 0.0;
    }

    // initialize signal buffer
    for(i=0; i<(x->x_window*2+x->x_n); i++)
        x->x_signalBuffer[i] = 0.0;

    // re-init window functions
    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);

    tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, x->x_numFilters, x->x_numFilters, x->x_window, x->x_sr);

    post("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void barkSpecFlux_tilde_overlap(t_barkSpecFlux_tilde *x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr/x->x_overlap;
    x->x_overlap = (o<1)?1:o;

    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void barkSpecFlux_tilde_x_windowFunction(t_barkSpecFlux_tilde *x, t_floatarg f)
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


static void barkSpecFlux_tilde_powerSpectrum(t_barkSpecFlux_tilde *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
    x->x_powerSpectrum = spec;

    if(x->x_powerSpectrum)
        post("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void barkSpecFlux_tilde_logSpectrum(t_barkSpecFlux_tilde *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
    x->x_logSpectrum = spec;

    if(x->x_logSpectrum)
        post("%s log spectrum enabled", x->x_objSymbol->s_name);
    else
        post("%s log spectrum disabled", x->x_objSymbol->s_name);
}


static void barkSpecFlux_tilde_mode(t_barkSpecFlux_tilde *x, t_symbol *m)
{
    if(!strcmp(m->s_name, "flux"))
        x->x_mode = mFlux;
    else if(!strcmp(m->s_name, "growth"))
        x->x_mode = mGrowth;
    else if(!strcmp(m->s_name, "decay"))
        x->x_mode = mDecay;
    else
        x->x_mode = mFlux;
}


static void barkSpecFlux_tilde_squaredDiff(t_barkSpecFlux_tilde *x, t_floatarg sd)
{
    sd = (sd<0)?0:sd;
    sd = (sd>1)?1:sd;
    x->x_squaredDiff = sd;

    if(x->x_squaredDiff)
        post("%s using squared difference", x->x_objSymbol->s_name);
    else
        post("%s using absolute value of difference", x->x_objSymbol->s_name);
}


static void barkSpecFlux_tilde_normalize(t_barkSpecFlux_tilde *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
    x->x_normalize = norm;

    if(x->x_normalize)
        post("%s spectrum normalization ON.", x->x_objSymbol->s_name);
    else
        post("%s spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void barkSpecFlux_tilde_separation(t_barkSpecFlux_tilde *x, t_floatarg s)
{
    if(s > x->x_window)
    {
        post("%s frame separation cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
        x->x_separation = x->x_windowHalf;
    }
    else if(s < 0)
    {
        post("%s frame separation must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
        x->x_separation = x->x_windowHalf;
    }
    else
        x->x_separation = s;

    post("%s frame separation: %i", x->x_objSymbol->s_name, x->x_separation);
}


static void *barkSpecFlux_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_barkSpecFlux_tilde *x = (t_barkSpecFlux_tilde *)pd_new(barkSpecFlux_tilde_class);
    t_float sepFloat;
    t_sampIdx i;

    x->x_flux = outlet_new(&x->x_obj, &s_float);
    x->x_fluxList = outlet_new(&x->x_obj, gensym("list"));

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch(argc)
    {
        case 3:
            x->x_window = atom_getfloat(argv);
            if(x->x_window<TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_barkSpacing = atom_getfloat(argv+1);
            if(x->x_barkSpacing<TID_MINBARKSPACING || x->x_barkSpacing>TID_MAXBARKSPACING)
            {
                x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
                post("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, TID_MINBARKSPACING, TID_MAXBARKSPACING, TID_BARKSPACINGDEFAULT);
            }

            sepFloat = atom_getfloat(argv+2);
            if(sepFloat > x->x_window)
            {
                post("%s WARNING: frame separation cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
                x->x_separation = x->x_window*0.5;
            }
            else if(sepFloat < 0)
            {
                post("%s WARNING: frame separation must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
                x->x_separation = x->x_window*0.5;
            }
            else
                x->x_separation = sepFloat;
            break;

        case 2:
            x->x_window = atom_getfloat(argv);
            if(x->x_window<TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_barkSpacing = atom_getfloat(argv+1);
            if(x->x_barkSpacing<TID_MINBARKSPACING || x->x_barkSpacing>TID_MAXBARKSPACING)
            {
                x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
                post("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, TID_MINBARKSPACING, TID_MAXBARKSPACING, TID_BARKSPACINGDEFAULT);
            }

            x->x_separation = x->x_window*0.5;
            break;

        case 1:
            x->x_window = atom_getfloat(argv);
            if(x->x_window<TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;

            x->x_separation = x->x_window*0.5;
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            x->x_separation = x->x_window*0.5;
            break;

        default:
            post("%s WARNING: Too many arguments supplied. Using default window size of %i, Bark spacing of %f, and frame separation of %i.", x->x_objSymbol->s_name, TID_WINDOWSIZEDEFAULT, TID_BARKSPACINGDEFAULT, TID_WINDOWSIZEDEFAULT*0.5);
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            x->x_separation = x->x_window*0.5;
            break;
    }

    x->x_windowHalf = x->x_window*0.5;
    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_windowFunction = blackman;
    x->x_normalize = false;
    x->x_powerSpectrum = false;
    x->x_logSpectrum = false;
    x->x_lastDspTime = clock_getlogicaltime();
    x->x_sizeFilterFreqs = 0;
    x->x_numFilters = 0; // this is just an init size that will be updated in createFilterbank anyway.
    x->x_squaredDiff = false; // absolute value by default
    x->x_specBandAvg = false;
    x->x_filterAvg = false;
    x->x_mode = mFlux;

    x->x_signalBuffer = (t_sample *)t_getbytes((x->x_window*2+x->x_n) * sizeof(t_sample));
    x->x_fftwInForwardWindow = (t_sample *)t_getbytes(x->x_window * sizeof(t_sample));
    x->x_fftwInBackWindow = (t_sample *)t_getbytes(x->x_window * sizeof(t_sample));

     for(i=0; i<(x->x_window*2+x->x_n); i++)
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

    // set up the FFTW output buffer. Is there no function to initialize it?
    x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
    x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

    // DFT plan
    x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
    x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for(i=0; i<x->x_window; i++)
     {
        x->x_fftwInForwardWindow[i] = 0.0;
        x->x_fftwInBackWindow[i] = 0.0;
    }

    // grab memory
    x->x_filterbank = (t_filter *)t_getbytes(0);
    x->x_filterFreqs = (t_float *)t_getbytes(0);

    x->x_sizeFilterFreqs = tIDLib_getBarkBoundFreqs(&x->x_filterFreqs, x->x_sizeFilterFreqs, x->x_barkSpacing, x->x_sr);

    // sizeFilterFreqs-2 is the correct number of filters, since we don't count the start point of the first filter, or the finish point of the last filter
    x->x_numFilters = x->x_sizeFilterFreqs-2;

    tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, 0, x->x_numFilters, x->x_window, x->x_sr);

    // create listOut memory
    x->x_listOut = (t_atom *)t_getbytes(x->x_numFilters*sizeof(t_atom));

    return (x);
}


static t_int *barkSpecFlux_tilde_perform(t_int *w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_barkSpecFlux_tilde *x = (t_barkSpecFlux_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];

     // shift signal buffer contents back.
    for(i=0; i<(x->x_window*2); i++)
        x->x_signalBuffer[i] = x->x_signalBuffer[i+n];

    // write new block to end of signal buffer.
    for(i=0; i<n; i++)
        x->x_signalBuffer[x->x_window*2+i] = in[i];

    x->x_lastDspTime = clock_getlogicaltime();

    return (w+4);
}


static void barkSpecFlux_tilde_dsp(t_barkSpecFlux_tilde *x, t_signal **sp)
{
    dsp_add(
        barkSpecFlux_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if( sp[0]->s_sr != (x->x_sr*x->x_overlap) )
    {
        x->x_sr = sp[0]->s_sr/x->x_overlap;

        tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, x->x_numFilters, x->x_numFilters, x->x_window, x->x_sr);
    };

// compare n to stored n and update/resize buffer if different
    if( sp[0]->s_n != x->x_n )
    {
        t_sampIdx i;

        x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window*2+x->x_n) * sizeof(t_sample), (x->x_window*2+sp[0]->s_n) * sizeof(t_sample));

        x->x_n = sp[0]->s_n;
        x->x_lastDspTime = clock_getlogicaltime();

        // init signal buffer
        for(i=0; i<(x->x_window*2+x->x_n); i++)
            x->x_signalBuffer[i] = 0.0;
    }
};

static void barkSpecFlux_tilde_free(t_barkSpecFlux_tilde *x)
{
    t_filterIdx i;

    // free the input buffer memory
    t_freebytes(x->x_signalBuffer, (x->x_window*2+x->x_n)*sizeof(t_sample));

    // free FFTW stuff
    t_freebytes(x->x_fftwInForwardWindow, (x->x_window)*sizeof(t_sample));
    t_freebytes(x->x_fftwInBackWindow, (x->x_window)*sizeof(t_sample));
    fftwf_free(x->x_fftwOutForwardWindow);
    fftwf_free(x->x_fftwOutBackWindow);
    fftwf_destroy_plan(x->x_fftwPlanForwardWindow);
    fftwf_destroy_plan(x->x_fftwPlanBackWindow);

    // free the window memory
    t_freebytes(x->x_blackman, x->x_window*sizeof(t_float));
    t_freebytes(x->x_cosine, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hamming, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hann, x->x_window*sizeof(t_float));

    // free filterFreqs memory
    t_freebytes(x->x_filterFreqs, x->x_sizeFilterFreqs*sizeof(t_float));

    // free the filterbank memory
    for(i=0; i<x->x_numFilters; i++)
        t_freebytes(x->x_filterbank[i].filter, x->x_filterbank[i].filterSize*sizeof(t_float));

    t_freebytes(x->x_filterbank, x->x_numFilters*sizeof(t_filter));

    t_freebytes(x->x_listOut, x->x_numFilters*sizeof(t_atom));
}

void barkSpecFlux_tilde_setup(void)
{
    barkSpecFlux_tilde_class =
    class_new(
        gensym("barkSpecFlux~"),
        (t_newmethod)barkSpecFlux_tilde_new,
        (t_method)barkSpecFlux_tilde_free,
        sizeof(t_barkSpecFlux_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)barkSpecFlux_tilde_new,
        gensym("timbreIDLib/barkSpecFlux~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN(barkSpecFlux_tilde_class, t_barkSpecFlux_tilde, x_f);

    class_addbang(barkSpecFlux_tilde_class, barkSpecFlux_tilde_bang);

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_createFilterbank,
        gensym("filterbank"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_print,
        gensym("print"),
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_separation,
        gensym("separation"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_overlap,
        gensym("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_mode,
        gensym("mode"),
        A_DEFSYMBOL,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_squaredDiff,
        gensym("squared_diff"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_normalize,
        gensym("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_x_windowFunction,
        gensym("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_powerSpectrum,
        gensym("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_logSpectrum,
        gensym("log_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_dsp,
        gensym("dsp"),
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_spec_band_avg,
        gensym("spec_band_avg"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_tilde_class,
        (t_method)barkSpecFlux_tilde_filter_avg,
        gensym("filter_avg"),
        A_DEFFLOAT,
        0
    );
}
