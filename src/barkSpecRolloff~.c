/*

barkSpecRolloff~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *barkSpecRolloff_tilde_class;

typedef struct _barkSpecRolloff_tilde
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
    double x_lastDspTime;
    t_sample *x_signalBuffer;
    t_sample *x_fftwIn;
    fftwf_complex *x_fftwOut;
    fftwf_plan x_fftwPlan;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_filterIdx x_sizeFilterFreqs;
    t_filterIdx x_numFilters;
    t_float *x_barkFreqList;
    t_float x_barkSpacing;
    t_float *x_filterFreqs;
    t_filter *x_filterbank;
    t_bool x_specBandAvg;
    t_bool x_filterAvg;
    t_float x_concentration;
    t_outlet *x_rolloff;
    t_float x_f;

} t_barkSpecRolloff_tilde;


/* ------------------------ barkSpecRolloff~ -------------------------------- */

static void barkSpecRolloff_tilde_bang (t_barkSpecRolloff_tilde *x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float energyTarget, energySum, rolloff, *windowFuncPtr;
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

    if (x->x_specBandAvg)
        tIDLib_specFilterBands (windowHalf + 1, x->x_numFilters, x->x_fftwIn, x->x_filterbank, false);
    else
        tIDLib_filterbankMultiply (x->x_fftwIn, false, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

    energyTarget = 0;
    for (i = 0; i < x->x_numFilters; i++)
        energyTarget += x->x_fftwIn[i];

    energyTarget *= x->x_concentration;

    energySum=0.0;
    i=0;

    while(energySum <= energyTarget)
    {
        energySum += x->x_fftwIn[i];
        i++;

        if (i>x->x_numFilters)
            break;
    }

    if (i>0)
        rolloff = x->x_barkFreqList[i-1]; // back up one because the last one went over...
    else
        rolloff = -1;

    outlet_float (x->x_rolloff, rolloff);
}

// STARTHERE
static void barkSpecRolloff_tilde_createFilterbank (t_barkSpecRolloff_tilde *x, t_floatarg bs)
{
    t_filterIdx i, oldNumFilters;

    x->x_barkSpacing = bs;

    if (x->x_barkSpacing < TID_MINBARKSPACING || x->x_barkSpacing > TID_MAXBARKSPACING)
    {
        x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
        post ("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, TID_MINBARKSPACING, TID_MAXBARKSPACING, TID_BARKSPACINGDEFAULT);
    }

    oldNumFilters = x->x_numFilters;

    x->x_sizeFilterFreqs = tIDLib_getBarkBoundFreqs (&x->x_filterFreqs, x->x_sizeFilterFreqs, x->x_barkSpacing, x->x_sr);

    x->x_numFilters = x->x_sizeFilterFreqs - 2;

    tIDLib_createFilterbank (x->x_filterFreqs, &x->x_filterbank, oldNumFilters, x->x_numFilters, x->x_window, x->x_sr);

    x->x_barkFreqList = (t_float *)t_resizebytes (x->x_barkFreqList, oldNumFilters * sizeof (t_float), x->x_numFilters * sizeof (t_float));

     for (i = 0; i < x->x_numFilters; i++)
        x->x_barkFreqList[i] = i*x->x_barkSpacing;
}


static void barkSpecRolloff_tilde_spec_band_avg (t_barkSpecRolloff_tilde *x, t_floatarg avg)
{
    avg = (avg < 0) ? 0 : avg;
    avg = (avg > 1) ? 1 : avg;
    x->x_specBandAvg = avg;

    if (x->x_specBandAvg)
        post ("%s: averaging energy in spectrum bands.", x->x_objSymbol->s_name);
    else
        post ("%s: using triangular filterbank.", x->x_objSymbol->s_name);
}


static void barkSpecRolloff_tilde_filter_avg (t_barkSpecRolloff_tilde *x, t_floatarg avg)
{
    avg = (avg < 0) ? 0 : avg;
    avg = (avg > 1) ? 1 : avg;
    x->x_filterAvg = avg;

    if (x->x_filterAvg)
        post ("%s: averaging energy in triangular filters.", x->x_objSymbol->s_name);
    else
        post ("%s: summing energy in triangular filters.", x->x_objSymbol->s_name);
}


static void barkSpecRolloff_tilde_print (t_barkSpecRolloff_tilde *x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post ("%s Bark spacing: %f", x->x_objSymbol->s_name, x->x_barkSpacing);
    post ("%s number of filters: %i", x->x_objSymbol->s_name, x->x_numFilters);
    post ("%s spectrum band averaging: %i", x->x_objSymbol->s_name, x->x_specBandAvg);
    post ("%s triangular filter averaging: %i", x->x_objSymbol->s_name, x->x_filterAvg);
    post ("%s concentration: %f", x->x_objSymbol->s_name, x->x_concentration);
}


static void barkSpecRolloff_tilde_concentration(t_barkSpecRolloff_tilde *x, t_floatarg c)
{
    if (c<0 || c>1.0)
        post ("%s concentration must be between 0.0 and 1.0.", x->x_objSymbol->s_name);
    else
    {
        x->x_concentration = c;
        post ("%s concentration: %f", x->x_objSymbol->s_name, x->x_concentration);
    }
}


static void barkSpecRolloff_tilde_window (t_barkSpecRolloff_tilde *x, t_floatarg w)
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

    // destroy old DFT plan, which depended on x->x_window
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

    // x_numFilters doesn't change with a change to x_window, so oldNumFilters and newNumFilters arguments are the same
    tIDLib_createFilterbank (x->x_filterFreqs, &x->x_filterbank, x->x_numFilters, x->x_numFilters, x->x_window, x->x_sr);

    post ("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void barkSpecRolloff_tilde_overlap (t_barkSpecRolloff_tilde *x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1) ? 1 : o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void barkSpecRolloff_tilde_windowFunction (t_barkSpecRolloff_tilde *x, t_floatarg f)
{
    f = (f < 0) ? 0 : f;
    f = (f > 4) ? 4 : f;
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


static void barkSpecRolloff_tilde_powerSpectrum (t_barkSpecRolloff_tilde *x, t_floatarg power)
{
    power = (power<0)?0:power;
    power = (power>1)?1:power;
    x->x_powerSpectrum = power;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void *barkSpecRolloff_tilde_new (t_symbol *s, int argc, t_atom *argv)
{
    t_barkSpecRolloff_tilde *x = (t_barkSpecRolloff_tilde *)pd_new (barkSpecRolloff_tilde_class);
    t_sampIdx i;

    x->x_rolloff = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 3:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_barkSpacing = atom_getfloat (argv + 1);
            if (x->x_barkSpacing < TID_MINBARKSPACING || x->x_barkSpacing > TID_MAXBARKSPACING)
            {
                x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
                post ("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, TID_MINBARKSPACING, TID_MAXBARKSPACING, TID_BARKSPACINGDEFAULT);
            }

            x->x_concentration = atom_getfloat (argv + 2);
            x->x_concentration = (x->x_concentration<0)?0:x->x_concentration;
            x->x_concentration = (x->x_concentration>1)?1:x->x_concentration;
            break;

        case 2:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_barkSpacing = atom_getfloat (argv + 1);
            if (x->x_barkSpacing < TID_MINBARKSPACING || x->x_barkSpacing > TID_MAXBARKSPACING)
            {
                x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
                post ("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, TID_MINBARKSPACING, TID_MAXBARKSPACING, TID_BARKSPACINGDEFAULT);
            }
            x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            break;

        case 1:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            break;

        default:
            post ("%s WARNING: Too many arguments supplied. Using default Bark spacing of %0.2f, and a concentration of %0.2f.", x->x_objSymbol->s_name, TID_BARKSPACINGDEFAULT, TID_SPECROLLOFF_DEFAULTCONCENTRATION);
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            break;
    }

    x->x_windowHalf = x->x_window * 0.5;
    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_windowFunction = blackman;
    x->x_powerSpectrum = false;
    x->x_lastDspTime = clock_getlogicaltime();
    x->x_sizeFilterFreqs = 0;
    x->x_numFilters = 0; // this is just an init size that will be updated in createFilterbank anyway.
    x->x_specBandAvg = false;
    x->x_filterAvg = false;

    x->x_signalBuffer = (t_sample *)t_getbytes ((x->x_window + x->x_n) * sizeof (t_sample));
    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));

    // initialize signal buffer
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

    // set up the FFTW output buffer
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);

    // DFT plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_window; i++)
        x->x_fftwIn[i] = 0.0;

    // grab memory
    x->x_filterbank = (t_filter *)t_getbytes (0);
    x->x_filterFreqs = (t_float *)t_getbytes (0);

    x->x_sizeFilterFreqs = tIDLib_getBarkBoundFreqs (&x->x_filterFreqs, x->x_sizeFilterFreqs, x->x_barkSpacing, x->x_sr);

    // sizeFilterFreqs - 2 is the correct number of filters, since we don't count the start point of the first filter, or the finish point of the last filter
    x->x_numFilters = x->x_sizeFilterFreqs - 2;

    tIDLib_createFilterbank (x->x_filterFreqs, &x->x_filterbank, 0, x->x_numFilters, x->x_window, x->x_sr);

    x->x_barkFreqList = (t_float *)t_getbytes (x->x_numFilters * sizeof (t_float));

     for (i = 0; i < x->x_numFilters; i++)
        x->x_barkFreqList[i] = i*x->x_barkSpacing;

    return (x);
}


static t_int *barkSpecRolloff_tilde_perform (t_int *w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_barkSpecRolloff_tilde *x = (t_barkSpecRolloff_tilde *)(w[1]);

    t_sample *in = (t_sample *)(w[2]);
    n = w[3];

     // shift signal buffer contents back.
    for (i = 0; i < x->x_window; i++)
        x->x_signalBuffer[i] = x->x_signalBuffer[i+n];

    // write new block to end of signal buffer.
    for (i = 0; i < n; i++)
        x->x_signalBuffer[x->x_window + i] = in[i];

    x->x_lastDspTime = clock_getlogicaltime();

    return (w + 4);
}


static void barkSpecRolloff_tilde_dsp (t_barkSpecRolloff_tilde *x, t_signal **sp)
{
    dsp_add (
        barkSpecRolloff_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if ( sp[0]->s_sr != x->x_sr * x->x_overlap )
    {
        x->x_sr = sp[0]->s_sr / x->x_overlap;

        tIDLib_createFilterbank (x->x_filterFreqs, &x->x_filterbank, x->x_numFilters, x->x_numFilters, x->x_window, x->x_sr);
    };

// compare n to stored n and update/resize buffer if different
    if ( sp[0]->s_n != x->x_n )
    {
        t_sampIdx i;

        x->x_signalBuffer = (t_sample *)t_resizebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample), (x->x_window + sp[0]->s_n) * sizeof (t_sample));

        x->x_n = sp[0]->s_n;
        x->x_lastDspTime = clock_getlogicaltime();

        // init signal buffer
        for (i = 0; i < x->x_window + x->x_n; i++)
            x->x_signalBuffer[i] = 0.0;
    }
};


static void barkSpecRolloff_tilde_free (t_barkSpecRolloff_tilde *x)
{
    t_filterIdx i;

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

    // free filterFreqs memory
    t_freebytes (x->x_filterFreqs, x->x_sizeFilterFreqs * sizeof (t_float));
    t_freebytes (x->x_barkFreqList, x->x_numFilters * sizeof (t_float));

    // free the filterbank memory
    for (i = 0; i < x->x_numFilters; i++)
        t_freebytes (x->x_filterbank[i].filter, x->x_filterbank[i].filterSize * sizeof (t_float));

    t_freebytes (x->x_filterbank, x->x_numFilters * sizeof (t_filter));
}


void barkSpecRolloff_tilde_setup (void)
{
    barkSpecRolloff_tilde_class =
    class_new (
        gensym ("barkSpecRolloff~"),
        (t_newmethod)barkSpecRolloff_tilde_new,
        (t_method)barkSpecRolloff_tilde_free,
        sizeof (t_barkSpecRolloff_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)barkSpecRolloff_tilde_new,
        gensym ("timbreIDLib/barkSpecRolloff~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (barkSpecRolloff_tilde_class, t_barkSpecRolloff_tilde, x_f);

    class_addbang (barkSpecRolloff_tilde_class, barkSpecRolloff_tilde_bang);

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_createFilterbank,
        gensym ("filterbank"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_concentration,
        gensym ("concentration"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_dsp,
        gensym ("dsp"),
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_spec_band_avg,
        gensym ("spec_band_avg"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        barkSpecRolloff_tilde_class,
        (t_method)barkSpecRolloff_tilde_filter_avg,
        gensym ("filter_avg"),
        A_DEFFLOAT,
        0
    );
}
