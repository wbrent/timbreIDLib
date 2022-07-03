/*

specBrightness~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* specBrightness_tilde_class;

typedef struct _specBrightness_tilde
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
    double x_lastDspTime;
    t_sample* x_signalBuffer;
    t_sample* x_fftwIn;
    fftwf_complex* x_fftwOut;
    fftwf_plan x_fftwPlan;
    t_float* x_blackman;
    t_float* x_cosine;
    t_float* x_hamming;
    t_float* x_hann;
    t_float x_freqBoundary;
    t_binIdx x_binBoundary;
    t_float* x_binFreqs;
    t_outlet* x_brightness;
    t_float x_f;
} t_specBrightness_tilde;


/* ------------------------ specBrightness~ -------------------------------- */

static void specBrightness_tilde_bang (t_specBrightness_tilde* x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float dividend, divisor, brightness;
    t_float* windowFuncPtr;
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

    dividend = divisor = brightness = 0.0;

    for (i = x->x_binBoundary; i <= windowHalf; i++)
        dividend += x->x_fftwIn[i];

    for (i = 0; i <= windowHalf; i++)
        divisor += x->x_fftwIn[i];

    if (divisor > 0.0)
        brightness = dividend / divisor;
    else
        brightness = -1.0;

    outlet_float (x->x_brightness, brightness);
}


static void specBrightness_tilde_boundary (t_specBrightness_tilde* x, t_floatarg b)
{
    if (b < 0 || b > x->x_sr * 0.5)
    {
        pd_error (x, "%s boundary frequency must be a positive real number and less than Nyquist.", x->x_objSymbol->s_name);
        return;
    }
    else
    {
        x->x_freqBoundary = b;
        x->x_binBoundary = tIDLib_nearestBinIndex (x->x_freqBoundary, x->x_binFreqs, x->x_windowHalf + 1);

        post ("%s boundary frequency: %0.2f", x->x_objSymbol->s_name, x->x_freqBoundary);
    }
}


static void specBrightness_tilde_print (t_specBrightness_tilde* x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post ("%s Frequency boundary: %0.2f", x->x_objSymbol->s_name, x->x_freqBoundary);
    post ("%s Bin boundary: %i", x->x_objSymbol->s_name, x->x_binBoundary);
}


static void specBrightness_tilde_window (t_specBrightness_tilde* x, t_floatarg w)
{
    t_sampIdx i, window, windowHalf;

    if (w < TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
    }
    else
        window = w;

    windowHalf = window * 0.5;

    x->x_signalBuffer = (t_sample *)t_resizebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample), (window + x->x_n) * sizeof (t_sample));
    x->x_fftwIn = (t_sample *)t_resizebytes (x->x_fftwIn, x->x_window * sizeof (t_sample), window * sizeof (t_sample));
    x->x_binFreqs = (t_float *)t_resizebytes (x->x_binFreqs, (x->x_windowHalf + 1) * sizeof (t_float), (windowHalf + 1) * sizeof (t_float));

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

    // freqs for each bin based on current window size and sample rate
    for (i = 0; i <= x->x_windowHalf; i++)
        x->x_binFreqs[i] = tIDLib_bin2freq (i, x->x_window, x->x_sr);

    x->x_binBoundary = tIDLib_nearestBinIndex (x->x_freqBoundary, x->x_binFreqs, x->x_windowHalf + 1);

    post ("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void specBrightness_tilde_overlap (t_specBrightness_tilde* x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1) ? 1 : o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void specBrightness_tilde_windowFunction (t_specBrightness_tilde* x, t_floatarg f)
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


static void specBrightness_tilde_powerSpectrum (t_specBrightness_tilde* x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void* specBrightness_tilde_new (t_symbol* s, int argc, t_atom* argv)
{
    t_specBrightness_tilde* x = (t_specBrightness_tilde *)pd_new (specBrightness_tilde_class);
    t_sampIdx i;

    x->x_brightness = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    // need samplerate in the switch() below, so assigning that here
    x->x_sr = TID_SAMPLERATEDEFAULT;

    switch (argc)
    {
        case 2:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_freqBoundary = atom_getfloat (argv + 1);
            if (x->x_freqBoundary < 0)
            {
                post ("%s boundary frequency must be a positive real number and less than Nyquist. Using default boundary of %0.2f Hz instead.", x->x_objSymbol->s_name, TID_SPECBRIGHTNESS_DEFAULTBOUND);
                x->x_freqBoundary = TID_SPECBRIGHTNESS_DEFAULTBOUND;
            }
            else if (x->x_freqBoundary > (x->x_sr * 0.5))
            {
                post ("%s boundary frequency must be a positive real number and less than Nyquist. Using default boundary of %0.2f Hz instead.", x->x_objSymbol->s_name, TID_SPECBRIGHTNESS_DEFAULTBOUND);
                x->x_freqBoundary = TID_SPECBRIGHTNESS_DEFAULTBOUND;
            };
            break;

        case 1:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_freqBoundary = TID_SPECBRIGHTNESS_DEFAULTBOUND;
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_freqBoundary = TID_SPECBRIGHTNESS_DEFAULTBOUND;
            break;

        default:
            post ("%s WARNING: Too many arguments supplied. Using default window size of %i and boundary of %0.2f Hz.", x->x_objSymbol->s_name, TID_WINDOWSIZEDEFAULT, TID_SPECBRIGHTNESS_DEFAULTBOUND);
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_freqBoundary = TID_SPECBRIGHTNESS_DEFAULTBOUND;
            break;
    }

    x->x_windowHalf = x->x_window * 0.5;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_windowFunction = blackman;
    x->x_powerSpectrum = false;
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

    x->x_binFreqs = (t_float *)t_getbytes ((x->x_windowHalf + 1) * sizeof (t_float));

    // freqs for each bin based on current window size and sample rate
    for (i = 0; i <= x->x_windowHalf; i++)
        x->x_binFreqs[i] = tIDLib_bin2freq (i, x->x_window, x->x_sr);

    x->x_binBoundary = tIDLib_nearestBinIndex (x->x_freqBoundary, x->x_binFreqs, x->x_windowHalf + 1);

    return (x);
}


static t_int* specBrightness_tilde_perform (t_int* w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_specBrightness_tilde* x = (t_specBrightness_tilde *)(w[1]);

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


static void specBrightness_tilde_dsp (t_specBrightness_tilde* x, t_signal** sp)
{
    dsp_add (
        specBrightness_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if (sp[0]->s_sr != x->x_sr * x->x_overlap)
    {
        t_sampIdx i;

        x->x_sr = sp[0]->s_sr / x->x_overlap;
        x->x_lastDspTime = clock_getlogicaltime();

        // freqs for each bin based on current window size and sample rate
        for (i = 0; i <= x->x_windowHalf; i++)
            x->x_binFreqs[i] = tIDLib_bin2freq (i, x->x_window, x->x_sr);

        x->x_binBoundary = tIDLib_nearestBinIndex (x->x_freqBoundary, x->x_binFreqs, x->x_windowHalf + 1);
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
    };
};


static void specBrightness_tilde_free (t_specBrightness_tilde* x)
{
    // free the input buffer memory
    t_freebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample));

    // free the binFreqs memory
    t_freebytes (x->x_binFreqs, (x->x_windowHalf + 1) * sizeof (t_float));

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

void specBrightness_tilde_setup (void)
{
    specBrightness_tilde_class =
    class_new (
        gensym ("specBrightness~"),
        (t_newmethod)specBrightness_tilde_new,
        (t_method)specBrightness_tilde_free,
        sizeof (t_specBrightness_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)specBrightness_tilde_new,
        gensym ("timbreIDLib/specBrightness~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (specBrightness_tilde_class, t_specBrightness_tilde, x_f);

    class_addbang (specBrightness_tilde_class, specBrightness_tilde_bang);

    class_addmethod (
        specBrightness_tilde_class,
        (t_method)specBrightness_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        specBrightness_tilde_class,
        (t_method)specBrightness_tilde_boundary,
        gensym ("boundary"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specBrightness_tilde_class,
        (t_method)specBrightness_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specBrightness_tilde_class,
        (t_method)specBrightness_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specBrightness_tilde_class,
        (t_method)specBrightness_tilde_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specBrightness_tilde_class,
        (t_method)specBrightness_tilde_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specBrightness_tilde_class,
        (t_method)specBrightness_tilde_dsp,
        gensym ("dsp"),
        0
    );
}
