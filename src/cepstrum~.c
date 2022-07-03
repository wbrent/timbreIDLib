/*

cepstrum~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* cepstrum_tilde_class;

typedef struct _cepstrum_tilde
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
    t_bool x_powerCepstrum;
    t_bool x_spectrumOffset;
    double x_lastDspTime;
    t_sample* x_signalBuffer;
    t_sample* x_fftwIn;
    fftwf_complex* x_fftwOut;
    fftwf_plan x_fftwForwardPlan;
    fftwf_plan x_fftwBackwardPlan;
    t_float* x_blackman;
    t_float* x_cosine;
    t_float* x_hamming;
    t_float* x_hann;
    t_atom* x_listOut;
    t_float x_f;
    t_outlet* x_ceps;

} t_cepstrum_tilde;


/* ------------------------ cepstrum~ -------------------------------- */

static void cepstrum_tilde_bang (t_cepstrum_tilde* x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float* windowFuncPtr, nRecip;
    double currentTime;

    window = x->x_window;
    windowHalf = x->x_windowHalf;

    nRecip = 1.0 / window;

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

    fftwf_execute (x->x_fftwForwardPlan);

    // put the result of power calc back in x_fftwIn
    tIDLib_power (windowHalf + 1, x->x_fftwOut, x->x_fftwIn);

    if ( !x->x_powerSpectrum)
        tIDLib_mag (windowHalf + 1, x->x_fftwIn);

    // add 1.0 to power or magnitude spectrum before taking the log and then IFT. Avoid large negative values from log (negativeNum)
    if (x->x_spectrumOffset)
        for (i = 0; i < windowHalf + 1; i++)
            x->x_fftwIn[i] += 1.0;

    tIDLib_log (windowHalf + 1, x->x_fftwIn);

    // copy forward DFT magnitude result into real part of backward DFT complex input buffer, and zero out the imaginary part. fftwOut is only N/2+1 points long, while fftwIn is N points long
    for (i = 0; i < windowHalf + 1; i++)
    {
        x->x_fftwOut[i][0] = x->x_fftwIn[i];
        x->x_fftwOut[i][1] = 0.0;
    }

    fftwf_execute (x->x_fftwBackwardPlan);

    for (i = 0; i < windowHalf + 1; i++)
        x->x_fftwIn[i] *= nRecip;

    // optionally square the cepstrum results for power cepstrum
    if (x->x_powerCepstrum)
        for (i = 0; i < windowHalf + 1; i++)
            x->x_fftwIn[i] = x->x_fftwIn[i] * x->x_fftwIn[i];

    for (i = 0; i < windowHalf + 1; i++)
        SETFLOAT (x->x_listOut + i, x->x_fftwIn[i]);

     outlet_list (x->x_ceps, 0, windowHalf + 1, x->x_listOut);
}


static void cepstrum_tilde_print (t_cepstrum_tilde* x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s power cepstrum: %i", x->x_objSymbol->s_name, x->x_powerCepstrum);
    post ("%s spectrum offset: %i", x->x_objSymbol->s_name, x->x_spectrumOffset);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
}


static void cepstrum_tilde_window (t_cepstrum_tilde* x, t_floatarg w)
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
    x->x_listOut = (t_atom *)t_resizebytes (x->x_listOut, (x->x_windowHalf + 1) * sizeof (t_atom), (windowHalf + 1) * sizeof (t_atom));

    x->x_blackman = (t_float *)t_resizebytes (x->x_blackman, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_cosine = (t_float *)t_resizebytes (x->x_cosine, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_hamming = (t_float *)t_resizebytes (x->x_hamming, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_hann = (t_float *)t_resizebytes (x->x_hann, x->x_window * sizeof (t_float), window * sizeof (t_float));

    x->x_window = window;
    x->x_windowHalf = windowHalf;

    // free the FFTW output buffer, and re-malloc according to new window
    fftwf_free (x->x_fftwOut);

    // destroy old plan, which depended on x->x_window
    fftwf_destroy_plan (x->x_fftwForwardPlan);
    fftwf_destroy_plan (x->x_fftwBackwardPlan);

    // allocate new fftwf_complex memory for the plan based on new window size
    x->x_fftwOut = (fftwf_complex *) fftwf_alloc_complex (windowHalf + 1);

    // create a new DFT plan based on new window size
    x->x_fftwForwardPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);
    x->x_fftwBackwardPlan = fftwf_plan_dft_c2r_1d (x->x_window, x->x_fftwOut, x->x_fftwIn, FFTWPLANNERFLAG);

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


static void cepstrum_tilde_overlap (t_cepstrum_tilde* x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1) ? 1 : o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void cepstrum_tilde_windowFunction (t_cepstrum_tilde* x, t_floatarg f)
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


// magnitude spectrum == 0, power spectrum == 1
static void cepstrum_tilde_powerSpectrum (t_cepstrum_tilde* x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void cepstrum_tilde_powerCepstrum (t_cepstrum_tilde* x, t_floatarg power)
{
    power = (power < 0) ? 0 : power;
    power = (power > 1) ? 1 : power;
    x->x_powerCepstrum = power;

    if (x->x_powerCepstrum)
        post ("%s reporting power cepstrum", x->x_objSymbol->s_name);
    else
        post ("%s reporting magnitude cepstrum", x->x_objSymbol->s_name);
}


static void cepstrum_tilde_spectrumOffset (t_cepstrum_tilde* x, t_floatarg offset)
{
    offset = (offset < 0) ? 0 : offset;
    offset = (offset > 1) ? 1 : offset;
    x->x_spectrumOffset = offset;

    if (x->x_spectrumOffset)
        post ("%s spectrum offset ON", x->x_objSymbol->s_name);
    else
        post ("%s spectrum offset OFF", x->x_objSymbol->s_name);
}


static void* cepstrum_tilde_new (t_symbol* s, int argc, t_atom* argv)
{
    t_cepstrum_tilde* x = (t_cepstrum_tilde *)pd_new (cepstrum_tilde_class);
    t_sampIdx i;

    x->x_ceps = outlet_new (&x->x_obj, gensym ("list"));

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
    x->x_powerCepstrum = false;
    x->x_lastDspTime = clock_getlogicaltime();

    x->x_signalBuffer = (t_sample *)t_getbytes ((x->x_window + x->x_n) * sizeof (t_sample));
    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));
    x->x_listOut = (t_atom *)t_getbytes ((x->x_windowHalf + 1) * sizeof (t_atom));

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

    // Forward DFT plan
    x->x_fftwForwardPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

    // Backward DFT plan
    x->x_fftwBackwardPlan = fftwf_plan_dft_c2r_1d (x->x_window, x->x_fftwOut, x->x_fftwIn, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
    for (i = 0; i < x->x_window; i++)
        x->x_fftwIn[i] = 0.0;

    return (x);
}


static t_int* cepstrum_tilde_perform (t_int* w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_cepstrum_tilde* x = (t_cepstrum_tilde *)(w[1]);

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


static void cepstrum_tilde_dsp (t_cepstrum_tilde* x, t_signal** sp)
{
    dsp_add (
        cepstrum_tilde_perform,
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

        post ("%s block size: %i", x->x_objSymbol->s_name, x->x_n);
    };
};

static void cepstrum_tilde_free (t_cepstrum_tilde* x)
{
    // free the input buffer memory
    t_freebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample));

    // free the list out memory
    t_freebytes (x->x_listOut, (x->x_windowHalf + 1) * sizeof (t_atom));

    // free FFTW stuff
    t_freebytes (x->x_fftwIn, x->x_window * sizeof (t_sample));
    fftwf_free (x->x_fftwOut);
    fftwf_destroy_plan (x->x_fftwForwardPlan);
    fftwf_destroy_plan (x->x_fftwBackwardPlan);

    // free the window memory
    t_freebytes (x->x_blackman, x->x_window * sizeof (t_float));
    t_freebytes (x->x_cosine, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hamming, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hann, x->x_window * sizeof (t_float));
}

void cepstrum_tilde_setup (void)
{
    cepstrum_tilde_class =
    class_new (
        gensym ("cepstrum~"),
        (t_newmethod)cepstrum_tilde_new,
        (t_method)cepstrum_tilde_free,
        sizeof (t_cepstrum_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)cepstrum_tilde_new,
        gensym ("timbreIDLib/cepstrum~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (cepstrum_tilde_class, t_cepstrum_tilde, x_f);

    class_addbang (cepstrum_tilde_class, cepstrum_tilde_bang);

    class_addmethod (
        cepstrum_tilde_class,
        (t_method)cepstrum_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        cepstrum_tilde_class,
        (t_method)cepstrum_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_tilde_class,
        (t_method)cepstrum_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_tilde_class,
        (t_method)cepstrum_tilde_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_tilde_class,
        (t_method)cepstrum_tilde_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_tilde_class,
        (t_method)cepstrum_tilde_powerCepstrum,
        gensym ("power_cepstrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_tilde_class,
        (t_method)cepstrum_tilde_spectrumOffset,
        gensym ("spectrum_offset"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_tilde_class,
        (t_method)cepstrum_tilde_dsp,
        gensym ("dsp"),
        0
    );
}
