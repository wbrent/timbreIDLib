/*

specIrregularity~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* specIrregularity_tilde_class;

typedef struct _specIrregularity_tilde
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_bool x_normalize;
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
    t_irregAlgoChoice x_algorithm;
    t_outlet* x_irregularity;
    t_float x_f;

} t_specIrregularity_tilde;


/* ------------------------ specIrregularity~ -------------------------------- */

static void specIrregularity_tilde_bang (t_specIrregularity_tilde* x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float divisor, irregularity;
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

    if (x->x_normalize)
        tIDLib_normal (windowHalf + 1, x->x_fftwIn);

    divisor = irregularity = 0.0;

    switch (x->x_algorithm)
    {
        case jensen:
            // Jensen
            for (i = 0; i <= windowHalf; i++)
            {
                if (i==windowHalf)
                    irregularity += x->x_fftwIn[i] * x->x_fftwIn[i];
                else
                    irregularity += powf (x->x_fftwIn[i] - x->x_fftwIn[i + 1], 2);

                divisor += x->x_fftwIn[i] * x->x_fftwIn[i];
            }

            if (divisor <= 0)
                irregularity = -1;
            else
                irregularity /= divisor;
            break;

        case krimphoff:
            // Krimphoff
            for (i = 1; i < windowHalf; i++)
            {
                t_float localAvg;
                localAvg = 0.0;

                for (j = 0; j < 3; j++)
                    localAvg += x->x_fftwIn[i - 1 + j];

                localAvg *= 0.333333333333;

                irregularity += fabs (x->x_fftwIn[i] - localAvg);
                //irregularity = log10(irregularity);
            }
            break;

        default:
            break;
    }

    outlet_float (x->x_irregularity, irregularity);
}


static void specIrregularity_tilde_print (t_specIrregularity_tilde* x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);

    switch (x->x_algorithm)
    {
        case jensen:
            post ("%s algorithm: Jensen", x->x_objSymbol->s_name);
            break;
        case krimphoff:
            post ("%s algorithm: Krimphoff", x->x_objSymbol->s_name);
            break;
        default:
            break;
    };
}


static void specIrregularity_tilde_algorithm (t_specIrregularity_tilde* x, t_floatarg a)
{
    a = (a < 0) ? 0 : a;
    a = (a > 1) ? 1 : a;
    x->x_algorithm = a;

    switch (x->x_algorithm)
    {
        case 0:
            post ("%s using Jensen irregularity.", x->x_objSymbol->s_name);
            break;
        case 1:
            post ("%s using Krimphoff irregularity.", x->x_objSymbol->s_name);
            break;
        default:
            break;
    };
}


static void specIrregularity_tilde_window (t_specIrregularity_tilde* x, t_floatarg w)
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

    post ("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void specIrregularity_tilde_overlap (t_specIrregularity_tilde* x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1) ? 1 : o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void specIrregularity_tilde_windowFunction (t_specIrregularity_tilde* x, t_floatarg f)
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


static void specIrregularity_tilde_normalize (t_specIrregularity_tilde* x, t_floatarg norm)
{
    norm = (norm < 0) ? 0 : norm;
    norm = (norm > 1) ? 1 : norm;
    x->x_normalize = norm;

    if (x->x_normalize)
        post ("%s spectrum normalization ON.", x->x_objSymbol->s_name);
    else
        post ("%s spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void specIrregularity_tilde_powerSpectrum (t_specIrregularity_tilde* x, t_floatarg power)
{
    power = (power < 0) ? 0 : power;
    power = (power > 1) ? 1 : power;
    x->x_powerSpectrum = power;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void* specIrregularity_tilde_new (t_symbol* s, int argc, t_atom* argv)
{
    t_specIrregularity_tilde* x = (t_specIrregularity_tilde *)pd_new (specIrregularity_tilde_class);
    t_sampIdx i;

    x->x_irregularity = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 2:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            if (atom_getfloat (argv + 1) >= 1)
                x->x_algorithm = krimphoff;
            else
                x->x_algorithm = jensen;
            break;

        case 1:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            x->x_algorithm = TID_SPECIRREGULARITY_DEFAULTALGO;
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_algorithm = TID_SPECIRREGULARITY_DEFAULTALGO;
            break;

        default:
            post ("%s WARNING: Too many arguments supplied. Using default window size of %i and Jensen algorithm", x->x_objSymbol->s_name, TID_WINDOWSIZEDEFAULT);
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_algorithm = TID_SPECIRREGULARITY_DEFAULTALGO;
            break;
    }

    x->x_windowHalf = x->x_window * 0.5;
    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_windowFunction = blackman;
    x->x_normalize = true;
    x->x_powerSpectrum = false;
    x->x_lastDspTime = clock_getlogicaltime();

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

    return (x);
}


static t_int* specIrregularity_tilde_perform (t_int* w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_specIrregularity_tilde* x = (t_specIrregularity_tilde *)(w[1]);

    t_sample* in = (t_sample *)(w[2]);
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


static void specIrregularity_tilde_dsp (t_specIrregularity_tilde* x, t_signal** sp)
{
    dsp_add (
        specIrregularity_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if (sp[0]->s_sr != x->x_sr * x->x_overlap)
    {
        x->x_sr = sp[0]->s_sr / x->x_overlap;
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
    }
};


static void specIrregularity_tilde_free (t_specIrregularity_tilde* x)
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


void specIrregularity_tilde_setup (void)
{
    specIrregularity_tilde_class =
    class_new (
        gensym ("specIrregularity~"),
        (t_newmethod)specIrregularity_tilde_new,
        (t_method)specIrregularity_tilde_free,
        sizeof (t_specIrregularity_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)specIrregularity_tilde_new,
        gensym ("timbreIDLib/specIrregularity~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (specIrregularity_tilde_class, t_specIrregularity_tilde, x_f);

    class_addbang (specIrregularity_tilde_class, specIrregularity_tilde_bang);

    class_addmethod (
        specIrregularity_tilde_class,
        (t_method)specIrregularity_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        specIrregularity_tilde_class,
        (t_method)specIrregularity_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specIrregularity_tilde_class,
        (t_method)specIrregularity_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specIrregularity_tilde_class,
        (t_method)specIrregularity_tilde_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specIrregularity_tilde_class,
        (t_method)specIrregularity_tilde_normalize,
        gensym ("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specIrregularity_tilde_class,
        (t_method)specIrregularity_tilde_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specIrregularity_tilde_class,
        (t_method)specIrregularity_tilde_algorithm,
        gensym ("algorithm"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specIrregularity_tilde_class,
        (t_method)specIrregularity_tilde_dsp,
        gensym ("dsp"),
        0
    );
}
