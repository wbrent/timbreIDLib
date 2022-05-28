/*

dct~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *dct_tilde_class;

typedef struct _dct_tilde
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_sampIdx x_window;
    t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_bool x_normalize;
    double x_lastDspTime;
    t_sample *x_signalBuffer;
    t_float *x_dctIn;
    t_float *x_dctOut;
    fftwf_plan x_fftwDctPlan;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_atom *x_listOut;
    t_outlet *x_dct;
    t_float x_f;
} t_dct_tilde;


/* ------------------------ dct~ -------------------------------- */

static void dct_tilde_bang (t_dct_tilde *x)
{
    t_sampIdx i, j, window, bangSample;
    t_float *windowFuncPtr;
    double currentTime;

    window = x->x_window;

    currentTime = clock_gettimesince (x->x_lastDspTime);
    bangSample = roundf ((currentTime / 1000.0) * x->x_sr);

    if (bangSample >= x->x_n)
        bangSample = x->x_n - 1;

    // construct analysis window using bangSample as the end of the window
    for (i = 0, j = bangSample; i < window; i++, j++)
        x->x_dctIn[i] = x->x_signalBuffer[j];

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
            x->x_dctIn[i] *= *windowFuncPtr;

    fftwf_execute (x->x_fftwDctPlan);

    if (x->x_normalize)
        tIDLib_normalPeak(x->x_window, x->x_dctOut);

    for (i = 0; i < window; i++)
        SETFLOAT (x->x_listOut + i, x->x_dctOut[i]);

     outlet_list (x->x_dct, 0, window, x->x_listOut);
}


static void dct_tilde_print (t_dct_tilde *x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_sampIdx)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
}


static void dct_tilde_window (t_dct_tilde *x, t_floatarg w)
{
    t_sampIdx i, window;

    window = w;

    if (window < TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
    }

    x->x_signalBuffer = (t_sample *)t_resizebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample), (window + x->x_n) * sizeof (t_sample));
    x->x_dctIn = (t_float *)t_resizebytes (x->x_dctIn, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_dctOut = (t_float *)t_resizebytes (x->x_dctOut, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_listOut = (t_atom *)t_resizebytes (x->x_listOut, x->x_window * sizeof (t_atom), window * sizeof (t_atom));

    x->x_blackman = (t_float *)t_resizebytes (x->x_blackman, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_cosine = (t_float *)t_resizebytes (x->x_cosine, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_hamming = (t_float *)t_resizebytes (x->x_hamming, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_hann = (t_float *)t_resizebytes (x->x_hann, x->x_window * sizeof (t_float), window * sizeof (t_float));

    x->x_window = window;

    // destroy old plan, which depended on x->x_window
    fftwf_destroy_plan (x->x_fftwDctPlan);

    // create a new DCT plan based on new window size
    x->x_fftwDctPlan = fftwf_plan_r2r_1d (x->x_window, x->x_dctIn, x->x_dctOut, FFTW_REDFT10, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_window; i++)
     {
        x->x_dctIn[i] = 0.0;
        x->x_dctOut[i] = 0.0;
    }

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


static void dct_tilde_overlap (t_dct_tilde *x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o<1.0)?1.0:o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void dct_tilde_windowFunction (t_dct_tilde *x, t_floatarg f)
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


static void dct_tilde_normalize (t_dct_tilde *x, t_floatarg norm)
{
    norm = (norm<0.0)?0.0:norm;
    norm = (norm>1.0)?1.0:norm;
    x->x_normalize = norm;

    if (x->x_normalize)
        post ("%s normalization ON.", x->x_objSymbol->s_name);
    else
        post ("%s normalization OFF.", x->x_objSymbol->s_name);
}


static void *dct_tilde_new (t_symbol *s, int argc, t_atom *argv)
{
    t_dct_tilde *x = (t_dct_tilde *)pd_new (dct_tilde_class);
    t_sampIdx i;

    x->x_dct = outlet_new (&x->x_obj, gensym ("list"));

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

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_windowFunction = blackman;
    x->x_normalize = true;
    x->x_lastDspTime = clock_getlogicaltime();

    x->x_signalBuffer = (t_sample *)t_getbytes ((x->x_window + x->x_n) * sizeof (t_sample));
    x->x_dctIn = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_dctOut = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_listOut = (t_atom *)t_getbytes (x->x_window * sizeof (t_atom));

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

    // DCT plan. FFTW_REDFT10 is the DCT-II
    x->x_fftwDctPlan = fftwf_plan_r2r_1d (x->x_window, x->x_dctIn, x->x_dctOut, FFTW_REDFT10, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_window; i++)
     {
        x->x_dctIn[i] = 0.0;
        x->x_dctOut[i] = 0.0;
    }

    return (x);
}


static t_int *dct_tilde_perform (t_int *w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_dct_tilde *x = (t_dct_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
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


static void dct_tilde_dsp (t_dct_tilde *x, t_signal **sp)
{
    dsp_add (
        dct_tilde_perform,
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


static void dct_tilde_free (t_dct_tilde *x)
{
    // free the input buffer memory
    t_freebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample));

    // free the list out memory
    t_freebytes (x->x_listOut, x->x_window * sizeof (t_atom));

    // free FFTW stuff
    t_freebytes (x->x_dctIn, x->x_window * sizeof (t_float));
    t_freebytes (x->x_dctOut, x->x_window * sizeof (t_float));
    fftwf_destroy_plan (x->x_fftwDctPlan);

    // free the window memory
    t_freebytes (x->x_blackman, x->x_window * sizeof (t_float));
    t_freebytes (x->x_cosine, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hamming, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hann, x->x_window * sizeof (t_float));
}

void dct_tilde_setup (void)
{
    dct_tilde_class =
    class_new (
        gensym ("dct~"),
        (t_newmethod)dct_tilde_new,
        (t_method)dct_tilde_free,
        sizeof (t_dct_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)dct_tilde_new,
        gensym ("timbreIDLib/dct~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (dct_tilde_class, t_dct_tilde, x_f);

    class_addbang (dct_tilde_class, dct_tilde_bang);

    class_addmethod (
        dct_tilde_class,
        (t_method)dct_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        dct_tilde_class,
        (t_method)dct_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        dct_tilde_class,
        (t_method)dct_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        dct_tilde_class,
        (t_method)dct_tilde_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        dct_tilde_class,
        (t_method)dct_tilde_normalize,
        gensym ("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        dct_tilde_class,
        (t_method)dct_tilde_dsp,
        gensym ("dsp"),
        0
    );
}
