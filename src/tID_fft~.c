/*

tID_fft~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* tID_fft_tilde_class;

typedef struct _tID_fft_tilde
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_sampIdx x_zeroPad;
    t_sampIdx x_zeroPadHalf;
    double x_lastDspTime;
    t_sample* x_signalBuffer;
    t_sample* x_fftwIn;
    t_sample* x_fftwInZeroPad;
    fftwf_complex* x_fftwOut;
    fftwf_complex* x_fftwOutZeroPad;
    fftwf_plan x_fftwPlan;
    fftwf_plan x_fftwPlanZeroPad;
    t_bool x_normalize;
    t_float* x_blackman;
    t_float* x_cosine;
    t_float* x_hamming;
    t_float* x_hann;
    t_atom* x_listOutReal;
    t_atom* x_listOutImag;
    t_atom* x_listOutRealZeroPad;
    t_atom* x_listOutImagZeroPad;
    t_outlet* x_realOut;
    t_outlet* x_imagOut;
    t_float x_f;
} t_tID_fft_tilde;


/* ------------------------ tID_fft~ -------------------------------- */

static void tID_fft_tilde_bang (t_tID_fft_tilde* x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float* windowFuncPtr, realMax, imagMax;
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

    if (x->x_zeroPad > 0)
    {
        // place windowed signal in x_fftwInZeroPad
        for (i = 0; i < x->x_window; i++)
            x->x_fftwInZeroPad[i] = x->x_fftwIn[i];

        // fill out the remaining space with zeros
        for (; i < x->x_zeroPad; i++)
            x->x_fftwInZeroPad[i] = 0.0;

        // execute the zero pad plan
        fftwf_execute (x->x_fftwPlanZeroPad);

        realMax = 0.0;
        imagMax = 0.0;

        if (x->x_normalize)
        {
            for (i = 0; i <= x->x_zeroPadHalf; i++)
            {
                t_float thisReal, thisImag;

                thisReal = fabsf (x->x_fftwOutZeroPad[i][0]);
                thisImag = fabsf (x->x_fftwOutZeroPad[i][1]);

                if (thisReal > realMax)
                    realMax = thisReal;

                if (thisImag > imagMax)
                    imagMax = thisImag;
            }

            realMax = (realMax == 0.0) ? 1.0 : realMax;
            imagMax = (imagMax == 0.0) ? 1.0 : imagMax;

            realMax = 1.0 / realMax;
            imagMax = 1.0 / imagMax;
        }

        for (i = 0; i <= x->x_zeroPadHalf; i++)
        {
            if (x->x_normalize)
            {
                SETFLOAT (x->x_listOutRealZeroPad + i, x->x_fftwOutZeroPad[i][0] * realMax);
                SETFLOAT (x->x_listOutImagZeroPad + i, x->x_fftwOutZeroPad[i][1] * imagMax);
            }
            else
            {
                SETFLOAT (x->x_listOutRealZeroPad + i, x->x_fftwOutZeroPad[i][0]);
                SETFLOAT (x->x_listOutImagZeroPad + i, x->x_fftwOutZeroPad[i][1]);
            }
        }

        outlet_list (x->x_imagOut, 0, x->x_zeroPadHalf + 1, x->x_listOutImagZeroPad);
        outlet_list (x->x_realOut, 0, x->x_zeroPadHalf + 1, x->x_listOutRealZeroPad);
    }
    else
    {
        fftwf_execute (x->x_fftwPlan);

        realMax = 0.0;
        imagMax = 0.0;

        if (x->x_normalize)
        {
            for (i = 0; i <= x->x_windowHalf; i++)
            {
                t_float thisReal, thisImag;

                thisReal = fabsf (x->x_fftwOut[i][0]);
                thisImag = fabsf (x->x_fftwOut[i][1]);

                if (thisReal > realMax)
                    realMax = thisReal;

                if (thisImag > imagMax)
                    imagMax = thisImag;
            }

            realMax = (realMax == 0.0) ? 1.0 : realMax;
            imagMax = (imagMax == 0.0) ? 1.0 : imagMax;

            realMax = 1.0 / realMax;
            imagMax = 1.0 / imagMax;
        }

        for (i = 0; i <= windowHalf; i++)
        {
            if (x->x_normalize)
            {
                SETFLOAT (x->x_listOutReal + i, x->x_fftwOut[i][0] * realMax);
                SETFLOAT (x->x_listOutImag + i, x->x_fftwOut[i][1] * imagMax);
            }
            else
            {
                SETFLOAT (x->x_listOutReal + i, x->x_fftwOut[i][0]);
                SETFLOAT (x->x_listOutImag + i, x->x_fftwOut[i][1]);
            }
        }

        outlet_list (x->x_imagOut, 0, windowHalf + 1, x->x_listOutImag);
        outlet_list (x->x_realOut, 0, windowHalf + 1, x->x_listOutReal);
    }
}


static void tID_fft_tilde_print (t_tID_fft_tilde* x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_sampIdx)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post ("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);

    if (x->x_zeroPad > x->x_window)
        post ("%s zero padding size: %i", x->x_objSymbol->s_name, x->x_zeroPad);
    else
        post ("%s zero padding OFF", x->x_objSymbol->s_name);
}


static void tID_fft_tilde_window (t_tID_fft_tilde* x, t_floatarg w)
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
    x->x_listOutReal = (t_atom *)t_resizebytes (x->x_listOutReal, (x->x_windowHalf + 1) * sizeof (t_atom), (windowHalf + 1) * sizeof (t_atom));
    x->x_listOutImag = (t_atom *)t_resizebytes (x->x_listOutImag, (x->x_windowHalf + 1) * sizeof (t_atom), (windowHalf + 1) * sizeof (t_atom));

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


static void tID_fft_tilde_overlap (t_tID_fft_tilde* x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1.0) ? 1.0 : o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void tID_fft_tilde_windowFunction (t_tID_fft_tilde* x, t_floatarg f)
{
    f = (f < 0.0) ? 0.0 : f;
    f = (f > 4.0) ? 4.0 : f;
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


static void tID_fft_tilde_zeroPad (t_tID_fft_tilde* x, t_floatarg z)
{
    t_sampIdx i, oldZeroPad, oldZeroPadHalf;

    oldZeroPad = x->x_zeroPad;
    oldZeroPadHalf = x->x_zeroPadHalf;

    // if the requested zero pad size is less or equal to current window size, zero pad should be 0 (OFF)
    z = (z <= x->x_window) ? 0.0 : z;
    x->x_zeroPad = z;
    x->x_zeroPadHalf = x->x_zeroPad * 0.5;

    // resize zeroPad FFTW input buffer
    x->x_fftwInZeroPad = (t_sample *)t_resizebytes (x->x_fftwInZeroPad, oldZeroPad * sizeof (t_sample), x->x_zeroPad * sizeof (t_sample));

    // free the zeroPad FFTW output buffer, and re-malloc according to new window
    fftwf_free (x->x_fftwOutZeroPad);

    // destroy old plan
    fftwf_destroy_plan (x->x_fftwPlanZeroPad);

    // allocate new fftwf_complex memory for the plan based on new zero pad size
    x->x_fftwOutZeroPad = (fftwf_complex *) fftwf_alloc_complex (x->x_zeroPadHalf + 1);

    // create a new DFT plan based on new window size
    x->x_fftwPlanZeroPad = fftwf_plan_dft_r2c_1d (x->x_zeroPad, x->x_fftwInZeroPad, x->x_fftwOutZeroPad, FFTWPLANNERFLAG);

    x->x_listOutRealZeroPad = (t_atom *)t_resizebytes (x->x_listOutRealZeroPad, (oldZeroPadHalf + 1) * sizeof (t_atom), (x->x_zeroPadHalf + 1) * sizeof (t_atom));
    x->x_listOutImagZeroPad = (t_atom *)t_resizebytes (x->x_listOutImagZeroPad, (oldZeroPadHalf + 1) * sizeof (t_atom), (x->x_zeroPadHalf + 1) * sizeof (t_atom));

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_zeroPad; i++)
        x->x_fftwInZeroPad[i] = 0.0;

    if (x->x_zeroPad>0)
        post ("%s zero padding size: %i", x->x_objSymbol->s_name, x->x_zeroPad);
    else
        post ("%s zero padding OFF", x->x_objSymbol->s_name);
}


static void tID_fft_tilde_normalize (t_tID_fft_tilde* x, t_floatarg norm)
{
    norm = (norm < 0) ? 0 : norm;
    norm = (norm > 1) ? 1 : norm;
    x->x_normalize = norm;

    if (x->x_normalize)
        post ("%s normalization ON.", x->x_objSymbol->s_name);
    else
        post ("%s normalization OFF.", x->x_objSymbol->s_name);
}


static void* tID_fft_tilde_new (t_symbol* s, int argc, t_atom* argv)
{
    t_tID_fft_tilde* x = (t_tID_fft_tilde *)pd_new (tID_fft_tilde_class);
    t_sampIdx i;

    x->x_realOut = outlet_new (&x->x_obj, gensym ("list"));
    x->x_imagOut = outlet_new (&x->x_obj, gensym ("list"));

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
    x->x_zeroPad = 0;
    x->x_zeroPadHalf = 0;
    x->x_normalize = false;
    x->x_lastDspTime = clock_getlogicaltime();

    x->x_signalBuffer = (t_sample *)t_getbytes ((x->x_window + x->x_n) * sizeof (t_sample));
    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));
    x->x_fftwInZeroPad = (t_sample *)t_getbytes (x->x_zeroPad * sizeof (t_sample));
    x->x_listOutReal = (t_atom *)t_getbytes ((x->x_windowHalf + 1) * sizeof (t_atom));
    x->x_listOutImag = (t_atom *)t_getbytes ((x->x_windowHalf + 1) * sizeof (t_atom));
    x->x_listOutRealZeroPad = (t_atom *)t_getbytes ((x->x_zeroPadHalf + 1) * sizeof (t_atom));
    x->x_listOutImagZeroPad = (t_atom *)t_getbytes ((x->x_zeroPadHalf + 1) * sizeof (t_atom));

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
    x->x_fftwOutZeroPad = (fftwf_complex *)fftwf_alloc_complex (x->x_zeroPadHalf + 1);

    // DFT plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);
    x->x_fftwPlanZeroPad = fftwf_plan_dft_r2c_1d (x->x_zeroPad, x->x_fftwInZeroPad, x->x_fftwOutZeroPad, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_window; i++)
        x->x_fftwIn[i] = 0.0;

    return (x);
}


static t_int* tID_fft_tilde_perform (t_int* w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_tID_fft_tilde* x = (t_tID_fft_tilde *)(w[1]);

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


static void tID_fft_tilde_dsp (t_tID_fft_tilde* x, t_signal** sp)
{
    dsp_add (
        tID_fft_tilde_perform,
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


static void tID_fft_tilde_free (t_tID_fft_tilde* x)
{
    // free the input buffer memory
    t_freebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample));

    // free the list out memory
    t_freebytes (x->x_listOutReal, (x->x_windowHalf + 1) * sizeof (t_atom));
    t_freebytes (x->x_listOutImag, (x->x_windowHalf + 1) * sizeof (t_atom));
    t_freebytes (x->x_listOutRealZeroPad, (x->x_zeroPadHalf + 1) * sizeof (t_atom));
    t_freebytes (x->x_listOutImagZeroPad, (x->x_zeroPadHalf + 1) * sizeof (t_atom));

    // free FFTW stuff
    t_freebytes (x->x_fftwIn, x->x_window * sizeof (t_sample));
    t_freebytes (x->x_fftwInZeroPad, (x->x_zeroPad) * sizeof (t_sample));
    fftwf_free (x->x_fftwOut);
    fftwf_free (x->x_fftwOutZeroPad);
    fftwf_destroy_plan (x->x_fftwPlan);
    fftwf_destroy_plan (x->x_fftwPlanZeroPad);

    // free the window memory
    t_freebytes (x->x_blackman, x->x_window * sizeof (t_float));
    t_freebytes (x->x_cosine, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hamming, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hann, x->x_window * sizeof (t_float));
}

void tID_fft_tilde_setup (void)
{
    tID_fft_tilde_class =
    class_new (
        gensym ("tID_fft~"),
        (t_newmethod)tID_fft_tilde_new,
        (t_method)tID_fft_tilde_free,
        sizeof (t_tID_fft_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)tID_fft_tilde_new,
        gensym ("timbreIDLib/tID_fft~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (tID_fft_tilde_class, t_tID_fft_tilde, x_f);

    class_addbang (tID_fft_tilde_class, tID_fft_tilde_bang);

    class_addmethod (
        tID_fft_tilde_class,
        (t_method)tID_fft_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        tID_fft_tilde_class,
        (t_method)tID_fft_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tID_fft_tilde_class,
        (t_method)tID_fft_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tID_fft_tilde_class,
        (t_method)tID_fft_tilde_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tID_fft_tilde_class,
        (t_method)tID_fft_tilde_zeroPad,
        gensym ("zero_pad"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tID_fft_tilde_class,
        (t_method)tID_fft_tilde_normalize,
        gensym ("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tID_fft_tilde_class,
        (t_method)tID_fft_tilde_dsp,
        gensym ("dsp"),
        0
    );
}
