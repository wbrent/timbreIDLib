/*

cepstrum

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* cepstrum_class;

typedef struct _cepstrum
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_bool x_powerSpectrum;
    t_bool x_powerCepstrum;
    t_bool x_spectrumOffset;
    t_sample* x_fftwIn;
    fftwf_complex* x_fftwOut;
    fftwf_plan x_fftwForwardPlan;
    fftwf_plan x_fftwBackwardPlan;
    t_float* x_blackman;
    t_float* x_cosine;
    t_float* x_hamming;
    t_float* x_hann;
    t_word* x_vec;
    t_symbol* x_arrayName;
    t_sampIdx x_arrayPoints;
    t_atom* x_listOut;
    t_outlet* x_featureList;
} t_cepstrum;


/* ------------------------ cepstrum -------------------------------- */
static void cepstrum_resizeWindow (t_cepstrum* x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx* endSamp)
{
    t_sampIdx oldWindowHalf, windowHalf;

    windowHalf = window * 0.5;
    oldWindowHalf = oldWindow * 0.5;

    if (window < TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        windowHalf = window * 0.5;
        post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);

        *endSamp = startSamp + window - 1;
        if (*endSamp >= x->x_arrayPoints)
            *endSamp = x->x_arrayPoints - 1;
    }

    // hang on to these values for next time
    x->x_window = window;
    x->x_windowHalf = windowHalf;

    x->x_fftwIn = (t_sample *)t_resizebytes (x->x_fftwIn, oldWindow * sizeof (t_sample), x->x_window * sizeof (t_sample));

    fftwf_free (x->x_fftwOut);
    fftwf_destroy_plan (x->x_fftwForwardPlan);
    fftwf_destroy_plan (x->x_fftwBackwardPlan);
    // set up a new FFTW output buffer
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);
    // FFTW plan
    x->x_fftwForwardPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);
    x->x_fftwBackwardPlan = fftwf_plan_dft_c2r_1d (x->x_window, x->x_fftwOut, x->x_fftwIn, FFTWPLANNERFLAG);

    x->x_blackman = (t_float *)t_resizebytes (x->x_blackman, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));
    x->x_cosine = (t_float *)t_resizebytes (x->x_cosine, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));
    x->x_hamming = (t_float *)t_resizebytes (x->x_hamming, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));
    x->x_hann = (t_float *)t_resizebytes (x->x_hann, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));

    tIDLib_blackmanWindow (x->x_blackman, x->x_window);
    tIDLib_cosineWindow (x->x_cosine, x->x_window);
    tIDLib_hammingWindow (x->x_hamming, x->x_window);
    tIDLib_hannWindow (x->x_hann, x->x_window);

    // resize x_listOut
    x->x_listOut = (t_atom *)t_resizebytes (x->x_listOut, (oldWindowHalf + 1) * sizeof (t_atom), (x->x_windowHalf + 1) * sizeof (t_atom));
}


static void cepstrum_analyze (t_cepstrum* x, t_floatarg start, t_floatarg n)
{
    t_sampIdx i, j, window, startSamp, endSamp;
    t_float nRecip;
    t_float* windowFuncPtr;
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        startSamp = (start < 0) ? 0 : start;

        if (n)
            endSamp = startSamp + n - 1;
        else
            endSamp = startSamp + x->x_window - 1;

        if (endSamp >= x->x_arrayPoints)
            endSamp = x->x_arrayPoints - 1;

        if (endSamp <= startSamp)
        {
            post ("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        window = endSamp - startSamp + 1;

        if (x->x_window != window)
            cepstrum_resizeWindow (x, x->x_window, window, startSamp, &endSamp);

        // construct analysis window
        for (i = 0, j = startSamp; j <= endSamp; i++, j++)
            x->x_fftwIn[i] = x->x_vec[j].w_float;

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
            for (i = 0; i < x->x_window; i++, windowFuncPtr++)
                x->x_fftwIn[i] *= *windowFuncPtr;

        fftwf_execute (x->x_fftwForwardPlan);

        tIDLib_power (x->x_windowHalf + 1, x->x_fftwOut, x->x_fftwIn);

        if ( !x->x_powerSpectrum)
            tIDLib_mag (x->x_windowHalf + 1, x->x_fftwIn);

        // add 1.0 to power or magnitude spectrum before taking the log and then IFT. Avoid large negative values from log (negativeNum)
        if (x->x_spectrumOffset)
            for (i = 0; i < x->x_windowHalf + 1; i++)
                x->x_fftwIn[i] += 1.0;

        tIDLib_log (x->x_windowHalf + 1, x->x_fftwIn);

        // copy forward DFT magnitude result into real part of backward DFT complex input buffer, and zero out the imaginary part. fftwOut is only N/2+1 points long, while fftwIn is N points long
        for (i = 0; i < x->x_windowHalf + 1; i++)
        {
            x->x_fftwOut[i][0] = x->x_fftwIn[i];
            x->x_fftwOut[i][1] = 0.0;
        }

        fftwf_execute (x->x_fftwBackwardPlan);

        // divide by N because FFTW's c2r transform doesn't normalize
        nRecip = 1.0 / x->x_window;

        for (i = 0; i < x->x_windowHalf + 1; i++)
            x->x_fftwIn[i] *= nRecip;

        // optionally square the cepstrum results for power cepstrum
        if (x->x_powerCepstrum)
            for (i = 0; i < x->x_windowHalf + 1; i++)
                x->x_fftwIn[i] = x->x_fftwIn[i] * x->x_fftwIn[i];

        for (i = 0; i < x->x_windowHalf + 1; i++)
            SETFLOAT (x->x_listOut + i, x->x_fftwIn[i]);

        outlet_list (x->x_featureList, 0, x->x_windowHalf + 1, x->x_listOut);
    }
}


static void cepstrum_chain_fftData (t_cepstrum* x, t_symbol* s, int argc, t_atom* argv)
{
    t_sampIdx i, windowHalf;
    t_float nRecip;

    // incoming fftData list should be 2*(N/2+1) elements long, so windowHalf is:
    windowHalf = argc - 2;
    windowHalf *= 0.5;

    // make sure that windowHalf == x->x_windowHalf in order to avoid an out of bounds memory read in the tIDLib_ functions below. we won't resize all memory based on an incoming chain_ command with a different window size. instead, just throw an error and exit
    if (windowHalf != x->x_windowHalf)
    {
        pd_error (x, "%s: window size of chain_ message (%lu) does not match current window size (%lu)", x->x_objSymbol->s_name, windowHalf * 2, x->x_window);
        return;
    }

    // fill the x_fftwOut buffer with the incoming fftData list, for both real and imag elements
    for (i = 0; i <= x->x_windowHalf; i++)
    {
        x->x_fftwOut[i][0] = atom_getfloat (argv + i);
        x->x_fftwOut[i][1] = atom_getfloat (argv + (x->x_windowHalf + 1) + i);
    }

    tIDLib_power (x->x_windowHalf + 1, x->x_fftwOut, x->x_fftwIn);

    if ( !x->x_powerSpectrum)
        tIDLib_mag (x->x_windowHalf + 1, x->x_fftwIn);

    // add 1.0 to power or magnitude spectrum before taking the log and then IFT. Avoid large negative values from log (negativeNum)
    if (x->x_spectrumOffset)
        for (i = 0; i < x->x_windowHalf + 1; i++)
            x->x_fftwIn[i] += 1.0;

    tIDLib_log (x->x_windowHalf + 1, x->x_fftwIn);

    // copy forward DFT magnitude result into real part of backward DFT complex input buffer, and zero out the imaginary part. fftwOut is only N/2+1 points long, while fftwIn is N points long
    for (i = 0; i < x->x_windowHalf + 1; i++)
    {
        x->x_fftwOut[i][0] = x->x_fftwIn[i];
        x->x_fftwOut[i][1] = 0.0;
    }

    fftwf_execute (x->x_fftwBackwardPlan);

    // divide by N because FFTW's c2r transform doesn't normalize
    nRecip = 1.0 / x->x_window;

    for (i = 0; i < x->x_windowHalf + 1; i++)
        x->x_fftwIn[i] *= nRecip;

    // optionally square the cepstrum results for power cepstrum
    if (x->x_powerCepstrum)
        for (i = 0; i < x->x_windowHalf + 1; i++)
            x->x_fftwIn[i] = x->x_fftwIn[i] * x->x_fftwIn[i];

    for (i = 0; i < x->x_windowHalf + 1; i++)
        SETFLOAT (x->x_listOut + i, x->x_fftwIn[i]);

    outlet_list (x->x_featureList, 0, x->x_windowHalf + 1, x->x_listOut);
}


static void cepstrum_chain_magSpec (t_cepstrum* x, t_symbol* s, int argc, t_atom* argv)
{
    t_sampIdx i, windowHalf;
    t_float nRecip;

    // incoming magSpec list should be N/2+1 elements long, so windowHalf is one less than this
    windowHalf = argc - 1;

    // make sure that windowHalf == x->x_windowHalf in order to avoid an out of bounds memory read in the tIDLib_ functions below. we won't resize all memory based on an incoming chain_ command with a different window size. instead, just throw an error and exit
    if (windowHalf != x->x_windowHalf)
    {
        pd_error (x, "%s: window size of chain_ message (%lu) does not match current window size (%lu)", x->x_objSymbol->s_name, windowHalf * 2, x->x_window);
        return;
    }

    // fill the x_fftwIn buffer with the incoming magSpec list
    for (i = 0; i <= x->x_windowHalf; i++)
        x->x_fftwIn[i] = atom_getfloat (argv + i);

    // add 1.0 to power or magnitude spectrum before taking the log and then IFT. Avoid large negative values from log (negativeNum)
    if (x->x_spectrumOffset)
        for (i = 0; i < x->x_windowHalf + 1; i++)
            x->x_fftwIn[i] += 1.0;

    tIDLib_log (x->x_windowHalf + 1, x->x_fftwIn);

    // copy forward DFT magnitude result into real part of backward DFT complex input buffer, and zero out the imaginary part. fftwOut is only N/2+1 points long, while fftwIn is N points long
    for (i = 0; i < x->x_windowHalf + 1; i++)
    {
        x->x_fftwOut[i][0] = x->x_fftwIn[i];
        x->x_fftwOut[i][1] = 0.0;
    }

    fftwf_execute (x->x_fftwBackwardPlan);

    // divide by N because FFTW's c2r transform doesn't normalize
    nRecip = 1.0 / x->x_window;

    for (i = 0; i < x->x_windowHalf + 1; i++)
        x->x_fftwIn[i] *= nRecip;

    // optionally square the cepstrum results for power cepstrum
    if (x->x_powerCepstrum)
        for (i = 0; i < x->x_windowHalf + 1; i++)
            x->x_fftwIn[i] = x->x_fftwIn[i] * x->x_fftwIn[i];

    for (i = 0; i < x->x_windowHalf + 1; i++)
        SETFLOAT (x->x_listOut + i, x->x_fftwIn[i]);

    outlet_list (x->x_featureList, 0, x->x_windowHalf + 1, x->x_listOut);
}


// analyze the whole damn array
static void cepstrum_bang (t_cepstrum* x)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx window, startSamp;
        startSamp = 0;
        window = x->x_arrayPoints;
        cepstrum_analyze (x, startSamp, window);
    }
}


static void cepstrum_set (t_cepstrum* x, t_symbol* s)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (s, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void cepstrum_print (t_cepstrum* x)
{
    post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s power cepstrum: %i", x->x_objSymbol->s_name, x->x_powerCepstrum);
    post ("%s spectrum offset: %i", x->x_objSymbol->s_name, x->x_spectrumOffset);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
}


static void cepstrum_samplerate (t_cepstrum* x, t_floatarg sr)
{
    if (sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void cepstrum_window (t_cepstrum* x, t_floatarg w)
{
    t_sampIdx endSamp;

    // catch negative window size arguments here while data type is still t_floatarg. once passed to _resizeWindow, it will be cast to type t_sampIdx and negative values will turn to garbage.
    w = (w < 0.0) ? TID_MINWINDOWSIZE : w;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow() requires that
    endSamp = 0;

    cepstrum_resizeWindow (x, x->x_window, w, 0, &endSamp);
}


static void cepstrum_windowFunction (t_cepstrum* x, t_floatarg f)
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


static void cepstrum_powerSpectrum (t_cepstrum* x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void cepstrum_powerCepstrum (t_cepstrum* x, t_floatarg power)
{
    power = (power < 0) ? 0 : power;
    power = (power > 1) ? 1 : power;
    x->x_powerCepstrum = power;

    if (x->x_powerCepstrum)
        post ("%s reporting power cepstrum", x->x_objSymbol->s_name);
    else
        post ("%s reporting magnitude cepstrum", x->x_objSymbol->s_name);
}


static void cepstrum_spectrumOffset (t_cepstrum* x, t_floatarg offset)
{
    offset = (offset < 0) ? 0 : offset;
    offset = (offset > 1) ? 1 : offset;
    x->x_spectrumOffset = offset;

    if (x->x_spectrumOffset)
        post ("%s spectrum offset ON", x->x_objSymbol->s_name);
    else
        post ("%s spectrum offset OFF", x->x_objSymbol->s_name);
}


static void* cepstrum_new (t_symbol* s, int argc, t_atom* argv)
{
    t_cepstrum* x = (t_cepstrum *)pd_new (cepstrum_class);
    t_sampIdx i;
//	t_garray *a;

    x->x_featureList = outlet_new (&x->x_obj, gensym ("list"));

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 1:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            break;

        case 0:
            post ("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze()
            x->x_arrayName = gensym ("NOARRAYSPECIFIED");
            break;

        default:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post ("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
            break;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_window = TID_WINDOWSIZEDEFAULT;
    x->x_windowHalf = x->x_window * 0.5;
    x->x_windowFunction = blackman;
    x->x_powerSpectrum = false;
    x->x_powerCepstrum = false;
    x->x_spectrumOffset = false;

    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));
    x->x_listOut = (t_atom *)t_getbytes ((x->x_windowHalf + 1) * sizeof (t_atom));

    // set up the FFTW output buffer. Is there no function to initialize it?
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);

    // Forward DFT plan
    x->x_fftwForwardPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

    // Backward DFT plan
    x->x_fftwBackwardPlan = fftwf_plan_dft_c2r_1d (x->x_window, x->x_fftwOut, x->x_fftwIn, FFTWPLANNERFLAG);

    for (i = 0; i < x->x_window; i++)
        x->x_fftwIn[i] = 0.0;

    x->x_blackman = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_cosine = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_hamming = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_hann = (t_float *)t_getbytes (x->x_window * sizeof (t_float));

     // initialize signal windowing functions
    tIDLib_blackmanWindow (x->x_blackman, x->x_window);
    tIDLib_cosineWindow (x->x_cosine, x->x_window);
    tIDLib_hammingWindow (x->x_hamming, x->x_window);
    tIDLib_hannWindow (x->x_hann, x->x_window);

    return (x);
}


static void cepstrum_free (t_cepstrum* x)
{
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


void cepstrum_setup (void)
{
    cepstrum_class =
    class_new (
        gensym ("cepstrum"),
        (t_newmethod)cepstrum_new,
        (t_method)cepstrum_free,
        sizeof (t_cepstrum),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)cepstrum_new,
        gensym ("timbreIDLib/cepstrum"),
        A_GIMME,
        0
    );

    class_addbang (cepstrum_class, cepstrum_bang);

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_chain_fftData,
        gensym ("chain_fftData"),
        A_GIMME,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_chain_magSpec,
        gensym ("chain_magSpec"),
        A_GIMME,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_set,
        gensym ("set"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_powerCepstrum,
        gensym ("power_cepstrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrum_class,
        (t_method)cepstrum_spectrumOffset,
        gensym ("spectrum_offset"),
        A_DEFFLOAT,
        0
    );
}
