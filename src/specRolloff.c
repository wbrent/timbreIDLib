/*

specRolloff

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *specRolloff_class;

typedef struct _specRolloff
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_bool x_powerSpectrum;
    t_sample *x_fftwIn;
    fftwf_complex *x_fftwOut;
    fftwf_plan x_fftwPlan;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_word *x_vec;
    t_symbol *x_arrayName;
    t_sampIdx x_arrayPoints;
    t_float x_concentration;
    t_float *x_binFreqs;
    t_outlet *x_rolloff;
} t_specRolloff;


/* ------------------------ specRolloff -------------------------------- */
static void specRolloff_resizeWindow (t_specRolloff *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
{
    t_sampIdx i, oldWindowHalf, windowHalf;

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
    x->x_binFreqs = (t_float *)t_resizebytes (x->x_binFreqs, (oldWindowHalf + 1) * sizeof (t_float), (x->x_windowHalf + 1) * sizeof (t_float));

    fftwf_free (x->x_fftwOut);
    fftwf_destroy_plan (x->x_fftwPlan);
    // set up a new FFTW output buffer
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);
    // FFTW plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

    x->x_blackman = (t_float *)t_resizebytes (x->x_blackman, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));
    x->x_cosine = (t_float *)t_resizebytes (x->x_cosine, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));
    x->x_hamming = (t_float *)t_resizebytes (x->x_hamming, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));
    x->x_hann = (t_float *)t_resizebytes (x->x_hann, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));

    tIDLib_blackmanWindow (x->x_blackman, x->x_window);
    tIDLib_cosineWindow (x->x_cosine, x->x_window);
    tIDLib_hammingWindow (x->x_hamming, x->x_window);
    tIDLib_hannWindow (x->x_hann, x->x_window);

     for (i = 0; i <= x->x_windowHalf; i++)
        x->x_binFreqs[i] = tIDLib_bin2freq(i, x->x_window, x->x_sr);
}


static void specRolloff_analyze (t_specRolloff *x, t_floatarg start, t_floatarg n)
{
    t_garray *a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, window, startSamp, endSamp;
        t_float energyTarget, energySum, rolloff, *windowFuncPtr;

        startSamp = (start < 0) ? 0 : start;

        if (n)
            endSamp = startSamp + n - 1;
        else
            endSamp = startSamp + x->x_window - 1;

        if (endSamp >= x->x_arrayPoints)
            endSamp = x->x_arrayPoints - 1;

        window = endSamp - startSamp + 1;

        if (endSamp <= startSamp)
        {
            post ("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        if (x->x_window != window)
            specRolloff_resizeWindow (x, x->x_window, window, startSamp, &endSamp);

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

        fftwf_execute (x->x_fftwPlan);

        tIDLib_power (x->x_windowHalf + 1, x->x_fftwOut, x->x_fftwIn);

        if ( !x->x_powerSpectrum)
            tIDLib_mag (x->x_windowHalf + 1, x->x_fftwIn);

        energyTarget = 0.0;
        for (i = 0; i <= x->x_windowHalf; i++)
            energyTarget += x->x_fftwIn[i];

        energyTarget *= x->x_concentration;

        energySum=0.0;
        i = 0;

        while(energySum <= energyTarget)
        {
            energySum += x->x_fftwIn[i];
            i++;

            if (i>x->x_windowHalf)
                break;
        }

        if (i>0)
            rolloff = x->x_binFreqs[i - 1]; // back up one because the last one went over...
        else
            rolloff =  - 1.0;

        outlet_float (x->x_rolloff, rolloff);
    }
}


static void specRolloff_chain_fftData (t_specRolloff *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;
    t_float energyTarget, energySum, rolloff;

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

    energyTarget = 0.0;
    for (i = 0; i <= x->x_windowHalf; i++)
        energyTarget += x->x_fftwIn[i];

    energyTarget *= x->x_concentration;

    energySum=0.0;
    i = 0;

    while(energySum <= energyTarget)
    {
        energySum += x->x_fftwIn[i];
        i++;

        if (i>x->x_windowHalf)
            break;
    }

    if (i>0)
        rolloff = x->x_binFreqs[i - 1]; // back up one because the last one went over...
    else
        rolloff =  - 1.0;

    outlet_float (x->x_rolloff, rolloff);
}


static void specRolloff_chain_magSpec (t_specRolloff *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;
    t_float energyTarget, energySum, rolloff;

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

    energyTarget = 0.0;
    for (i = 0; i <= x->x_windowHalf; i++)
        energyTarget += x->x_fftwIn[i];

    energyTarget *= x->x_concentration;

    energySum=0.0;
    i = 0;

    while(energySum <= energyTarget)
    {
        energySum += x->x_fftwIn[i];
        i++;

        if (i>x->x_windowHalf)
            break;
    }

    if (i>0)
        rolloff = x->x_binFreqs[i - 1]; // back up one because the last one went over...
    else
        rolloff =  - 1.0;

    outlet_float (x->x_rolloff, rolloff);
}


// analyze the whole damn array
static void specRolloff_bang (t_specRolloff *x)
{
    t_garray *a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx window, startSamp;
        startSamp = 0;
        window = x->x_arrayPoints;
        specRolloff_analyze (x, startSamp, window);
    }
}


static void specRolloff_concentration(t_specRolloff *x, t_floatarg c)
{
    if (c<0 || c>1.0)
        post ("%s concentration must be between 0.0 and 1.0.", x->x_objSymbol->s_name);
    else
        x->x_concentration = c;

    post ("%s concentration: %f", x->x_objSymbol->s_name, x->x_concentration);
}


static void specRolloff_set (t_specRolloff *x, t_symbol *s)
{
    t_garray *a;

    if ( !(a = (t_garray *)pd_findbyclass (s, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void specRolloff_print (t_specRolloff *x)
{
    post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post ("%s concentration: %f", x->x_objSymbol->s_name, x->x_concentration);
}


static void specRolloff_samplerate (t_specRolloff *x, t_floatarg sr)
{
    t_sampIdx i;

    if (sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;

     for (i = 0; i <= x->x_windowHalf; i++)
        x->x_binFreqs[i] = tIDLib_bin2freq(i, x->x_window, x->x_sr);
}


static void specRolloff_window (t_specRolloff *x, t_floatarg w)
{
    t_sampIdx endSamp;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow () requires that
    endSamp = 0;

    specRolloff_resizeWindow (x, x->x_window, w, 0, &endSamp);
}


static void specRolloff_windowFunction (t_specRolloff *x, t_floatarg f)
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


static void specRolloff_powerSpectrum (t_specRolloff *x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void *specRolloff_new (t_symbol *s, int argc, t_atom *argv)
{
    t_specRolloff *x = (t_specRolloff *)pd_new (specRolloff_class);
    t_sampIdx i;
//	t_garray *a;

    x->x_rolloff = outlet_new (&x->x_obj, gensym ("list"));

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    // x->x_sr is used in argument check, so must be assigned before switch
    x->x_sr = TID_SAMPLERATEDEFAULT;

    switch (argc)
    {
        case 2:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_concentration = atom_getfloat (argv + 1);
            if (x->x_concentration<0)
            {
                post ("%s concentration must be between 0.0 and 1.0. Using default concentration of %0.2f instead.", x->x_objSymbol->s_name, TID_SPECROLLOFF_DEFAULTCONCENTRATION);
                x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            }
            else if (x->x_concentration>1.0)
            {
                post ("%s concentration must be between 0.0 and 1.0. Using default concentration of %0.2f instead.", x->x_objSymbol->s_name, TID_SPECROLLOFF_DEFAULTCONCENTRATION);
                x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            };
            break;

        case 1:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            break;

        case 0:
            post ("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze ()
            x->x_arrayName = gensym ("NOARRAYSPECIFIED");
            x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            break;

        default:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post ("%s WARNING: extra arguments ignored. Using default concentration of %0.2f.", x->x_objSymbol->s_name, TID_SPECROLLOFF_DEFAULTCONCENTRATION);
            x->x_concentration = TID_SPECROLLOFF_DEFAULTCONCENTRATION;
            break;
    }

    x->x_window = TID_WINDOWSIZEDEFAULT;
    x->x_windowHalf = x->x_window * 0.5;
    x->x_windowFunction = blackman;
    x->x_powerSpectrum = false;

    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);

    // FFTW plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

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

    x->x_binFreqs = (t_float *)t_getbytes ((x->x_windowHalf + 1) * sizeof (t_float));

     for (i = 0; i <= x->x_windowHalf; i++)
        x->x_binFreqs[i] = tIDLib_bin2freq(i, x->x_window, x->x_sr);

    return (x);
}


static void specRolloff_free (t_specRolloff *x)
{
    // free FFTW stuff
    t_freebytes (x->x_fftwIn, x->x_window * sizeof (t_sample));
    fftwf_free (x->x_fftwOut);
    fftwf_destroy_plan (x->x_fftwPlan);

    // free the window memory
    t_freebytes (x->x_blackman, x->x_window * sizeof (t_float));
    t_freebytes (x->x_cosine, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hamming, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hann, x->x_window * sizeof (t_float));

    // free x_binFreqs memory
    t_freebytes (x->x_binFreqs, (x->x_windowHalf + 1) * sizeof (t_float));
}


void specRolloff_setup (void)
{
    specRolloff_class =
    class_new (
        gensym ("specRolloff"),
        (t_newmethod)specRolloff_new,
        (t_method)specRolloff_free,
        sizeof (t_specRolloff),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)specRolloff_new,
        gensym ("timbreIDLib/specRolloff"),
        A_GIMME,
        0
    );

    class_addbang (specRolloff_class, specRolloff_bang);

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_chain_fftData,
        gensym ("chain_fftData"),
        A_GIMME,
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_chain_magSpec,
        gensym ("chain_magSpec"),
        A_GIMME,
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_set,
        gensym ("set"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        specRolloff_class,
        (t_method)specRolloff_concentration,
        gensym ("concentration"),
        A_DEFFLOAT,
        0
    );
}
