 /*

chroma

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* chroma_class;

typedef struct _chroma
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_bool x_normalize;
    t_bool x_powerSpectrum;
    t_binIdx* x_binRanges;
    t_float x_loFreq;
    t_float x_hiFreq;
    t_float x_pitchTolerance;
    t_float x_energyThresh;
    t_uChar x_numChroma;
    t_float x_resolution;
    t_float x_microtune;
    t_float* x_pitchClasses;
    t_sample* x_fftwIn;
    fftwf_complex* x_fftwOut;
    fftwf_plan x_fftwPlan;
    t_float* x_blackman;
    t_float* x_cosine;
    t_float* x_hamming;
    t_float* x_hann;
    t_word* x_vec;
    t_symbol* x_arrayName;
    t_sampIdx x_arrayPoints;
    t_atom* x_listOut;
    t_outlet* x_chroma;
} t_chroma;


/* ------------------------ chroma -------------------------------- */
static void chroma_resizeWindow (t_chroma* x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx* endSamp)
{
    t_sampIdx windowHalf;

    windowHalf = window * 0.5;

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
}


static void chroma_analyze (t_chroma* x, t_floatarg start, t_floatarg n)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, window, startSamp, endSamp;
        t_float* windowFuncPtr, maxEnergySum, chromaSums[x->x_numChroma];

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
            chroma_resizeWindow (x, x->x_window, window, startSamp, &endSamp);

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

        maxEnergySum = 0.0;

        for (i = 0; i < x->x_numChroma; i++)
        {
            t_uInt cardinality;
            chromaSums[i] = 0.0;

            cardinality = tIDLib_getPitchBinRanges (x->x_binRanges, x->x_pitchClasses[i], x->x_loFreq, x->x_hiFreq, x->x_pitchTolerance, x->x_window, x->x_sr);

            //startpost ("chroma %lu bin ranges ", i);

            j = 0;

            while (x->x_binRanges[j] != ULONG_MAX)
            {
                t_binIdx k, numBins;

                //startpost ("%lu ", x->x_binRanges[j]);
                //startpost ("%lu ", x->x_binRanges[j+1]);

                // just in case j+1 is somehow ULONG_MAX, abort
                if (x->x_binRanges[j + 1] == ULONG_MAX)
                    break;

                numBins = x->x_binRanges[j + 1] - x->x_binRanges[j] + 1;

                // sum all the energy in the binRange for thisPitch
                for (k = 0; k < numBins; k++)
                {
                    t_float thisEnergy;
                    thisEnergy = x->x_fftwIn[x->x_binRanges[j] + k];

                    if (thisEnergy >= x->x_energyThresh)
                        chromaSums[i] += thisEnergy;
                }

                j += 2;
            }
            //endpost();

            // divide by the cardinality to account for the fact that different pitch classes will have energy in a different number of bins
            chromaSums[i] /= cardinality;

            if (chromaSums[i] > maxEnergySum)
                maxEnergySum = chromaSums[i];
        }

        // safety to make sure we don't get a division by zero
        maxEnergySum = (maxEnergySum <= 0.0) ? 1.0 : maxEnergySum;

        // neutralize maxEnergySum if not normalizing
        if ( !x->x_normalize)
            maxEnergySum = 1.0;

        for (i = 0; i < x->x_numChroma; i++)
            SETFLOAT (x->x_listOut + i, chromaSums[i] / maxEnergySum);

        outlet_list (x->x_chroma, 0, x->x_numChroma, x->x_listOut);
    }
}


static void chroma_chain_fftData (t_chroma* x, t_symbol* s, int argc, t_atom* argv)
{
    t_sampIdx i, j, windowHalf;
    t_float maxEnergySum, chromaSums[x->x_numChroma];

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

    maxEnergySum = 0.0;

    for (i = 0; i < x->x_numChroma; i++)
    {
        t_uInt cardinality;
        chromaSums[i] = 0.0;

        cardinality = tIDLib_getPitchBinRanges (x->x_binRanges, x->x_pitchClasses[i], x->x_loFreq, x->x_hiFreq, x->x_pitchTolerance, x->x_window, x->x_sr);
        //startpost ("chroma %lu bin ranges ", i);

        j = 0;

        while (x->x_binRanges[j] != ULONG_MAX)
        {
            t_binIdx k, numBins;

            //startpost ("%lu ", x->x_binRanges[j]);
            //startpost ("%lu ", x->x_binRanges[j+1]);

            // just in case j+1 is somehow ULONG_MAX, abort
            if (x->x_binRanges[j + 1] == ULONG_MAX)
                break;

            numBins = x->x_binRanges[j + 1] - x->x_binRanges[j] + 1;

            // sum all the energy in the binRange for thisPitch
            for (k = 0; k < numBins; k++)
            {
                t_float thisEnergy;
                thisEnergy = x->x_fftwIn[x->x_binRanges[j] + k];

                if (thisEnergy >= x->x_energyThresh)
                    chromaSums[i] += thisEnergy;
            }

            j += 2;
        }
        //endpost();

        // divide by the cardinality to account for the fact that different pitch classes will have energy in a different number of bins
        chromaSums[i] /= cardinality;

        if (chromaSums[i] > maxEnergySum)
            maxEnergySum = chromaSums[i];
    }

    // safety to make sure we don't get a division by zero
    maxEnergySum = (maxEnergySum <= 0.0) ? 1.0 : maxEnergySum;

    // neutralize maxEnergySum if not normalizing
    if ( !x->x_normalize)
        maxEnergySum = 1.0;

    for (i = 0; i < x->x_numChroma; i++)
        SETFLOAT (x->x_listOut + i, chromaSums[i] / maxEnergySum);

    outlet_list (x->x_chroma, 0, x->x_numChroma, x->x_listOut);
}


static void chroma_chain_magSpec (t_chroma* x, t_symbol* s, int argc, t_atom* argv)
{
    t_sampIdx i, j, windowHalf;
    t_float maxEnergySum, chromaSums[x->x_numChroma];

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

    maxEnergySum = 0.0;

    for (i = 0; i < x->x_numChroma; i++)
    {
        t_uInt cardinality;
        chromaSums[i] = 0.0;

        cardinality = tIDLib_getPitchBinRanges (x->x_binRanges, x->x_pitchClasses[i], x->x_loFreq, x->x_hiFreq, x->x_pitchTolerance, x->x_window, x->x_sr);
        //startpost ("chroma %lu bin ranges ", i);

        j = 0;

        while (x->x_binRanges[j] != ULONG_MAX)
        {
            t_binIdx k, numBins;

            //startpost ("%lu ", x->x_binRanges[j]);
            //startpost ("%lu ", x->x_binRanges[j+1]);

            // just in case j+1 is somehow ULONG_MAX, abort
            if (x->x_binRanges[j + 1] == ULONG_MAX)
                break;

            numBins = x->x_binRanges[j + 1] - x->x_binRanges[j] + 1;

            // sum all the energy in the binRange for thisPitch
            for (k = 0; k < numBins; k++)
            {
                t_float thisEnergy;
                thisEnergy = x->x_fftwIn[x->x_binRanges[j] + k];

                if (thisEnergy >= x->x_energyThresh)
                    chromaSums[i] += thisEnergy;
            }

            j += 2;
        }
        //endpost();

        // divide by the cardinality to account for the fact that different pitch classes will have energy in a different number of bins
        chromaSums[i] /= cardinality;

        if (chromaSums[i] > maxEnergySum)
            maxEnergySum = chromaSums[i];
    }

    // safety to make sure we don't get a division by zero
    maxEnergySum = (maxEnergySum <= 0.0) ? 1.0 : maxEnergySum;

    // neutralize maxEnergySum if not normalizing
    if ( !x->x_normalize)
        maxEnergySum = 1.0;

    for (i = 0; i < x->x_numChroma; i++)
        SETFLOAT (x->x_listOut + i, chromaSums[i] / maxEnergySum);

    outlet_list (x->x_chroma, 0, x->x_numChroma, x->x_listOut);
}


// analyze the whole damn array
static void chroma_bang (t_chroma* x)
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
        chroma_analyze (x, startSamp, window);
    }
}


static void chroma_set (t_chroma* x, t_symbol* s)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (s, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void chroma_print (t_chroma* x)
{
    t_uChar i;

    post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post ("%s frequency range: %0.2f to %0.2fHz", x->x_objSymbol->s_name, x->x_loFreq, x->x_hiFreq);
    post ("%s pitch tolerance: %f", x->x_objSymbol->s_name, x->x_pitchTolerance);
    post ("%s spectral energy threshold: %f", x->x_objSymbol->s_name, x->x_energyThresh);
    post ("%s resolution: %i", x->x_objSymbol->s_name, x->x_numChroma);

    startpost ("%s Bass MIDI pitches: ", x->x_objSymbol->s_name);

    for (i = 0; i < x->x_numChroma; i++)
        startpost ("%0.2f ", x->x_pitchClasses[i]);

    endpost();
}


static void chroma_samplerate (t_chroma* x, t_floatarg sr)
{
    if (sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void chroma_window (t_chroma* x, t_floatarg w)
{
    t_sampIdx endSamp;

    // catch negative window size arguments here while data type is still t_floatarg. once passed to _resizeWindow, it will be cast to type t_sampIdx and negative values will turn to garbage.
    w = (w < 0.0) ? TID_MINWINDOWSIZE : w;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow() requires that
    endSamp = 0;

    chroma_resizeWindow (x, x->x_window, w, 0, &endSamp);
}


static void chroma_windowFunction (t_chroma* x, t_floatarg f)
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


static void chroma_powerSpectrum (t_chroma* x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void chroma_normalize (t_chroma* x, t_floatarg norm)
{
    norm = (norm < 0) ? 0 : norm;
    norm = (norm > 1) ? 1 : norm;
    x->x_normalize = norm;

    if (x->x_normalize)
        post ("%s normalization ON.", x->x_objSymbol->s_name);
    else
        post ("%s normalization OFF.", x->x_objSymbol->s_name);
}


static void chroma_pitchTol (t_chroma* x, t_floatarg tol)
{
    tol = (tol < 0) ? 0 : tol;
    tol = (tol > 0.5) ? 0.5 : tol;
    x->x_pitchTolerance = tol;
}


static void chroma_freqRange (t_chroma* x, t_floatarg loFreq, t_floatarg hiFreq)
{
    t_float nyquist;

    nyquist = x->x_sr * 0.5;

    loFreq = (loFreq < 0) ? 0 : loFreq;
    loFreq = (loFreq > nyquist) ? nyquist : loFreq;

    hiFreq = (hiFreq < 0) ? 0 : hiFreq;
    hiFreq = (hiFreq > nyquist) ? nyquist : hiFreq;

    if (hiFreq < loFreq)
    {
        t_float tmpFreq;
        tmpFreq = hiFreq;
        hiFreq = loFreq;
        loFreq = tmpFreq;
    }

    x->x_loFreq = loFreq;
    x->x_hiFreq = hiFreq;
}


static void chroma_energyThresh (t_chroma* x, t_floatarg thresh)
{
    thresh = (thresh < 0) ? 0 : thresh;
    x->x_energyThresh = thresh;
}


static void chroma_microtune (t_chroma* x, t_floatarg cents)
{
    t_uChar i;

    cents = (cents < -100.0) ? -100.0 : cents;
    cents = (cents > 100.0) ? 100.0 : cents;

    x->x_microtune = cents / 100.0;

    for (i = 0; i < x->x_numChroma; i++)
        x->x_pitchClasses[i] += x->x_microtune;
}


static void chroma_resolution (t_chroma* x, t_symbol* r)
{
    t_uChar i, oldNumChroma;
    t_float basePitch;

    oldNumChroma = x->x_numChroma;

    if ( !strcmp (r->s_name, "full"))
    {
        x->x_numChroma = 12;
        x->x_resolution = 1.0;
    }
    else if ( !strcmp (r->s_name, "half"))
    {
        x->x_numChroma = 24;
        x->x_resolution = 0.5;
    }
    else if ( !strcmp (r->s_name, "third"))
    {
        x->x_numChroma = 36;
        x->x_resolution = 0.333333;
    }
    else
    {
        x->x_numChroma = 12;
        x->x_resolution = 1.0;
        post ("%s resolution must be either \"full\" (12), \"half\" (24), or \"third\" (36). Using default of %i", x->x_objSymbol->s_name, x->x_numChroma);
    }

    // resize the pitch class and list out memory since x_numChroma changed
    x->x_pitchClasses = (t_float *)t_resizebytes (x->x_pitchClasses, oldNumChroma * sizeof (t_float), x->x_numChroma * sizeof (t_float));
    x->x_listOut = (t_atom *)t_resizebytes (x->x_listOut, oldNumChroma * sizeof (t_atom), x->x_numChroma * sizeof (t_atom));

    // beginning of pitch class array starts at lowest C on piano
    basePitch = 24 + x->x_microtune;

    for (i = 0; i < x->x_numChroma - (3 / x->x_resolution); i++)
    {
        x->x_pitchClasses[i] = basePitch;
        basePitch += x->x_resolution;
    }

    // fill out the last part starting with lowest A on piano
    basePitch = 21 + x->x_microtune;

    for (; i < x->x_numChroma; i++)
    {
        x->x_pitchClasses[i] = basePitch;
        basePitch += x->x_resolution;
    }
}


static void* chroma_new (t_symbol* s, int argc, t_atom* argv)
{
    t_chroma* x = (t_chroma *)pd_new (chroma_class);
    t_sampIdx i;

    x->x_chroma = outlet_new (&x->x_obj, gensym ("list"));

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    x->x_numChroma = 12;
    x->x_resolution = 1.0;
    x->x_microtune = 0.0;

      x->x_pitchClasses = (t_float *)t_getbytes (x->x_numChroma * sizeof (t_float));

    // these are the lowest MIDI pitches being considered for C through B
    x->x_pitchClasses[0] = 24.0;
    x->x_pitchClasses[1] = 25.0;
    x->x_pitchClasses[2] = 26.0;
    x->x_pitchClasses[3] = 27.0;
    x->x_pitchClasses[4] = 28.0;
    x->x_pitchClasses[5] = 29.0;
    x->x_pitchClasses[6] = 30.0;
    x->x_pitchClasses[7] = 31.0;
    x->x_pitchClasses[8] = 32.0;
    x->x_pitchClasses[9] = 21.0;
    x->x_pitchClasses[10] = 22.0;
    x->x_pitchClasses[11] = 23.0;

    switch (argc)
    {
        case 4:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_loFreq = atom_getfloat (argv + 1);
            x->x_hiFreq = atom_getfloat (argv + 2);
            x->x_pitchTolerance = atom_getfloat (argv + 3);
            break;

        case 3:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_loFreq = atom_getfloat (argv + 1);
            x->x_hiFreq = atom_getfloat (argv + 2);
            x->x_pitchTolerance = 0.1;
            break;

        case 2:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_loFreq = atom_getfloat (argv + 1);
            x->x_hiFreq = 5000.0;
            x->x_pitchTolerance = 0.1;
            break;

        case 1:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_loFreq = 50.0;
            x->x_hiFreq = 5000.0;
            x->x_pitchTolerance = 0.1;
            break;

        case 0:
            post ("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze()
            x->x_arrayName = gensym ("NOARRAYSPECIFIED");
            x->x_loFreq = 50.0;
            x->x_hiFreq = 5000.0;
            x->x_pitchTolerance = 0.1;
            break;

        default:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_loFreq = atom_getfloat (argv + 1);
            x->x_hiFreq = atom_getfloat (argv + 2);
            x->x_pitchTolerance = atom_getfloat (argv + 3);
            post ("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
            break;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;

    chroma_freqRange (x, x->x_loFreq, x->x_hiFreq);

    x->x_pitchTolerance = (x->x_pitchTolerance < 0) ? 0 : x->x_pitchTolerance;
    x->x_pitchTolerance = (x->x_pitchTolerance > 0.5) ? 0.5 : x->x_pitchTolerance;
    x->x_energyThresh = 0.0;

    x->x_window = TID_WINDOWSIZEDEFAULT;
    x->x_windowHalf = x->x_window * 0.5;
    x->x_windowFunction = blackman;
    x->x_normalize = false;
    x->x_powerSpectrum = true;

    // make the bin ranges memory 20x the window size so that we can find an upper and lower bin for more than 8 octaves of a given pitch class
    x->x_binRanges = (t_binIdx *)t_getbytes (TID_PBINRANGEBUFSIZE * sizeof (t_binIdx));
    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));
    x->x_listOut = (t_atom *)t_getbytes (x->x_numChroma * sizeof (t_atom));

    // set up the FFTW output buffer. Is there no function to initialize it?
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);

    // FFTW plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG); // FFTWPLANNERFLAG may be slower than FFTWPLANNERFLAG but more efficient after the first run?

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


static void chroma_free (t_chroma* x)
{
    // free the list out and pitch class memory
    t_freebytes (x->x_listOut, x->x_numChroma * sizeof (t_atom));
    t_freebytes (x->x_pitchClasses, x->x_numChroma * sizeof (t_float));

    // free the bin ranges memory
    t_freebytes (x->x_binRanges, TID_PBINRANGEBUFSIZE * sizeof (t_binIdx));

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


void chroma_setup (void)
{
    chroma_class =
    class_new (
        gensym ("chroma"),
        (t_newmethod)chroma_new,
        (t_method)chroma_free,
        sizeof (t_chroma),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)chroma_new,
        gensym ("timbreIDLib/chroma"),
        A_GIMME,
        0
    );

    class_addbang (chroma_class, chroma_bang);

    class_addmethod (
        chroma_class,
        (t_method)chroma_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_chain_fftData,
        gensym ("chain_fftData"),
        A_GIMME,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_chain_magSpec,
        gensym ("chain_magSpec"),
        A_GIMME,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_set,
        gensym ("set"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_normalize,
        gensym ("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_pitchTol,
        gensym ("pitch_tolerance"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_freqRange,
        gensym ("freq_range"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_energyThresh,
        gensym ("thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_microtune,
        gensym ("microtune"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        chroma_class,
        (t_method)chroma_resolution,
        gensym ("resolution"),
        A_DEFSYMBOL,
        0
    );
}
