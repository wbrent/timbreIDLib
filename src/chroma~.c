/*

chroma~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *chroma_tilde_class;

typedef struct _chroma_tilde
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
    t_binIdx *x_binRanges;
    t_float x_loFreq;
    t_float x_hiFreq;
    t_float x_pitchTolerance;
    t_float x_energyThresh;
    t_uChar x_numChroma;
    t_float x_resolution;
    t_float x_microtune;
    t_float *x_pitchClasses;
    double x_lastDspTime;
    t_sample *x_signalBuffer;
    t_sample *x_fftwIn;
    fftwf_complex *x_fftwOut;
    fftwf_plan x_fftwPlan;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_atom *x_listOut;
    t_outlet *x_chroma;
    t_float x_f;
} t_chroma_tilde;


/* ------------------------ chroma~ -------------------------------- */

static void chroma_tilde_bang (t_chroma_tilde *x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float *windowFuncPtr, maxEnergySum, chromaSums[x->x_numChroma];
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

    if (!x->x_powerSpectrum)
        tIDLib_mag (windowHalf + 1, x->x_fftwIn);

    maxEnergySum = 0.0;

    for (i = 0; i < x->x_numChroma; i++)
    {
        t_uInt cardinality;
        chromaSums[i] = 0.0;

        cardinality = tIDLib_getPitchBinRanges(x->x_binRanges, x->x_pitchClasses[i], x->x_loFreq, x->x_hiFreq, x->x_pitchTolerance, x->x_window, x->x_sr);

        //post ("PC[%lu] cardinality: %i", i, cardinality);
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
        //endpost ();

        // divide by the cardinality to account for the fact that different pitch classes will have energy in a different number of bins
        chromaSums[i] /= cardinality;

        if (chromaSums[i] > maxEnergySum)
            maxEnergySum = chromaSums[i];
    }

    // safety to make sure we don't get a division by zero
    maxEnergySum = (maxEnergySum <= 0.0) ? 1.0 : maxEnergySum;

    // neutralize maxEnergySum if not normalizing
    if (!x->x_normalize)
        maxEnergySum = 1.0;

    for (i = 0; i < x->x_numChroma; i++)
        SETFLOAT (x->x_listOut + i, chromaSums[i] / maxEnergySum);

    outlet_list (x->x_chroma, 0, x->x_numChroma, x->x_listOut);
}


static void chroma_tilde_print(t_chroma_tilde *x)
{
    t_uChar i;

    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_sampIdx)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
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

    endpost ();
}


static void chroma_tilde_window (t_chroma_tilde *x, t_floatarg w)
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


static void chroma_tilde_overlap (t_chroma_tilde *x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1.0) ? 1.0 : o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void chroma_tilde_windowFunction (t_chroma_tilde *x, t_floatarg f)
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


static void chroma_tilde_normalize (t_chroma_tilde *x, t_floatarg norm)
{
    norm = (norm < 0.0) ? 0.0 : norm;
    norm = (norm > 1.0) ? 1.0 : norm;
    x->x_normalize = norm;

    if (x->x_normalize)
        post ("%s normalization ON.", x->x_objSymbol->s_name);
    else
        post ("%s normalization OFF.", x->x_objSymbol->s_name);
}


static void chroma_tilde_powerSpectrum (t_chroma_tilde *x, t_floatarg spec)
{
    spec = (spec < 0.0) ? 0.0 : spec;
    spec = (spec > 1.0) ? 1.0 : spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void chroma_tilde_pitchTol (t_chroma_tilde *x, t_floatarg tol)
{
    tol = (tol < 0) ? 0 : tol;
    tol = (tol > 0.5) ? 0.5 : tol;
    x->x_pitchTolerance = tol;
}


static void chroma_tilde_freqRange (t_chroma_tilde *x, t_floatarg loFreq, t_floatarg hiFreq)
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


static void chroma_tilde_energyThresh (t_chroma_tilde *x, t_floatarg thresh)
{
    thresh = (thresh < 0) ? 0 : thresh;
    x->x_energyThresh = thresh;
}


static void chroma_tilde_microtune (t_chroma_tilde *x, t_floatarg cents)
{
    t_uChar i;

    cents = (cents < -100.0) ? -100.0 : cents;
    cents = (cents > 100.0) ? 100.0 : cents;

    x->x_microtune = cents / 100.0;

    for (i = 0; i < x->x_numChroma; i++)
        x->x_pitchClasses[i] += x->x_microtune;
}


static void chroma_tilde_resolution (t_chroma_tilde *x, t_symbol *r)
{
    t_uChar i, oldNumChroma;
    t_float basePitch;

    oldNumChroma = x->x_numChroma;

    if (!strcmp (r->s_name, "full"))
    {
        x->x_numChroma = 12;
        x->x_resolution = 1.0;
    }
    else if (!strcmp (r->s_name, "half"))
    {
        x->x_numChroma = 24;
        x->x_resolution = 0.5;
    }
    else if (!strcmp (r->s_name, "third"))
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


static void *chroma_tilde_new (t_symbol *s, int argc, t_atom *argv)
{
    t_chroma_tilde *x = (t_chroma_tilde *)pd_new (chroma_tilde_class);
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
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            x->x_loFreq = atom_getfloat(argv + 1);
            x->x_hiFreq = atom_getfloat(argv + 2);
            x->x_pitchTolerance = atom_getfloat(argv + 3);
            break;

        case 3:
            x->x_window = atom_getfloat(argv);
            if(x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            x->x_loFreq = atom_getfloat(argv + 1);
            x->x_hiFreq = atom_getfloat(argv + 2);
            x->x_pitchTolerance = 0.1;
            break;

        case 2:
            x->x_window = atom_getfloat(argv);
            if(x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            x->x_loFreq = atom_getfloat(argv + 1);
            x->x_hiFreq = 5000.0;
            x->x_pitchTolerance = 0.1;
            break;

        case 1:
            x->x_window = atom_getfloat(argv);
            if(x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            x->x_loFreq = 50.0;
            x->x_hiFreq = 5000.0;
            x->x_pitchTolerance = 0.1;
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_loFreq = 50.0;
            x->x_hiFreq = 5000.0;
            x->x_pitchTolerance = 0.1;
            break;

        default:
            x->x_window = atom_getfloat(argv);
            if(x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            x->x_loFreq = atom_getfloat(argv + 1);
            x->x_hiFreq = atom_getfloat(argv + 2);
            x->x_pitchTolerance = atom_getfloat(argv + 3);
            post ("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
            break;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;

    chroma_tilde_freqRange(x, x->x_loFreq, x->x_hiFreq);

    x->x_pitchTolerance = (x->x_pitchTolerance<0)?0:x->x_pitchTolerance;
    x->x_pitchTolerance = (x->x_pitchTolerance>0.5)?0.5:x->x_pitchTolerance;
    x->x_energyThresh = 0.0;

    x->x_windowHalf = x->x_window * 0.5;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_windowFunction = blackman;
    x->x_normalize = false;
    x->x_powerSpectrum = true;
    x->x_lastDspTime = clock_getlogicaltime();

    x->x_signalBuffer = (t_sample *)t_getbytes ((x->x_window + x->x_n) * sizeof (t_sample));
    x->x_binRanges = (t_binIdx *)t_getbytes (TID_PBINRANGEBUFSIZE * sizeof (t_binIdx));
    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));
    x->x_listOut = (t_atom *)t_getbytes (x->x_numChroma * sizeof (t_atom));

    for (i = 0; i < (x->x_window + x->x_n); i++)
        x->x_signalBuffer[i] = 0.0;

    x->x_blackman = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_cosine = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_hamming = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_hann = (t_float *)t_getbytes (x->x_window * sizeof (t_float));

     // initialize signal windowing functions
    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);

    // set up the FFTW output buffer.
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf + 1);

    // DFT plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_window; i++)
        x->x_fftwIn[i] = 0.0;

    return (x);
}


static t_int *chroma_tilde_perform (t_int *w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_chroma_tilde *x = (t_chroma_tilde *)(w[1]);

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


static void chroma_tilde_dsp (t_chroma_tilde *x, t_signal **sp)
{
    dsp_add (
        chroma_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if(sp[0]->s_sr != x->x_sr * x->x_overlap)
    {
        x->x_sr = sp[0]->s_sr / x->x_overlap;
        x->x_lastDspTime = clock_getlogicaltime();
    };

// compare n to stored n and update/resize buffer if different
    if(sp[0]->s_n != x->x_n)
    {
        t_sampIdx i;

        x->x_signalBuffer = (t_sample *)t_resizebytes (x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample), (x->x_window + sp[0]->s_n) * sizeof (t_sample));

        x->x_n = sp[0]->s_n;
        x->x_lastDspTime = clock_getlogicaltime();

        // init signal buffer
        for (i = 0; i < (x->x_window + x->x_n); i++)
            x->x_signalBuffer[i] = 0.0;

        post ("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    };
};


static void chroma_tilde_free (t_chroma_tilde *x)
{
    // free the input buffer memory
    t_freebytes(x->x_signalBuffer, (x->x_window + x->x_n) * sizeof (t_sample));

    // free the list out and pitch class memory
    t_freebytes(x->x_listOut, x->x_numChroma * sizeof (t_atom));
    t_freebytes(x->x_pitchClasses, x->x_numChroma * sizeof (t_float));

    // free the bin ranges memory
    t_freebytes(x->x_binRanges, TID_PBINRANGEBUFSIZE * sizeof (t_binIdx));

    // free FFTW stuff
    t_freebytes(x->x_fftwIn, (x->x_window) * sizeof (t_sample));
    fftwf_free(x->x_fftwOut);
    fftwf_destroy_plan(x->x_fftwPlan);

    // free the window memory
    t_freebytes(x->x_blackman, x->x_window * sizeof (t_float));
    t_freebytes(x->x_cosine, x->x_window * sizeof (t_float));
    t_freebytes(x->x_hamming, x->x_window * sizeof (t_float));
    t_freebytes(x->x_hann, x->x_window * sizeof (t_float));
}

void chroma_tilde_setup (void)
{
    chroma_tilde_class =
    class_new(
        gensym("chroma~"),
        (t_newmethod)chroma_tilde_new,
        (t_method)chroma_tilde_free,
        sizeof (t_chroma_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)chroma_tilde_new,
        gensym("timbreIDLib/chroma~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN(chroma_tilde_class, t_chroma_tilde, x_f);

    class_addbang(chroma_tilde_class, chroma_tilde_bang);

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_print,
        gensym("print"),
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_overlap,
        gensym("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_windowFunction,
        gensym("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_normalize,
        gensym("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_powerSpectrum,
        gensym("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_dsp,
        gensym("dsp"),
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_pitchTol,
        gensym("pitch_tolerance"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_freqRange,
        gensym("freq_range"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_energyThresh,
        gensym("thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_microtune,
        gensym("microtune"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        chroma_tilde_class,
        (t_method)chroma_tilde_resolution,
        gensym("resolution"),
        A_DEFSYMBOL,
        0
    );
}
