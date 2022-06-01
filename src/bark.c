/*

bark

Copyright 2010 William Brent

bark is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

bark is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

t_float bark_weights_dB[] = {-69.9, -60.4, -51.4, -43.3, -36.6, -30.3, -24.3, -19.5, -14.8, -10.7, -7.5, -4.8, -2.6, -0.8, 0.0, 0.6, 0.5, 0.0, -0.1, 0.5, 1.5, 3.6, 5.9, 6.5, 4.2, -2.6, -10.2, -10.0, -2.8};

t_float bark_weights_freqs[] = {20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500};

static t_class* bark_class;

typedef struct _bark
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_bool x_debug;
    t_float x_sr;

    t_word* x_vec;
    t_symbol* x_arrayName;
    t_sampIdx x_arrayPoints;

    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_sampIdx x_hop;

    t_windowFunction x_windowFunction;
    t_bool x_powerSpectrum;
    t_bool x_normalize;
    t_bool x_useWeights;
    t_bool x_spew;

    t_float x_barkSpacing;
    t_filterIdx x_numFilters;
    t_filterIdx x_sizeFilterFreqs;
    t_float* x_filterFreqs;
    t_filter* x_filterbank;
    t_bool x_specBandAvg;
    t_bool x_filterAvg;
    t_float* x_loudWeights;

    t_float x_prevTotalGrowth;
    t_float x_loThresh;
    t_float x_hiThresh;
    t_float x_minvel;
    t_filterIdx x_loBin; // these aren't really bins, but filters
    t_filterIdx x_hiBin;
    t_float x_debounceTime;
    t_uInt x_debounceSamp;
    t_float x_maskDecay;
    t_uInt x_maskPeriods;
    t_filterIdx* x_numPeriods; // t_filterIdx type because this buffer is used to check making per filter band

    // this should remain a large integer of some sort, because it is used to keep track of the running sum of samples processed. set to UINT_MAX as off flag
    t_uInt x_debounceActive;
    t_bool x_haveHit;

    t_float* x_blackman;
    t_float* x_cosine;
    t_float* x_hamming;
    t_float* x_hann;

    t_sample* x_fftwIn;
    fftwf_complex* x_fftwOut;
    fftwf_plan x_fftwPlan;

    t_float* x_mask;
    t_float* x_growth;
    t_atom* x_growthList;

    t_outlet* x_outputList;
    t_outlet* x_growthOut;
    t_outlet* x_timeOut;

} t_bark;


/* ---------------- utility functions ---------------------- */

static void bark_create_loudness_weighting (t_bark* x)
{
    t_filterIdx i;
    t_float barkSum, *barkFreqs;

    barkFreqs = (t_float *)t_getbytes (x->x_numFilters * sizeof (t_float));

    barkSum = x->x_barkSpacing;

    for (i = 0; i < x->x_numFilters; i++)
    {
        barkFreqs[i] = tIDLib_bark2freq (barkSum);
        barkSum += x->x_barkSpacing;
    }

    for (i = 0; i < x->x_numFilters; i++)
    {
        t_binIdx nearIdx;
        t_float nearFreq, diffFreq, diffdB, dBint;

        nearIdx = tIDLib_nearestBinIndex (barkFreqs[i], bark_weights_freqs, TID_NUMWEIGHTPOINTS);
        nearFreq = bark_weights_freqs[nearIdx];
        diffdB = 0;

        // this doesn't have to be if/else'd into a greater/less situation.  later on i should write a more general interpolation solution, and maybe move it up to 4 points instead.
        if (barkFreqs[i]>nearFreq)
        {
            if (nearIdx <= TID_NUMWEIGHTPOINTS-2)
            {
                diffFreq = (barkFreqs[i] - nearFreq)/(bark_weights_freqs[nearIdx + 1] - nearFreq);
                diffdB = diffFreq * (bark_weights_dB[nearIdx + 1] - bark_weights_dB[nearIdx]);
            }

            dBint = bark_weights_dB[nearIdx] + diffdB;
        }
        else
        {
            if (nearIdx>0)
            {
                diffFreq = (barkFreqs[i] - bark_weights_freqs[nearIdx - 1])/(nearFreq - bark_weights_freqs[nearIdx - 1]);
                diffdB = diffFreq * (bark_weights_dB[nearIdx] - bark_weights_dB[nearIdx - 1]);
            }

            dBint = bark_weights_dB[nearIdx - 1] + diffdB;
        }

        if (x->x_powerSpectrum)
            x->x_loudWeights[i] = pow (10.0, dBint*0.1);
        else
            x->x_loudWeights[i] = pow (10.0, dBint*0.05);

    }
}

/* ---------------- END utility functions ---------------------- */


/* ------------------------ bark -------------------------------- */

static void bark_windowFunction (t_bark* x, t_floatarg f)
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

static void bark_thresh (t_bark* x, t_floatarg lo, t_floatarg hi)
{
    if (hi < lo)
    {
        post ("%s WARNING: high threshold less than low threshold.", x->x_objSymbol->s_name);
        x->x_hiThresh = lo;
        x->x_loThresh = hi;

        x->x_loThresh = (x->x_loThresh < 0) ? -1 : x->x_loThresh;

//		post ("bark: low thresh: %f, high thresh: %f", x->x_loThresh, x->x_hiThresh);
    }
    else
    {
        x->x_hiThresh = hi;
        x->x_loThresh = lo;

        x->x_loThresh = (x->x_loThresh < 0) ? -1 : x->x_loThresh;

//		post ("bark: low thresh: %0.2f, high thresh: %0.2f", x->x_loThresh, x->x_hiThresh);
    }
}

static void bark_minvel (t_bark* x, t_floatarg mv)
{
    x->x_minvel = (mv < 0) ? 0 : mv;
}

static void bark_filter_range (t_bark* x, t_floatarg lo, t_floatarg hi)
{
    if (hi < lo)
    {
        t_float tmp;

        tmp = hi;
        hi = lo;
        lo = tmp;

        x->x_loBin = (lo < 0) ? 0 : lo;
        x->x_hiBin = (hi >= x->x_numFilters) ? x->x_numFilters - 1 : hi;

        post ("%s WARNING: high band less than low band. Reversing order.", x->x_objSymbol->s_name);
    }
    else
    {
        x->x_loBin = (lo < 0) ? 0 : lo;
        x->x_hiBin = (hi >= x->x_numFilters) ? x->x_numFilters - 1 : hi;
    }
}

static void bark_mask (t_bark* x, t_floatarg per, t_floatarg dec)
{
    x->x_maskPeriods = per;
    x->x_maskDecay = dec;

    x->x_maskDecay = (x->x_maskDecay < 0.05) ? 0.05 : x->x_maskDecay;
    x->x_maskDecay = (x->x_maskDecay > 0.95) ? 0.95 : x->x_maskDecay;

//	post ("bark: masking for %i periods and decaying by %i%% thereafter.", x->x_maskPeriods, (int)(x->x_maskDecay*100));
}

static void bark_debounce (t_bark* x, t_floatarg db)
{
    if (x->x_debounceTime >= 0)
    {
        x->x_debounceTime = db;
        x->x_debounceSamp = x->x_debounceTime * 0.001 * x->x_sr;
        post ("%s debounce time: %0.2f ms, %i samples", x->x_objSymbol->s_name, x->x_debounceTime, x->x_debounceSamp);
    }
    else
        pd_error (x, "%s debounce time must be >= 0 ms", x->x_objSymbol->s_name);
}

static void bark_spew (t_bark* x, t_floatarg spew)
{
    spew = (spew < 0) ? 0 : spew;
    spew = (spew > 1) ? 1 : spew;
    x->x_spew = spew;

    post ("%s spew mode: %i", x->x_objSymbol->s_name, x->x_spew);
}

static void bark_print (t_bark* x)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
        post ("%s total frames in array: %i", x->x_objSymbol->s_name, (int)floor ((x->x_arrayPoints - x->x_window) / x->x_hop));
    }

    post ("%s window size: %i", x->x_objSymbol->s_name, (t_sampIdx)x->x_window);
    post ("%s hop: %i", x->x_objSymbol->s_name, x->x_hop);
    post ("%s Bark spacing: %0.2f", x->x_objSymbol->s_name, x->x_barkSpacing);
    post ("%s number of filters: %i", x->x_objSymbol->s_name, x->x_numFilters);
    post ("%s bin range: %i through %i (inclusive)", x->x_objSymbol->s_name, x->x_loBin, x->x_hiBin);
    post ("%s low thresh: %0.2f, high thresh: %0.2f", x->x_objSymbol->s_name, x->x_loThresh, x->x_hiThresh);
    post ("%s minvel: %f", x->x_objSymbol->s_name, x->x_minvel);
    post ("%s mask periods: %i, mask decay: %0.2f", x->x_objSymbol->s_name, x->x_maskPeriods, x->x_maskDecay);
    post ("%s debounce time: %0.2f", x->x_objSymbol->s_name, x->x_debounceTime);
    post ("%s normalization: %i", x->x_objSymbol->s_name, x->x_normalize);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s spectrum band averaging: %i", x->x_objSymbol->s_name, x->x_specBandAvg);
    post ("%s triangular filter averaging: %i", x->x_objSymbol->s_name, x->x_filterAvg);
    post ("%s loudness weights: %i", x->x_objSymbol->s_name, x->x_useWeights);
    post ("%s spew mode: %i", x->x_objSymbol->s_name, x->x_spew);
    post ("%s debug mode: %i", x->x_objSymbol->s_name, x->x_debug);
}


static void bark_debug (t_bark* x, t_floatarg debug)
{
    debug = (debug < 0) ? 0 : debug;
    debug = (debug > 1) ? 1 : debug;
    x->x_debug = debug;

    if (x->x_debug)
        post ("%s debug mode ON", x->x_objSymbol->s_name);
    else
        post ("%s debug mode OFF", x->x_objSymbol->s_name);
}

static void bark_use_weights (t_bark* x, t_floatarg w)
{
    w = (w < 0) ? 0 : w;
    w = (w > 1) ? 1 : w;
    x->x_useWeights = w;

    if (x->x_useWeights)
        post ("%s using loudness weighting", x->x_objSymbol->s_name);
    else
        post ("%s using unweighted spectrum", x->x_objSymbol->s_name);
}


static void bark_spec_band_avg (t_bark* x, t_floatarg avg)
{
    avg = (avg < 0) ? 0 : avg;
    avg = (avg > 1) ? 1 : avg;
    x->x_specBandAvg = avg;

    if (x->x_specBandAvg)
        post ("%s: averaging energy in spectrum bands.", x->x_objSymbol->s_name);
    else
        post ("%s: using triangular filterbank.", x->x_objSymbol->s_name);
}


static void bark_filter_avg (t_bark* x, t_floatarg avg)
{
    avg = (avg < 0) ? 0 : avg;
    avg = (avg > 1) ? 1 : avg;
    x->x_filterAvg = avg;

    if (x->x_filterAvg)
        post ("%s: averaging energy in triangular filters.", x->x_objSymbol->s_name);
    else
        post ("%s: summing energy in triangular filters.", x->x_objSymbol->s_name);
}


static void bark_powerSpectrum (t_bark* x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    bark_create_loudness_weighting (x);

    if (x->x_powerSpectrum)
        post ("%s: using power spectrum.", x->x_objSymbol->s_name);
    else
        post ("%s: using magnitude spectrum.", x->x_objSymbol->s_name);
}


static void bark_normalize (t_bark* x, t_floatarg norm)
{
    norm = (norm < 0) ? 0 : norm;
    norm = (norm > 1) ? 1 : norm;
    x->x_normalize = norm;

    if (x->x_normalize)
        post ("%s: spectrum normalization ON.", x->x_objSymbol->s_name);
    else
        post ("%s: spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void bark_samplerate (t_bark* x, t_floatarg sr)
{
    t_filterIdx i, oldNumFilters;

    oldNumFilters = x->x_numFilters;

    if (sr < TID_MINSAMPLERATE)
    {
        pd_error (x, "%s: samplerate must be at least %i. default value of %i used instead.", x->x_objSymbol->s_name, TID_MINSAMPLERATE, TID_SAMPLERATEDEFAULT);
        x->x_sr = TID_SAMPLERATEDEFAULT;
    }
    else
        x->x_sr = sr;

    x->x_debounceSamp = x->x_debounceTime*0.001 * x->x_sr;

    x->x_sizeFilterFreqs = tIDLib_getBarkBoundFreqs (&x->x_filterFreqs, x->x_sizeFilterFreqs, x->x_barkSpacing, x->x_sr);

    // sizeFilterFreqs - 2 is the correct number of filters, since we don't count the start point of the first filter, or the finish point of the last filter
    x->x_numFilters = x->x_sizeFilterFreqs - 2;

    tIDLib_createFilterbank (x->x_filterFreqs, &x->x_filterbank, oldNumFilters, x->x_numFilters, x->x_window, x->x_sr);

    x->x_loBin = 0;
    x->x_hiBin = x->x_numFilters - 1;

    x->x_mask = (t_float *)t_resizebytes (x->x_mask, oldNumFilters * sizeof (t_float), x->x_numFilters * sizeof (t_float));
    x->x_growth = (t_float *)t_resizebytes (x->x_growth, oldNumFilters * sizeof (t_float), x->x_numFilters * sizeof (t_float));
    x->x_numPeriods = (t_filterIdx *)t_resizebytes (x->x_numPeriods, oldNumFilters * sizeof (t_filterIdx), x->x_numFilters * sizeof (t_filterIdx));
    x->x_growthList = (t_atom *)t_resizebytes (x->x_growthList, oldNumFilters * sizeof (t_atom), x->x_numFilters * sizeof (t_atom));
    x->x_loudWeights = (t_float *)t_resizebytes (x->x_loudWeights, oldNumFilters * sizeof (t_float), x->x_numFilters * sizeof (t_float));

    for (i = 0; i < x->x_numFilters; i++)
    {
        x->x_mask[i] = 0.0;
        x->x_growth[i] = 0.0;
        x->x_numPeriods[i] = 0.0;
        SETFLOAT (x->x_growthList + i, 0.0);
        x->x_loudWeights[i] = 0.0;
    }
}


static void* bark_new (t_symbol* s, int argc, t_atom* argv)
{
    t_bark* x = (t_bark *)pd_new (bark_class);
    t_sampIdx i;
//	t_garray *a;

    x->x_timeOut = outlet_new (&x->x_obj, &s_float);
    x->x_growthOut = outlet_new (&x->x_obj, &s_float);
    x->x_outputList = outlet_new (&x->x_obj, gensym ("list"));

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 4:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_window = atom_getfloat (argv + 1);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_hop = atom_getfloat (argv + 2);

            if (x->x_hop<1)
            {
                post ("%s WARNING: requested hop is less than 1 sample. Quarter of window size (%i samples) used instead.", x->x_objSymbol->s_name, x->x_hop);
                x->x_hop = x->x_window*0.25;
            }

            x->x_barkSpacing = atom_getfloat (argv + 3);
            if (x->x_barkSpacing < TID_MINBARKSPACING || x->x_barkSpacing > TID_MAXBARKSPACING)
            {
                x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
                post ("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, TID_MINBARKSPACING, TID_MAXBARKSPACING, TID_BARKSPACINGDEFAULT);
            }
            break;

        case 3:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_window = atom_getfloat (argv + 1);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_hop = atom_getfloat (argv + 2);

            if (x->x_hop<1)
            {
                x->x_hop = x->x_window*0.25;
                post ("%s WARNING: requested hop is less than 1 sample. Quarter of window size (%i samples) used instead.", x->x_objSymbol->s_name, x->x_hop);
            }

            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            break;

        case 2:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_window = atom_getfloat (argv + 1);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_hop = x->x_window*0.25;
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            break;

        case 1:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_hop = TID_WINDOWSIZEDEFAULT*0.25;
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            break;

        case 0:
            pd_error (x, "%s: no array specified.", x->x_objSymbol->s_name);
            x->x_arrayName = gensym ("NOARRAYSPECIFIED");
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_hop = TID_WINDOWSIZEDEFAULT*0.25;
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            break;

        default:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post ("%s WARNING: Too many arguments supplied. Using default window size of %i, hop of %i, and Bark spacing of %0.2f.", x->x_objSymbol->s_name, TID_WINDOWSIZEDEFAULT, (int)(TID_WINDOWSIZEDEFAULT*0.25), TID_BARKSPACINGDEFAULT);
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_hop = TID_WINDOWSIZEDEFAULT*0.25;
            x->x_barkSpacing = TID_BARKSPACINGDEFAULT;
            break;
    }

    x->x_windowHalf = x->x_window * 0.5;
    x->x_debug = false;
    x->x_spew = false;
    x->x_sr = 44100.0;
    x->x_windowFunction = blackman;
    x->x_powerSpectrum = true;
    x->x_normalize = false;
    x->x_loThresh = 3;
    x->x_hiThresh = 7;
    x->x_minvel = 1.0;
    x->x_haveHit = false;
    x->x_debounceTime = 100;
    x->x_debounceSamp = x->x_debounceTime*0.001 * x->x_sr;
    x->x_debounceActive = 0;
    x->x_maskDecay = 0.7;
    x->x_maskPeriods = 4;
    x->x_numFilters = 0;
    x->x_prevTotalGrowth = 0.0;
    x->x_useWeights = false;
    x->x_specBandAvg = false;
    x->x_filterAvg = false;

    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));

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

    x->x_loBin = 0;
    x->x_hiBin = x->x_numFilters - 1;

    x->x_mask = (t_float *)t_getbytes (x->x_numFilters * sizeof (t_float));
    x->x_growth = (t_float *)t_getbytes (x->x_numFilters * sizeof (t_float));
    x->x_numPeriods = (t_filterIdx *)t_getbytes (x->x_numFilters * sizeof (t_filterIdx));
    x->x_growthList = (t_atom *)t_getbytes (x->x_numFilters * sizeof (t_atom));
    x->x_loudWeights = (t_float *)t_getbytes (x->x_numFilters * sizeof (t_float));

    for (i = 0; i < x->x_numFilters; i++)
    {
        x->x_mask[i] = 0.0;
        x->x_growth[i] = 0.0;
        x->x_numPeriods[i] = 0.0;
        SETFLOAT (x->x_growthList + i, 0.0);
        x->x_loudWeights[i] = 0.0;
    }

    bark_create_loudness_weighting (x);

    return (x);
}


static void bark_analyze (t_bark* x, t_floatarg startTime, t_floatarg endTime)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, window, windowHalf, hop, nFrames, frame, sampRange, startSamp, endSamp;
        t_float totalGrowth, totalVel, *windowFuncPtr;

        nFrames = 0;

        window = x->x_window;
        windowHalf = x->x_windowHalf;
        hop = x->x_hop;

        if (startTime>endTime)
        {
            pd_error (x, "%s: invalid time range", x->x_objSymbol->s_name);
            return;
        }
        else if (startTime<endTime)
        {
            startTime = (startTime < 0.0) ? 0.0 : startTime;
            endTime = (endTime < 0.0) ? 0.0 : endTime;

            startSamp = floor (startTime * x->x_sr);
            endSamp = floor (endTime * x->x_sr);

            if (endSamp < x->x_arrayPoints)
                sampRange = endSamp - startSamp + 1;
            else
            {
                pd_error (x, "%s: invalid time range", x->x_objSymbol->s_name);
                return;
            }
        }
        else
        {
            sampRange = x->x_arrayPoints;
            startSamp = 0;
            endSamp = x->x_arrayPoints - 1;
        }

        nFrames = floor ((sampRange-window)/hop);

        // init mask to zero
        for (i = 0; i < x->x_numFilters; i++)
            x->x_mask[i] = 0.0;

        for (frame=0; frame<nFrames; frame++)
        {
            // fill buffer with <window> samples
            for (i = 0, j = frame * hop + startSamp; i < window; i++, j++)
                x->x_fftwIn[i] = x->x_vec[j].w_float;

            totalGrowth = 0.0;
            totalVel = 0.0;

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
                tIDLib_specFilterBands (windowHalf + 1, x->x_numFilters, x->x_fftwIn, x->x_filterbank, x->x_normalize);
            else
                tIDLib_filterbankMultiply (x->x_fftwIn, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

            // optional loudness weighting
            if (x->x_useWeights)
                for (i = 0; i < x->x_numFilters; i++)
                    x->x_fftwIn[i] *= x->x_loudWeights[i];

            for (i = 0; i < x->x_numFilters; i++)
                totalVel += x->x_fftwIn[i];

            // init growth list to zero
            for (i = 0; i < x->x_numFilters; i++)
                x->x_growth[i] = 0.0;

            for (i = 0; i < x->x_numFilters; i++)
            {
                // from p.3 of Puckette/Apel/Zicarelli, 1998
                // salt divisor with + 1.0e-15 in case previous power was zero
                if (x->x_fftwIn[i] > x->x_mask[i])
                    x->x_growth[i] = x->x_fftwIn[i] / (x->x_mask[i] + 1.0e-15) - 1.0;

                if (i>=x->x_loBin && i <= x->x_hiBin && x->x_growth[i]>0)
                    totalGrowth += x->x_growth[i];

                SETFLOAT (x->x_growthList + i, x->x_growth[i]);
            }

            if (frame*hop+startSamp >= x->x_debounceActive)
                x->x_debounceActive = UINT_MAX;

            if (totalVel >= x->x_minvel && totalGrowth > x->x_hiThresh && !x->x_haveHit && x->x_debounceActive==UINT_MAX)
            {
                if (x->x_debug)
                    post ("%s peak: %f", x->x_objSymbol->s_name, totalGrowth);

                x->x_haveHit = true;
                x->x_debounceActive = frame*hop+startSamp + x->x_debounceSamp;
            }
            else if (x->x_haveHit && x->x_loThresh>0 && totalGrowth < x->x_loThresh) // if loThresh is an actual value (not -1), then wait until growth drops below that value before reporting attack
            {
                if (x->x_debug)
                    post ("%s drop: %f", x->x_objSymbol->s_name, totalGrowth);

                x->x_haveHit = false;

                // don't output data if spew will do it anyway below
                if ( !x->x_spew)
                {
                    outlet_list (x->x_outputList, 0, x->x_numFilters, x->x_growthList);
                    outlet_float (x->x_growthOut, totalGrowth);
                }

                // add half a window of samples as a fudge factor. note that since this NRT and we can look into the future, all attack reports will be roughly a half window EARLY.  in RT, everything is a half window LATE because the point of reference is the END of the window.  here, it's the BEGINNING of the window.
                outlet_float (x->x_timeOut, (frame*hop+startSamp + windowHalf) / x->x_sr);
            }
            else if (x->x_haveHit && x->x_loThresh<0 && totalGrowth < x->x_prevTotalGrowth) // if loThresh == -1, report attack as soon as growth shows any decay at all
            {
                if (x->x_debug)
                    post ("%s drop: %f", x->x_objSymbol->s_name, totalGrowth);

                x->x_haveHit = false;

                // don't output data if spew will do it anyway below
                if ( !x->x_spew)
                {
                    outlet_list (x->x_outputList, 0, x->x_numFilters, x->x_growthList);
                    outlet_float (x->x_growthOut, totalGrowth);
                }

                // add half a window of samples as a fudge factor. note that since this NRT and we can look into the future, all attack reports will be roughly a half window EARLY.  in RT, everything is a half window LATE because the point of reference is the END of the window.  here, it's the BEGINNING of the window.
                outlet_float (x->x_timeOut, (frame*hop+startSamp + windowHalf) / x->x_sr);
            }

            if (x->x_spew)
            {
                outlet_list (x->x_outputList, 0, x->x_numFilters, x->x_growthList);
                outlet_float (x->x_growthOut, totalGrowth);
            }

            // update mask
            for (i = 0; i < x->x_numFilters; i++)
            {
                if (x->x_fftwIn[i] > x->x_mask[i])
                {
                    x->x_mask[i] = x->x_fftwIn[i];
                    x->x_numPeriods[i] = 0;
                }
                else
                    if (++x->x_numPeriods[i] >= x->x_maskPeriods)
                        x->x_mask[i] *= x->x_maskDecay;
            }

            x->x_prevTotalGrowth = totalGrowth;
        }

        post ("%s: analyzed %i frames", x->x_objSymbol->s_name, nFrames);
    }
}


static void bark_free (t_bark* x)
{
    t_filterIdx i;

    // free FFTW stuff
    t_freebytes (x->x_fftwIn, x->x_window * sizeof (t_sample));
    fftwf_free (x->x_fftwOut);
    fftwf_destroy_plan (x->x_fftwPlan);

    // free the output list
    t_freebytes (x->x_growthList, x->x_numFilters * sizeof (t_atom));

    // free the window memory
    t_freebytes (x->x_blackman, x->x_window * sizeof (t_float));
    t_freebytes (x->x_cosine, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hamming, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hann, x->x_window * sizeof (t_float));

    // free the mask memory
    t_freebytes (x->x_mask, x->x_numFilters * sizeof (t_float));

    // free the growth record memory
    t_freebytes (x->x_growth, x->x_numFilters * sizeof (t_float));

    // free the mask counter memory
    t_freebytes (x->x_numPeriods, x->x_numFilters * sizeof (t_filterIdx));

    // free the loudness weights memory
    t_freebytes (x->x_loudWeights, x->x_numFilters * sizeof (t_float));

    // free the filterFreqs memory
    t_freebytes (x->x_filterFreqs, x->x_sizeFilterFreqs * sizeof (t_float));

    // free the filterbank memory
    for (i = 0; i < x->x_numFilters; i++)
        t_freebytes (x->x_filterbank[i].filter, x->x_filterbank[i].filterSize * sizeof (t_float));

    t_freebytes (x->x_filterbank, x->x_numFilters * sizeof (t_filter));
}


void bark_setup (void)
{
    bark_class =
    class_new (
        gensym ("bark"),
        (t_newmethod)bark_new,
        (t_method)bark_free,
        sizeof (t_bark),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)bark_new,
        gensym ("timbreIDLib/bark"),
        A_GIMME,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_thresh,
        gensym ("thresh"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_minvel,
        gensym ("minvel"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_filter_range,
        gensym ("filter_range"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_mask,
        gensym ("mask"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_debounce,
        gensym ("debounce"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_spew,
        gensym ("spew"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_debug,
        gensym ("debug"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_use_weights,
        gensym ("loudness"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_normalize,
        gensym ("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_spec_band_avg,
        gensym ("spec_band_avg"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        bark_class,
        (t_method)bark_filter_avg,
        gensym ("filter_avg"),
        A_DEFFLOAT,
        0
    );
}
