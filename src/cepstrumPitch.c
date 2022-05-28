/*

cepstrumPitch

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *cepstrumPitch_class;

typedef struct _cepstrumPitch
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_bool x_powerSpectrum;
    t_bool x_powerCepstrum;
    t_bool x_spectrumOffset;
    t_float x_loFreq;
    t_float x_hiFreq;
    t_float x_thresh;
    t_sample *x_fftwIn;
    fftwf_complex *x_fftwOut;
    fftwf_plan x_fftwForwardPlan;
    fftwf_plan x_fftwBackwardPlan;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_word *x_vec;
    t_symbol *x_arrayName;
    t_sampIdx x_arrayPoints;
    t_outlet *x_pitch;
} t_cepstrumPitch;


/* ------------------------ cepstrumPitch -------------------------------- */
static void cepstrumPitch_resizeWindow (t_cepstrumPitch *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
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
}


static void cepstrumPitch_analyze (t_cepstrumPitch *x, t_floatarg start, t_floatarg n)
{
    t_sampIdx i, j, binCount, window, startSamp, endSamp;
    t_binIdx loFreqBin, hiFreqBin, maxValIdx;
    t_float nRecip, *windowFuncPtr, maxVal, pitch, loFreqBinFloat, hiFreqBinFloat, mean, std, sum;
    t_garray *a;

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

        window = endSamp - startSamp + 1;
        nRecip = 1.0/window;

        if (endSamp <= startSamp)
        {
            post ("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        if (x->x_window != window)
            cepstrumPitch_resizeWindow (x, x->x_window, window, startSamp, &endSamp);

        loFreqBinFloat = roundf (x->x_sr/x->x_loFreq);
        hiFreqBinFloat = roundf (x->x_sr/x->x_hiFreq);

        // these must be after the potential window resize, which would change x_windowHalf and therefore valid bounds for hiFreqBin and loFreqBin
        hiFreqBin = (hiFreqBinFloat<0)?0:hiFreqBinFloat;
        hiFreqBin = (hiFreqBin>x->x_windowHalf)?x->x_windowHalf:hiFreqBin;

        loFreqBin = (loFreqBinFloat<0)?0:loFreqBinFloat;
        loFreqBin = (loFreqBin>x->x_windowHalf)?x->x_windowHalf:loFreqBin;

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

        // add 1.0 to power or magnitude spectrum before taking the log and then IFT. Avoid large negative values from log(negativeNum). MPM (McCleod Pitch Method)
        if (x->x_spectrumOffset)
            for (i = 0; i < x->x_windowHalf + 1; i++)
                x->x_fftwIn[i] += 1.0;

        tIDLib_log(x->x_windowHalf + 1, x->x_fftwIn);

        // copy forward DFT magnitude result into real part of backward DFT complex input buffer, and zero out the imaginary part. fftwOut is only N/2+1 points long, while fftwIn is N points long
        for (i = 0; i < x->x_windowHalf + 1; i++)
        {
            x->x_fftwOut[i][0] = x->x_fftwIn[i];
            x->x_fftwOut[i][1] = 0.0;
        }

        fftwf_execute (x->x_fftwBackwardPlan);

        for (i = 0; i < x->x_windowHalf + 1; i++)
            x->x_fftwIn[i] *= nRecip;

        // optionally square the cepstrum results for power cepstrum
        if (x->x_powerCepstrum)
            for (i = 0; i < x->x_windowHalf + 1; i++)
                x->x_fftwIn[i] = x->x_fftwIn[i]*x->x_fftwIn[i];

        maxVal = 0.0;
        maxValIdx = 0;

        sum=mean=std=0.0;
        binCount=0;

        // traverse from hiFreq to loFreq because the high frequency cepstrum bin is lower than the low frequency cepstrum bin
        for (i=hiFreqBin; i < =loFreqBin; i++, binCount++)
        {
            // check that loFreqBin doesn't go above Nyquist bin
            if (i>=x->x_windowHalf)
                break;

            // accumulate a sum to get the mean below
            sum += x->x_fftwIn[i];

            if (x->x_fftwIn[i]>maxVal)
            {
                maxVal = x->x_fftwIn[i];
                maxValIdx = i;
            }
        }

        mean = sum/binCount;
        sum = 0.0;
        binCount = 0;

        // center & square the data
        for (i=hiFreqBin; i < =loFreqBin; i++, binCount++)
        {
            x->x_fftwIn[i] -= mean;
            x->x_fftwIn[i] *= x->x_fftwIn[i];
            sum += x->x_fftwIn[i];
        }

        // get standard deviation
        std = sum/(binCount-1);
        std = sqrt(std);

        // see if maxVal is above the mean by more than x_thresh standard deviations
        if ( fabs (maxVal-mean) > (x->x_thresh*std) )
            pitch = ftom(x->x_sr/((t_float)maxValIdx));
        else
            pitch = -1500.0;

        outlet_float (x->x_pitch, pitch);
    }
}


// analyze the whole damn array
static void cepstrumPitch_bang (t_cepstrumPitch *x)
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
        cepstrumPitch_analyze (x, startSamp, window);
    }
}


static void cepstrumPitch_set (t_cepstrumPitch *x, t_symbol *s)
{
    t_garray *a;

    if ( !(a = (t_garray *)pd_findbyclass (s, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void cepstrumPitch_print (t_cepstrumPitch *x)
{
    post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s power cepstrum: %i", x->x_objSymbol->s_name, x->x_powerCepstrum);
    post ("%s spectrum offset: %i", x->x_objSymbol->s_name, x->x_spectrumOffset);
    post ("%s pitch range: %0.2f to %0.2f", x->x_objSymbol->s_name, ftom(x->x_loFreq), ftom(x->x_hiFreq));
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
}


static void cepstrumPitch_samplerate (t_cepstrumPitch *x, t_floatarg sr)
{
    if (sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void cepstrumPitch_window (t_cepstrumPitch *x, t_floatarg w)
{
    t_sampIdx endSamp;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow () requires that
    endSamp = 0;

    cepstrumPitch_resizeWindow (x, x->x_window, w, 0, &endSamp);
}


static void cepstrumPitch_windowFunction (t_cepstrumPitch *x, t_floatarg f)
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


static void cepstrumPitch_powerSpectrum (t_cepstrumPitch *x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void cepstrumPitch_powerCepstrum(t_cepstrumPitch *x, t_floatarg power)
{
    power = (power<0)?0:power;
    power = (power>1)?1:power;
    x->x_powerCepstrum = power;

    if (x->x_powerCepstrum)
        post ("%s using power cepstrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude cepstrum", x->x_objSymbol->s_name);
}


static void cepstrumPitch_spectrumOffset (t_cepstrumPitch *x, t_floatarg offset)
{
    offset = (offset<0)?0:offset;
    offset = (offset>1)?1:offset;
    x->x_spectrumOffset = offset;

    if (x->x_spectrumOffset)
        post ("%s spectrum offset ON", x->x_objSymbol->s_name);
    else
        post ("%s spectrum offset OFF", x->x_objSymbol->s_name);
}


static void cepstrumPitch_pitchRange(t_cepstrumPitch *x, t_floatarg low, t_floatarg hi)
{
    low = (low<0)?0:low;
    low = (low>20000)?20000:low;

    hi = (hi < 0)?0:hi;
    hi = (hi>20000)?20000:hi;

    if (low>hi)
    {
        t_float tmp;
        tmp = hi;
        hi = low;
        low = tmp;
    }

    x->x_loFreq = mtof(low);
    x->x_hiFreq = mtof(hi);

    post ("%s pitch range: %0.2f to %0.2f", x->x_objSymbol->s_name, ftom(x->x_loFreq), ftom(x->x_hiFreq));
}


static void cepstrumPitch_threshold(t_cepstrumPitch *x, t_floatarg thresh)
{
    x->x_thresh = (thresh<0)?0:thresh;

    post ("%s peak threshold: %0.4f standard deviations above mean", x->x_objSymbol->s_name, x->x_thresh);
}


static void *cepstrumPitch_new (t_symbol *s, int argc, t_atom *argv)
{
    t_cepstrumPitch *x = (t_cepstrumPitch *)pd_new (cepstrumPitch_class);
    t_sampIdx i;
//	t_garray *a;

    x->x_pitch = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 3:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_loFreq = atom_getfloat (argv + 1);
            x->x_hiFreq = atom_getfloat (argv + 2);
            x->x_loFreq = mtof(x->x_loFreq);
            x->x_hiFreq = mtof(x->x_hiFreq);
            x->x_loFreq = (x->x_loFreq<0)?0:x->x_loFreq;
            x->x_hiFreq = (x->x_hiFreq<0)?0:x->x_hiFreq;
            x->x_loFreq = (x->x_loFreq>20000)?20000:x->x_loFreq;
            x->x_hiFreq = (x->x_hiFreq>20000)?20000:x->x_hiFreq;
            break;

        case 2:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_loFreq = atom_getfloat (argv + 1);
            x->x_hiFreq = x->x_loFreq+12;
            x->x_loFreq = mtof(x->x_loFreq);
            x->x_hiFreq = mtof(x->x_hiFreq);
            x->x_loFreq = (x->x_loFreq<0)?0:x->x_loFreq;
            x->x_hiFreq = (x->x_hiFreq<0)?0:x->x_hiFreq;
            x->x_loFreq = (x->x_loFreq>20000)?20000:x->x_loFreq;
            x->x_hiFreq = (x->x_hiFreq>20000)?20000:x->x_hiFreq;
            break;

        case 1:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_loFreq = 50;
            x->x_hiFreq = 500;
            break;

        case 0:
            post ("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze ()
            x->x_arrayName = gensym ("NOARRAYSPECIFIED");
            x->x_loFreq = 50;
            x->x_hiFreq = 500;
            break;

        default:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_loFreq = 50;
            x->x_hiFreq = 500;
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post ("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
            break;
    }

    if (x->x_loFreq>x->x_hiFreq)
    {
        t_float tmp;
        tmp = x->x_hiFreq;
        x->x_hiFreq = x->x_loFreq;
        x->x_loFreq = tmp;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_window = TID_WINDOWSIZEDEFAULT;
    x->x_windowHalf = x->x_window * 0.5;
    x->x_windowFunction = rectangular;
    x->x_powerSpectrum = true;
    x->x_powerCepstrum = false;
    x->x_spectrumOffset = true;
    x->x_thresh = 0.0;

    x->x_fftwIn = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));

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


static void cepstrumPitch_free (t_cepstrumPitch *x)
{
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


void cepstrumPitch_setup (void)
{
    cepstrumPitch_class =
    class_new (
        gensym ("cepstrumPitch"),
        (t_newmethod)cepstrumPitch_new,
        (t_method)cepstrumPitch_free,
        sizeof (t_cepstrumPitch),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)cepstrumPitch_new,
        gensym ("timbreIDLib/cepstrumPitch"),
        A_GIMME,
        0
    );

    class_addbang (cepstrumPitch_class, cepstrumPitch_bang);

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_set,
        gensym ("set"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_powerCepstrum,
        gensym ("power_cepstrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_spectrumOffset,
        gensym ("spectrum_offset"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_pitchRange,
        gensym ("pitch_range"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        cepstrumPitch_class,
        (t_method)cepstrumPitch_threshold,
        gensym ("threshold"),
        A_DEFFLOAT,
        0
    );
}
