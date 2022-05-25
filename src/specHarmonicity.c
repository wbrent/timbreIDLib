/*

specHarmonicity

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *specHarmonicity_class;

typedef struct _specHarmonicity
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_bool x_powerSpectrum;
    t_bool x_inputFund;
    t_float x_fundFreq;
    t_float x_minFund;
    t_float x_maxFund;
    t_float x_threshPct;
    t_uShortInt x_maxPeaks;
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
    t_outlet *x_harm;
    t_outlet *x_inHarm;
} t_specHarmonicity;


/* ------------------------ specHarmonicity -------------------------------- */
static void specHarmonicity_resizeWindow(t_specHarmonicity *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
{
    t_sampIdx windowHalf;

    windowHalf = window * 0.5;

    if(window<TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        windowHalf = window * 0.5;
        post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);

        *endSamp = startSamp + window-1;
        if(*endSamp >= x->x_arrayPoints)
            *endSamp = x->x_arrayPoints-1;
    }

    // hang on to these values for next time
    x->x_window = window;
    x->x_windowHalf = windowHalf;

    x->x_fftwIn = (t_sample *)t_resizebytes(x->x_fftwIn, oldWindow*sizeof(t_sample), x->x_window*sizeof(t_sample));

    fftwf_free(x->x_fftwOut);
    fftwf_destroy_plan(x->x_fftwPlan);
    // set up a new FFTW output buffer
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
    // FFTW plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG);

    x->x_blackman = (t_float *)t_resizebytes(x->x_blackman, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_cosine = (t_float *)t_resizebytes(x->x_cosine, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_hamming = (t_float *)t_resizebytes(x->x_hamming, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_hann = (t_float *)t_resizebytes(x->x_hann, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);
}


static void specHarmonicity_analyze(t_specHarmonicity *x, t_floatarg start, t_floatarg n)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, window, startSamp, endSamp;
        t_uShortInt numPeaks;
        t_float *windowFuncPtr, fund, harmSpacing, halfHarmSpacing, harm, inHarm, harmDividend, inHarmDividend, divisor, *flagsBuf, minPeakVal, maxPeakVal, thresh, *peakFreqs, *peakAmps;

        startSamp = (start<0)?0:start;

        if(n)
            endSamp = startSamp + n-1;
        else
            endSamp = startSamp + x->x_window-1;

        if(endSamp >= x->x_arrayPoints)
            endSamp = x->x_arrayPoints-1;

        window = endSamp-startSamp+1;

        if(endSamp <= startSamp)
        {
            post("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        if(x->x_window != window)
            specHarmonicity_resizeWindow(x, x->x_window, window, startSamp, &endSamp);

        // construct analysis window
        for(i=0, j=startSamp; j<=endSamp; i++, j++)
            x->x_fftwIn[i] = x->x_vec[j].w_float;

        windowFuncPtr = x->x_blackman;

        switch(x->x_windowFunction)
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
        if(x->x_windowFunction!=rectangular)
            for(i=0; i<x->x_window; i++, windowFuncPtr++)
                x->x_fftwIn[i] *= *windowFuncPtr;

        fftwf_execute(x->x_fftwPlan);

        tIDLib_power(x->x_windowHalf+1, x->x_fftwOut, x->x_fftwIn);

        if(!x->x_powerSpectrum)
            tIDLib_mag(x->x_windowHalf+1, x->x_fftwIn);

        // find all significant peaks
        flagsBuf = (t_float *)t_getbytes((x->x_windowHalf+1)*sizeof(t_float));
        minPeakVal = FLT_MAX;
        maxPeakVal = -FLT_MAX;
        numPeaks = 0;

        tIDLib_peaksValleys(x->x_windowHalf+1, x->x_fftwIn, flagsBuf, &minPeakVal, &maxPeakVal);

        thresh = maxPeakVal * (x->x_threshPct/100.0);
        peakFreqs = (t_float *)t_getbytes(0);
        peakAmps = (t_float *)t_getbytes(0);

        for(i=0; i<=x->x_windowHalf; i++)
        {
            // 0.5 in the flagsBuf means a half peak, which we'll ignore
            if(flagsBuf[i]>0.5)
            {
                t_float thisAmp;

                thisAmp = x->x_fftwIn[i];

                if(thisAmp>=thresh)
                {
                    peakFreqs = (t_float *)t_resizebytes(peakFreqs, numPeaks*sizeof(t_float), (numPeaks+1)*sizeof(t_float));
                    peakAmps = (t_float *)t_resizebytes(peakAmps, numPeaks*sizeof(t_float), (numPeaks+1)*sizeof(t_float));

                    peakAmps[numPeaks] = thisAmp;
                    peakFreqs[numPeaks] = tIDLib_bin2freq(i, x->x_window, x->x_sr);
                    numPeaks++;

                    if(numPeaks>=x->x_maxPeaks)
                        break;
                }
            }
        }

        t_freebytes(flagsBuf, (x->x_windowHalf+1)*sizeof(t_float));

        harm = inHarm = harmDividend = inHarmDividend = divisor = 0.0;

        if(x->x_inputFund)
            fund = x->x_fundFreq;
        else
        {
            if(peakFreqs[0]==0.0)
                fund = peakFreqs[1];
            else
                fund = peakFreqs[0];
        }

        if(fund<x->x_minFund || fund>x->x_maxFund)
        {
            harmDividend = -numPeaks; // to make harm value -1.0
            goto earlyExit;
        }

        harmSpacing = fund;
        halfHarmSpacing = harmSpacing*0.5;

        for(i=0; i<numPeaks; i++)
        {
            t_float thisAmp;

            thisAmp = peakAmps[i];

            if(thisAmp>0.0)
            {
                t_float thisFreq, deviation;
                t_uShortInt roundedHarm;

                thisFreq = peakFreqs[i];
                roundedHarm = roundf(thisFreq/fund);
                deviation = fabs(thisFreq - (roundedHarm*fund));
                inHarmDividend += deviation * thisAmp;
                divisor += thisAmp;

                harmDividend += (halfHarmSpacing-deviation)/halfHarmSpacing;
            }
        }

        earlyExit:

        t_freebytes(peakAmps, numPeaks*sizeof(t_float));
        t_freebytes(peakFreqs, numPeaks*sizeof(t_float));

        if(divisor<=0.0 || fund<=0.0)
            inHarm = -1.0;
        else
            inHarm = (2*inHarmDividend)/(divisor*fund);

        if(numPeaks<=0)
            harm = -1.0;
        else
            harm = harmDividend/numPeaks;

        outlet_float(x->x_inHarm, inHarm);
        outlet_float(x->x_harm, harm);
    }
}


static void specHarmonicity_chain_fftData(t_specHarmonicity *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;
    t_uShortInt numPeaks;
    t_float fund, harmSpacing, halfHarmSpacing, harm, inHarm, harmDividend, inHarmDividend, divisor, *flagsBuf, minPeakVal, maxPeakVal, thresh, *peakFreqs, *peakAmps;

    // incoming fftData list should be 2*(N/2+1) elements long, so windowHalf is:
    windowHalf = argc-2;
    windowHalf *= 0.5;

    // make sure that windowHalf == x->x_windowHalf in order to avoid an out of bounds memory read in the tIDLib_ functions below. we won't resize all memory based on an incoming chain_ command with a different window size. instead, just throw an error and exit
    if(windowHalf!=x->x_windowHalf)
    {
        pd_error(x, "%s: window size of chain_ message (%lu) does not match current window size (%lu)", x->x_objSymbol->s_name, windowHalf*2, x->x_window);
        return;
    }

    // fill the x_fftwOut buffer with the incoming fftData list, for both real and imag elements
    for(i=0; i<=x->x_windowHalf; i++)
    {
        x->x_fftwOut[i][0] = atom_getfloat(argv+i);
        x->x_fftwOut[i][1] = atom_getfloat(argv+(x->x_windowHalf+1)+i);
    }

    tIDLib_power(x->x_windowHalf+1, x->x_fftwOut, x->x_fftwIn);

    if(!x->x_powerSpectrum)
        tIDLib_mag(x->x_windowHalf+1, x->x_fftwIn);

    // find all significant peaks
    flagsBuf = (t_float *)t_getbytes((x->x_windowHalf+1)*sizeof(t_float));
    minPeakVal = FLT_MAX;
    maxPeakVal = -FLT_MAX;
    numPeaks = 0;

    tIDLib_peaksValleys(x->x_windowHalf+1, x->x_fftwIn, flagsBuf, &minPeakVal, &maxPeakVal);

    thresh = maxPeakVal * (x->x_threshPct/100.0);
    peakFreqs = (t_float *)t_getbytes(0);
    peakAmps = (t_float *)t_getbytes(0);

    for(i=0; i<=x->x_windowHalf; i++)
    {
        // 0.5 in the flagsBuf means a half peak, which we'll ignore
        if(flagsBuf[i]>0.5)
        {
            t_float thisAmp;

            thisAmp = x->x_fftwIn[i];

            if(thisAmp>=thresh)
            {
                peakFreqs = (t_float *)t_resizebytes(peakFreqs, numPeaks*sizeof(t_float), (numPeaks+1)*sizeof(t_float));
                peakAmps = (t_float *)t_resizebytes(peakAmps, numPeaks*sizeof(t_float), (numPeaks+1)*sizeof(t_float));

                peakAmps[numPeaks] = thisAmp;
                peakFreqs[numPeaks] = tIDLib_bin2freq(i, x->x_window, x->x_sr);
                numPeaks++;

                if(numPeaks>=x->x_maxPeaks)
                    break;
            }
        }
    }

    t_freebytes(flagsBuf, (x->x_windowHalf+1)*sizeof(t_float));

    harm = inHarm = harmDividend = inHarmDividend = divisor = 0.0;

    if(x->x_inputFund)
        fund = x->x_fundFreq;
    else
    {
        if(peakFreqs[0]==0.0)
            fund = peakFreqs[1];
        else
            fund = peakFreqs[0];
    }

    if(fund<x->x_minFund || fund>x->x_maxFund)
    {
        harmDividend = -numPeaks; // to make harm value -1.0
        goto earlyExit;
    }

    harmSpacing = fund;
    halfHarmSpacing = harmSpacing*0.5;

    for(i=0; i<numPeaks; i++)
    {
        t_float thisAmp;

        thisAmp = peakAmps[i];

        if(thisAmp>0.0)
        {
            t_float thisFreq, deviation;
            t_uShortInt roundedHarm;

            thisFreq = peakFreqs[i];
            roundedHarm = roundf(thisFreq/fund);
            deviation = fabs(thisFreq - (roundedHarm*fund));
            inHarmDividend += deviation * thisAmp;
            divisor += thisAmp;

            harmDividend += (halfHarmSpacing-deviation)/halfHarmSpacing;
        }
    }

    earlyExit:

    t_freebytes(peakAmps, numPeaks*sizeof(t_float));
    t_freebytes(peakFreqs, numPeaks*sizeof(t_float));

    if(divisor<=0.0 || fund<=0.0)
        inHarm = -1.0;
    else
        inHarm = (2*inHarmDividend)/(divisor*fund);

    if(numPeaks<=0)
        harm = -1.0;
    else
        harm = harmDividend/numPeaks;

    outlet_float(x->x_inHarm, inHarm);
    outlet_float(x->x_harm, harm);
}


static void specHarmonicity_chain_magSpec(t_specHarmonicity *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;
    t_uShortInt numPeaks;
    t_float fund, harmSpacing, halfHarmSpacing, harm, inHarm, harmDividend, inHarmDividend, divisor, *flagsBuf, minPeakVal, maxPeakVal, thresh, *peakFreqs, *peakAmps;

    // incoming magSpec list should be N/2+1 elements long, so windowHalf is one less than this
    windowHalf = argc-1;

    // make sure that windowHalf == x->x_windowHalf in order to avoid an out of bounds memory read in the tIDLib_ functions below. we won't resize all memory based on an incoming chain_ command with a different window size. instead, just throw an error and exit
    if(windowHalf!=x->x_windowHalf)
    {
        pd_error(x, "%s: window size of chain_ message (%lu) does not match current window size (%lu)", x->x_objSymbol->s_name, windowHalf*2, x->x_window);
        return;
    }

    // fill the x_fftwIn buffer with the incoming magSpec list
    for(i=0; i<=x->x_windowHalf; i++)
        x->x_fftwIn[i] = atom_getfloat(argv+i);

    // find all significant peaks
    flagsBuf = (t_float *)t_getbytes((x->x_windowHalf+1)*sizeof(t_float));
    minPeakVal = FLT_MAX;
    maxPeakVal = -FLT_MAX;
    numPeaks = 0;

    tIDLib_peaksValleys(x->x_windowHalf+1, x->x_fftwIn, flagsBuf, &minPeakVal, &maxPeakVal);

    thresh = maxPeakVal * (x->x_threshPct/100.0);
    peakFreqs = (t_float *)t_getbytes(0);
    peakAmps = (t_float *)t_getbytes(0);

    for(i=0; i<=x->x_windowHalf; i++)
    {
        // 0.5 in the flagsBuf means a half peak, which we'll ignore
        if(flagsBuf[i]>0.5)
        {
            t_float thisAmp;

            thisAmp = x->x_fftwIn[i];

            if(thisAmp>=thresh)
            {
                peakFreqs = (t_float *)t_resizebytes(peakFreqs, numPeaks*sizeof(t_float), (numPeaks+1)*sizeof(t_float));
                peakAmps = (t_float *)t_resizebytes(peakAmps, numPeaks*sizeof(t_float), (numPeaks+1)*sizeof(t_float));

                peakAmps[numPeaks] = thisAmp;
                peakFreqs[numPeaks] = tIDLib_bin2freq(i, x->x_window, x->x_sr);
                numPeaks++;

                if(numPeaks>=x->x_maxPeaks)
                    break;
            }
        }
    }

    t_freebytes(flagsBuf, (x->x_windowHalf+1)*sizeof(t_float));

    harm = inHarm = harmDividend = inHarmDividend = divisor = 0.0;

    if(x->x_inputFund)
        fund = x->x_fundFreq;
    else
    {
        if(peakFreqs[0]==0.0)
            fund = peakFreqs[1];
        else
            fund = peakFreqs[0];
    }

    if(fund<x->x_minFund || fund>x->x_maxFund)
    {
        harmDividend = -numPeaks; // to make harm value -1.0
        goto earlyExit;
    }

    harmSpacing = fund;
    halfHarmSpacing = harmSpacing*0.5;

    for(i=0; i<numPeaks; i++)
    {
        t_float thisAmp;

        thisAmp = peakAmps[i];

        if(thisAmp>0.0)
        {
            t_float thisFreq, deviation;
            t_uShortInt roundedHarm;

            thisFreq = peakFreqs[i];
            roundedHarm = roundf(thisFreq/fund);
            deviation = fabs(thisFreq - (roundedHarm*fund));
            inHarmDividend += deviation * thisAmp;
            divisor += thisAmp;

            harmDividend += (halfHarmSpacing-deviation)/halfHarmSpacing;
        }
    }

    earlyExit:

    t_freebytes(peakAmps, numPeaks*sizeof(t_float));
    t_freebytes(peakFreqs, numPeaks*sizeof(t_float));

    if(divisor<=0.0 || fund<=0.0)
        inHarm = -1.0;
    else
        inHarm = (2*inHarmDividend)/(divisor*fund);

    if(numPeaks<=0)
        harm = -1.0;
    else
        harm = harmDividend/numPeaks;

    outlet_float(x->x_inHarm, inHarm);
    outlet_float(x->x_harm, harm);
}


// analyze the whole damn array
static void specHarmonicity_bang(t_specHarmonicity *x)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx window, startSamp;
        startSamp = 0;
        window = x->x_arrayPoints;
        specHarmonicity_analyze(x, startSamp, window);
    }
}


static void specHarmonicity_set(t_specHarmonicity *x, t_symbol *s)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void specHarmonicity_print(t_specHarmonicity *x)
{
    post("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post("%s use input fundamental: %i", x->x_objSymbol->s_name, x->x_inputFund);
    post("%s current fundamental frequency: %f", x->x_objSymbol->s_name, x->x_fundFreq);
    post("%s minimum fundamental frequency: %f", x->x_objSymbol->s_name, x->x_minFund);
    post("%s maximum fundamental frequency: %f", x->x_objSymbol->s_name, x->x_maxFund);
    post("%s spectral peak threshold: %f", x->x_objSymbol->s_name, x->x_threshPct);
    post("%s maximum spectral peaks to consider: %i", x->x_objSymbol->s_name, x->x_maxPeaks);
}


static void specHarmonicity_samplerate(t_specHarmonicity *x, t_floatarg sr)
{
    if(sr<TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void specHarmonicity_window(t_specHarmonicity *x, t_floatarg w)
{
    t_sampIdx endSamp;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow() requires that
    endSamp = 0;

    specHarmonicity_resizeWindow(x, x->x_window, w, 0, &endSamp);
}


static void specHarmonicity_windowFunction(t_specHarmonicity *x, t_floatarg f)
{
    f = (f<0)?0:f;
    f = (f>4)?4:f;
    x->x_windowFunction = f;

    switch(x->x_windowFunction)
    {
        case rectangular:
            post("%s window function: rectangular.", x->x_objSymbol->s_name);
            break;
        case blackman:
            post("%s window function: blackman.", x->x_objSymbol->s_name);
            break;
        case cosine:
            post("%s window function: cosine.", x->x_objSymbol->s_name);
            break;
        case hamming:
            post("%s window function: hamming.", x->x_objSymbol->s_name);
            break;
        case hann:
            post("%s window function: hann.", x->x_objSymbol->s_name);
            break;
        default:
            break;
    };
}


static void specHarmonicity_powerSpectrum(t_specHarmonicity *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
    x->x_powerSpectrum = spec;

    if(x->x_powerSpectrum)
        post("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void specHarmonicity_maxPeaks(t_specHarmonicity *x, t_floatarg max)
{
    max = (max<1.0)?1.0:max;
    max = (max>(x->x_window/4.0))?(x->x_window/4.0):max;
    x->x_maxPeaks = max;

    post("%s maximum spectral peaks to consider: %i.", x->x_objSymbol->s_name, x->x_maxPeaks);
}


static void specHarmonicity_inputFund(t_specHarmonicity *x, t_floatarg useFund)
{
    useFund = (useFund<0)?0:useFund;
    useFund = (useFund>1)?1:useFund;
    x->x_inputFund = useFund;

    if(x->x_inputFund)
        post("%s using incoming fundamental.", x->x_objSymbol->s_name);
    else
        post("%s using first spectral peak as fundamental.", x->x_objSymbol->s_name);
}


static void specHarmonicity_peakThresh(t_specHarmonicity *x, t_floatarg thresh)
{
    thresh = (thresh<0.0)?0.0:thresh;
    thresh = (thresh>100.0)?100.0:thresh;
    x->x_threshPct = thresh;

    post("%s spectral peak thresh: %0.2f%% of maximum peak amplitude.", x->x_objSymbol->s_name, x->x_threshPct);
}


static void specHarmonicity_fundFreq(t_specHarmonicity *x, t_floatarg fund)
{
    if(fund <= 0.0)
        x->x_fundFreq = 0.0;
    else
        x->x_fundFreq = fund;
}


static void specHarmonicity_minFund(t_specHarmonicity *x, t_floatarg min)
{
    if (min < 0.0 || min > 20000.0)
        pd_error(x, "%s: minimum fundamental frequency must be between 0 and 20kHz.", x->x_objSymbol->s_name);
    else
        x->x_minFund = min;
}


static void specHarmonicity_maxFund(t_specHarmonicity *x, t_floatarg max)
{
    if (max < 0.0 || max > 20000.0)
        pd_error(x, "%s: maximum fundamental frequency must be between 0 and 20kHz.", x->x_objSymbol->s_name);
    else
        x->x_maxFund = max;
}


static void *specHarmonicity_new(t_symbol *s, int argc, t_atom *argv)
{
    t_specHarmonicity *x = (t_specHarmonicity *)pd_new(specHarmonicity_class);
    t_sampIdx i;
//	t_garray *a;

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("fund"));
    x->x_harm = outlet_new(&x->x_obj, &s_float);
    x->x_inHarm = outlet_new(&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch(argc)
    {
        case 1:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            break;

        case 0:
            post("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze()
            x->x_arrayName = gensym("NOARRAYSPECIFIED");
            break;

        default:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
            break;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_window = TID_WINDOWSIZEDEFAULT;
    x->x_windowHalf = x->x_window*0.5;
    x->x_windowFunction = blackman;
    x->x_powerSpectrum = false;
    x->x_inputFund = false;
    x->x_fundFreq = 100.0;
    x->x_minFund = 30.0;
    x->x_maxFund = 4000.0;
    x->x_threshPct = 5.0;
    x->x_maxPeaks = 24;

    x->x_fftwIn = (t_sample *)t_getbytes(x->x_window*sizeof(t_sample));

    // set up the FFTW output buffer. Is there no function to initialize it?
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

    // FFTW plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG); // FFTWPLANNERFLAG may be slower than FFTWPLANNERFLAG but more efficient after the first run?

    for(i=0; i<x->x_window; i++)
        x->x_fftwIn[i] = 0.0;

      x->x_blackman = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
      x->x_cosine = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
      x->x_hamming = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
      x->x_hann = (t_float *)t_getbytes(x->x_window*sizeof(t_float));

     // initialize signal windowing functions
    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);

    return (x);
}


static void specHarmonicity_free(t_specHarmonicity *x)
{
    // free FFTW stuff
    t_freebytes(x->x_fftwIn, (x->x_window)*sizeof(t_sample));
    fftwf_free(x->x_fftwOut);
    fftwf_destroy_plan(x->x_fftwPlan);

    // free the window memory
    t_freebytes(x->x_blackman, x->x_window*sizeof(t_float));
    t_freebytes(x->x_cosine, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hamming, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hann, x->x_window*sizeof(t_float));
}


void specHarmonicity_setup(void)
{
    specHarmonicity_class =
    class_new(
        gensym("specHarmonicity"),
        (t_newmethod)specHarmonicity_new,
        (t_method)specHarmonicity_free,
        sizeof(t_specHarmonicity),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)specHarmonicity_new,
        gensym("timbreIDLib/specHarmonicity"),
        A_GIMME,
        0
    );

    class_addbang(specHarmonicity_class, specHarmonicity_bang);

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_analyze,
        gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_chain_fftData,
        gensym("chain_fftData"),
        A_GIMME,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_chain_magSpec,
        gensym("chain_magSpec"),
        A_GIMME,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_set,
        gensym("set"),
        A_SYMBOL,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_print,
        gensym("print"),
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_samplerate,
        gensym("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_windowFunction,
        gensym("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_powerSpectrum,
        gensym("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_maxPeaks,
        gensym("max_peaks"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_peakThresh,
        gensym("peak_thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_inputFund,
        gensym("input_fund"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_fundFreq,
        gensym("fund"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_minFund,
        gensym("min_fund"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specHarmonicity_class,
        (t_method)specHarmonicity_maxFund,
        gensym("max_fund"),
        A_DEFFLOAT,
        0
    );
}
