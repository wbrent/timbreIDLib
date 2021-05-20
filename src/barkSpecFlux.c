/*

barkSpecFlux

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *barkSpecFlux_class;

typedef struct _barkSpecFlux
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_bool x_normalize;
    t_bool x_powerSpectrum;
    t_bool x_logSpectrum;
    t_float *x_fftwInForwardWindow;
    t_float *x_fftwInBackWindow;
    fftwf_complex *x_fftwOutForwardWindow;
    fftwf_complex *x_fftwOutBackWindow;
    fftwf_plan x_fftwPlanForwardWindow;
    fftwf_plan x_fftwPlanBackWindow;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_word *x_vec;
    t_symbol *x_arrayName;
    t_sampIdx x_arrayPoints;
    t_filterIdx x_sizeFilterFreqs;
    t_filterIdx x_numFilters;
    t_float x_barkSpacing;
    t_float *x_filterFreqs;
    t_filter *x_filterbank;
    t_bool x_specBandAvg;
    t_bool x_filterAvg;
    t_fluxMode x_mode;
    t_bool x_squaredDiff;
    t_uInt x_separation;
    t_atom *x_listOut;
    t_outlet *x_fluxList;
    t_outlet *x_flux;
} t_barkSpecFlux;


/* ------------------------ barkSpecFlux -------------------------------- */
static void barkSpecFlux_resizeWindow(t_barkSpecFlux *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
{
    t_sampIdx i, windowHalf;

    windowHalf = window * 0.5;

    if(window<MINWINDOWSIZE)
    {
        window = WINDOWSIZEDEFAULT;
        windowHalf = window * 0.5;
        post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);

        *endSamp = startSamp + window-1;
        if(*endSamp >= x->x_arrayPoints)
            *endSamp = x->x_arrayPoints-1;
    }

    // hang on to these values for next time
    x->x_window = window;
    x->x_windowHalf = windowHalf;

    if(x->x_separation > x->x_window)
    {
        post("%s WARNING: window size change resulted in a separation greater than the current window size. Using default separation of a half window size.", x->x_objSymbol->s_name);
        x->x_separation = x->x_window*0.5;
    }

    x->x_fftwInForwardWindow = (t_float *)t_resizebytes(x->x_fftwInForwardWindow, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_fftwInBackWindow = (t_float *)t_resizebytes(x->x_fftwInBackWindow, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

    fftwf_free(x->x_fftwOutForwardWindow);
    fftwf_free(x->x_fftwOutBackWindow);
    fftwf_destroy_plan(x->x_fftwPlanForwardWindow);
    fftwf_destroy_plan(x->x_fftwPlanBackWindow);

    // set up a new FFTW output buffer
    x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
    x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

    // FFTW plan
    x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
    x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

    x->x_blackman = (t_float *)t_resizebytes(x->x_blackman, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_cosine = (t_float *)t_resizebytes(x->x_cosine, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_hamming = (t_float *)t_resizebytes(x->x_hamming, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_hann = (t_float *)t_resizebytes(x->x_hann, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);

    for(i=0; i<x->x_window; i++)
    {
        x->x_fftwInForwardWindow[i] = 0.0;
        x->x_fftwInBackWindow[i] = 0.0;
    }

    // x_numFilters doesn't change with a change to x_window, so oldNumFilters and newNumFilters arguments are the same
    tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, x->x_numFilters, x->x_numFilters, x->x_window, x->x_sr);
}


static void barkSpecFlux_analyze(t_barkSpecFlux *x, t_floatarg start, t_floatarg n)
{
    t_sampIdx i, j, window, startSamp, endSamp, startSampBack, endSampBack;
    t_float flux, *windowFuncPtr;
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        startSamp = (start<0)?0:start;

        if(n)
            endSamp = startSamp + n-1;
        else
            endSamp = startSamp + x->x_window-1;

        if(endSamp >= x->x_arrayPoints-1)
            endSamp = x->x_arrayPoints-1;

        window = endSamp-startSamp+1;

        if(endSamp <= startSamp)
        {
            post("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        if(x->x_window != window)
            barkSpecFlux_resizeWindow(x, x->x_window, window, startSamp, &endSamp);

        // construct forward analysis window
        for(i=0, j=startSamp; j<=endSamp; i++, j++)
            x->x_fftwInForwardWindow[i] = x->x_vec[j].w_float;

        // do these sample start/end location calculations AFTER the potential call to resizeWindow(), as x->x_window may have changed
        if(startSamp>=x->x_separation)
            startSampBack = startSamp - x->x_separation;
        else
            startSampBack = 0;

        endSampBack = startSampBack + x->x_window-1;

        // construct back analysis window x->x_separation frames earlier
        for(i=0, j=startSampBack; j<=endSampBack; i++, j++)
            x->x_fftwInBackWindow[i] = x->x_vec[j].w_float;

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

        // if x_windowFunction == 0, skip the windowing (rectangular)
        if(x->x_windowFunction!=rectangular)
            for(i=0; i<x->x_window; i++, windowFuncPtr++)
                x->x_fftwInForwardWindow[i] *= *windowFuncPtr;

        fftwf_execute(x->x_fftwPlanForwardWindow);

        // put the result of power calc back in x_fftwIn
        tIDLib_power(x->x_windowHalf+1, x->x_fftwOutForwardWindow, x->x_fftwInForwardWindow);

        if(!x->x_powerSpectrum)
            tIDLib_mag(x->x_windowHalf+1, x->x_fftwInForwardWindow);

        if(x->x_specBandAvg)
            tIDLib_specFilterBands(x->x_windowHalf+1, x->x_numFilters, x->x_fftwInForwardWindow, x->x_filterbank, x->x_normalize);
        else
            tIDLib_filterbankMultiply(x->x_fftwInForwardWindow, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

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

        // if x_windowFunction == 0, skip the windowing (rectangular)
        if(x->x_windowFunction!=rectangular)
            for(i=0; i<x->x_window; i++, windowFuncPtr++)
                x->x_fftwInBackWindow[i] *= *windowFuncPtr;

        fftwf_execute(x->x_fftwPlanBackWindow);

        // put the result of power calc back in x_fftwIn
        tIDLib_power(x->x_windowHalf+1, x->x_fftwOutBackWindow, x->x_fftwInBackWindow);

        if(!x->x_powerSpectrum)
            tIDLib_mag(x->x_windowHalf+1, x->x_fftwInBackWindow);

        if(x->x_specBandAvg)
            tIDLib_specFilterBands(x->x_windowHalf+1, x->x_numFilters, x->x_fftwInBackWindow, x->x_filterbank, x->x_normalize);
        else
            tIDLib_filterbankMultiply(x->x_fftwInBackWindow, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

        flux=0;

        for(i=0; i<x->x_numFilters; i++)
        {
            t_float diff, val;

        if (x->x_logSpectrum)
        {
            t_float logForward, logBack;

            if (x->x_fftwInForwardWindow[i] == 0.0)
                logForward = 0.0;
            else
                logForward = log (x->x_fftwInForwardWindow[i]);

            if (x->x_fftwInBackWindow[i] == 0.0)
                logBack = 0.0;
            else
                logBack = log (x->x_fftwInBackWindow[i]);

            diff = logForward - logBack;
        }
        else
            diff = x->x_fftwInForwardWindow[i] - x->x_fftwInBackWindow[i];

            switch(x->x_mode)
            {
                case mGrowth:
                    diff = (diff<0)?0:diff;
                    break;
                case mDecay:
                    diff = (diff>0)?0:diff;
                    break;
                default:
                    break;
            }

            if(x->x_squaredDiff)
                val = diff*diff;
            else
                val = fabs(diff);

            SETFLOAT(x->x_listOut+i, diff);
            flux += val;
        }

        outlet_list(x->x_fluxList, 0, x->x_numFilters, x->x_listOut);
        outlet_float(x->x_flux, flux);
    }
}


static void barkSpecFlux_chain_fftData(t_barkSpecFlux *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;
    t_float flux;

    // for barkSpecFlux fftData in particular:
    // incoming fftData list should be 4*(N/2+1) elements long, so windowHalf is:
    windowHalf = argc-4;
    windowHalf *= 0.25;

    // make sure that windowHalf == x->x_windowHalf in order to avoid an out of bounds memory read in the tIDLib_ functions below. we won't resize all memory based on an incoming chain_ command with a different window size. instead, just throw an error and exit
    if(windowHalf!=x->x_windowHalf)
    {
        pd_error(x, "%s: window size of chain_ message (%lu) does not match current window size (%lu)", x->x_objSymbol->s_name, windowHalf*2, x->x_window);
        return;
    }

    // fill the x_fftwOut buffers with the incoming fftData list, for both real and imag elements
    // for specFlux in particular, the first 2*(N/2+1) elements in the atom list are for the FORWARD window complex results. The second set of 2*(N/2+1) elements are for the BACK window complex results.
    for(i=0; i<=x->x_windowHalf; i++)
    {
        x->x_fftwOutForwardWindow[i][0] = atom_getfloat(argv+i);
        x->x_fftwOutForwardWindow[i][1] = atom_getfloat(argv+(x->x_windowHalf+1)+i);
        x->x_fftwOutBackWindow[i][0] = atom_getfloat(argv+(x->x_window+2)+i);
        x->x_fftwOutBackWindow[i][1] = atom_getfloat(argv+(x->x_window+x->x_windowHalf+3)+i);
    }

    // put the result of power calc back in x_fftwIn
    tIDLib_power(x->x_windowHalf+1, x->x_fftwOutForwardWindow, x->x_fftwInForwardWindow);
    tIDLib_power(x->x_windowHalf+1, x->x_fftwOutBackWindow, x->x_fftwInBackWindow);

    if(!x->x_powerSpectrum)
    {
        tIDLib_mag(x->x_windowHalf+1, x->x_fftwInForwardWindow);
        tIDLib_mag(x->x_windowHalf+1, x->x_fftwInBackWindow);
    }

    if(x->x_specBandAvg)
    {
        tIDLib_specFilterBands(windowHalf+1, x->x_numFilters, x->x_fftwInForwardWindow, x->x_filterbank, x->x_normalize);
        tIDLib_specFilterBands(windowHalf+1, x->x_numFilters, x->x_fftwInBackWindow, x->x_filterbank, x->x_normalize);
    }
    else
    {
        tIDLib_filterbankMultiply(x->x_fftwInForwardWindow, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);
        tIDLib_filterbankMultiply(x->x_fftwInBackWindow, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);
    }

    flux=0;

    for(i=0; i<x->x_numFilters; i++)
    {
        t_float diff, val;

        if (x->x_logSpectrum)
        {
            t_float logForward, logBack;

            if (x->x_fftwInForwardWindow[i] == 0.0)
                logForward = 0.0;
            else
                logForward = log (x->x_fftwInForwardWindow[i]);

            if (x->x_fftwInBackWindow[i] == 0.0)
                logBack = 0.0;
            else
                logBack = log (x->x_fftwInBackWindow[i]);

            diff = logForward - logBack;
        }
        else
            diff = x->x_fftwInForwardWindow[i] - x->x_fftwInBackWindow[i];

        switch(x->x_mode)
        {
            case mGrowth:
                diff = (diff<0)?0:diff;
                break;
            case mDecay:
                diff = (diff>0)?0:diff;
                break;
            default:
                break;
        }

        if(x->x_squaredDiff)
            val = diff*diff;
        else
            val = fabs(diff);

        SETFLOAT(x->x_listOut+i, diff);
        flux += val;
    }

    outlet_list(x->x_fluxList, 0, x->x_numFilters, x->x_listOut);
    outlet_float(x->x_flux, flux);
}


static void barkSpecFlux_chain_magSpec(t_barkSpecFlux *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;
    t_float flux;

    // for barkSpecFlux magSpec in particular:
    // incoming magSpec list should be 2*(N/2+1) elements long, so windowHalf is:
    windowHalf = argc-2;
    windowHalf *= 0.5;

    // make sure that windowHalf == x->x_windowHalf in order to avoid an out of bounds memory read in the tIDLib_ functions below. we won't resize all memory based on an incoming chain_ command with a different window size. instead, just throw an error and exit
    if(windowHalf!=x->x_windowHalf)
    {
        pd_error(x, "%s: window size of chain_ message (%lu) does not match current window size (%lu)", x->x_objSymbol->s_name, windowHalf*2, x->x_window);
        return;
    }

    // fill the x_fftwIn buffers with the incoming magSpec lists
    // for barkSpecFlux in particular, the first N/2+1 elements in the atom list are for the FORWARD window magnitudes. The second set of N/2+1 elements are for the BACK window magnitudes.
    for(i=0; i<=x->x_windowHalf; i++)
    {
        x->x_fftwInForwardWindow[i] = atom_getfloat(argv+i);
        x->x_fftwInBackWindow[i] = atom_getfloat(argv+(x->x_windowHalf+1)+i);
    }

    if(x->x_specBandAvg)
    {
        tIDLib_specFilterBands(windowHalf+1, x->x_numFilters, x->x_fftwInForwardWindow, x->x_filterbank, x->x_normalize);
        tIDLib_specFilterBands(windowHalf+1, x->x_numFilters, x->x_fftwInBackWindow, x->x_filterbank, x->x_normalize);
    }
    else
    {
        tIDLib_filterbankMultiply(x->x_fftwInForwardWindow, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);
        tIDLib_filterbankMultiply(x->x_fftwInBackWindow, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);
    }

    flux=0.0;

    for(i=0; i<x->x_numFilters; i++)
    {
        t_float diff, val;

        if (x->x_logSpectrum)
        {
            t_float logForward, logBack;

            if (x->x_fftwInForwardWindow[i] == 0.0)
                logForward = 0.0;
            else
                logForward = log (x->x_fftwInForwardWindow[i]);

            if (x->x_fftwInBackWindow[i] == 0.0)
                logBack = 0.0;
            else
                logBack = log (x->x_fftwInBackWindow[i]);

            diff = logForward - logBack;
        }
        else
            diff = x->x_fftwInForwardWindow[i] - x->x_fftwInBackWindow[i];

        switch(x->x_mode)
        {
            case mGrowth:
                diff = (diff<0)?0:diff;
                break;
            case mDecay:
                diff = (diff>0)?0:diff;
                break;
            default:
                break;
        }

        if(x->x_squaredDiff)
            val = diff*diff;
        else
            val = fabs(diff);

        SETFLOAT(x->x_listOut+i, diff);
        flux += val;
    }

    outlet_list(x->x_fluxList, 0, x->x_numFilters, x->x_listOut);
    outlet_float(x->x_flux, flux);
}


static void barkSpecFlux_chain_barkSpec(t_barkSpecFlux *x, t_symbol *s, int argc, t_atom *argv)
{
    t_filterIdx i;
    t_float flux;

    // make sure that argc == 2*(x->x_numFilters) in order to avoid an out of bounds memory read below. we won't resize all memory based on an incoming chain_ command with a different size. instead, just throw an error and exit
    if(argc!=(x->x_numFilters)*2)
    {
        pd_error(x, "%s: length of chain_ message (%i) does not match current number of Bark filters (%i)", x->x_objSymbol->s_name, argc, x->x_numFilters);
        return;
    }

    // fill the x_fftwIn buffer with the incoming magSpec list
    for(i=0; i<x->x_numFilters; i++)
    {
        x->x_fftwInForwardWindow[i] = atom_getfloat(argv+i);
        x->x_fftwInBackWindow[i] = atom_getfloat(argv+x->x_numFilters+i);
    }

    flux=0.0;

    for(i=0; i<x->x_numFilters; i++)
    {
        t_float diff, val;

        if (x->x_logSpectrum)
        {
            t_float logForward, logBack;

            if (x->x_fftwInForwardWindow[i] == 0.0)
                logForward = 0.0;
            else
                logForward = log (x->x_fftwInForwardWindow[i]);

            if (x->x_fftwInBackWindow[i] == 0.0)
                logBack = 0.0;
            else
                logBack = log (x->x_fftwInBackWindow[i]);

            diff = logForward - logBack;
        }
        else
            diff = x->x_fftwInForwardWindow[i] - x->x_fftwInBackWindow[i];

        switch(x->x_mode)
        {
            case mGrowth:
                diff = (diff<0)?0:diff;
                break;
            case mDecay:
                diff = (diff>0)?0:diff;
                break;
            default:
                break;
        }

        if(x->x_squaredDiff)
            val = diff*diff;
        else
            val = fabs(diff);

        SETFLOAT(x->x_listOut+i, diff);
        flux += val;
    }

    outlet_list(x->x_fluxList, 0, x->x_numFilters, x->x_listOut);
    outlet_float(x->x_flux, flux);
}


// analyze the whole damn array
static void barkSpecFlux_bang(t_barkSpecFlux *x)
{
    t_sampIdx window, startSamp;
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        startSamp = 0;
        window = x->x_arrayPoints;
        barkSpecFlux_analyze(x, startSamp, window);
    }
}


static void barkSpecFlux_createFilterbank(t_barkSpecFlux *x, t_floatarg bs)
{
    t_filterIdx oldNumFilters;

    x->x_barkSpacing = bs;

    if(x->x_barkSpacing<MINBARKSPACING || x->x_barkSpacing>MAXBARKSPACING)
    {
        x->x_barkSpacing = BARKSPACINGDEFAULT;
        post("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, MINBARKSPACING, MAXBARKSPACING, BARKSPACINGDEFAULT);
    }

    oldNumFilters = x->x_numFilters;

    x->x_sizeFilterFreqs = tIDLib_getBarkBoundFreqs(&x->x_filterFreqs, x->x_sizeFilterFreqs, x->x_barkSpacing, x->x_sr);

    x->x_numFilters = x->x_sizeFilterFreqs-2;

    tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, oldNumFilters, x->x_numFilters, x->x_window, x->x_sr);

    // resize listOut memory
    x->x_listOut = (t_atom *)t_resizebytes(x->x_listOut, oldNumFilters*sizeof(t_atom), x->x_numFilters*sizeof(t_atom));
}


static void barkSpecFlux_spec_band_avg(t_barkSpecFlux *x, t_floatarg avg)
{
    avg = (avg<0)?0:avg;
    avg = (avg>1)?1:avg;
    x->x_specBandAvg = avg;

    if(x->x_specBandAvg)
        post("%s: averaging energy in spectrum bands.", x->x_objSymbol->s_name);
    else
        post("%s: using triangular filterbank.", x->x_objSymbol->s_name);
}


static void barkSpecFlux_filter_avg(t_barkSpecFlux *x, t_floatarg avg)
{
    avg = (avg<0)?0:avg;
    avg = (avg>1)?1:avg;
    x->x_filterAvg = avg;

    if(x->x_filterAvg)
        post("%s: averaging energy in triangular filters.", x->x_objSymbol->s_name);
    else
        post("%s: summing energy in triangular filters.", x->x_objSymbol->s_name);
}


static void barkSpecFlux_set(t_barkSpecFlux *x, t_symbol *s)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void barkSpecFlux_print(t_barkSpecFlux *x)
{
    post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)x->x_sr);
    post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
    post("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post("%s log spectrum: %i", x->x_objSymbol->s_name, x->x_logSpectrum);
    post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post("%s separation: %i", x->x_objSymbol->s_name, x->x_separation);
    post("%s squared difference: %i", x->x_objSymbol->s_name, x->x_squaredDiff);

    switch(x->x_mode)
    {
        case mFlux:
            post("%s mode: flux", x->x_objSymbol->s_name);
            break;
        case mGrowth:
            post("%s mode: growth", x->x_objSymbol->s_name);
            break;
        case mDecay:
            post("%s mode: decay", x->x_objSymbol->s_name);
            break;
        default:
            post("%s mode: flux", x->x_objSymbol->s_name);
            break;
    }
}


static void barkSpecFlux_samplerate(t_barkSpecFlux *x, t_floatarg sr)
{
    if(sr<MINSAMPLERATE)
        x->x_sr = MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void barkSpecFlux_window(t_barkSpecFlux *x, t_floatarg w)
{
    t_sampIdx endSamp;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow() requires that
    endSamp = 0;

    barkSpecFlux_resizeWindow(x, x->x_window, w, 0, &endSamp);
}


static void barkSpecFlux_windowFunction(t_barkSpecFlux *x, t_floatarg f)
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


static void barkSpecFlux_powerSpectrum(t_barkSpecFlux *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
    x->x_powerSpectrum = spec;

    if(x->x_powerSpectrum)
        post("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void barkSpecFlux_logSpectrum(t_barkSpecFlux *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
    x->x_logSpectrum = spec;

    if(x->x_logSpectrum)
        post("%s log spectrum enabled", x->x_objSymbol->s_name);
    else
        post("%s log spectrum disabled", x->x_objSymbol->s_name);
}


static void barkSpecFlux_normalize(t_barkSpecFlux *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
    x->x_normalize = norm;

    if(x->x_normalize)
        post("%s spectrum normalization ON.", x->x_objSymbol->s_name);
    else
        post("%s spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void barkSpecFlux_separation(t_barkSpecFlux *x, t_floatarg s)
{
    if(s > x->x_window)
    {
        post("%s frame separation cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
        x->x_separation = x->x_windowHalf;
    }
    else if(s < 0)
    {
        post("%s frame separation must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
        x->x_separation = x->x_windowHalf;
    }
    else
    {
        x->x_separation = s;
    }

    post("%s frame separation: %i", x->x_objSymbol->s_name, x->x_separation);

}


static void barkSpecFlux_squaredDiff(t_barkSpecFlux *x, t_floatarg sd)
{
    sd = (sd<0)?0:sd;
    sd = (sd>1)?1:sd;
    x->x_squaredDiff = sd;

    if(x->x_squaredDiff)
        post("%s using squared difference", x->x_objSymbol->s_name);
    else
        post("%s using absolute value of difference", x->x_objSymbol->s_name);
}


static void barkSpecFlux_mode(t_barkSpecFlux *x, t_symbol *m)
{
    if(!strcmp(m->s_name, "flux"))
        x->x_mode = mFlux;
    else if(!strcmp(m->s_name, "growth"))
        x->x_mode = mGrowth;
    else if(!strcmp(m->s_name, "decay"))
        x->x_mode = mDecay;
    else
        x->x_mode = mFlux;
}


static void *barkSpecFlux_new(t_symbol *s, int argc, t_atom *argv)
{
    t_barkSpecFlux *x = (t_barkSpecFlux *)pd_new(barkSpecFlux_class);
    t_sampIdx i;
    t_float sepFloat;
//	t_garray *a;

    x->x_flux = outlet_new(&x->x_obj, &s_float);
    x->x_fluxList = outlet_new(&x->x_obj, gensym("list"));

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch(argc)
    {
        case 3:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_barkSpacing = atom_getfloat(argv+1);
            if(x->x_barkSpacing<MINBARKSPACING || x->x_barkSpacing>MAXBARKSPACING)
            {
                x->x_barkSpacing = BARKSPACINGDEFAULT;
                post("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, MINBARKSPACING, MAXBARKSPACING, BARKSPACINGDEFAULT);
            }

            sepFloat = atom_getfloat(argv+2);
            if(sepFloat > WINDOWSIZEDEFAULT)
            {
                post("%s WARNING: frame separation cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
                x->x_separation = WINDOWSIZEDEFAULT*0.5;
            }
            else if(sepFloat < 0)
            {
                post("%s WARNING: frame separation must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
                x->x_separation = WINDOWSIZEDEFAULT*0.5;
            }
            else
                x->x_separation = sepFloat;
            break;

        case 2:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_barkSpacing = atom_getfloat(argv+1);
            if(x->x_barkSpacing<MINBARKSPACING || x->x_barkSpacing>MAXBARKSPACING)
            {
                x->x_barkSpacing = BARKSPACINGDEFAULT;
                post("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, MINBARKSPACING, MAXBARKSPACING, BARKSPACINGDEFAULT);
            }

            x->x_separation = WINDOWSIZEDEFAULT*0.5;
            break;

        case 1:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_barkSpacing = BARKSPACINGDEFAULT;
            x->x_separation = WINDOWSIZEDEFAULT*0.5;
            break;

        case 0:
            post("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze()
            x->x_arrayName = gensym("NOARRAYSPECIFIED");
            x->x_separation = WINDOWSIZEDEFAULT*0.5;
            break;

        default:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post("%s WARNING: Too many arguments supplied. Using default Bark spacing of %0.2f Barks, and frame separation of %i samples.", x->x_objSymbol->s_name, BARKSPACINGDEFAULT, (t_uInt)(WINDOWSIZEDEFAULT*0.5));
            x->x_barkSpacing = BARKSPACINGDEFAULT;
            x->x_separation = WINDOWSIZEDEFAULT*0.5;
            break;
    }

    x->x_sr = SAMPLERATEDEFAULT;
    x->x_window = WINDOWSIZEDEFAULT;
    x->x_windowHalf = x->x_window*0.5;
    x->x_windowFunction = blackman;
    x->x_normalize = false;
    x->x_powerSpectrum = false;
    x->x_logSpectrum = false;
    x->x_specBandAvg = false;
    x->x_filterAvg = false;
    x->x_squaredDiff = false; // absolute value by default
    x->x_mode = mFlux;

    x->x_fftwInForwardWindow = (t_float *)t_getbytes(x->x_window * sizeof(t_float));
    x->x_fftwInBackWindow = (t_float *)t_getbytes(x->x_window * sizeof(t_float));

    x->x_listOut = (t_atom *)t_getbytes(x->x_numFilters*sizeof(t_atom));

    // set up the FFTW output buffer. Is there no function to initialize it?
    x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
    x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

    // DFT plan
    x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
    x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
    for(i=0; i<x->x_window; i++)
    {
        x->x_fftwInForwardWindow[i] = 0.0;
        x->x_fftwInBackWindow[i] = 0.0;
    }

      x->x_blackman = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
      x->x_cosine = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
      x->x_hamming = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
      x->x_hann = (t_float *)t_getbytes(x->x_window*sizeof(t_float));

     // initialize signal windowing functions
    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);

    // grab memory
    x->x_filterbank = (t_filter *)t_getbytes(0);
    x->x_filterFreqs = (t_float *)t_getbytes(0);

    x->x_sizeFilterFreqs = tIDLib_getBarkBoundFreqs(&x->x_filterFreqs, x->x_sizeFilterFreqs, x->x_barkSpacing, x->x_sr);

    // sizeFilterFreqs-2 is the correct number of filters, since we don't count the start point of the first filter, or the finish point of the last filter
    x->x_numFilters = x->x_sizeFilterFreqs-2;

    tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, 0, x->x_numFilters, x->x_window, x->x_sr);

    // create listOut memory
    x->x_listOut = (t_atom *)t_getbytes(x->x_numFilters*sizeof(t_atom));

    return (x);
}


static void barkSpecFlux_free(t_barkSpecFlux *x)
{
    t_filterIdx i;

    // free the list out memory
    t_freebytes(x->x_listOut, x->x_numFilters*sizeof(t_atom));

    // free FFTW stuff
    t_freebytes(x->x_fftwInForwardWindow, (x->x_window)*sizeof(t_float));
    t_freebytes(x->x_fftwInBackWindow, (x->x_window)*sizeof(t_float));
    fftwf_free(x->x_fftwOutForwardWindow);
    fftwf_free(x->x_fftwOutBackWindow);
    fftwf_destroy_plan(x->x_fftwPlanForwardWindow);
    fftwf_destroy_plan(x->x_fftwPlanBackWindow);

    // free the window memory
    t_freebytes(x->x_blackman, x->x_window*sizeof(t_float));
    t_freebytes(x->x_cosine, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hamming, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hann, x->x_window*sizeof(t_float));

    // free filterFreqs memory
    t_freebytes(x->x_filterFreqs, x->x_sizeFilterFreqs*sizeof(t_float));

    // free the filterbank memory
    for(i=0; i<x->x_numFilters; i++)
        t_freebytes(x->x_filterbank[i].filter, x->x_filterbank[i].size*sizeof(t_float));

    t_freebytes(x->x_filterbank, x->x_numFilters*sizeof(t_filter));
}


void barkSpecFlux_setup(void)
{
    barkSpecFlux_class =
    class_new(
        gensym("barkSpecFlux"),
        (t_newmethod)barkSpecFlux_new,
        (t_method)barkSpecFlux_free,
        sizeof(t_barkSpecFlux),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)barkSpecFlux_new,
        gensym("timbreIDLib/barkSpecFlux"),
        A_GIMME,
        0
    );

    class_addbang(barkSpecFlux_class, barkSpecFlux_bang);

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_analyze,
        gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_chain_fftData,
        gensym("chain_fftData"),
        A_GIMME,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_chain_magSpec,
        gensym("chain_magSpec"),
        A_GIMME,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_chain_barkSpec,
        gensym("chain_barkSpec"),
        A_GIMME,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_set,
        gensym("set"),
        A_SYMBOL,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_separation,
        gensym("separation"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_mode,
        gensym("mode"),
        A_DEFSYMBOL,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_squaredDiff,
        gensym("squared_diff"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_print,
        gensym("print"),
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_samplerate,
        gensym("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_windowFunction,
        gensym("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_powerSpectrum,
        gensym("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_logSpectrum,
        gensym("log_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_normalize,
        gensym("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_createFilterbank,
        gensym("filterbank"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_spec_band_avg,
        gensym("spec_band_avg"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpecFlux_class,
        (t_method)barkSpecFlux_filter_avg,
        gensym("filter_avg"),
        A_DEFFLOAT,
        0
    );
}
