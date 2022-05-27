/*

specFlux

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *specFlux_class;

typedef struct _specFlux
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
    t_sample *x_fftwInForwardWindow;
    t_sample *x_fftwInBackWindow;
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
    t_fluxMode x_mode;
    t_bool x_squaredDiff;
    t_uInt x_separation;
    t_atom *x_listOut;
    t_outlet *x_fluxList;
    t_outlet *x_flux;
} t_specFlux;


/* ------------------------ specFlux -------------------------------- */
static void specFlux_resizeWindow(t_specFlux *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
{
    t_sampIdx i, oldWindowHalf, windowHalf;

    windowHalf = window * 0.5;
    oldWindowHalf = oldWindow * 0.5;

    if(window < TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        windowHalf = window * 0.5;
        post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);

        *endSamp = startSamp + window - 1;
        if(*endSamp >= x->x_arrayPoints)
            *endSamp = x->x_arrayPoints - 1;
    }

    // hang on to these values for next time
    x->x_window = window;
    x->x_windowHalf = windowHalf;

    if(x->x_separation > x->x_window)
    {
        post("%s WARNING: window size change resulted in a separation greater than the current window size. Using default separation of a half window size.", x->x_objSymbol->s_name);
        x->x_separation = x->x_window * 0.5;
    }

    x->x_fftwInForwardWindow = (t_sample *)t_resizebytes(x->x_fftwInForwardWindow, oldWindow * sizeof(t_sample), x->x_window * sizeof(t_sample));
    x->x_fftwInBackWindow = (t_sample *)t_resizebytes(x->x_fftwInBackWindow, oldWindow * sizeof(t_sample), x->x_window * sizeof(t_sample));

    fftwf_free(x->x_fftwOutForwardWindow);
    fftwf_free(x->x_fftwOutBackWindow);
    fftwf_destroy_plan(x->x_fftwPlanForwardWindow);
    fftwf_destroy_plan(x->x_fftwPlanBackWindow);

    // set up a new FFTW output buffer
    x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf + 1);
    x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf + 1);

    // FFTW plan
    x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
    x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

    x->x_blackman = (t_float *)t_resizebytes(x->x_blackman, oldWindow * sizeof(t_float), x->x_window * sizeof(t_float));
    x->x_cosine = (t_float *)t_resizebytes(x->x_cosine, oldWindow * sizeof(t_float), x->x_window * sizeof(t_float));
    x->x_hamming = (t_float *)t_resizebytes(x->x_hamming, oldWindow * sizeof(t_float), x->x_window * sizeof(t_float));
    x->x_hann = (t_float *)t_resizebytes(x->x_hann, oldWindow * sizeof(t_float), x->x_window * sizeof(t_float));

    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);

    for(i = 0; i < x->x_window; i++)
    {
        x->x_fftwInForwardWindow[i] = 0.0;
        x->x_fftwInBackWindow[i] = 0.0;
    }

    // resize x_listOut
    x->x_listOut = (t_atom *)t_resizebytes(x->x_listOut, (oldWindowHalf + 1) * sizeof(t_atom), (x->x_windowHalf + 1) * sizeof(t_atom));
}


static void specFlux_analyze(t_specFlux *x, t_floatarg start, t_floatarg n)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, window, startSamp, endSamp, startSampBack, endSampBack;
        t_float flux, *windowFuncPtr;

        startSamp = (start < 0) ? 0 : start;

        if(n)
            endSamp = startSamp + n - 1;
        else
            endSamp = startSamp + x->x_window - 1;

        if(endSamp >= x->x_arrayPoints - 1)
            endSamp = x->x_arrayPoints - 1;

        window = endSamp - startSamp + 1;

        if(endSamp <= startSamp)
        {
            post("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        if(x->x_window != window)
            specFlux_resizeWindow(x, x->x_window, window, startSamp, &endSamp);

        // construct forward analysis window
        for(i = 0, j = startSamp; j <= endSamp; i++, j++)
            x->x_fftwInForwardWindow[i] = x->x_vec[j].w_float;

        // do these sample start/end location calculations AFTER the potential call to resizeWindow(), as x->x_window may have changed
        if(startSamp>=x->x_separation)
            startSampBack = startSamp - x->x_separation;
        else
            startSampBack = 0;

        endSampBack = startSampBack + x->x_window - 1;

        // construct back analysis window x->x_separation frames earlier
        for(i=0, j=startSampBack; j<=endSampBack; i++, j++)
            x->x_fftwInBackWindow[i] = x->x_vec[j].w_float;

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

        // if x_windowFunction == 0, skip the windowing (rectangular)
        if(x->x_windowFunction != rectangular)
            for(i = 0; i < x->x_window; i++, windowFuncPtr++)
                x->x_fftwInForwardWindow[i] *= *windowFuncPtr;

        fftwf_execute(x->x_fftwPlanForwardWindow);

        // put the result of power calc back in x_fftwIn
        tIDLib_power(x->x_windowHalf + 1, x->x_fftwOutForwardWindow, x->x_fftwInForwardWindow);

        if(!x->x_powerSpectrum)
            tIDLib_mag(x->x_windowHalf + 1, x->x_fftwInForwardWindow);

        if(x->x_normalize)
            tIDLib_normal(x->x_windowHalf + 1, x->x_fftwInForwardWindow);

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
        if(x->x_windowFunction != rectangular)
            for(i = 0; i < x->x_window; i++, windowFuncPtr++)
                x->x_fftwInBackWindow[i] *= *windowFuncPtr;

        fftwf_execute(x->x_fftwPlanBackWindow);

        // put the result of power calc back in x_fftwIn
        tIDLib_power(x->x_windowHalf + 1, x->x_fftwOutBackWindow, x->x_fftwInBackWindow);

        if(!x->x_powerSpectrum)
            tIDLib_mag(x->x_windowHalf + 1, x->x_fftwInBackWindow);

        if(x->x_normalize)
            tIDLib_normal(x->x_windowHalf + 1, x->x_fftwInBackWindow);

        flux=0.0;

        for(i = 0; i <= x->x_windowHalf; i++)
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

        outlet_list(x->x_fluxList, 0, x->x_windowHalf + 1, x->x_listOut);
        outlet_float(x->x_flux, flux);
    }
}


static void specFlux_chain_fftData(t_specFlux *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;
    t_float flux;

    // for specFlux fftData in particular:
    // incoming fftData list should be 4*(N/2+1) elements long, so windowHalf is:
    windowHalf = argc-4;
    windowHalf *= 0.25;

    // make sure that windowHalf == x->x_windowHalf in order to avoid an out of bounds memory read in the tIDLib_ functions below. we won't resize all memory based on an incoming chain_ command with a different window size. instead, just throw an error and exit
    if(windowHalf != x->x_windowHalf)
    {
        pd_error(x, "%s: window size of chain_ message (%lu) does not match current window size (%lu)", x->x_objSymbol->s_name, windowHalf * 2, x->x_window);
        return;
    }

    // fill the x_fftwOut buffers with the incoming fftData list, for both real and imag elements
    // for specFlux in particular, the first 2*(N/2+1) elements in the atom list are for the FORWARD window complex results. The second set of 2*(N/2+1) elements are for the BACK window complex results.
    for(i = 0; i <= x->x_windowHalf; i++)
    {
        x->x_fftwOutForwardWindow[i][0] = atom_getfloat(argv + i);
        x->x_fftwOutForwardWindow[i][1] = atom_getfloat(argv + (x->x_windowHalf + 1) + i);
        x->x_fftwOutBackWindow[i][0] = atom_getfloat(argv+(x->x_window+2)+i);
        x->x_fftwOutBackWindow[i][1] = atom_getfloat(argv+(x->x_window+x->x_windowHalf+3)+i);
    }

    // put the result of power calc back in x_fftwIn
    tIDLib_power(x->x_windowHalf + 1, x->x_fftwOutForwardWindow, x->x_fftwInForwardWindow);
    tIDLib_power(x->x_windowHalf + 1, x->x_fftwOutBackWindow, x->x_fftwInBackWindow);

    if(!x->x_powerSpectrum)
    {
        tIDLib_mag(x->x_windowHalf + 1, x->x_fftwInForwardWindow);
        tIDLib_mag(x->x_windowHalf + 1, x->x_fftwInBackWindow);
    }

    if(x->x_normalize)
    {
        tIDLib_normal(x->x_windowHalf + 1, x->x_fftwInForwardWindow);
        tIDLib_normal(x->x_windowHalf + 1, x->x_fftwInBackWindow);
    }

    flux=0.0;

    for(i = 0; i <= x->x_windowHalf; i++)
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

    outlet_list(x->x_fluxList, 0, x->x_windowHalf + 1, x->x_listOut);
    outlet_float(x->x_flux, flux);
}


static void specFlux_chain_magSpec(t_specFlux *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;
    t_float flux;

    // for specFlux magSpec in particular:
    // incoming magSpec list should be 2*(N/2+1) elements long, so windowHalf is:
    windowHalf = argc - 2;
    windowHalf *= 0.5;

    // make sure that windowHalf == x->x_windowHalf in order to avoid an out of bounds memory read in the tIDLib_ functions below. we won't resize all memory based on an incoming chain_ command with a different window size. instead, just throw an error and exit
    if(windowHalf != x->x_windowHalf)
    {
        pd_error(x, "%s: window size of chain_ message (%lu) does not match current window size (%lu)", x->x_objSymbol->s_name, windowHalf * 2, x->x_window);
        return;
    }

    // fill the x_fftwIn buffers with the incoming magSpec lists
    // for specFlux in particular, the first N/2+1 elements in the atom list are for the FORWARD window magnitudes. The second set of N/2+1 elements are for the BACK window magnitudes.
    for(i = 0; i <= x->x_windowHalf; i++)
    {
        x->x_fftwInForwardWindow[i] = atom_getfloat(argv + i);
        x->x_fftwInBackWindow[i] = atom_getfloat(argv + (x->x_windowHalf + 1) + i);
    }

    flux=0.0;

    for(i = 0; i <= x->x_windowHalf; i++)
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

    outlet_list(x->x_fluxList, 0, x->x_windowHalf + 1, x->x_listOut);
    outlet_float(x->x_flux, flux);
}


// analyze the whole damn array
static void specFlux_bang(t_specFlux *x)
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
        specFlux_analyze(x, startSamp, window);
    }
}


static void specFlux_set(t_specFlux *x, t_symbol *s)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void specFlux_print(t_specFlux *x)
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


static void specFlux_samplerate(t_specFlux *x, t_floatarg sr)
{
    if(sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void specFlux_window(t_specFlux *x, t_floatarg w)
{
    t_sampIdx endSamp;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow() requires that
    endSamp = 0;

    specFlux_resizeWindow(x, x->x_window, w, 0, &endSamp);
}


static void specFlux_windowFunction(t_specFlux *x, t_floatarg f)
{
    f = (f < 0) ? 0 : f;
    f = (f > 4) ? 4 : f;
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


static void specFlux_powerSpectrum(t_specFlux *x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    if(x->x_powerSpectrum)
        post("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void specFlux_logSpectrum(t_specFlux *x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_logSpectrum = spec;

    if(x->x_logSpectrum)
        post("%s log spectrum enabled", x->x_objSymbol->s_name);
    else
        post("%s log spectrum disabled", x->x_objSymbol->s_name);
}


static void specFlux_normalize(t_specFlux *x, t_floatarg norm)
{
    norm = (norm < 0) ? 0 : norm;
    norm = (norm > 1) ? 1 : norm;
    x->x_normalize = norm;

    if(x->x_normalize)
        post("%s spectrum normalization ON.", x->x_objSymbol->s_name);
    else
        post("%s spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void specFlux_separation(t_specFlux *x, t_floatarg s)
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
        x->x_separation = s;

    post("%s frame separation: %i", x->x_objSymbol->s_name, x->x_separation);
}


static void specFlux_squaredDiff(t_specFlux *x, t_floatarg sd)
{
    sd = (sd<0)?0:sd;
    sd = (sd>1)?1:sd;
    x->x_squaredDiff = sd;

    if(x->x_squaredDiff)
        post("%s using squared difference", x->x_objSymbol->s_name);
    else
        post("%s using absolute value of difference", x->x_objSymbol->s_name);
}


static void specFlux_mode(t_specFlux *x, t_symbol *m)
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


static void *specFlux_new(t_symbol *s, int argc, t_atom *argv)
{
    t_specFlux *x = (t_specFlux *)pd_new(specFlux_class);
    t_float sepFloat;
    t_sampIdx i;
//	t_garray *a;

    x->x_flux = outlet_new(&x->x_obj, &s_float);
    x->x_fluxList = outlet_new(&x->x_obj, gensym("list"));

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch(argc)
    {
        case 2:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            sepFloat = atom_getfloat(argv + 1);
            if(sepFloat > TID_WINDOWSIZEDEFAULT)
            {
                post("%s WARNING: frame separation cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
                x->x_separation = TID_WINDOWSIZEDEFAULT*0.5;
            }
            else if(sepFloat < 0)
            {
                post("%s WARNING: frame separation must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
                x->x_separation = TID_WINDOWSIZEDEFAULT*0.5;
            }
            else
                x->x_separation = sepFloat;
            break;

        case 1:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            x->x_separation = TID_WINDOWSIZEDEFAULT*0.5;
            break;

        case 0:
            post("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze()
            x->x_arrayName = gensym("NOARRAYSPECIFIED");
            x->x_separation = TID_WINDOWSIZEDEFAULT*0.5;
            break;

        default:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post("%s WARNING: Too many arguments supplied. Using default frame separation of %i samples.", x->x_objSymbol->s_name, (t_sampIdx)(TID_WINDOWSIZEDEFAULT*0.5));
            x->x_separation = TID_WINDOWSIZEDEFAULT*0.5;
            break;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_window = TID_WINDOWSIZEDEFAULT;
    x->x_windowHalf = x->x_window * 0.5;
    x->x_windowFunction = blackman;
    x->x_normalize = false;
    x->x_powerSpectrum = false;
    x->x_logSpectrum = false;
    x->x_squaredDiff = false;
    x->x_mode = mFlux;

    x->x_fftwInForwardWindow = (t_sample *)t_getbytes(x->x_window * sizeof(t_sample));
    x->x_fftwInBackWindow = (t_sample *)t_getbytes(x->x_window * sizeof(t_sample));

    x->x_listOut = (t_atom *)t_getbytes((x->x_windowHalf + 1) * sizeof(t_atom));

    // set up the FFTW output buffer. Is there no function to initialize it?
    x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf + 1);
    x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf + 1);

    // DFT plan
    x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
    x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for(i = 0; i < x->x_window; i++)
     {
        x->x_fftwInForwardWindow[i] = 0.0;
        x->x_fftwInBackWindow[i] = 0.0;
    }

      x->x_blackman = (t_float *)t_getbytes(x->x_window * sizeof(t_float));
      x->x_cosine = (t_float *)t_getbytes(x->x_window * sizeof(t_float));
      x->x_hamming = (t_float *)t_getbytes(x->x_window * sizeof(t_float));
      x->x_hann = (t_float *)t_getbytes(x->x_window * sizeof(t_float));

     // initialize signal windowing functions
    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);

    return (x);
}


static void specFlux_free(t_specFlux *x)
{
    // free the list out memory
    t_freebytes(x->x_listOut, (x->x_windowHalf + 1) * sizeof(t_atom));

    // free FFTW stuff
    t_freebytes(x->x_fftwInForwardWindow, (x->x_window) * sizeof(t_sample));
    t_freebytes(x->x_fftwInBackWindow, (x->x_window) * sizeof(t_sample));
    fftwf_free(x->x_fftwOutForwardWindow);
    fftwf_free(x->x_fftwOutBackWindow);
    fftwf_destroy_plan(x->x_fftwPlanForwardWindow);
    fftwf_destroy_plan(x->x_fftwPlanBackWindow);

    // free the window memory
    t_freebytes(x->x_blackman, x->x_window * sizeof(t_float));
    t_freebytes(x->x_cosine, x->x_window * sizeof(t_float));
    t_freebytes(x->x_hamming, x->x_window * sizeof(t_float));
    t_freebytes(x->x_hann, x->x_window * sizeof(t_float));
}


void specFlux_setup(void)
{
    specFlux_class =
    class_new(
        gensym("specFlux"),
        (t_newmethod)specFlux_new,
        (t_method)specFlux_free,
        sizeof(t_specFlux),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)specFlux_new,
        gensym("timbreIDLib/specFlux"),
        A_GIMME,
        0
    );

    class_addbang(specFlux_class, specFlux_bang);

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_analyze,
        gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_chain_fftData,
        gensym("chain_fftData"),
        A_GIMME,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_chain_magSpec,
        gensym("chain_magSpec"),
        A_GIMME,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_set,
        gensym("set"),
        A_SYMBOL,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_print,
        gensym("print"),
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_samplerate,
        gensym("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_windowFunction,
        gensym("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_powerSpectrum,
        gensym("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_logSpectrum,
        gensym("log_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_normalize,
        gensym("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_separation,
        gensym("separation"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_squaredDiff,
        gensym("squared_diff"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        specFlux_class,
        (t_method)specFlux_mode,
        gensym("mode"),
        A_DEFSYMBOL,
        0
    );
}
