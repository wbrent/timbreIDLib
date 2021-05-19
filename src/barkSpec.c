/*

barkSpec

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *barkSpec_class;

typedef struct _barkSpec
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_bool x_normalize;
    t_bool x_powerSpectrum;
    t_float *x_fftwIn;
    fftwf_complex *x_fftwOut;
    fftwf_plan x_fftwPlan;
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
    t_atom *x_listOut;
    t_outlet *x_featureList;
} t_barkSpec;


/* ------------------------ barkSpec -------------------------------- */
static void barkSpec_resizeWindow(t_barkSpec *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
{
    t_sampIdx windowHalf;

    windowHalf = window * 0.5;

    // FFT must be at least MINWINDOWSIZE points long
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

    x->x_fftwIn = (t_float *)t_resizebytes(x->x_fftwIn, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

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

    // x_numFilters doesn't change with a change to x_window, so oldNumFilters and newNumFilters arguments are the same
    tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, x->x_numFilters, x->x_numFilters, x->x_window, x->x_sr);
}


static void barkSpec_analyze(t_barkSpec *x, t_floatarg start, t_floatarg n)
{
    t_sampIdx i, j, window, startSamp, endSamp;
    t_float *windowFuncPtr;
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        startSamp = (start<0.0)?0:start;

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
            barkSpec_resizeWindow(x, x->x_window, window, startSamp, &endSamp);

        // construct analysis window
        for(i=0, j=startSamp; j<=endSamp; i++, j++)
            x->x_fftwIn[i] = x->x_vec[j].w_float;

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

        if(x->x_specBandAvg)
            tIDLib_specFilterBands(x->x_windowHalf+1, x->x_numFilters, x->x_fftwIn, x->x_filterbank, x->x_normalize);
        else
            tIDLib_filterbankMultiply(x->x_fftwIn, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

        for(i=0; i<x->x_numFilters; i++)
            SETFLOAT(x->x_listOut+i, x->x_fftwIn[i]);

        outlet_list(x->x_featureList, 0, x->x_numFilters, x->x_listOut);
    }
}


static void barkSpec_chain_fftData(t_barkSpec *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;

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

    if(x->x_specBandAvg)
        tIDLib_specFilterBands(windowHalf+1, x->x_numFilters, x->x_fftwIn, x->x_filterbank, x->x_normalize);
    else
        tIDLib_filterbankMultiply(x->x_fftwIn, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

    for(i=0; i<x->x_numFilters; i++)
        SETFLOAT(x->x_listOut+i, x->x_fftwIn[i]);

    outlet_list(x->x_featureList, 0, x->x_numFilters, x->x_listOut);
}


static void barkSpec_chain_magSpec(t_barkSpec *x, t_symbol *s, int argc, t_atom *argv)
{
    t_sampIdx i, windowHalf;

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

    if(x->x_specBandAvg)
        tIDLib_specFilterBands(x->x_windowHalf+1, x->x_numFilters, x->x_fftwIn, x->x_filterbank, x->x_normalize);
    else
        tIDLib_filterbankMultiply(x->x_fftwIn, x->x_normalize, x->x_filterAvg, x->x_filterbank, x->x_numFilters);

    for(i=0; i<x->x_numFilters; i++)
        SETFLOAT(x->x_listOut+i, x->x_fftwIn[i]);

    outlet_list(x->x_featureList, 0, x->x_numFilters, x->x_listOut);
}


// analyze the whole damn array
static void barkSpec_bang(t_barkSpec *x)
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
        barkSpec_analyze(x, startSamp, window);
    }
}


static void barkSpec_createFilterbank(t_barkSpec *x, t_floatarg bs)
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


static void barkSpec_spec_band_avg(t_barkSpec *x, t_floatarg avg)
{
    avg = (avg<0)?0:avg;
    avg = (avg>1)?1:avg;
    x->x_specBandAvg = avg;

    if(x->x_specBandAvg)
        post("%s: averaging energy in spectrum bands.", x->x_objSymbol->s_name);
    else
        post("%s: using triangular filterbank.", x->x_objSymbol->s_name);
}


static void barkSpec_filter_avg(t_barkSpec *x, t_floatarg avg)
{
    avg = (avg<0)?0:avg;
    avg = (avg>1)?1:avg;
    x->x_filterAvg = avg;

    if(x->x_filterAvg)
        post("%s: averaging energy in triangular filters.", x->x_objSymbol->s_name);
    else
        post("%s: summing energy in triangular filters.", x->x_objSymbol->s_name);
}


static void barkSpec_set(t_barkSpec *x, t_symbol *s)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void barkSpec_print(t_barkSpec *x)
{
    post("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
    post("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post("%s Bark spacing: %f", x->x_objSymbol->s_name, x->x_barkSpacing);
    post("%s number of filters: %i", x->x_objSymbol->s_name, x->x_numFilters);
    post("%s spectrum band averaging: %i", x->x_objSymbol->s_name, x->x_specBandAvg);
    post("%s triangular filter averaging: %i", x->x_objSymbol->s_name, x->x_filterAvg);
}


static void barkSpec_samplerate(t_barkSpec *x, t_floatarg sr)
{
    if(sr<MINSAMPLERATE)
        x->x_sr = MINSAMPLERATE;
    else
        x->x_sr = sr;

    tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, x->x_numFilters, x->x_numFilters, x->x_window, x->x_sr);

}


static void barkSpec_window(t_barkSpec *x, t_floatarg w)
{
    t_sampIdx endSamp;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow() requires that
    endSamp = 0;

    barkSpec_resizeWindow(x, x->x_window, w, 0, &endSamp);
}


static void barkSpec_windowFunction(t_barkSpec *x, t_floatarg f)
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


static void barkSpec_powerSpectrum(t_barkSpec *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
    x->x_powerSpectrum = spec;

    if(x->x_powerSpectrum)
        post("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void barkSpec_normalize(t_barkSpec *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
    x->x_normalize = norm;

    if(x->x_normalize)
        post("%s spectrum normalization ON.", x->x_objSymbol->s_name);
    else
        post("%s spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void *barkSpec_new(t_symbol *s, int argc, t_atom *argv)
{
    t_barkSpec *x = (t_barkSpec *)pd_new(barkSpec_class);
    t_sampIdx i;
//	t_garray *a;

    x->x_featureList = outlet_new(&x->x_obj, gensym("list"));

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
            x->x_barkSpacing = atom_getfloat(argv+1);
            if(x->x_barkSpacing<MINBARKSPACING || x->x_barkSpacing>MAXBARKSPACING)
            {
                x->x_barkSpacing = BARKSPACINGDEFAULT;
                post("%s WARNING: Bark spacing must be between %f and %f Barks. Using default spacing of %f instead.", x->x_objSymbol->s_name, MINBARKSPACING, MAXBARKSPACING, BARKSPACINGDEFAULT);
            }
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
            break;

        case 0:
            post("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze()
            x->x_arrayName = gensym("NOARRAYSPECIFIED");
            x->x_barkSpacing = BARKSPACINGDEFAULT;
            break;

        default:
            x->x_arrayName = atom_getsymbol(argv);
            /*
            if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
                pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post("%s WARNING: Too many arguments supplied. Using default Bark spacing of %f.", x->x_objSymbol->s_name, BARKSPACINGDEFAULT);
            x->x_barkSpacing = BARKSPACINGDEFAULT;
            break;
    }

    x->x_sr = SAMPLERATEDEFAULT;
    x->x_window = WINDOWSIZEDEFAULT;
    x->x_windowHalf = x->x_window*0.5;
    x->x_windowFunction = blackman;
    x->x_normalize = true;
    x->x_powerSpectrum = false;
    x->x_sizeFilterFreqs = 0;
    x->x_numFilters = 0; // this is just an init size that will be updated in createFilterbank anyway.
    x->x_specBandAvg = false;
    x->x_filterAvg = false;

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


static void barkSpec_free(t_barkSpec *x)
{
    t_filterIdx i;

    // free the listOut memory
    t_freebytes(x->x_listOut, x->x_numFilters*sizeof(t_atom));

    // free FFTW stuff
    t_freebytes(x->x_fftwIn, (x->x_window)*sizeof(t_float));
    fftwf_free(x->x_fftwOut);
    fftwf_destroy_plan(x->x_fftwPlan);

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


void barkSpec_setup(void)
{
    barkSpec_class =
    class_new(
        gensym("barkSpec"),
        (t_newmethod)barkSpec_new,
        (t_method)barkSpec_free,
        sizeof(t_barkSpec),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)barkSpec_new,
        gensym("timbreIDLib/barkSpec"),
        A_GIMME,
        0
    );

    class_addbang(barkSpec_class, barkSpec_bang);

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_analyze,
        gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_chain_fftData,
        gensym("chain_fftData"),
        A_GIMME,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_chain_magSpec,
        gensym("chain_magSpec"),
        A_GIMME,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_set,
        gensym("set"),
        A_SYMBOL,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_print,
        gensym("print"),
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_samplerate,
        gensym("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_windowFunction,
        gensym("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_powerSpectrum,
        gensym("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_normalize,
        gensym("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_createFilterbank,
        gensym("filterbank"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_spec_band_avg,
        gensym("spec_band_avg"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        barkSpec_class,
        (t_method)barkSpec_filter_avg,
        gensym("filter_avg"),
        A_DEFFLOAT,
        0
    );
}
