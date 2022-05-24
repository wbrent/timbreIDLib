/*

dct

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *dct_class;

typedef struct _dct
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_windowFunction x_windowFunction;
    t_bool x_normalize;
    t_float *x_dctIn;
    t_float *x_dctOut;
    fftwf_plan x_fftwDctPlan;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_word *x_vec;
    t_symbol *x_arrayName;
    t_sampIdx x_arrayPoints;
    t_atom *x_listOut;
    t_outlet *x_dct;
} t_dct;


/* ------------------------ dct -------------------------------- */
static void dct_resizeWindow(t_dct *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
{
    if(window<MINWINDOWSIZE)
    {
        window = WINDOWSIZEDEFAULT;
        post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);

        *endSamp = startSamp + window-1;
        if(*endSamp >= x->x_arrayPoints)
            *endSamp = x->x_arrayPoints-1;
    }

    // hang on to these values for next time
    x->x_window = window;

    x->x_dctIn = (t_float *)t_resizebytes(x->x_dctIn, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_dctOut = (t_float *)t_resizebytes(x->x_dctOut, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

    fftwf_destroy_plan(x->x_fftwDctPlan);
    // FFTW plan
    x->x_fftwDctPlan = fftwf_plan_r2r_1d(x->x_window, x->x_dctIn, x->x_dctOut, FFTW_REDFT10, FFTWPLANNERFLAG);

    x->x_blackman = (t_float *)t_resizebytes(x->x_blackman, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_cosine = (t_float *)t_resizebytes(x->x_cosine, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_hamming = (t_float *)t_resizebytes(x->x_hamming, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
    x->x_hann = (t_float *)t_resizebytes(x->x_hann, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

    tIDLib_blackmanWindow(x->x_blackman, x->x_window);
    tIDLib_cosineWindow(x->x_cosine, x->x_window);
    tIDLib_hammingWindow(x->x_hamming, x->x_window);
    tIDLib_hannWindow(x->x_hann, x->x_window);

    // resize x_listOut
    x->x_listOut = (t_atom *)t_resizebytes(x->x_listOut, oldWindow*sizeof(t_atom), x->x_window*sizeof(t_atom));
}


static void dct_analyze(t_dct *x, t_floatarg start, t_floatarg n)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, window, startSamp, endSamp;
        t_float *windowFuncPtr;

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
            dct_resizeWindow(x, x->x_window, window, startSamp, &endSamp);

        // construct analysis window
        for(i=0, j=startSamp; j<=endSamp; i++, j++)
            x->x_dctIn[i] = x->x_vec[j].w_float;

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
                x->x_dctIn[i] *= *windowFuncPtr;

        fftwf_execute(x->x_fftwDctPlan);

        if(x->x_normalize)
            tIDLib_normalPeak(x->x_window, x->x_dctOut);

        for(i=0; i<x->x_window; i++)
            SETFLOAT(x->x_listOut+i, x->x_dctOut[i]);

        outlet_list(x->x_dct, 0, x->x_window, x->x_listOut);
    }
}


// analyze the whole damn array
static void dct_bang(t_dct *x)
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
        dct_analyze(x, startSamp, window);
    }
}


static void dct_set(t_dct *x, t_symbol *s)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void dct_print(t_dct *x)
{
    post("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
    post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
}


static void dct_samplerate(t_dct *x, t_floatarg sr)
{
    if(sr<MINSAMPLERATE)
        x->x_sr = MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void dct_window(t_dct *x, t_floatarg w)
{
    t_sampIdx endSamp;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow() requires that
    endSamp = 0;

    dct_resizeWindow(x, x->x_window, w, 0, &endSamp);
}


static void dct_windowFunction(t_dct *x, t_floatarg f)
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


static void dct_normalize(t_dct *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
    x->x_normalize = norm;

    if(x->x_normalize)
        post("%s normalization ON.", x->x_objSymbol->s_name);
    else
        post("%s normalization OFF.", x->x_objSymbol->s_name);
}


static void *dct_new(t_symbol *s, int argc, t_atom *argv)
{
    t_dct *x = (t_dct *)pd_new(dct_class);
    t_sampIdx i;
//	t_garray *a;

    x->x_dct = outlet_new(&x->x_obj, gensym("list"));

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

    x->x_sr = SAMPLERATEDEFAULT;
    x->x_window = WINDOWSIZEDEFAULT;
    x->x_windowFunction = blackman;
    x->x_normalize = true;

    x->x_dctIn = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
    x->x_dctOut = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
    x->x_listOut = (t_atom *)t_getbytes(x->x_window*sizeof(t_atom));

    // DCT plan. FFTW_REDFT10 is the DCT-II
    x->x_fftwDctPlan = fftwf_plan_r2r_1d(x->x_window, x->x_dctIn, x->x_dctOut, FFTW_REDFT10, FFTWPLANNERFLAG);

    for(i=0; i<x->x_window; i++)
    {
        x->x_dctIn[i] = 0.0;
        x->x_dctOut[i] = 0.0;
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

    return (x);
}


static void dct_free(t_dct *x)
{
    // free the list out memory
    t_freebytes(x->x_listOut, x->x_window*sizeof(t_atom));

    // free FFTW stuff
    t_freebytes(x->x_dctIn, x->x_window*sizeof(t_float));
    t_freebytes(x->x_dctOut, x->x_window*sizeof(t_float));
    fftwf_destroy_plan(x->x_fftwDctPlan);

    // free the window memory
    t_freebytes(x->x_blackman, x->x_window*sizeof(t_float));
    t_freebytes(x->x_cosine, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hamming, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hann, x->x_window*sizeof(t_float));
}


void dct_setup(void)
{
    dct_class =
    class_new(
        gensym("dct"),
        (t_newmethod)dct_new,
        (t_method)dct_free,
        sizeof(t_dct),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)dct_new,
        gensym("timbreIDLib/dct"),
        A_GIMME,
        0
    );

    class_addbang(dct_class, dct_bang);

    class_addmethod(
        dct_class,
        (t_method)dct_analyze,
        gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        dct_class,
        (t_method)dct_set,
        gensym("set"),
        A_SYMBOL,
        0
    );

    class_addmethod(
        dct_class,
        (t_method)dct_print,
        gensym("print"),
        0
    );

    class_addmethod(
        dct_class,
        (t_method)dct_samplerate,
        gensym("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        dct_class,
        (t_method)dct_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        dct_class,
        (t_method)dct_windowFunction,
        gensym("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        dct_class,
        (t_method)dct_normalize,
        gensym("normalize"),
        A_DEFFLOAT,
        0
    );
}
