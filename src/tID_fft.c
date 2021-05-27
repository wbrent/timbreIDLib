/*

tID_fft

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *tID_fft_class;

typedef struct _tID_fft
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_sampIdx x_zeroPad;
    t_sampIdx x_zeroPadHalf;
    t_float *x_fftwIn;
    t_float *x_fftwInZeroPad;
    fftwf_complex *x_fftwOut;
    fftwf_complex *x_fftwOutZeroPad;
    fftwf_plan x_fftwPlan;
    fftwf_plan x_fftwPlanZeroPad;
    t_bool x_normalize;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_word *x_vec;
    t_symbol *x_arrayName;
    t_sampIdx x_arrayPoints;
    t_atom *x_listOutReal;
    t_atom *x_listOutImag;
    t_atom *x_listOutRealZeroPad;
    t_atom *x_listOutImagZeroPad;
    t_outlet *x_realOut;
    t_outlet *x_imagOut;
} t_tID_fft;


/* ------------------------ tID_fft -------------------------------- */
static void tID_fft_resizeWindow(t_tID_fft *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
{
    t_sampIdx oldWindowHalf, windowHalf;

    windowHalf = window * 0.5;
    oldWindowHalf = oldWindow*0.5;

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

    // resize x_listOutReal and x_listOutImag
    x->x_listOutReal = (t_atom *)t_resizebytes(x->x_listOutReal, (oldWindowHalf+1)*sizeof(t_atom), (x->x_windowHalf+1)*sizeof(t_atom));
    x->x_listOutImag = (t_atom *)t_resizebytes(x->x_listOutImag, (oldWindowHalf+1)*sizeof(t_atom), (x->x_windowHalf+1)*sizeof(t_atom));
}


static void tID_fft_analyze(t_tID_fft *x, t_floatarg start, t_floatarg n)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, window, startSamp, endSamp;
        t_float *windowFuncPtr, realMax, imagMax;

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
            tID_fft_resizeWindow(x, x->x_window, window, startSamp, &endSamp);

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

        if(x->x_zeroPad>0)
        {
            // place windowed signal in x_fftwInZeroPad
            for(i=0; i<x->x_window; i++)
                x->x_fftwInZeroPad[i] = x->x_fftwIn[i];

            // fill out the remaining space with zeros
            for(; i<x->x_zeroPad; i++)
                x->x_fftwInZeroPad[i] = 0.0;

            // execute the zero pad plan
            fftwf_execute(x->x_fftwPlanZeroPad);

            realMax = 0.0;
            imagMax = 0.0;

            if(x->x_normalize)
            {
                for(i=0; i<=x->x_zeroPadHalf; i++)
                {
                    t_float thisReal, thisImag;

                    thisReal = fabsf(x->x_fftwOutZeroPad[i][0]);
                    thisImag = fabsf(x->x_fftwOutZeroPad[i][1]);

                    if(thisReal>realMax)
                        realMax = thisReal;

                    if(thisImag>imagMax)
                        imagMax = thisImag;
                }

                realMax = (realMax == 0.0) ? 1.0 : realMax;
                imagMax = (imagMax == 0.0) ? 1.0 : imagMax;

                realMax = 1.0/realMax;
                imagMax = 1.0/imagMax;
            }

            for(i=0; i<=x->x_zeroPadHalf; i++)
            {
                if(x->x_normalize)
                {
                    SETFLOAT(x->x_listOutRealZeroPad+i, x->x_fftwOutZeroPad[i][0]*realMax);
                    SETFLOAT(x->x_listOutImagZeroPad+i, x->x_fftwOutZeroPad[i][1]*imagMax);
                }
                else
                {
                    SETFLOAT(x->x_listOutRealZeroPad+i, x->x_fftwOutZeroPad[i][0]);
                    SETFLOAT(x->x_listOutImagZeroPad+i, x->x_fftwOutZeroPad[i][1]);
                }
            }

            outlet_list(x->x_imagOut, 0, x->x_zeroPadHalf+1, x->x_listOutImagZeroPad);
            outlet_list(x->x_realOut, 0, x->x_zeroPadHalf+1, x->x_listOutRealZeroPad);
        }
        else
        {
            fftwf_execute(x->x_fftwPlan);

            realMax = 0.0;
            imagMax = 0.0;

            if(x->x_normalize)
            {
                for(i=0; i<=x->x_windowHalf; i++)
                {
                    t_float thisReal, thisImag;

                    thisReal = fabsf(x->x_fftwOut[i][0]);
                    thisImag = fabsf(x->x_fftwOut[i][1]);

                    if(thisReal>realMax)
                        realMax = thisReal;

                    if(thisImag>imagMax)
                        imagMax = thisImag;
                }

                realMax = (realMax == 0.0) ? 1.0 : realMax;
                imagMax = (imagMax == 0.0) ? 1.0 : imagMax;

                realMax = 1.0/realMax;
                imagMax = 1.0/imagMax;
            }

            for(i=0; i<=x->x_windowHalf; i++)
            {
                if(x->x_normalize)
                {
                    SETFLOAT(x->x_listOutReal+i, x->x_fftwOut[i][0]*realMax);
                    SETFLOAT(x->x_listOutImag+i, x->x_fftwOut[i][1]*imagMax);
                }
                else
                {
                    SETFLOAT(x->x_listOutReal+i, x->x_fftwOut[i][0]);
                    SETFLOAT(x->x_listOutImag+i, x->x_fftwOut[i][1]);
                }
            }

            outlet_list(x->x_imagOut, 0, x->x_windowHalf+1, x->x_listOutImag);
            outlet_list(x->x_realOut, 0, x->x_windowHalf+1, x->x_listOutReal);
        }
    }
}


// analyze the whole damn array
static void tID_fft_bang(t_tID_fft *x)
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
        tID_fft_analyze(x, startSamp, window);
    }
}


static void tID_fft_set(t_tID_fft *x, t_symbol *s)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void tID_fft_print(t_tID_fft *x)
{
    post("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);

    if(x->x_zeroPad>x->x_window)
        post("%s zero padding size: %i", x->x_objSymbol->s_name, x->x_zeroPad);
    else
        post("%s zero padding OFF", x->x_objSymbol->s_name);
}


static void tID_fft_samplerate(t_tID_fft *x, t_floatarg sr)
{
    if(sr<MINSAMPLERATE)
        x->x_sr = MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void tID_fft_window(t_tID_fft *x, t_floatarg w)
{
    t_sampIdx endSamp;

    // have to pass in an address to a dummy t_sampIdx value since _resizeWindow() requires that
    endSamp = 0;

    tID_fft_resizeWindow(x, x->x_window, w, 0, &endSamp);
}


static void tID_fft_windowFunction(t_tID_fft *x, t_floatarg f)
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


static void tID_fft_zeroPad(t_tID_fft *x, t_floatarg z)
{
    t_sampIdx i, oldZeroPad, oldZeroPadHalf;

    oldZeroPad = x->x_zeroPad;
    oldZeroPadHalf = x->x_zeroPadHalf;

    // if the requested zero pad size is less or equal to current window size, zero pad should be 0 (OFF)
    z = (z<=x->x_window)?0.0:z;
    x->x_zeroPad = z;
    x->x_zeroPadHalf = x->x_zeroPad * 0.5;

    // resize zeroPad FFTW input buffer
    x->x_fftwInZeroPad = (t_float *)t_resizebytes(x->x_fftwInZeroPad, oldZeroPad * sizeof(t_float), x->x_zeroPad * sizeof(t_float));

    // free the zeroPad FFTW output buffer, and re-malloc according to new window
    fftwf_free(x->x_fftwOutZeroPad);

    // destroy old plan
    fftwf_destroy_plan(x->x_fftwPlanZeroPad);

    // allocate new fftwf_complex memory for the plan based on new zero pad size
    x->x_fftwOutZeroPad = (fftwf_complex *) fftwf_alloc_complex(x->x_zeroPadHalf+1);

    // create a new DFT plan based on new window size
    x->x_fftwPlanZeroPad = fftwf_plan_dft_r2c_1d(x->x_zeroPad, x->x_fftwInZeroPad, x->x_fftwOutZeroPad, FFTWPLANNERFLAG);

    x->x_listOutRealZeroPad = (t_atom *)t_resizebytes(x->x_listOutRealZeroPad, (oldZeroPadHalf+1) * sizeof(t_atom), (x->x_zeroPadHalf+1) * sizeof(t_atom));
    x->x_listOutImagZeroPad = (t_atom *)t_resizebytes(x->x_listOutImagZeroPad, (oldZeroPadHalf+1) * sizeof(t_atom), (x->x_zeroPadHalf+1) * sizeof(t_atom));

    // we're supposed to initialize the input array after we create the plan
     for(i=0; i<x->x_zeroPad; i++)
        x->x_fftwInZeroPad[i] = 0.0;

    if(x->x_zeroPad>0)
        post("%s zero padding size: %i", x->x_objSymbol->s_name, x->x_zeroPad);
    else
        post("%s zero padding OFF", x->x_objSymbol->s_name);
}


static void tID_fft_normalize(t_tID_fft *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
    x->x_normalize = norm;

    if(x->x_normalize)
        post("%s normalization ON.", x->x_objSymbol->s_name);
    else
        post("%s normalization OFF.", x->x_objSymbol->s_name);
}


static void *tID_fft_new(t_symbol *s, int argc, t_atom *argv)
{
    t_tID_fft *x = (t_tID_fft *)pd_new(tID_fft_class);
    t_sampIdx i;
//	t_garray *a;

    x->x_realOut = outlet_new(&x->x_obj, gensym("list"));
    x->x_imagOut = outlet_new(&x->x_obj, gensym("list"));

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
    x->x_windowHalf = x->x_window*0.5;
    x->x_windowFunction = blackman;
    x->x_zeroPad = 0;
    x->x_zeroPadHalf = 0;
    x->x_normalize = false;

    x->x_fftwIn = (t_sample *)t_getbytes(x->x_window*sizeof(t_sample));
    x->x_fftwInZeroPad = (t_float *)t_getbytes(x->x_zeroPad * sizeof(t_float));
    x->x_listOutReal = (t_atom *)t_getbytes((x->x_windowHalf+1)*sizeof(t_atom));
    x->x_listOutImag = (t_atom *)t_getbytes((x->x_windowHalf+1)*sizeof(t_atom));
    x->x_listOutRealZeroPad = (t_atom *)t_getbytes((x->x_zeroPadHalf+1)*sizeof(t_atom));
    x->x_listOutImagZeroPad = (t_atom *)t_getbytes((x->x_zeroPadHalf+1)*sizeof(t_atom));

    // set up the FFTW output buffer. Is there no function to initialize it?
    x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
    x->x_fftwOutZeroPad = (fftwf_complex *)fftwf_alloc_complex(x->x_zeroPadHalf+1);

    // FFTW plan
    x->x_fftwPlan = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwIn, x->x_fftwOut, FFTWPLANNERFLAG); // FFTWPLANNERFLAG may be slower than FFTWPLANNERFLAG but more efficient after the first run?
    x->x_fftwPlanZeroPad = fftwf_plan_dft_r2c_1d(x->x_zeroPad, x->x_fftwInZeroPad, x->x_fftwOutZeroPad, FFTWPLANNERFLAG);

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


static void tID_fft_free(t_tID_fft *x)
{
    // free the list out memory
    t_freebytes(x->x_listOutReal, (x->x_windowHalf+1)*sizeof(t_atom));
    t_freebytes(x->x_listOutImag, (x->x_windowHalf+1)*sizeof(t_atom));
    t_freebytes(x->x_listOutRealZeroPad, (x->x_zeroPadHalf+1)*sizeof(t_atom));
    t_freebytes(x->x_listOutImagZeroPad, (x->x_zeroPadHalf+1)*sizeof(t_atom));

    // free FFTW stuff
    t_freebytes(x->x_fftwIn, (x->x_window)*sizeof(t_float));
    t_freebytes(x->x_fftwInZeroPad, (x->x_zeroPad)*sizeof(t_float));
    fftwf_free(x->x_fftwOut);
    fftwf_free(x->x_fftwOutZeroPad);
    fftwf_destroy_plan(x->x_fftwPlan);
    fftwf_destroy_plan(x->x_fftwPlanZeroPad);

    // free the window memory
    t_freebytes(x->x_blackman, x->x_window*sizeof(t_float));
    t_freebytes(x->x_cosine, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hamming, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hann, x->x_window*sizeof(t_float));
}


void tID_fft_setup(void)
{
    tID_fft_class =
    class_new(
        gensym("tID_fft"),
        (t_newmethod)tID_fft_new,
        (t_method)tID_fft_free,
        sizeof(t_tID_fft),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)tID_fft_new,
        gensym("timbreIDLib/tID_fft"),
        A_GIMME,
        0
    );

    class_addbang(tID_fft_class, tID_fft_bang);

    class_addmethod(
        tID_fft_class,
        (t_method)tID_fft_analyze,
        gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        tID_fft_class,
        (t_method)tID_fft_set,
        gensym("set"),
        A_SYMBOL,
        0
    );

    class_addmethod(
        tID_fft_class,
        (t_method)tID_fft_print,
        gensym("print"),
        0
    );

    class_addmethod(
        tID_fft_class,
        (t_method)tID_fft_samplerate,
        gensym("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        tID_fft_class,
        (t_method)tID_fft_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        tID_fft_class,
        (t_method)tID_fft_windowFunction,
        gensym("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        tID_fft_class,
        (t_method)tID_fft_zeroPad,
        gensym("zero_pad"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        tID_fft_class,
        (t_method)tID_fft_normalize,
        gensym("normalize"),
        A_DEFFLOAT,
        0
    );
}
