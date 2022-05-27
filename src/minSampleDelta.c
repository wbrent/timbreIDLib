/*

minSampleDelta

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *minSampleDelta_class;

typedef struct _minSampleDelta
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_float *x_analysisBuffer;
    t_word *x_vec;
    t_symbol *x_arrayName;
    t_sampIdx x_arrayPoints;
    t_outlet *x_minSampleDeltaIdx;
    t_outlet *x_minSampleDelta;
} t_minSampleDelta;


/* ------------------------ minSampleDelta -------------------------------- */

static void minSampleDelta_analyze(t_minSampleDelta *x, t_floatarg start, t_floatarg n)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, oldWindow, window, startSamp, endSamp, minIdx;
        t_float min;

        startSamp = (start < 0) ? 0 : start;

        if(n)
            endSamp = startSamp + n - 1;
        else
            endSamp = startSamp + x->x_window - 1;

        if(endSamp >= x->x_arrayPoints)
            endSamp = x->x_arrayPoints - 1;

        window = endSamp - startSamp + 1;

        if(endSamp <= startSamp)
        {
            post("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        if(x->x_window != window)
        {
            oldWindow = x->x_window;

            // window must be at least 4 points long
            if(window < TID_MINWINDOWSIZE)
            {
                window = TID_WINDOWSIZEDEFAULT;
                post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            // hang on to these values for next time
            x->x_window = window;

            endSamp = startSamp + x->x_window - 1;
            if(endSamp > x->x_arrayPoints)
                endSamp = x->x_arrayPoints - 1;

            x->x_analysisBuffer = (t_float *)t_resizebytes(x->x_analysisBuffer, oldWindow * sizeof(t_float), x->x_window * sizeof(t_float));

        }

        // construct analysis window
        for(i = 0, j = startSamp; j <= endSamp; i++, j++)
            x->x_analysisBuffer[i] = x->x_vec[j].w_float;

        if(startSamp>0)
        {
            min = fabs(x->x_analysisBuffer[0] - x->x_vec[startSamp-1].w_float);
            minIdx = 0;
        }
        else
        {
            min = FLT_MAX;
            minIdx = ULONG_MAX;
        }

        for(i=1; i<x->x_window; i++)
        {
            t_float thisDiff;

            thisDiff = fabs(x->x_analysisBuffer[i] - x->x_analysisBuffer[i-1]);

            if(thisDiff<min)
            {
                min = thisDiff;
                minIdx = i;
            }
        }

        outlet_float(x->x_minSampleDeltaIdx, minIdx);
        outlet_float(x->x_minSampleDelta, min);
    }
}


// analyze the whole damn array
static void minSampleDelta_bang(t_minSampleDelta *x)
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
        minSampleDelta_analyze(x, startSamp, window);
    }
}


static void minSampleDelta_set(t_minSampleDelta *x, t_symbol *s)
{
    t_garray *a;

    if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void minSampleDelta_print(t_minSampleDelta *x)
{
    post("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
}


static void minSampleDelta_samplerate(t_minSampleDelta *x, t_floatarg sr)
{
    if(sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void *minSampleDelta_new(t_symbol *s, int argc, t_atom *argv)
{
    t_minSampleDelta *x = (t_minSampleDelta *)pd_new(minSampleDelta_class);
//	t_garray *a;

    x->x_minSampleDelta = outlet_new(&x->x_obj, &s_float);
    x->x_minSampleDeltaIdx = outlet_new(&x->x_obj, &s_float);

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

    x->x_analysisBuffer = (t_sample *)t_getbytes(x->x_window * sizeof(t_sample));

    return (x);
}


static void minSampleDelta_free(t_minSampleDelta *x)
{
    // free the input buffer memory
    t_freebytes(x->x_analysisBuffer, x->x_window * sizeof(t_sample));
}


void minSampleDelta_setup(void)
{
    minSampleDelta_class =
    class_new(
        gensym("minSampleDelta"),
        (t_newmethod)minSampleDelta_new,
        (t_method)minSampleDelta_free,
        sizeof(t_minSampleDelta),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)minSampleDelta_new,
        gensym("timbreIDLib/minSampleDelta"),
        A_GIMME,
        0
    );

    class_addbang(minSampleDelta_class, minSampleDelta_bang);

    class_addmethod(
        minSampleDelta_class,
        (t_method)minSampleDelta_analyze,
        gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        minSampleDelta_class,
        (t_method)minSampleDelta_set,
        gensym("set"),
        A_SYMBOL,
        0
    );

    class_addmethod(
        minSampleDelta_class,
        (t_method)minSampleDelta_print,
        gensym("print"),
        0
    );

    class_addmethod(
        minSampleDelta_class,
        (t_method)minSampleDelta_samplerate,
        gensym("samplerate"),
        A_DEFFLOAT,
        0
    );
}
