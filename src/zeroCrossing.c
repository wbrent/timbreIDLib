/*

zeroCrossing

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *zeroCrossing_class;

typedef struct _zeroCrossing
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_bool x_normalize;
    t_float *x_analysisBuffer;
    t_word *x_vec;
    t_symbol *x_arrayName;
    t_sampIdx x_arrayPoints;
    t_outlet *x_crossings;
} t_zeroCrossing;


/* ------------------------ zeroCrossing -------------------------------- */

static void zeroCrossing_analyze (t_zeroCrossing *x, t_floatarg start, t_floatarg n)
{
    t_garray *a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, oldWindow, window, startSamp, endSamp;
        t_float crossings;

        startSamp = (start < 0) ? 0 : start;

        if (n)
            endSamp = startSamp + n - 1;
        else
            endSamp = startSamp + x->x_window - 1;

        if (endSamp >= x->x_arrayPoints)
            endSamp = x->x_arrayPoints - 1;

        window = endSamp - startSamp + 1;

        if (endSamp <= startSamp)
        {
            post ("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        if (x->x_window != window)
        {
            oldWindow = x->x_window;

            // window must be at least 4 points long
            if (window < TID_MINWINDOWSIZE)
            {
                window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            // hang on to these values for next time
            x->x_window = window;

            endSamp = startSamp + x->x_window - 1;
            if (endSamp > x->x_arrayPoints)
                endSamp = x->x_arrayPoints - 1;

            x->x_analysisBuffer = (t_float *)t_resizebytes (x->x_analysisBuffer, oldWindow * sizeof (t_float), x->x_window * sizeof (t_float));

        }

        // construct analysis window
        for (i = 0, j = startSamp; j <= endSamp; i++, j++)
            x->x_analysisBuffer[i] = x->x_vec[j].w_float;

        crossings=0.0;
        crossings = tIDLib_zeroCrossingRate(window, x->x_analysisBuffer, x->x_normalize);
        outlet_float (x->x_crossings, crossings);
    }
}


// analyze the whole damn array
static void zeroCrossing_bang (t_zeroCrossing *x)
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
        zeroCrossing_analyze (x, startSamp, window);
    }
}


static void zeroCrossing_set (t_zeroCrossing *x, t_symbol *s)
{
    t_garray *a;

    if ( !(a = (t_garray *)pd_findbyclass (s, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}

static void zeroCrossing_normalize (t_zeroCrossing *x, t_floatarg n)
{
    n = (n < 0) ? 0 : n;
    x->x_normalize = (n > 1) ? 1 : n;

    post ("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
}

static void zeroCrossing_print (t_zeroCrossing *x)
{
    post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
}


static void zeroCrossing_samplerate (t_zeroCrossing *x, t_floatarg sr)
{
    if (sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void *zeroCrossing_new (t_symbol *s, int argc, t_atom *argv)
{
    t_zeroCrossing *x = (t_zeroCrossing *)pd_new (zeroCrossing_class);
//	t_garray *a;

    x->x_crossings = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 1:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            break;

        case 0:
            post ("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze ()
            x->x_arrayName = gensym ("NOARRAYSPECIFIED");
            break;

        default:
            x->x_arrayName = atom_getsymbol (argv);
            /*
            if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
                pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
            else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
                pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
            */
            post ("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
            break;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_window = TID_WINDOWSIZEDEFAULT;
    x->x_normalize = false;

    x->x_analysisBuffer = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));

    return (x);
}


static void zeroCrossing_free (t_zeroCrossing *x)
{
    // free the input buffer memory
    t_freebytes (x->x_analysisBuffer, x->x_window * sizeof (t_sample));
}


void zeroCrossing_setup (void)
{
    zeroCrossing_class =
    class_new (
        gensym ("zeroCrossing"),
        (t_newmethod)zeroCrossing_new,
        (t_method)zeroCrossing_free,
        sizeof (t_zeroCrossing),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)zeroCrossing_new,
        gensym ("timbreIDLib/zeroCrossing"),
        A_GIMME,
        0
    );

    class_addbang (zeroCrossing_class, zeroCrossing_bang);

    class_addmethod (
        zeroCrossing_class,
        (t_method)zeroCrossing_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        zeroCrossing_class,
        (t_method)zeroCrossing_set,
        gensym ("set"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        zeroCrossing_class,
        (t_method)zeroCrossing_normalize,
        gensym ("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        zeroCrossing_class,
        (t_method)zeroCrossing_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        zeroCrossing_class,
        (t_method)zeroCrossing_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );
}
