/*

autoCorrPitch

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* autoCorrPitch_class;

typedef struct _autoCorrPitch
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_float* x_analysisBuffer;
    t_word* x_vec;
    t_symbol* x_arrayName;
    t_sampIdx x_arrayPoints;
    t_float x_thresh;
    t_outlet* x_pitchOut;
} t_autoCorrPitch;


/* ------------------------ autoCorrPitch -------------------------------- */

static void autoCorrPitch_analyze (t_autoCorrPitch* x, t_floatarg start, t_floatarg n)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, rhoBufferSize, oldWindow, window, startSamp, endSamp;
        t_float pitch;
        t_float* rhoBuffer;

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

        rhoBufferSize = (x->x_window * 2) - 1;

        rhoBuffer = (t_float *)t_getbytes (rhoBufferSize * sizeof (t_float));

        // init rhoBuffer with zeros
        for (i = 0; i < rhoBufferSize; i++)
            rhoBuffer[i] = 0.0;

        tIDLib_autoCorr (x->x_window, x->x_analysisBuffer, rhoBufferSize, rhoBuffer, false);

        // pass the rhoBuffer to tIDLib_autoCorrPitch to get the estimated period in samples
        pitch = tIDLib_autoCorrPeriod (rhoBufferSize, rhoBuffer, x->x_window, x->x_thresh);

        // if tIDLib_autoCorrPitch() returns a period of 0 samples, we output the failure value of -1500
        if ( !pitch)
        {
            pitch = -1500.0;
        }
        else
        {
            // divide the sampling rate by the period in samples to get the frequency estimate in Hz
            pitch = x->x_sr / pitch;

            // convert Hz to MIDI units
            pitch = ftom (pitch);
        }

        // free local memory
        t_freebytes (rhoBuffer, rhoBufferSize * sizeof (t_float));

        outlet_float (x->x_pitchOut, pitch);
    }
}


// analyze the whole damn array
static void autoCorrPitch_bang (t_autoCorrPitch* x)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx window, startSamp;
        startSamp = 0;
        window = x->x_arrayPoints;
        autoCorrPitch_analyze (x, startSamp, window);
    }
}


static void autoCorrPitch_set (t_autoCorrPitch* x, t_symbol* s)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (s, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void autoCorrPitch_print (t_autoCorrPitch* x)
{
    post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s threshold: %0.2f", x->x_objSymbol->s_name, x->x_thresh);
}


static void autoCorrPitch_samplerate (t_autoCorrPitch* x, t_floatarg sr)
{
    if (sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void autoCorrPitch_threshold (t_autoCorrPitch* x, t_floatarg thresh)
{
    thresh = (thresh < 0.0) ? 0.0 : thresh;
    thresh = (thresh > 100.0) ? 100.0 : thresh;

    x->x_thresh = thresh;
}


static void* autoCorrPitch_new (t_symbol* s, int argc, t_atom* argv)
{
    t_autoCorrPitch* x = (t_autoCorrPitch *)pd_new (autoCorrPitch_class);
//	t_garray *a;

    x->x_pitchOut = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 1:
            x->x_arrayName = atom_getsymbol (argv);
            break;

        case 0:
            post ("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze()
            x->x_arrayName = gensym ("NOARRAYSPECIFIED");
            break;

        default:
            x->x_arrayName = atom_getsymbol (argv);
            post ("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
            break;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_window = TID_WINDOWSIZEDEFAULT;
    x->x_thresh = 60.0;

    x->x_analysisBuffer = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));

    return (x);
}


static void autoCorrPitch_free (t_autoCorrPitch* x)
{
    // free the input buffer memory
    t_freebytes (x->x_analysisBuffer, x->x_window * sizeof (t_sample));
}


void autoCorrPitch_setup (void)
{
    autoCorrPitch_class =
    class_new (
        gensym ("autoCorrPitch"),
        (t_newmethod)autoCorrPitch_new,
        (t_method)autoCorrPitch_free,
        sizeof (t_autoCorrPitch),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)autoCorrPitch_new,
        gensym ("timbreIDLib/autoCorrPitch"),
        A_GIMME,
        0
    );

    class_addbang (autoCorrPitch_class, autoCorrPitch_bang);

    class_addmethod (
        autoCorrPitch_class,
        (t_method)autoCorrPitch_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        autoCorrPitch_class,
        (t_method)autoCorrPitch_set,
        gensym ("set"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        autoCorrPitch_class,
        (t_method)autoCorrPitch_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        autoCorrPitch_class,
        (t_method)autoCorrPitch_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        autoCorrPitch_class,
        (t_method)autoCorrPitch_threshold,
        gensym ("threshold"),
        A_DEFFLOAT,
        0
    );
}
