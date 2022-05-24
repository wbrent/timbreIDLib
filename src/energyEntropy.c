/*

energyEntropy

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"
#define SUBWINDOWSPERMIDWINDOWDEFAULT 16
#define MINSUBWINDOWSPERMIDWINDOW 2

static t_class *energyEntropy_class;

typedef struct _energyEntropy
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_sampIdx x_subWindowSize;
    t_sampIdx x_subWindowsPerMidTermWindow;
    t_float* x_analysisBuffer;
    t_word* x_vec;
    t_symbol* x_arrayName;
    t_sampIdx x_arrayPoints;
    t_outlet* x_entropyOutlet;
} t_energyEntropy;


/* ------------------------ energyEntropy -------------------------------- */

static void energyEntropy_analyze (t_energyEntropy* x, t_floatarg start, t_floatarg n)
{
    t_garray *a;

    if (!(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if (!garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, bufLen, oldSubWindow, subWindow, startSamp, endSamp;
        t_float energyEntropyResult;

        startSamp = (start < 0) ? 0 : start;

        if (n)
            endSamp = (startSamp + (n * x->x_subWindowsPerMidTermWindow)) - 1;
        else
            endSamp = (startSamp + (x->x_subWindowSize * x->x_subWindowsPerMidTermWindow)) - 1;

        if (endSamp >= x->x_arrayPoints)
            endSamp = x->x_arrayPoints - 1;

        // TODO: maybe make a function to get startSamp/endSamp in tIDLib since all NRT objects use this code
        bufLen = endSamp - startSamp + 1;
        subWindow = bufLen / x->x_subWindowsPerMidTermWindow;

        if (endSamp <= startSamp)
        {
            post ("%s: bad range of samples.", x->x_objSymbol->s_name);
            return;
        }

        if (x->x_subWindowSize != subWindow)
        {
            oldSubWindow = x->x_subWindowSize;

            // window must be at least 4 points long
            if (subWindow < MINWINDOWSIZE)
            {
                subWindow = WINDOWSIZEDEFAULT;
                post ("%s WARNING: sub-window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
            }

            // hang on to this value for next time
            x->x_subWindowSize = subWindow;

            endSamp = (startSamp + (x->x_subWindowSize * x->x_subWindowsPerMidTermWindow)) - 1;

            if (endSamp > x->x_arrayPoints)
                endSamp = x->x_arrayPoints - 1;

            x->x_analysisBuffer = (t_float *)t_resizebytes (x->x_analysisBuffer, oldSubWindow * sizeof (t_float), x->x_subWindowSize * sizeof (t_float));
        }

        // construct analysis window
        for (i = 0, j = startSamp; j <= endSamp; i++, j++)
            x->x_analysisBuffer[i] = x->x_vec[j].w_float;

        energyEntropyResult = tIDLib_sigEnergyEntropy (x->x_subWindowSize, x->x_subWindowsPerMidTermWindow, x->x_analysisBuffer);

        outlet_float (x->x_entropyOutlet, energyEntropyResult);
    }
}


// analyze the whole damn array
static void energyEntropy_bang (t_energyEntropy* x)
{
    t_garray *a;

    if (!(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if (!garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx window, startSamp;
        startSamp = 0;
        window = x->x_arrayPoints / x->x_subWindowsPerMidTermWindow;
        energyEntropy_analyze (x, startSamp, window);
    }
}


static void energyEntropy_set (t_energyEntropy* x, t_symbol* s)
{
    t_garray *a;

    if (!(a = (t_garray *)pd_findbyclass (s, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if (!garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void energyEntropy_print (t_energyEntropy* x)
{
    post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post ("%s sub-window size: %i", x->x_objSymbol->s_name, x->x_subWindowSize);
    post ("%s sub-windows per mid-term window: %i", x->x_objSymbol->s_name, x->x_subWindowsPerMidTermWindow);
}


static void energyEntropy_samplerate (t_energyEntropy* x, t_floatarg sr)
{
    if (sr < MINSAMPLERATE)
        x->x_sr = MINSAMPLERATE;
    else
        x->x_sr = sr;
}


static void *energyEntropy_new (t_symbol* s, int argc, t_atom* argv)
{
    t_energyEntropy* x = (t_energyEntropy *)pd_new (energyEntropy_class);
    t_sampIdx i, bufLen;

    x->x_entropyOutlet = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 3:
            x->x_arrayName = atom_getsymbol (argv);

            x->x_subWindowSize = atom_getfloat (argv + 1);
            if (x->x_subWindowSize < MINWINDOWSIZE)
            {
                x->x_subWindowSize = WINDOWSIZEDEFAULT;
                post ("%s WARNING: sub-window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
            }

            x->x_subWindowsPerMidTermWindow = atom_getfloat (argv + 2);
            if (x->x_subWindowsPerMidTermWindow < MINSUBWINDOWSPERMIDWINDOW)
            {
                post ("%s WARNING: cannot use %i sub-windows per mid-term window. Using default of %i sub-windows per mid-term window.", x->x_objSymbol->s_name, x->x_subWindowsPerMidTermWindow, SUBWINDOWSPERMIDWINDOWDEFAULT);
                x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            }
            break;

        case 2:
            x->x_arrayName = atom_getsymbol (argv);

            x->x_subWindowSize = atom_getfloat (argv + 1);
            if (x->x_subWindowSize < MINWINDOWSIZE)
            {
                x->x_subWindowSize = WINDOWSIZEDEFAULT;
                post ("%s WARNING: sub-window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
            }

            x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            break;

        case 1:
            x->x_arrayName = atom_getsymbol (argv);
            x->x_subWindowSize = WINDOWSIZEDEFAULT;
            x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            break;

        case 0:
            post ("%s: no array specified.", x->x_objSymbol->s_name);
            // a bogus array name to trigger the safety check in _analyze()
            x->x_arrayName = gensym ("NOARRAYSPECIFIED");
            x->x_subWindowSize = WINDOWSIZEDEFAULT;
            x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            break;

        default:
            x->x_arrayName = atom_getsymbol (argv);
            post ("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
            x->x_subWindowSize = WINDOWSIZEDEFAULT;
            x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            break;
    }

    x->x_sr = SAMPLERATEDEFAULT;

    // analysisBuffer should be subWindowSize * subWindowsPerMidTermWindow
    bufLen = x->x_subWindowSize * x->x_subWindowsPerMidTermWindow;

    x->x_analysisBuffer = (t_sample *)t_getbytes (bufLen * sizeof (t_sample));

    for (i = 0; i < bufLen; i++)
        x->x_analysisBuffer[i] = 0.0;

    return (x);
}


static void energyEntropy_free (t_energyEntropy* x)
{
  t_sampIdx bufLen;

  bufLen = x->x_subWindowSize * x->x_subWindowsPerMidTermWindow;

  // free the analysis buffer memory
  t_freebytes (x->x_analysisBuffer, bufLen * sizeof (t_float));
}


void energyEntropy_setup (void)
{
    energyEntropy_class =
    class_new (
        gensym ("energyEntropy"),
        (t_newmethod)energyEntropy_new,
        (t_method)energyEntropy_free,
        sizeof (t_energyEntropy),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)energyEntropy_new,
        gensym ("timbreIDLib/energyEntropy"),
        A_GIMME,
        0
    );

    class_addbang (energyEntropy_class, energyEntropy_bang);

    class_addmethod (
        energyEntropy_class,
        (t_method)energyEntropy_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        energyEntropy_class,
        (t_method)energyEntropy_set,
        gensym ("set"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        energyEntropy_class,
        (t_method)energyEntropy_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        energyEntropy_class,
        (t_method)energyEntropy_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );
}
