/*

attackTime

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* attackTime_class;

typedef struct _attackTime
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_sampIdx x_window;
    t_float* x_analysisBuffer;
    t_word* x_vec;
    t_symbol* x_arrayName;
    t_sampIdx x_arrayPoints;
    t_uShortInt x_numSampsThresh;
    t_float x_sampMagThresh;
    t_sampIdx x_maxSearchRange;
    t_float* x_searchBuffer;
    t_outlet* x_peakSampIdx;
    t_outlet* x_attackStartIdx;
    t_outlet* x_attackTime;
} t_attackTime;


/* ------------------------ attackTime -------------------------------- */

static void attackTime_analyze (t_attackTime* x, t_floatarg start, t_floatarg n)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (x->x_arrayName, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
    else
    {
        t_sampIdx i, j, oldWindow, window, startSamp, endSamp, peakSampIdx, attackStartIdx;
        t_float attackTime, peakSampVal;

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

        tIDLib_peakSample (x->x_window, x->x_analysisBuffer, &peakSampIdx, &peakSampVal);
        peakSampIdx += startSamp; // add startSamp back so we can find the peak sample index relative to the whole array

        i = x->x_maxSearchRange;
        j = peakSampIdx;

        while (i--)
        {
            if (j == 0)
                x->x_searchBuffer[i] = x->x_vec[j].w_float;
            else
            {
                x->x_searchBuffer[i] = x->x_vec[j].w_float;
                j--;
            }
        }

        attackTime = 0.0;

        // send searchBuffer to routine to find the point where sample magnitude is below x_sampMagThresh for at least x_numSampsThresh samples
        attackStartIdx = tIDLib_findAttackStartSamp (x->x_maxSearchRange, x->x_searchBuffer, x->x_sampMagThresh, x->x_numSampsThresh);

        // if the index returned is ULONG_MAX, the search failed
        if (attackStartIdx==ULONG_MAX)
            attackTime = -1.0;
        else
        {
            // attack duration in samples is the end of buffer index (where the peak sample was) minus the start index
            attackTime = (x->x_maxSearchRange - attackStartIdx)/x->x_sr;
            attackTime *= 1000.0; // convert seconds to milliseconds
            // overwrite attackStartIdx to be the index relative to the entire table
            attackStartIdx = peakSampIdx - (x->x_maxSearchRange - attackStartIdx);
        }

        outlet_float (x->x_peakSampIdx, peakSampIdx);
        outlet_float (x->x_attackStartIdx, attackStartIdx);
        outlet_float (x->x_attackTime, attackTime);
    }
}


// analyze the whole damn array
static void attackTime_bang (t_attackTime* x)
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
        attackTime_analyze (x, startSamp, window);
    }
}


static void attackTime_set (t_attackTime* x, t_symbol* s)
{
    t_garray* a;

    if ( !(a = (t_garray *)pd_findbyclass (s, garray_class)))
        pd_error (x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
    else if ( !garray_getfloatwords (a, (int *)&x->x_arrayPoints, &x->x_vec))
        pd_error (x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
    else
        x->x_arrayName = s;
}


static void attackTime_print (t_attackTime* x)
{
    post ("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s maximum search range: %0.2f ms, %i samples", x->x_objSymbol->s_name, (x->x_maxSearchRange/x->x_sr)*1000.0, x->x_maxSearchRange);
    post ("%s sample magnitude threshold for finding onset: %f", x->x_objSymbol->s_name, x->x_sampMagThresh);
    post ("%s minimum sample threshold for finding onset: %i", x->x_objSymbol->s_name, x->x_numSampsThresh);
}


static void attackTime_maxSearchRange (t_attackTime* x, t_floatarg range)
{
    t_sampIdx i, newRange;

    range = (range<5.0)?5.0:range;
    range /= 1000.0;
    newRange = roundf (range*x->x_sr);

    x->x_searchBuffer = (t_float *)t_resizebytes (x->x_searchBuffer, x->x_maxSearchRange * sizeof (t_float), newRange * sizeof (t_float));

    x->x_maxSearchRange = newRange;

    for (i = 0; i < x->x_maxSearchRange; i++)
        x->x_searchBuffer[i] = 0.0;

    post ("%s maximum search range: %0.2f ms, %i samples", x->x_objSymbol->s_name, (x->x_maxSearchRange/x->x_sr)*1000.0, x->x_maxSearchRange);

}


static void attackTime_sampMagThresh (t_attackTime* x, t_floatarg thresh)
{
    x->x_sampMagThresh = (thresh<0.0)?0.0:thresh;

    post ("%s sample magnitude threshold for finding onset: %f", x->x_objSymbol->s_name, x->x_sampMagThresh);
}


static void attackTime_numSampsThresh (t_attackTime* x, t_floatarg thresh)
{
    thresh = (thresh<0.0)?0.0:thresh;
    thresh = (thresh>x->x_maxSearchRange)?x->x_maxSearchRange:thresh;
    x->x_numSampsThresh = thresh;

    post ("%s minimum sample threshold for finding onset: %i", x->x_objSymbol->s_name, x->x_numSampsThresh);
}


static void attackTime_samplerate (t_attackTime* x, t_floatarg sr)
{
    t_float rangeMs;

    // get the old range in milliseconds since we'll need to change the search buffer size with a call to _maxSearchRange() below
    rangeMs = (x->x_maxSearchRange/x->x_sr)*1000.0;

    if (sr < TID_MINSAMPLERATE)
        x->x_sr = TID_MINSAMPLERATE;
    else
        x->x_sr = sr;

    attackTime_maxSearchRange (x, rangeMs);
}


static void* attackTime_new (t_symbol* s, int argc, t_atom* argv)
{
    t_attackTime* x = (t_attackTime *)pd_new (attackTime_class);
//	t_garray *a;

    x->x_attackTime = outlet_new (&x->x_obj, &s_float);
    x->x_attackStartIdx = outlet_new (&x->x_obj, &s_float);
    x->x_peakSampIdx = outlet_new (&x->x_obj, &s_float);

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
            // a bogus array name to trigger the safety check in _analyze()
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

    x->x_numSampsThresh = 10;
    x->x_sampMagThresh = 0.005;
    x->x_maxSearchRange = x->x_sr * 2.0; // two seconds

    x->x_searchBuffer = (t_float *)t_getbytes (x->x_maxSearchRange * sizeof (t_float));
    x->x_analysisBuffer = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));

    return (x);
}


static void attackTime_free (t_attackTime* x)
{
    // free the input buffer memory
    t_freebytes (x->x_searchBuffer, x->x_maxSearchRange * sizeof (t_float));
    t_freebytes (x->x_analysisBuffer, x->x_window * sizeof (t_sample));
}


void attackTime_setup (void)
{
    attackTime_class =
    class_new (
        gensym ("attackTime"),
        (t_newmethod)attackTime_new,
        (t_method)attackTime_free,
        sizeof (t_attackTime),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)attackTime_new,
        gensym ("timbreIDLib/attackTime"),
        A_GIMME,
        0
    );

    class_addbang (attackTime_class, attackTime_bang);

    class_addmethod (
        attackTime_class,
        (t_method)attackTime_analyze,
        gensym ("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_class,
        (t_method)attackTime_set,
        gensym ("set"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        attackTime_class,
        (t_method)attackTime_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        attackTime_class,
        (t_method)attackTime_maxSearchRange,
        gensym ("max_search_range"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_class,
        (t_method)attackTime_sampMagThresh,
        gensym ("onset_samp_mag_thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_class,
        (t_method)attackTime_numSampsThresh,
        gensym ("onset_num_samps_thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_class,
        (t_method)attackTime_samplerate,
        gensym ("samplerate"),
        A_DEFFLOAT,
        0
    );
}
