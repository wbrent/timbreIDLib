/*

attackTime~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* attackTime_tilde_class;

typedef struct _attackTime_tilde
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_sampIdx x_window;
    t_uShortInt x_overlap;
    double x_lastDspTime;
    t_sample* x_signalBuffer;
    t_float* x_analysisBuffer;
    t_float* x_searchBuffer;
    t_uShortInt x_numSampsThresh;
    t_float x_sampMagThresh;
    t_sampIdx x_maxSearchRange;
    t_outlet* x_peakSampIdx;
    t_outlet* x_attackStartIdx;
    t_outlet* x_attackTime;
    t_float x_f;
} t_attackTime_tilde;


/* ------------------------ attackTime~ -------------------------------- */

static void attackTime_tilde_bang (t_attackTime_tilde* x)
{
    t_sampIdx i, j, window, bangSample, startSample, peakSampIdx, attackStartIdx;
    t_float attackTime, peakSampVal;
    double currentTime;

    window = x->x_window;

    currentTime = clock_gettimesince (x->x_lastDspTime);
    bangSample = roundf ((currentTime / 1000.0) * x->x_sr);

    if (bangSample >= x->x_n)
        bangSample = x->x_n - 1;

    // took a while to get this calculation right, but it seems correct now. remember that bangSample is always between 0 and 63 (or x_n - 1), and finding startSample within x_signalBuffer involves a few other steps.
    startSample = (x->x_maxSearchRange + x->x_n) - bangSample - window - 1;

    // construct analysis window
    for (i = 0, j = startSample; i < window; i++, j++)
        x->x_analysisBuffer[i] = x->x_signalBuffer[j];

    tIDLib_peakSample (window, x->x_analysisBuffer, &peakSampIdx, &peakSampVal);

    peakSampIdx += startSample; // add startSample back so we can find the peak sample index relative to x_signalBuffer

    i = x->x_maxSearchRange;
    j = peakSampIdx;

    while (i--)
    {
        if (j == 0)
            x->x_searchBuffer[i] = x->x_signalBuffer[j];
        else
        {
            x->x_searchBuffer[i] = x->x_signalBuffer[j];
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
        attackTime = (x->x_maxSearchRange - attackStartIdx) / x->x_sr;
        attackTime *= 1000.0; // convert seconds to milliseconds
        // overwrite attackStartIdx to be the index relative to the entire table
        attackStartIdx = peakSampIdx - (x->x_maxSearchRange - attackStartIdx);
    }

    outlet_float (x->x_peakSampIdx, peakSampIdx);
    outlet_float (x->x_attackStartIdx, attackStartIdx);
    outlet_float (x->x_attackTime, attackTime);
}


static void attackTime_tilde_print (t_attackTime_tilde* x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s maximum search range: %0.2f ms, %i samples", x->x_objSymbol->s_name, (x->x_maxSearchRange / x->x_sr) * 1000.0, x->x_maxSearchRange);
    post ("%s sample magnitude threshold for finding onset: %f", x->x_objSymbol->s_name, x->x_sampMagThresh);
    post ("%s minimum sample threshold for finding onset: %i", x->x_objSymbol->s_name, x->x_numSampsThresh);
}


static void attackTime_tilde_maxSearchRange (t_attackTime_tilde* x, t_floatarg range)
{
    t_sampIdx i, newRange;

    range = (range<5.0)?5.0:range;
    newRange = roundf ((range / 1000.0) * x->x_sr);

    x->x_searchBuffer = (t_float *)t_resizebytes (x->x_searchBuffer, x->x_maxSearchRange * sizeof (t_float), newRange * sizeof (t_float));
    x->x_signalBuffer = (t_float *)t_resizebytes (x->x_signalBuffer, (x->x_maxSearchRange + x->x_n) * sizeof (t_float), (newRange + x->x_n) * sizeof (t_float));

    x->x_maxSearchRange = newRange;

    for (i = 0; i < x->x_maxSearchRange; i++)
        x->x_searchBuffer[i] = 0.0;

    for (i = 0; i < x->x_maxSearchRange + x->x_n; i++)
        x->x_signalBuffer[i] = 0.0;

    post ("%s maximum search range: %0.2f ms, %i samples", x->x_objSymbol->s_name, (x->x_maxSearchRange / x->x_sr) * 1000.0, x->x_maxSearchRange);

}


static void attackTime_tilde_sampMagThresh (t_attackTime_tilde* x, t_floatarg thresh)
{
    x->x_sampMagThresh = (thresh < 0.0) ? 0.0 : thresh;

    post ("%s sample magnitude threshold for finding onset: %f", x->x_objSymbol->s_name, x->x_sampMagThresh);
}


static void attackTime_tilde_numSampsThresh (t_attackTime_tilde* x, t_floatarg thresh)
{
    thresh = (thresh < 0.0) ? 0.0 : thresh;
    thresh = (thresh > x->x_maxSearchRange)?x->x_maxSearchRange:thresh;
    x->x_numSampsThresh = thresh;

    post ("%s minimum sample threshold for finding onset: %i", x->x_objSymbol->s_name, x->x_numSampsThresh);
}


static void attackTime_tilde_window (t_attackTime_tilde* x, t_floatarg w)
{
    t_sampIdx i, window;

    if (w < TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
    }
    else
        window = w;

    x->x_analysisBuffer = (t_float *)t_resizebytes (x->x_analysisBuffer, x->x_window * sizeof (t_float), window * sizeof (t_float));

    x->x_window = window;

     for (i = 0; i < x->x_window; i++)
        x->x_analysisBuffer[i] = 0.0;

    post ("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void attackTime_tilde_overlap (t_attackTime_tilde* x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1) ? 1 : o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void* attackTime_tilde_new (t_symbol* s, int argc, t_atom* argv)
{
    t_attackTime_tilde* x = (t_attackTime_tilde *)pd_new (attackTime_tilde_class);
    t_sampIdx i;

    x->x_attackTime = outlet_new (&x->x_obj, &s_float);
    x->x_attackStartIdx = outlet_new (&x->x_obj, &s_float);
    x->x_peakSampIdx = outlet_new (&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 1:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            break;

        default:
            post ("%s WARNING: Too many arguments supplied. Using default window size of %i.", x->x_objSymbol->s_name, TID_WINDOWSIZEDEFAULT);
            x->x_window = TID_WINDOWSIZEDEFAULT;
            break;
    }

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_lastDspTime = clock_getlogicaltime();

    x->x_numSampsThresh = 10;
    x->x_sampMagThresh = 0.005;
    x->x_maxSearchRange = x->x_sr * 2.0; // two seconds

    x->x_signalBuffer = (t_sample *)t_getbytes ((x->x_maxSearchRange + x->x_n) * sizeof (t_sample));
    x->x_analysisBuffer = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
    x->x_searchBuffer = (t_float *)t_getbytes (x->x_maxSearchRange * sizeof (t_float));

    // initialize signal buffer
    for (i = 0; i < x->x_maxSearchRange + x->x_n; i++)
        x->x_signalBuffer[i] = 0.0;

    // initialize analysis buffer
    for (i = 0; i < x->x_window; i++)
        x->x_analysisBuffer[i] = 0.0;

     for (i = 0; i < x->x_maxSearchRange; i++)
        x->x_searchBuffer[i] = 0.0;

    return (x);
}


static t_int* attackTime_tilde_perform (t_int* w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_attackTime_tilde* x = (t_attackTime_tilde *)(w[1]);

    t_sample* in = (t_sample *)(w[2]);
    n = w[3];

     // shift signal buffer contents back.
    for (i = 0; i < x->x_maxSearchRange; i++)
        x->x_signalBuffer[i] = x->x_signalBuffer[i + n];

    // write new block to end of signal buffer.
    for (i = 0; i < n; i++)
        x->x_signalBuffer[x->x_maxSearchRange + i] = in[i];

    x->x_lastDspTime = clock_getlogicaltime();

    return (w + 4);
}


static void attackTime_tilde_dsp (t_attackTime_tilde* x, t_signal** sp)
{
    dsp_add (
        attackTime_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

    // compare sr to stored sr and update if different
    if (sp[0]->s_sr != x->x_sr * x->x_overlap)
    {
        t_sampIdx i, newRange;
        t_float range;

        range = (x->x_maxSearchRange / x->x_sr) * 1000.0;

        x->x_sr = sp[0]->s_sr / x->x_overlap;

        newRange = roundf ((range / 1000.0) * x->x_sr);

        x->x_searchBuffer = (t_float *)t_resizebytes (x->x_searchBuffer, x->x_maxSearchRange * sizeof (t_float), newRange * sizeof (t_float));
        x->x_signalBuffer = (t_float *)t_resizebytes (x->x_signalBuffer, (x->x_maxSearchRange + x->x_n) * sizeof (t_float), (newRange + x->x_n) * sizeof (t_float));

        x->x_maxSearchRange = newRange;

        for (i = 0; i < x->x_maxSearchRange; i++)
            x->x_searchBuffer[i] = 0.0;

        for (i = 0; i < x->x_maxSearchRange + x->x_n; i++)
            x->x_signalBuffer[i] = 0.0;

        post ("%s maximum search range: %0.2f ms, %i samples", x->x_objSymbol->s_name, (x->x_maxSearchRange / x->x_sr) * 1000.0, x->x_maxSearchRange);

        x->x_lastDspTime = clock_getlogicaltime();
    };

// compare n to stored n and update/resize buffer if different
    if (sp[0]->s_n != x->x_n)
    {
        t_sampIdx i;

        x->x_signalBuffer = (t_sample *)t_resizebytes (x->x_signalBuffer, (x->x_maxSearchRange + x->x_n) * sizeof (t_sample), (x->x_maxSearchRange + sp[0]->s_n) * sizeof (t_sample));

        x->x_n = sp[0]->s_n;

        // init signal buffer
        for (i = 0; i < (x->x_maxSearchRange + x->x_n); i++)
            x->x_signalBuffer[i] = 0.0;

        x->x_lastDspTime = clock_getlogicaltime();
    }
};


static void attackTime_tilde_free (t_attackTime_tilde* x)
{
    t_freebytes (x->x_analysisBuffer, x->x_window * sizeof (t_float));
    t_freebytes (x->x_signalBuffer, (x->x_maxSearchRange + x->x_n) * sizeof (t_sample));
    t_freebytes (x->x_searchBuffer, x->x_maxSearchRange * sizeof (t_float));
}


void attackTime_tilde_setup (void)
{
    attackTime_tilde_class =
    class_new (
        gensym ("attackTime~"),
        (t_newmethod)attackTime_tilde_new,
        (t_method)attackTime_tilde_free,
        sizeof (t_attackTime_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)attackTime_tilde_new,
        gensym ("timbreIDLib/attackTime~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (attackTime_tilde_class, t_attackTime_tilde, x_f);

    class_addbang (attackTime_tilde_class, attackTime_tilde_bang);

    class_addmethod (
        attackTime_tilde_class,
        (t_method)attackTime_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        attackTime_tilde_class,
        (t_method)attackTime_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_tilde_class,
        (t_method)attackTime_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_tilde_class,
        (t_method)attackTime_tilde_maxSearchRange,
        gensym ("max_search_range"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_tilde_class,
        (t_method)attackTime_tilde_sampMagThresh,
        gensym ("onset_samp_mag_thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_tilde_class,
        (t_method)attackTime_tilde_numSampsThresh,
        gensym ("onset_num_samps_thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        attackTime_tilde_class,
        (t_method)attackTime_tilde_dsp,
        gensym ("dsp"),
        0
    );
}
