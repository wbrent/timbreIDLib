/*

waveSlope~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *waveSlope_tilde_class;

typedef struct _waveSlope_tilde
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_uShortInt x_overlap;
    t_sampIdx x_window;
    t_bool x_normalize;
    double x_lastDspTime;
    t_sample *x_signalBuffer;
    t_float *x_analysisBuffer;
    t_outlet *x_slope;
    t_float x_f;

} t_waveSlope_tilde;


/* ------------------------ waveSlope~ -------------------------------- */

static void waveSlope_tilde_bang(t_waveSlope_tilde *x)
{
    t_sampIdx i, j, window, bangSample;
    t_float slope;
    double currentTime;

    window = x->x_window;

    currentTime = clock_gettimesince(x->x_lastDspTime);
    bangSample = roundf((currentTime / 1000.0) * x->x_sr);

    if(bangSample >= x->x_n)
        bangSample = x->x_n - 1;

    // construct analysis window using bangSample as the end of the window
    for(i = 0, j = bangSample; i < window; i++, j++)
        x->x_analysisBuffer[i] = x->x_signalBuffer[j];

    if(x->x_normalize)
        tIDLib_normalPeak(window, x->x_analysisBuffer);

    slope=0.0;

    slope = tIDLib_fitLineSlope(window, x->x_analysisBuffer);

    outlet_float(x->x_slope, slope);
}


static void waveSlope_tilde_window(t_waveSlope_tilde *x, t_floatarg w)
{
    t_sampIdx i, window;

    window = w;

    if(window < TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
    }

    x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window + x->x_n) * sizeof(t_sample), (window + x->x_n) * sizeof(t_sample));
    x->x_analysisBuffer = (t_sample *)t_resizebytes(x->x_analysisBuffer, x->x_window * sizeof(t_sample), window * sizeof(t_sample));

    x->x_window = window;

    // init signal buffer
    for(i = 0; i < x->x_window + x->x_n; i++)
        x->x_signalBuffer[i] = 0.0;

    // init analysis buffer
    for(i = 0; i < x->x_window; i++)
        x->x_analysisBuffer[i] = 0.0;

    post("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void waveSlope_tilde_overlap(t_waveSlope_tilde *x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1) ? 1 : o;

    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void waveSlope_tilde_normalize(t_waveSlope_tilde *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
    x->x_normalize = norm;

    if(x->x_normalize)
        post("%s waveform normalization ON.", x->x_objSymbol->s_name);
    else
        post("%s waveform normalization OFF.", x->x_objSymbol->s_name);
}


static void waveSlope_tilde_print(t_waveSlope_tilde *x)
{
    post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
}


static void *waveSlope_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_waveSlope_tilde *x = (t_waveSlope_tilde *)pd_new(waveSlope_tilde_class);
    t_sampIdx i;

    x->x_slope = outlet_new(&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch(argc)
    {
        case 1:
            x->x_window = atom_getfloat(argv);
            if(x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            break;

        default:
            post("%s WARNING: Too many arguments supplied. Using default window size of %i.", x->x_objSymbol->s_name, TID_WINDOWSIZEDEFAULT);
            x->x_window = TID_WINDOWSIZEDEFAULT;
            break;
    }


    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_normalize = false;
    x->x_lastDspTime = clock_getlogicaltime();

    x->x_signalBuffer = (t_sample *)t_getbytes((x->x_window + x->x_n) * sizeof(t_sample));
    x->x_analysisBuffer = (t_float *)t_getbytes(x->x_window * sizeof(t_float));

     for(i = 0; i < x->x_window + x->x_n; i++)
        x->x_signalBuffer[i] = 0.0;

     for(i = 0; i < x->x_window; i++)
        x->x_analysisBuffer[i] = 0.0;

    return (x);
}


static t_int *waveSlope_tilde_perform(t_int *w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_waveSlope_tilde *x = (t_waveSlope_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];

     // shift signal buffer contents back.
    for(i = 0; i < x->x_window; i++)
        x->x_signalBuffer[i] = x->x_signalBuffer[i+n];

    // write new block to end of signal buffer.
    for(i = 0; i < n; i++)
        x->x_signalBuffer[x->x_window + i] = in[i];

    x->x_lastDspTime = clock_getlogicaltime();

    return (w + 4);
}


static void waveSlope_tilde_dsp(t_waveSlope_tilde *x, t_signal **sp)
{
    dsp_add(
        waveSlope_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if(sp[0]->s_sr != x->x_sr * x->x_overlap)
    {
        x->x_sr = sp[0]->s_sr / x->x_overlap;
        x->x_lastDspTime = clock_getlogicaltime();
    };

// compare n to stored n and update/resize buffer if different
    if(sp[0]->s_n != x->x_n)
    {
        t_sampIdx i;

        x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window + x->x_n) * sizeof(t_sample), (x->x_window + sp[0]->s_n) * sizeof(t_sample));

        x->x_n = sp[0]->s_n;
        x->x_lastDspTime = clock_getlogicaltime();

        // init signal buffer
        for(i = 0; i < x->x_window + x->x_n; i++)
            x->x_signalBuffer[i] = 0.0;
    };
};

static void waveSlope_tilde_free(t_waveSlope_tilde *x)
{
    // free the input buffer memory
    t_freebytes(x->x_signalBuffer, (x->x_window + x->x_n) * sizeof(t_sample));

    // free the analysis buffer memory
    t_freebytes(x->x_analysisBuffer, x->x_window * sizeof(t_float));

}

void waveSlope_tilde_setup(void)
{
    waveSlope_tilde_class =
    class_new(
        gensym("waveSlope~"),
        (t_newmethod)waveSlope_tilde_new,
        (t_method)waveSlope_tilde_free,
        sizeof(t_waveSlope_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)waveSlope_tilde_new,
        gensym("timbreIDLib/waveSlope~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN(waveSlope_tilde_class, t_waveSlope_tilde, x_f);

    class_addbang(waveSlope_tilde_class, waveSlope_tilde_bang);

    class_addmethod(
        waveSlope_tilde_class,
        (t_method)waveSlope_tilde_print,
        gensym("print"),
        0
    );

    class_addmethod(
        waveSlope_tilde_class,
        (t_method)waveSlope_tilde_normalize,
        gensym("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        waveSlope_tilde_class,
        (t_method)waveSlope_tilde_window,
        gensym("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        waveSlope_tilde_class,
        (t_method)waveSlope_tilde_overlap,
        gensym("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        waveSlope_tilde_class,
        (t_method)waveSlope_tilde_dsp,
        gensym("dsp"),
        0
    );
}
