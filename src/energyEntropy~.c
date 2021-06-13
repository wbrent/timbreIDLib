/*

energyEntropy~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"
#define SUBWINDOWSPERMIDWINDOWDEFAULT 16
#define MINSUBWINDOWSPERMIDWINDOW 2

static t_class *energyEntropy_tilde_class;

typedef struct _energyEntropy_tilde
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_float x_n;
    t_uShortInt x_overlap;
    t_sampIdx x_subWindowSize;
    t_sampIdx x_subWindowsPerMidTermWindow;
    double x_lastDspTime;
    t_sample *x_signalBuffer;
    t_float *x_analysisBuffer;
    t_outlet *x_entropyOutlet;
    t_float x_f;

} t_energyEntropy_tilde;


/* ------------------------ energyEntropy~ -------------------------------- */

static void energyEntropy_tilde_bang(t_energyEntropy_tilde *x)
{
    t_sampIdx i, j, window, bangSample;
    t_float energyEntropyResult;
    double currentTime;

    window = x->x_subWindowSize;

    currentTime = clock_gettimesince(x->x_lastDspTime);
    bangSample = roundf((currentTime/1000.0)*x->x_sr);

    if(bangSample >= x->x_n)
        bangSample = x->x_n-1;

    // construct analysis window using bangSample as the end of the window
    for(i=0, j=bangSample; i<window; i++, j++)
        x->x_analysisBuffer[i] = x->x_signalBuffer[j];

    energyEntropyResult = 0.0;

    // energyEntropyResult = tIDLib_sigenergyEntropy (window, x->x_analysisBuffer, x->x_normalize, !x->x_power, x->x_db);

    outlet_float(x->x_entropyOutlet, energyEntropyResult);
}


static void energyEntropy_tilde_subWindow(t_energyEntropy_tilde *x, t_floatarg w)
{
    t_sampIdx i, window;

    window = w;

    if(window<MINWINDOWSIZE)
    {
        window = WINDOWSIZEDEFAULT;
        post("%s WARNING: sub-window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
    }

    x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_subWindowSize+x->x_n) * sizeof(t_sample), (window+x->x_n) * sizeof(t_sample));
    x->x_analysisBuffer = (t_sample *)t_resizebytes(x->x_analysisBuffer, x->x_subWindowSize * sizeof(t_sample), window * sizeof(t_sample));

    x->x_subWindowSize = window;

    // init signal buffer
    for(i=0; i<(x->x_subWindowSize+x->x_n); i++)
        x->x_signalBuffer[i] = 0.0;

    // init analysis buffer
    for(i=0; i<x->x_subWindowSize; i++)
        x->x_analysisBuffer[i] = 0.0;

    post("%s sub-window size: %i", x->x_objSymbol->s_name, x->x_subWindowSize);
}


static void energyEntropy_tilde_midTermWindow(t_energyEntropy_tilde *x, t_floatarg w)
{
    t_sampIdx i, subsPerMid;

    subsPerMid = w;

    if(subsPerMid < MINSUBWINDOWSPERMIDWINDOW)
    {
        post("%s WARNING: cannot use %i sub-windows per mid-term window. Using default of %i sub-windows per mid-term window.", x->x_objSymbol->s_name, subsPerMid, SUBWINDOWSPERMIDWINDOWDEFAULT);

        subsPerMid = SUBWINDOWSPERMIDWINDOWDEFAULT;
    }

    x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_subWindowSize+x->x_n) * sizeof(t_sample), (subsPerMid+x->x_n) * sizeof(t_sample));
    x->x_analysisBuffer = (t_sample *)t_resizebytes(x->x_analysisBuffer, x->x_subWindowSize * sizeof(t_sample), subsPerMid * sizeof(t_sample));

    x->x_subWindowsPerMidTermWindow = subsPerMid;

    // init signal buffer
    for(i=0; i<(x->x_subWindowSize+x->x_n); i++)
        x->x_signalBuffer[i] = 0.0;

    // init analysis buffer
    for(i=0; i<x->x_subWindowSize; i++)
        x->x_analysisBuffer[i] = 0.0;

    post("%s sub-windows per mid-term window: %i", x->x_objSymbol->s_name, x->x_subWindowsPerMidTermWindow);
}


static void energyEntropy_tilde_overlap(t_energyEntropy_tilde *x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr/x->x_overlap;
    x->x_overlap = (o<1)?1:o;

    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void energyEntropy_tilde_print(t_energyEntropy_tilde *x)
{
    post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr/x->x_overlap));
    post("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post("%s sub-window size: %i", x->x_objSymbol->s_name, x->x_subWindowSize);
    post("%s sub-windows per mid-term window: %i", x->x_objSymbol->s_name, x->x_subWindowsPerMidTermWindow);
}


static void *energyEntropy_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_energyEntropy_tilde *x = (t_energyEntropy_tilde *)pd_new(energyEntropy_tilde_class);
    t_sampIdx i;

    x->x_entropyOutlet = outlet_new(&x->x_obj, &s_float);

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch(argc)
    {
        case 2:
            x->x_subWindowSize = atom_getfloat(argv);
            if (x->x_subWindowSize < MINWINDOWSIZE)
            {
                x->x_subWindowSize = WINDOWSIZEDEFAULT;
                post("%s WARNING: sub-window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
            }

            x->x_subWindowsPerMidTermWindow = atom_getfloat(argv + 1);
            if (x->x_subWindowsPerMidTermWindow < MINSUBWINDOWSPERMIDWINDOW)
            {
                post("%s WARNING: cannot use %i sub-windows per mid-term window. Using default of %i sub-windows per mid-term window.", x->x_objSymbol->s_name, x->x_subWindowsPerMidTermWindow, SUBWINDOWSPERMIDWINDOWDEFAULT);
                x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            }
            break;

        case 1:
            x->x_subWindowSize = atom_getfloat(argv);
            if (x->x_subWindowSize < MINWINDOWSIZE)
            {
                x->x_subWindowSize = WINDOWSIZEDEFAULT;
                post("%s WARNING: sub-window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
            }

            x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            break;

        case 0:
            x->x_subWindowSize = WINDOWSIZEDEFAULT;
            x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            break;

        default:
            post("%s WARNING: Too many arguments supplied. Using default sub-window size of %i and %i sub-windows per mid-term window.", x->x_objSymbol->s_name, WINDOWSIZEDEFAULT, SUBWINDOWSPERMIDWINDOWDEFAULT);
            x->x_subWindowSize = WINDOWSIZEDEFAULT;
            x->x_subWindowsPerMidTermWindow = SUBWINDOWSPERMIDWINDOWDEFAULT;
            break;
    }


    x->x_sr = SAMPLERATEDEFAULT;
    x->x_n = BLOCKSIZEDEFAULT;
    x->x_overlap = 1;
    x->x_lastDspTime = clock_getlogicaltime();

    x->x_signalBuffer = (t_sample *)t_getbytes((x->x_subWindowSize+x->x_n) * sizeof(t_sample));
    x->x_analysisBuffer = (t_float *)t_getbytes(x->x_subWindowSize * sizeof(t_float));

     for(i=0; i<(x->x_subWindowSize+x->x_n); i++)
        x->x_signalBuffer[i] = 0.0;

     for(i=0; i<x->x_subWindowSize; i++)
        x->x_analysisBuffer[i] = 0.0;

    return (x);
}


static t_int *energyEntropy_tilde_perform(t_int *w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_energyEntropy_tilde *x = (t_energyEntropy_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];

     // shift signal buffer contents back.
    for(i=0; i<x->x_subWindowSize; i++)
        x->x_signalBuffer[i] = x->x_signalBuffer[i+n];

    // write new block to end of signal buffer.
    for(i=0; i<n; i++)
        x->x_signalBuffer[x->x_subWindowSize+i] = in[i];

    x->x_lastDspTime = clock_getlogicaltime();

    return (w+4);
}


static void energyEntropy_tilde_dsp(t_energyEntropy_tilde *x, t_signal **sp)
{
    dsp_add(
        energyEntropy_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if(sp[0]->s_sr != (x->x_sr*x->x_overlap))
    {
        x->x_sr = sp[0]->s_sr/x->x_overlap;
        x->x_lastDspTime = clock_getlogicaltime();
    };

// compare n to stored n and update/resize buffer if different
    if(sp[0]->s_n != x->x_n)
    {
        t_sampIdx i;

        x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_subWindowSize+x->x_n) * sizeof(t_sample), (x->x_subWindowSize+sp[0]->s_n) * sizeof(t_sample));

        x->x_n = sp[0]->s_n;
        x->x_lastDspTime = clock_getlogicaltime();

        // init signal buffer
        for(i=0; i<(x->x_subWindowSize+x->x_n); i++)
            x->x_signalBuffer[i] = 0.0;
    };
};

static void energyEntropy_tilde_free(t_energyEntropy_tilde *x)
{
    // free the input buffer memory
    t_freebytes(x->x_signalBuffer, (x->x_subWindowSize+x->x_n)*sizeof(t_sample));

    // free the analysis buffer memory
    t_freebytes(x->x_analysisBuffer, x->x_subWindowSize*sizeof(t_float));

}

void energyEntropy_tilde_setup(void)
{
    energyEntropy_tilde_class =
    class_new(
        gensym("energyEntropy~"),
        (t_newmethod)energyEntropy_tilde_new,
        (t_method)energyEntropy_tilde_free,
        sizeof(t_energyEntropy_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)energyEntropy_tilde_new,
        gensym("timbreIDLib/energyEntropy~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN(energyEntropy_tilde_class, t_energyEntropy_tilde, x_f);

    class_addbang(energyEntropy_tilde_class, energyEntropy_tilde_bang);

    class_addmethod(
        energyEntropy_tilde_class,
        (t_method)energyEntropy_tilde_print,
        gensym("print"),
        0
    );

    class_addmethod(
        energyEntropy_tilde_class,
        (t_method)energyEntropy_tilde_subWindow,
        gensym("sub_window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        energyEntropy_tilde_class,
        (t_method)energyEntropy_tilde_midTermWindow,
        gensym("mid_term_window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        energyEntropy_tilde_class,
        (t_method)energyEntropy_tilde_overlap,
        gensym("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        energyEntropy_tilde_class,
        (t_method)energyEntropy_tilde_dsp,
        gensym("dsp"),
        0
    );
}
