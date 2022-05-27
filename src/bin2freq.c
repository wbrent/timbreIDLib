/*

bin2freq - Calculate the frequency of a bin number for a given window size and sample rate.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *bin2freq_class;

typedef struct _bin2freq
{
    t_object x_obj;
    t_float x_n;
    t_float x_sr;
    t_outlet *x_freq;

} t_bin2freq;


/* ------------------------ bin2freq -------------------------------- */
static void bin2freq_calculate(t_bin2freq *x, t_float bin)
{
    if(bin>=0.0 && bin<x->x_n)
    {
        t_float freq;

        freq = tIDLib_bin2freq(bin, x->x_n, x->x_sr);

        outlet_float(x->x_freq, freq);
    }
    else
    {
        pd_error(x, "bin2freq: bin number must be between 0 and %lu", (t_binIdx)(x->x_n - 1));
    }
}


static void bin2freq_setWinSampRate(t_bin2freq *x, t_float n, t_float sr)
{
    if(sr)
    {
        x->x_n = n;
        x->x_sr = sr;
    }
    else if(n)
    {
        x->x_n = n;
          x->x_sr = TID_SAMPLERATEDEFAULT;
    }
    else
    {
        x->x_n = TID_WINDOWSIZEDEFAULT;
        x->x_sr = TID_SAMPLERATEDEFAULT;
    }
}


static void bin2freq_print(t_bin2freq *x)
{
    post("bin2freq window size: %i", (t_sampIdx)x->x_n);
    post("bin2freq samplerate: %i", (t_sampIdx)x->x_sr);
}


static void *bin2freq_new(t_float n, t_float sr)
{
    t_bin2freq *x = (t_bin2freq *)pd_new(bin2freq_class);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("setWinSampRate"));
    x->x_freq = outlet_new(&x->x_obj, &s_float);

    bin2freq_setWinSampRate(x, n, sr);

    return (x);
}


void bin2freq_setup(void)
{
    bin2freq_class =
    class_new(
        gensym("bin2freq"),
        (t_newmethod)bin2freq_new,
        0,
        sizeof(t_bin2freq),
        CLASS_DEFAULT,
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addcreator(
        (t_newmethod)bin2freq_new,
        gensym("timbreIDLib/bin2freq"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        bin2freq_class,
        (t_method)bin2freq_print,
        gensym("print"),
        0
    );

    class_addfloat(
        bin2freq_class,
        (t_method)bin2freq_calculate
    );

    class_addmethod(
        bin2freq_class,
        (t_method)bin2freq_setWinSampRate,
        gensym("setWinSampRate"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_sethelpsymbol(bin2freq_class, gensym("tID-conversion"));
}
