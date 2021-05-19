/*

freq2bin - Calculate the nearest bin number of a frequency for a given window size and sample rate.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *freq2bin_class;

typedef struct _freq2bin
{
    t_object x_obj;
    t_float x_n;
    t_float x_sr;
    t_outlet *x_bin;

} t_freq2bin;


/* ------------------------ freq2bin -------------------------------- */
static void freq2bin_calculate(t_freq2bin *x, t_float f)
{
    t_float freq;

    freq = f;

    if(freq>=0.0 && freq<x->x_sr)
    {
        t_float bin;

        bin = tIDLib_freq2bin(freq, x->x_n, x->x_sr);

        outlet_float(x->x_bin, bin);
    }
    else
    {
        pd_error(x, "freq2bin: frequency must be between 0 and %i", (int)x->x_sr);
    }
}

static void freq2bin_setWinSampRate(t_freq2bin *x, t_float n, t_float sr)
{
    if(sr)
    {
        x->x_n = n;
        x->x_sr = sr;
    }
    else if(n)
    {
        x->x_n = n;
          x->x_sr = SAMPLERATEDEFAULT;
    }
    else
    {
        x->x_n = WINDOWSIZEDEFAULT;
        x->x_sr = SAMPLERATEDEFAULT;
    }
}

static void freq2bin_print(t_freq2bin *x)
{
    post("freq2bin window size: %i", (t_sampIdx)x->x_n);
    post("freq2bin samplerate: %i", (t_sampIdx)x->x_sr);
}

static void *freq2bin_new(t_float n, t_float sr)
{
    t_freq2bin *x = (t_freq2bin *)pd_new(freq2bin_class);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("setWinSampRate"));
    x->x_bin = outlet_new(&x->x_obj, &s_float);

    freq2bin_setWinSampRate(x, n, sr);

    return (x);
}

void freq2bin_setup(void)
{
    freq2bin_class =
    class_new(
        gensym("freq2bin"),
        (t_newmethod)freq2bin_new,
        0,
        sizeof(t_freq2bin),
        CLASS_DEFAULT,
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addcreator(
        (t_newmethod)freq2bin_new,
        gensym("timbreIDLib/freq2bin"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        freq2bin_class,
        (t_method)freq2bin_print,
        gensym("print"),
        0
    );

    class_addfloat(
        freq2bin_class,
        (t_method)freq2bin_calculate
    );

    class_addmethod(
        freq2bin_class,
        (t_method)freq2bin_setWinSampRate,
        gensym("setWinSampRate"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_sethelpsymbol(freq2bin_class, gensym("tID-conversion"));
}
