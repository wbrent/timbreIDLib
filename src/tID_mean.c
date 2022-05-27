/*

tID_mean - Calculate mean of a list of numbers.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *tID_mean_class;

typedef struct _tID_mean
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_outlet *x_mean;

} t_tID_mean;


/* ------------------------ tID_mean -------------------------------- */

static void tID_mean_calculate(t_tID_mean *x, t_symbol *s, int argc, t_atom *argv)
{
    t_float n, sum, mean, *input;
    t_sampIdx i;

    n = argc;

    if(n < 1)
        pd_error(x, "%s: too few elements in list.", x->x_objSymbol->s_name);
    else
    {
        // create local memory
        input = (t_float *)t_getbytes(n * sizeof(t_float));

        for(i = 0; i < n; i++)
            input[i] = atom_getfloat(argv + i);

        sum = 0.0;

        for(i = 0; i < n; i++)
            sum += input[i];

        mean = sum/n;

        outlet_float(x->x_mean, mean);

        // free local memory
        t_freebytes(input, n * sizeof(t_float));
    }
}

static void *tID_mean_new(void)
{
    t_tID_mean *x = (t_tID_mean *)pd_new(tID_mean_class);
    x->x_mean = outlet_new(&x->x_obj, &s_float);

    x->x_objSymbol = gensym("tID_mean");

    return (x);
}

void tID_mean_setup(void)
{
    tID_mean_class =
    class_new(
        gensym("tID_mean"),
        (t_newmethod)tID_mean_new,
        0,
        sizeof(t_tID_mean),
        CLASS_DEFAULT,
        0
    );

    class_addcreator(
        (t_newmethod)tID_mean_new,
        gensym("timbreIDLib/tID_mean"),
        A_GIMME,
        0
    );

    class_addlist(
        tID_mean_class,
        (t_method)tID_mean_calculate
    );
}
