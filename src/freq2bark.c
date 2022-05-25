/*

freq2bark - Calculate the Bark frequency of a given frequency in Hertz

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *freq2bark_class;

typedef struct _freq2bark
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_freq2barkFormula x_formula;
    t_outlet *x_barkFreq;

} t_freq2bark;


// TODO: could add an option for the other formula: 13*arctan(0.00076*freq) + 3.5*arctan((freq/7500)^2)
/* ------------------------ freq2bark -------------------------------- */
static void freq2bark_calculate(t_freq2bark *x, t_float f)
{
    t_float freq;

    freq = f;

    if(freq>=0.0 && freq<=TID_MAXBARKFREQ)
    {
        t_float barkFreq;

        barkFreq = tIDLib_freq2bark(freq);

        outlet_float(x->x_barkFreq, barkFreq);
    }
    else
        pd_error(x, "%s: frequency must be between 0 and %f Hz", x->x_objSymbol->s_name, TID_MAXBARKFREQ);
}

static void *freq2bark_new(t_symbol *s, int argc, t_atom *argv)
{
    t_freq2bark *x = (t_freq2bark *)pd_new(freq2bark_class);
    x->x_barkFreq = outlet_new(&x->x_obj, &s_float);

    x->x_objSymbol = s;

    // will use x->x_formula in future
    switch(argc)
    {
        case 1:
            x->x_formula = atom_getfloat(argv);
            break;
        default:
            x->x_formula = freq2barkFormula0;
            break;
    }

    return (x);
}

void freq2bark_setup(void)
{
    freq2bark_class =
    class_new(
        gensym("freq2bark"),
        (t_newmethod)freq2bark_new,
        0,
        sizeof(t_freq2bark),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)freq2bark_new,
        gensym("timbreIDLib/freq2bark"),
        A_GIMME,
        0
    );

    class_addfloat(
        freq2bark_class,
        (t_method)freq2bark_calculate
    );

    class_sethelpsymbol(freq2bark_class, gensym("tID-conversion"));
}
