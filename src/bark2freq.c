/*

bark2freq

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *bark2freq_class;

typedef struct _bark2freq
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_bark2freqFormula x_formula;
    t_outlet *x_freq;

} t_bark2freq;


// TODO: could add an option for the other formula: 13*arctan(0.00076*freq) + 3.5*arctan((freq/7500)^2)
/* ------------------------ bark2freq -------------------------------- */
static void bark2freq_calculate(t_bark2freq *x, t_float b)
{
    t_float bark;

    bark = b;

    if (bark>=0.0 && bark<=TID_MAXBARKS)
    {
        t_float freq;

        freq = tIDLib_bark2freq(bark);
        outlet_float (x->x_freq, freq);
    }
    else
        pd_error (x, "%s: Bark frequency must be between 0 and %f Barks", x->x_objSymbol->s_name, TID_MAXBARKS);
}


static void *bark2freq_new (t_symbol *s, int argc, t_atom *argv)
{
    t_bark2freq *x = (t_bark2freq *)pd_new (bark2freq_class);
    x->x_freq = outlet_new (&x->x_obj, &s_float);

    x->x_objSymbol = s;

    // will use x->x_formula in future
    switch (argc)
    {
        case 1:
            x->x_formula = atom_getfloat (argv);
            break;
        default:
            x->x_formula = bark2freqFormula0;
            break;
    }

    return (x);
}


void bark2freq_setup (void)
{
    bark2freq_class =
    class_new (
        gensym ("bark2freq"),
        (t_newmethod)bark2freq_new,
        0,
        sizeof (t_bark2freq),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)bark2freq_new,
        gensym ("timbreIDLib/bark2freq"),
        A_GIMME,
        0
    );

    class_addfloat (
        bark2freq_class,
        (t_method)bark2freq_calculate
    );

    class_sethelpsymbol (bark2freq_class, gensym ("tID-conversion"));
}
