/*

mel2freq - Calculate the frequency in Hertz of a given mel frequency

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* mel2freq_class;

typedef struct _mel2freq
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_mel2freqFormula x_formula;
    t_outlet* x_freq;

} t_mel2freq;


/* ------------------------ mel2freq -------------------------------- */
static void mel2freq_calculate (t_mel2freq* x, t_float m)
{
    t_float mel, freq;

    mel = m;

    if (mel>=0.0 && mel<=TID_MAXMELS)
    {
        freq = tIDLib_mel2freq (mel);
        outlet_float (x->x_freq, freq);
    }
    else
        pd_error (x, "%s: mel frequency must be between 0 and %f mels", x->x_objSymbol->s_name, TID_MAXMELS);
}

static void* mel2freq_new (t_symbol* s, int argc, t_atom* argv)
{
    t_mel2freq* x = (t_mel2freq *)pd_new (mel2freq_class);
    x->x_freq = outlet_new (&x->x_obj, &s_float);

    x->x_objSymbol = s;

    // will use x->x_formula in future
    switch (argc)
    {
        case 1:
            x->x_formula = atom_getfloat (argv);
            break;
        default:
            x->x_formula = mel2freqFormula0;
            break;
    }

    return (x);
}

void mel2freq_setup (void)
{
    mel2freq_class =
    class_new (
        gensym ("mel2freq"),
        (t_newmethod)mel2freq_new,
        0,
        sizeof (t_mel2freq),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)mel2freq_new,
        gensym ("timbreIDLib/mel2freq"),
        A_GIMME,
        0
    );

    class_addfloat (
        mel2freq_class,
        (t_method)mel2freq_calculate
    );

    class_sethelpsymbol (mel2freq_class, gensym ("tID-conversion"));
}
