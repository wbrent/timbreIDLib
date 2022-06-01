/*

freq2mel - Calculate the mel frequency of a given frequency in Hertz

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* freq2mel_class;

typedef struct _freq2mel
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_freq2melFormula x_formula;
    t_outlet* x_melFreq;

} t_freq2mel;


/* ------------------------ freq2mel -------------------------------- */
static void freq2mel_calculate (t_freq2mel* x, t_float f)
{
    t_float freq;

    freq = f;

    if (freq>=0.0 && freq<=TID_MAXMELFREQ)
    {
        t_float melFreq;

        melFreq = tIDLib_freq2mel (freq);

        outlet_float (x->x_melFreq, melFreq);
    }
    else
        pd_error (x, "%s: frequency must be between 0 and %f Hz", x->x_objSymbol->s_name, TID_MAXMELFREQ);
}

static void* freq2mel_new (t_symbol* s, int argc, t_atom* argv)
{
    t_freq2mel* x = (t_freq2mel *)pd_new (freq2mel_class);
    x->x_melFreq = outlet_new (&x->x_obj, &s_float);

    x->x_objSymbol = s;

    // will use x->x_formula in future
    switch (argc)
    {
        case 1:
            x->x_formula = atom_getfloat (argv);
            break;
        default:
            x->x_formula = freq2melFormula0;
            break;
    }

    return (x);
}

void freq2mel_setup (void)
{
    freq2mel_class =
    class_new (
        gensym ("freq2mel"),
        (t_newmethod)freq2mel_new,
        0,
        sizeof (t_freq2mel),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)freq2mel_new,
        gensym ("timbreIDLib/freq2mel"),
        A_GIMME,
        0
    );

    class_addfloat (
        freq2mel_class,
        (t_method)freq2mel_calculate
    );

    class_sethelpsymbol (freq2mel_class, gensym ("tID-conversion"));
}
