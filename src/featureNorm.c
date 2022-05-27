/*

featureNorm - an external for normalizing features based on given min/max value attribute lists.

Copyright 2010 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *featureNorm_class;

typedef struct _featureNorm
{
    t_object  x_obj;
    t_symbol *x_objSymbol;
    t_float *x_minValues;
    t_float *x_maxValues;
    t_attributeIdx x_featureLength;
    t_bool x_mode;
    t_atom *x_listOut;
    t_outlet *x_featureList;
} t_featureNorm;


static void featureNorm_allocMem(t_featureNorm *x)
{
    // grab atom list memory
    x->x_listOut = (t_atom *)t_getbytes((x->x_featureLength) * sizeof(t_atom));

    // grab new memory
    x->x_minValues = (t_float *)t_getbytes((x->x_featureLength) * sizeof(t_float));
    x->x_maxValues = (t_float *)t_getbytes((x->x_featureLength) * sizeof(t_float));
}

static void featureNorm_initMem(t_featureNorm *x)
{
    t_attributeIdx i;

    // clear the atom list
    for(i = 0; i < x->x_featureLength; i++)
    {
        x->x_minValues[i] = 0.0;
        x->x_maxValues[i] = 0.0;
        SETFLOAT(x->x_listOut+i, 0.0);
    }
}

static void featureNorm_free(t_featureNorm *x)
{
    // free listOut memory
    t_freebytes(x->x_listOut, (x->x_featureLength) * sizeof(t_atom));

    // free the min/max values memory
    t_freebytes(x->x_minValues, (x->x_featureLength) * sizeof(t_float));
    t_freebytes(x->x_maxValues, (x->x_featureLength) * sizeof(t_float));
}

static void featureNorm_normalize(t_featureNorm *x, t_symbol *s, int argc, t_atom *argv)
{
    t_attributeIdx i;

    if(x->x_featureLength != (t_attributeIdx)argc)
    {
        pd_error(x, "%s: input length does not match current length setting. input ignored.", x->x_objSymbol->s_name);
        return;
    }
    else
    {
        switch(x->x_mode)
        {
            case 0:
                for(i = 0; i < x->x_featureLength; i++)
                {
                    t_float thisAtt, numer, denom;

                    thisAtt = atom_getfloat(argv + i);

                    thisAtt = (thisAtt<x->x_minValues[i])?x->x_minValues[i]:thisAtt;
                    thisAtt = (thisAtt>x->x_maxValues[i])?x->x_maxValues[i]:thisAtt;

                    numer = thisAtt - x->x_minValues[i];
                    denom = x->x_maxValues[i] - x->x_minValues[i];

                    if(denom<=0.0)
                    {
                        pd_error(x, "%s: invalid max/min values. input ignored.", x->x_objSymbol->s_name);
                        return;
                    }

                    SETFLOAT(x->x_listOut+i, numer/denom);
                }
                break;
            case 1:
                for(i = 0; i < x->x_featureLength; i++)
                {
                    t_float thisAtt, numer, denom;

                    thisAtt = atom_getfloat(argv + i);

                    thisAtt = (thisAtt<x->x_minValues[i])?x->x_minValues[i]:thisAtt;
                    thisAtt = (thisAtt>x->x_maxValues[i])?x->x_maxValues[i]:thisAtt;

                    numer = thisAtt - x->x_minValues[i];
                    denom = x->x_maxValues[i] - x->x_minValues[i];

                    if(denom<=0.0)
                    {
                        pd_error(x, "%s: invalid max/min values. input ignored.", x->x_objSymbol->s_name);
                        return;
                    }

                    SETFLOAT(x->x_listOut+i, (numer/denom)*2.0-1.0);
                }
                break;
            default:
                break;
        }
    }

    outlet_list(x->x_featureList, 0, x->x_featureLength, x->x_listOut);
}

static void featureNorm_minValues(t_featureNorm *x, t_symbol *s, int argc, t_atom *argv)
{
    t_attributeIdx i;

    if(x->x_featureLength != (t_attributeIdx)argc)
    {
        pd_error(x, "%s: min values list length does not match current length setting. value list ignored.", x->x_objSymbol->s_name);
        return;
    }
    else
    {
        for(i = 0; i < x->x_featureLength; i++)
            x->x_minValues[i] = atom_getfloat(argv + i);
    }
}

static void featureNorm_maxValues(t_featureNorm *x, t_symbol *s, int argc, t_atom *argv)
{
    t_attributeIdx i;

    if(x->x_featureLength != (t_attributeIdx)argc)
    {
        pd_error(x, "%s: max values list length does not match current length setting. value list ignored.", x->x_objSymbol->s_name);
        return;
    }
    else
    {
        for(i = 0; i < x->x_featureLength; i++)
            x->x_maxValues[i] = atom_getfloat(argv + i);
    }
}

static void featureNorm_print(t_featureNorm *x)
{
    t_attributeIdx i;

    post("%s feature length: %i", x->x_objSymbol->s_name, x->x_featureLength);
    switch(x->x_mode)
    {
        case 0:
            post("%s mode: 0 to 1 normalization", x->x_objSymbol->s_name);
            break;
        case 1:
            post("%s mode: -1 to 1 normalization", x->x_objSymbol->s_name);
            break;
        default:
            post("%s mode: 0 to 1 normalization", x->x_objSymbol->s_name);
            break;
    }

    startpost("minimum values: ");
    for(i = 0; i < x->x_featureLength; i++)
        startpost("%f ", x->x_minValues[i]);

    endpost();

    startpost("maximum values: ");
    for(i = 0; i < x->x_featureLength; i++)
        startpost("%f ", x->x_maxValues[i]);

    endpost();
}

static void featureNorm_clear(t_featureNorm *x)
{
    featureNorm_free(x);
    featureNorm_allocMem(x);
    featureNorm_initMem(x);
}

static void featureNorm_length(t_featureNorm *x, t_floatarg len)
{
    len = (len<1)?1:len;

    // free memory first
    featureNorm_free(x);

    x->x_featureLength = len;

    featureNorm_allocMem(x);
    featureNorm_initMem(x);
}

static void featureNorm_mode(t_featureNorm *x, t_floatarg m)
{
    m = (m<0)?0:m;
    m = (m>1)?1:m;

    x->x_mode = m;
}

static void *featureNorm_new(t_symbol *s, int argc, t_atom *argv)
{
    t_featureNorm *x = (t_featureNorm *)pd_new(featureNorm_class);
    t_float featureLength;

    x->x_featureList = outlet_new(&x->x_obj, gensym("list"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("minValues"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("maxValues"));

    x->x_objSymbol = gensym("featureNorm");

    switch(argc)
    {
        case 0:
            x->x_featureLength = 50;
            x->x_mode = 0;
            break;
        case 1:
            featureLength = atom_getfloat(argv);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;
            x->x_mode = 0;
            break;
        default:
            featureLength = atom_getfloat(argv);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;
            x->x_mode = atom_getfloat(argv + 1);
            break;
    }

     featureNorm_allocMem(x);
     featureNorm_initMem(x);

    return (void *)x;
}

void featureNorm_setup(void) {

    featureNorm_class = class_new(gensym("featureNorm"),
        (t_newmethod)featureNorm_new,
        (t_method)featureNorm_free,
        sizeof(t_featureNorm),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)featureNorm_new,
        gensym("timbreIDLib/featureNorm"),
        A_GIMME,
        0
    );

    class_addlist(featureNorm_class, (t_method)featureNorm_normalize);

    class_addmethod(
        featureNorm_class,
        (t_method)featureNorm_normalize,
        gensym("normalize"),
        A_GIMME,
        0
    );

    class_addmethod(
        featureNorm_class,
        (t_method)featureNorm_minValues,
        gensym("minValues"),
        A_GIMME,
        0
    );

    class_addmethod(
        featureNorm_class,
        (t_method)featureNorm_maxValues,
        gensym("maxValues"),
        A_GIMME,
        0
    );

    class_addmethod(
        featureNorm_class,
        (t_method)featureNorm_print,
        gensym("print"),
        0
    );

    class_addmethod(
        featureNorm_class,
        (t_method)featureNorm_clear,
        gensym("clear"),
        0
    );

    class_addmethod(
        featureNorm_class,
        (t_method)featureNorm_length,
        gensym("length"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        featureNorm_class,
        (t_method)featureNorm_mode,
        gensym("mode"),
        A_DEFFLOAT,
        0
    );
}
