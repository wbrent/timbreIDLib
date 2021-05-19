/*

featureDelta - an external for taking the inter-frame delta of features based on incoming feature lists.

Copyright 2010 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *featureDelta_class;

typedef enum
{
    deltaDiff = 0,
    deltaAbs,
    deltaSquared
} t_deltaMode;

typedef enum
{
    deltaPos = 0,
    deltaNeg,
    deltaBoth
} t_deltaDirection;

typedef struct _featureDelta
{
    t_object  x_obj;
    t_symbol *x_objSymbol;
    t_float *x_prevFeature;
    t_attributeIdx x_featureLength;
    t_deltaMode x_mode;
    t_deltaDirection x_direction;
    t_atom *x_listOut;
    t_outlet *x_featureList;
} t_featureDelta;


static void featureDelta_allocMem(t_featureDelta *x)
{
    // grab atom list memory
    x->x_listOut = (t_atom *)t_getbytes((x->x_featureLength)*sizeof(t_atom));

    // grab new memory
    x->x_prevFeature = (t_float *)t_getbytes((x->x_featureLength)*sizeof(t_float));
}

static void featureDelta_initMem(t_featureDelta *x)
{
    t_attributeIdx i;

    // clear the atom list
    for(i=0; i<x->x_featureLength; i++)
    {
        x->x_prevFeature[i] = 0.0;
        SETFLOAT(x->x_listOut+i, 0.0);
    }
}

static void featureDelta_free(t_featureDelta *x)
{
    // free listOut memory
    t_freebytes(x->x_listOut, (x->x_featureLength)*sizeof(t_atom));

    // free the previous feature memory
    t_freebytes(x->x_prevFeature, (x->x_featureLength)*sizeof(t_float));
}

static void featureDelta_delta(t_featureDelta *x, t_symbol *s, int argc, t_atom *argv)
{
    t_attributeIdx i;

    if(x->x_featureLength != (t_attributeIdx)argc)
    {
        pd_error(x, "%s: input length does not match current length setting. input ignored.", x->x_objSymbol->s_name);
        return;
    }
    else
    {
        for(i=0; i<x->x_featureLength; i++)
        {
            t_float thisDiff;

            thisDiff = atom_getfloat(argv+i) - x->x_prevFeature[i];

            switch(x->x_direction)
            {
                case deltaPos:
                    thisDiff = (thisDiff<0)?0:thisDiff;
                    break;
                case deltaNeg:
                    thisDiff = (thisDiff>0)?0:thisDiff;
                    break;
                default:
                    break;
            }

            switch(x->x_mode)
            {
                case deltaAbs:
                    thisDiff = fabs(thisDiff);
                    break;
                case deltaSquared:
                    thisDiff = thisDiff*thisDiff;
                    break;
                default:
                    break;
            }

            SETFLOAT(x->x_listOut+i, thisDiff);
        }
    }

    outlet_list(x->x_featureList, 0, x->x_featureLength, x->x_listOut);
}


static void featureDelta_prevFeature(t_featureDelta *x, t_symbol *s, int argc, t_atom *argv)
{
    t_attributeIdx i;

    if(x->x_featureLength != (t_attributeIdx)argc)
    {
        pd_error(x, "%s: previous feature list length does not match current length setting. input ignored.", x->x_objSymbol->s_name);
        return;
    }
    else
    {
        for(i=0; i<x->x_featureLength; i++)
            x->x_prevFeature[i] = atom_getfloat(argv+i);
    }
}


static void featureDelta_print(t_featureDelta *x)
{
    post("%s feature length: %i", x->x_objSymbol->s_name, x->x_featureLength);
    switch(x->x_mode)
    {
        case deltaDiff:
            post("%s mode: difference", x->x_objSymbol->s_name);
            break;
        case deltaAbs:
            post("%s mode: absolute value of difference", x->x_objSymbol->s_name);
            break;
        case deltaSquared:
            post("%s mode: squared difference", x->x_objSymbol->s_name);
            break;
        default:
            post("%s mode: difference", x->x_objSymbol->s_name);
            break;
    }

    switch(x->x_direction)
    {
        case deltaPos:
            post("%s direction: positive delta only", x->x_objSymbol->s_name);
            break;
        case deltaNeg:
            post("%s direction: negative delta only", x->x_objSymbol->s_name);
            break;
        case deltaBoth:
            post("%s direction: all delta", x->x_objSymbol->s_name);
            break;
        default:
            post("%s direction: all delta", x->x_objSymbol->s_name);
            break;
    }
}

static void featureDelta_clear(t_featureDelta *x)
{
    featureDelta_free(x);
    featureDelta_allocMem(x);
    featureDelta_initMem(x);
}

static void featureDelta_length(t_featureDelta *x, t_floatarg len)
{
    len = (len<1)?1:len;

    // free memory first
    featureDelta_free(x);

    x->x_featureLength = len;

    featureDelta_allocMem(x);
    featureDelta_initMem(x);
}

static void featureDelta_mode(t_featureDelta *x, t_symbol *m)
{
    if(!strcmp(m->s_name, "diff"))
        x->x_mode = deltaDiff;
    else if(!strcmp(m->s_name, "abs"))
        x->x_mode = deltaAbs;
    else if(!strcmp(m->s_name, "squared"))
        x->x_mode = deltaSquared;
    else
        x->x_mode = deltaDiff;
}

static void featureDelta_direction(t_featureDelta *x, t_symbol *d)
{
    if(!strcmp(d->s_name, "pos"))
        x->x_direction = deltaPos;
    else if(!strcmp(d->s_name, "neg"))
        x->x_direction = deltaNeg;
    else if(!strcmp(d->s_name, "both"))
        x->x_direction = deltaBoth;
    else
        x->x_direction = deltaBoth;
}

static void *featureDelta_new(t_symbol *s, int argc, t_atom *argv)
{
    t_featureDelta *x = (t_featureDelta *)pd_new(featureDelta_class);
    t_float featureLength;
    t_symbol *mode;
    t_symbol *direction;

    x->x_featureList = outlet_new(&x->x_obj, gensym("list"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("prevFeature"));

    x->x_objSymbol = gensym("featureDelta");

    switch(argc)
    {
        case 0:
            x->x_featureLength = 50;
            x->x_mode = deltaDiff;
            x->x_direction = deltaBoth;
            break;
        case 1:
            featureLength = atom_getfloat(argv);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;
            x->x_mode = deltaDiff;
            x->x_direction = deltaBoth;
            break;
        case 2:
            featureLength = atom_getfloat(argv);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;

            mode = atom_getsymbol(argv+1);

            if(!strcmp(mode->s_name, "diff"))
                x->x_mode = deltaDiff;
            else if(!strcmp(mode->s_name, "abs"))
                x->x_mode = deltaAbs;
            else if(!strcmp(mode->s_name, "squared"))
                x->x_mode = deltaSquared;
            else
                x->x_mode = deltaDiff;

            x->x_direction = deltaBoth;
            break;
        case 3:
            featureLength = atom_getfloat(argv);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;

            mode = atom_getsymbol(argv+1);

            if(!strcmp(mode->s_name, "diff"))
                x->x_mode = deltaDiff;
            else if(!strcmp(mode->s_name, "abs"))
                x->x_mode = deltaAbs;
            else if(!strcmp(mode->s_name, "squared"))
                x->x_mode = deltaSquared;
            else
                x->x_mode = deltaDiff;

            direction = atom_getsymbol(argv+2);

            if(!strcmp(direction->s_name, "pos"))
                x->x_direction = deltaPos;
            else if(!strcmp(direction->s_name, "neg"))
                x->x_direction = deltaNeg;
            else if(!strcmp(direction->s_name, "both"))
                x->x_direction = deltaBoth;
            else
                x->x_direction = deltaBoth;

            break;
        default:
            featureLength = atom_getfloat(argv);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;

            mode = atom_getsymbol(argv+1);

            if(!strcmp(mode->s_name, "diff"))
                x->x_mode = deltaDiff;
            else if(!strcmp(mode->s_name, "abs"))
                x->x_mode = deltaAbs;
            else if(!strcmp(mode->s_name, "squared"))
                x->x_mode = deltaSquared;
            else
                x->x_mode = deltaDiff;

            direction = atom_getsymbol(argv+2);

            if(!strcmp(direction->s_name, "pos"))
                x->x_direction = deltaPos;
            else if(!strcmp(direction->s_name, "neg"))
                x->x_direction = deltaNeg;
            else if(!strcmp(direction->s_name, "both"))
                x->x_direction = deltaBoth;
            else
                x->x_direction = deltaBoth;

            break;
    }

     featureDelta_allocMem(x);
     featureDelta_initMem(x);

    return (void *)x;
}

void featureDelta_setup(void) {

    featureDelta_class = class_new(gensym("featureDelta"),
        (t_newmethod)featureDelta_new,
        (t_method)featureDelta_free,
        sizeof(t_featureDelta),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)featureDelta_new,
        gensym("timbreIDLib/featureDelta"),
        A_GIMME,
        0
    );

    class_addlist(featureDelta_class, (t_method)featureDelta_delta);

    class_addmethod(
        featureDelta_class,
        (t_method)featureDelta_prevFeature,
        gensym("prevFeature"),
        A_GIMME,
        0
    );


    class_addmethod(
        featureDelta_class,
        (t_method)featureDelta_delta,
        gensym("delta"),
        A_GIMME,
        0
    );

    class_addmethod(
        featureDelta_class,
        (t_method)featureDelta_print,
        gensym("print"),
        0
    );

    class_addmethod(
        featureDelta_class,
        (t_method)featureDelta_clear,
        gensym("clear"),
        0
    );

    class_addmethod(
        featureDelta_class,
        (t_method)featureDelta_length,
        gensym("length"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        featureDelta_class,
        (t_method)featureDelta_mode,
        gensym("mode"),
        A_DEFSYMBOL,
        0
    );

    class_addmethod(
        featureDelta_class,
        (t_method)featureDelta_direction,
        gensym("direction"),
        A_DEFSYMBOL,
        0
    );
}
