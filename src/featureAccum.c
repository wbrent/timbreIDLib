/*

featureAccum - an external for accumulating multiple feature frames over time.

Copyright 2010 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *featureAccum_class;

typedef enum
{
    concat = 0,
    sum,
    mean,
    sma
} t_accumMode;

typedef struct _featureAccum
{
    t_object  x_obj;
    t_symbol *x_objSymbol;
    t_instance *x_instances;
    t_bool x_spew;
    t_accumMode x_mode;
    t_attributeIdx x_featureLength;
    t_attributeIdx x_numFrames;
    t_attributeIdx x_concatCurrentFrame;
    t_attributeIdx x_meanFrameCount;
    t_atom *x_listOut;
    t_outlet *x_featureList;
} t_featureAccum;


static void featureAccum_allocMem(t_featureAccum *x)
{
    t_attributeIdx i;

    // grab atom list memory
    x->x_listOut = (t_atom *)t_getbytes((x->x_featureLength*x->x_numFrames)*sizeof(t_atom));

    // grab new database memory
    // need numFrames+1 for an extra row for sma mode to sum to
    x->x_instances = (t_instance *)t_getbytes((x->x_numFrames+1)*sizeof(t_instance));

    // grab featureLength floats per row
    for(i=0; i<x->x_numFrames+1; i++)
        x->x_instances[i].data = (float *)t_getbytes(x->x_featureLength*sizeof(float));
}

static void featureAccum_initMem(t_featureAccum *x)
{
    t_attributeIdx i, j;

    // clear the atom list
    for(i=0; i<x->x_featureLength*x->x_numFrames; i++)
        SETFLOAT(x->x_listOut+i, 0.0);

    // initialize each row of database to zero
    for(i=0; i<x->x_numFrames+1; i++)
        for(j=0; j<x->x_featureLength; j++)
            x->x_instances[i].data[j] = 0.0;
}

static void featureAccum_free(t_featureAccum *x)
{
    t_attributeIdx i;

    // free listOut memory
    t_freebytes(x->x_listOut, (x->x_featureLength*x->x_numFrames)*sizeof(t_atom));

    // free the database memory rows
    for(i=0; i<x->x_numFrames+1; i++)
        t_freebytes(x->x_instances[i].data, x->x_featureLength*sizeof(float));

    // free the database memory itself
    t_freebytes(x->x_instances, (x->x_numFrames+1)*sizeof(t_instance));

}

static void featureAccum_sum(t_featureAccum *x)
{
    t_attributeIdx i;

    for(i=0; i<x->x_featureLength; i++)
        SETFLOAT(x->x_listOut+i, x->x_instances[0].data[i]);
        // only deal with first row
    outlet_list(x->x_featureList, 0, x->x_featureLength, x->x_listOut);
}

static void featureAccum_mean(t_featureAccum *x)
{
    t_attributeIdx i;

    for(i=0; i<x->x_featureLength; i++)
        SETFLOAT(x->x_listOut+i, x->x_instances[0].data[i]/x->x_meanFrameCount);
        // only deal with first row
    outlet_list(x->x_featureList, 0, x->x_featureLength, x->x_listOut);
}

static void featureAccum_accum(t_featureAccum *x, t_symbol *s, int argc, t_atom *argv)
{
    t_attributeIdx i, j, count, totalFeat;

    if(x->x_featureLength != (t_attributeIdx)argc)
    {
        pd_error(x, "%s: input length does not match current length setting. input ignored.", x->x_objSymbol->s_name);
        return;
    }
    else
    {
        switch(x->x_mode)
        {
            case concat:
            case sma:
                for(i=0; i<x->x_featureLength; i++)
                    x->x_instances[x->x_concatCurrentFrame].data[i] = atom_getfloat(argv+i);
                break;
            default: // for modes 1 and 2, add to previous contents. first row only
                for(i=0; i<x->x_featureLength; i++)
                    x->x_instances[0].data[i] += atom_getfloat(argv+i);
                break;
        }
    }

    // advance the concat frame count
    x->x_concatCurrentFrame++;
    // advance the frame counter for mean
    x->x_meanFrameCount++;

    switch(x->x_mode)
    {
        case concat:
            if((x->x_concatCurrentFrame==x->x_numFrames) || x->x_spew)
            {
                totalFeat = x->x_featureLength * x->x_numFrames;

                for(count=0, i=x->x_numFrames-x->x_concatCurrentFrame; count<x->x_numFrames; count++, i++)
                    for(j=0; j<x->x_featureLength; j++)
                        SETFLOAT(x->x_listOut+((i%x->x_numFrames)*x->x_featureLength)+j, x->x_instances[count].data[j]);

                outlet_list(x->x_featureList, 0, totalFeat, x->x_listOut);
                // reset currentFrame if we've hit numFrames only
                // do nothing if we're in spew mode and haven't hit numFrames
                x->x_concatCurrentFrame = (x->x_concatCurrentFrame>=x->x_numFrames)?0:x->x_concatCurrentFrame;
            }
            break;

        case sum:
            if(x->x_spew)
                featureAccum_sum(x);
            break;

        case mean:
            if(x->x_spew)
                featureAccum_mean(x);
            break;

        case sma:
            // clear the extra row of x_instances
            for(i=0; i<x->x_featureLength; i++)
                x->x_instances[x->x_numFrames].data[i] = 0.0;

            // sum all rows into the extra row
            for(i=0; i<x->x_numFrames; i++)
                for(j=0; j<x->x_featureLength; j++)
                    x->x_instances[x->x_numFrames].data[j] += x->x_instances[i].data[j];

            // divide the sum by x->x_numFrames
            for(i=0; i<x->x_featureLength; i++)
                SETFLOAT(x->x_listOut+i, x->x_instances[x->x_numFrames].data[i]/x->x_numFrames);

            outlet_list(x->x_featureList, 0, x->x_featureLength, x->x_listOut);

            // still have to reset concatCurrentFrame so that database row writing position wraps back to beginning
            x->x_concatCurrentFrame = (x->x_concatCurrentFrame>=x->x_numFrames)?0:x->x_concatCurrentFrame;
            break;

        default:
            break;
    }
}

static void featureAccum_print(t_featureAccum *x)
{
    post("%s num_frames: %i", x->x_objSymbol->s_name, x->x_numFrames);
    post("%s feature length: %i", x->x_objSymbol->s_name, x->x_featureLength);
    post("%s spew: %i", x->x_objSymbol->s_name, x->x_spew);
    switch(x->x_mode)
    {
        case concat:
            post("%s mode: concat", x->x_objSymbol->s_name);
            break;
        case sum:
            post("%s mode: sum", x->x_objSymbol->s_name);
            break;
        case mean:
            post("%s mode: mean", x->x_objSymbol->s_name);
            break;
        case sma:
            post("%s mode: simple moving average", x->x_objSymbol->s_name);
            break;
        default:
            post("%s mode: concat", x->x_objSymbol->s_name);
            break;
    }
}

static void featureAccum_clear(t_featureAccum *x)
{
    featureAccum_free(x);
    featureAccum_allocMem(x);
    featureAccum_initMem(x);

    x->x_concatCurrentFrame = 0;
    x->x_meanFrameCount = 0;
}

static void featureAccum_numFrames(t_featureAccum *x, t_floatarg num)
{
    num = (num<1)?1:num;

    // free memory first
    featureAccum_free(x);

    // then update numFrames
    x->x_numFrames = num;
    x->x_concatCurrentFrame = 0;

    featureAccum_allocMem(x);
    featureAccum_initMem(x);
}

static void featureAccum_length(t_featureAccum *x, t_floatarg len)
{
    len = (len<1)?1:len;

    // free memory first
    featureAccum_free(x);

    // then update numFrames
    x->x_featureLength = len;
    x->x_concatCurrentFrame = 0;

    featureAccum_allocMem(x);
    featureAccum_initMem(x);
}

static void featureAccum_spew(t_featureAccum *x, t_floatarg s)
{
    s = (s<=0)?0:s;
    s = (s>=1)?1:s;
    x->x_spew = s;
}

static void featureAccum_mode(t_featureAccum *x, t_symbol *m)
{
    if(!strcmp(m->s_name, "concat"))
        x->x_mode = concat;
    else if(!strcmp(m->s_name, "sum"))
        x->x_mode = sum;
    else if(!strcmp(m->s_name, "mean"))
        x->x_mode = mean;
    else if(!strcmp(m->s_name, "sma"))
        x->x_mode = sma;
    else
        x->x_mode = concat;

    featureAccum_clear(x);
}


static void *featureAccum_new(t_symbol *s, int argc, t_atom *argv)
{
    t_featureAccum *x = (t_featureAccum *)pd_new(featureAccum_class);
    t_symbol *mode;
    t_float numFrames, featureLength, spew;

    x->x_featureList = outlet_new(&x->x_obj, gensym("list"));

    x->x_objSymbol = gensym("featureAccum");

    switch(argc)
    {
        case 0:
            x->x_numFrames = 5;
            x->x_featureLength = 50;
            x->x_spew = false;
            x->x_mode = concat;
            break;
        case 1:
            numFrames = atom_getfloat(argv);
            numFrames = (numFrames<1)?1:numFrames;
            x->x_numFrames = numFrames;
            x->x_featureLength = 50;
            x->x_spew = 0;
            x->x_mode = concat;
            break;
        case 2:
            numFrames = atom_getfloat(argv);
            numFrames = (numFrames<1)?1:numFrames;
            x->x_numFrames = numFrames;
            featureLength = atom_getfloat(argv+1);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;
            x->x_spew = false;
            x->x_mode = concat;
            break;
        case 3:
            numFrames = atom_getfloat(argv);
            numFrames = (numFrames<1)?1:numFrames;
            x->x_numFrames = numFrames;
            featureLength = atom_getfloat(argv+1);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;
            spew = atom_getfloat(argv+2);
            spew = (spew>1)?1:spew;
            spew = (spew<0)?0:spew;
            x->x_spew = spew;
            x->x_mode = concat;
            break;
        case 4:
            numFrames = atom_getfloat(argv);
            numFrames = (numFrames<1)?1:numFrames;
            x->x_numFrames = numFrames;
            featureLength = atom_getfloat(argv+1);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;
            spew = atom_getfloat(argv+2);
            spew = (spew>1)?1:spew;
            spew = (spew<0)?0:spew;
            x->x_spew = spew;

            mode = atom_getsymbol(argv+3);

            if(!strcmp(mode->s_name, "concat"))
                x->x_mode = concat;
            else if(!strcmp(mode->s_name, "sum"))
                x->x_mode = sum;
            else if(!strcmp(mode->s_name, "mean"))
                x->x_mode = mean;
            else if(!strcmp(mode->s_name, "sma"))
                x->x_mode = sma;
            else
                x->x_mode = concat;

            break;
        default:
            numFrames = atom_getfloat(argv);
            numFrames = (numFrames<1)?1:numFrames;
            x->x_numFrames = numFrames;
            featureLength = atom_getfloat(argv+1);
            featureLength = (featureLength<1)?1:featureLength;
            x->x_featureLength = featureLength;
            spew = atom_getfloat(argv+2);
            spew = (spew>1)?1:spew;
            spew = (spew<0)?0:spew;
            x->x_spew = spew;

            mode = atom_getsymbol(argv+3);

            if(!strcmp(mode->s_name, "concat"))
                x->x_mode = concat;
            else if(!strcmp(mode->s_name, "sum"))
                x->x_mode = sum;
            else if(!strcmp(mode->s_name, "mean"))
                x->x_mode = mean;
            else if(!strcmp(mode->s_name, "sma"))
                x->x_mode = sma;
            else
                x->x_mode = concat;

            break;
    }

    x->x_concatCurrentFrame = 0;
    x->x_meanFrameCount = 0;

     featureAccum_allocMem(x);
     featureAccum_initMem(x);

    return (void *)x;
}


void featureAccum_setup(void) {

    featureAccum_class = class_new(gensym("featureAccum"),
        (t_newmethod)featureAccum_new,
        (t_method)featureAccum_free,
        sizeof(t_featureAccum),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator(
        (t_newmethod)featureAccum_new,
        gensym("timbreIDLib/featureAccum"),
        A_GIMME,
        0
    );

    class_addlist(featureAccum_class, (t_method)featureAccum_accum);

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_accum,
        gensym("accum"),
        A_GIMME,
        0
    );

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_print,
        gensym("print"),
        0
    );

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_clear,
        gensym("clear"),
        0
    );

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_numFrames,
        gensym("num_frames"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_length,
        gensym("length"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_sum,
        gensym("sum"),
        0
    );

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_mean,
        gensym("mean"),
        0
    );

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_spew,
        gensym("spew"),
        A_DEFFLOAT,
        0
    );

    class_addmethod(
        featureAccum_class,
        (t_method)featureAccum_mode,
        gensym("mode"),
        A_DEFSYMBOL,
        0
    );
}
