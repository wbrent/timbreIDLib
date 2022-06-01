/*

binWrangler - an external for organizing multiple frame feature vectors into component tracks.

Copyright 2010 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include "tIDLib.h"

static t_class* binWrangler_class;

typedef struct _binWrangler
{
    t_object  x_obj;
    t_symbol* x_objSymbol;
    t_instance* x_instances;
    t_bool x_spew;
    t_attributeIdx x_featureLength;
    t_attributeIdx x_numFrames;
    t_attributeIdx x_currentFrame;
    t_atom* x_listOut;
    t_outlet* x_binList;
} t_binWrangler;


static void binWrangler_accum (t_binWrangler* x, t_symbol* s, int argc, t_atom* argv)
{
    t_attributeIdx i, j, count, totalbins;

    if (x->x_featureLength != (t_attributeIdx)argc)
    {
        pd_error (x, "%s: input length does not match current length setting. input ignored.", x->x_objSymbol->s_name);
        return;
    }
    else
        for (i = 0; i < x->x_featureLength; i++)
            x->x_instances[x->x_currentFrame].data[i] = atom_getfloat (argv + i);

    x->x_currentFrame++;

    if ((x->x_currentFrame == x->x_numFrames) || x->x_spew)
    {
        totalbins = x->x_numFrames*x->x_featureLength;

        for (i = 0; i < x->x_featureLength; i++)
            for (count = 0, j = x->x_numFrames - x->x_currentFrame; count < x->x_numFrames; count++, j++)
                SETFLOAT (x->x_listOut+(i*x->x_numFrames)+(j%x->x_numFrames), x->x_instances[count].data[i]);

        outlet_list (x->x_binList, 0, totalbins, x->x_listOut);

        x->x_currentFrame = (x->x_currentFrame == x->x_numFrames)?0:x->x_currentFrame;
    }
}

static void binWrangler_print (t_binWrangler* x)
{
    post ("%s num_frames: %i", x->x_objSymbol->s_name, x->x_numFrames);
    post ("%s feature length: %i", x->x_objSymbol->s_name, x->x_featureLength);
    post ("%s spew mode: %i", x->x_objSymbol->s_name, x->x_spew);
}

static void binWrangler_clear (t_binWrangler* x)
{
    t_attributeIdx i, j;

    // free the database memory
    for (i = 0; i < x->x_numFrames; i++)
        t_freebytes (x->x_instances[i].data, x->x_featureLength * sizeof (float));

    t_freebytes (x->x_instances, x->x_numFrames * sizeof (t_instance));

    x->x_currentFrame = 0;

    for (i = 0; i < x->x_featureLength * x->x_numFrames; i++)
        SETFLOAT (x->x_listOut + i, 0.0);

    x->x_instances = (t_instance *)t_getbytes (x->x_numFrames * sizeof (t_instance));

    for (i = 0; i < x->x_numFrames; i++)
        x->x_instances[i].data = (float *)t_getbytes (x->x_featureLength * sizeof (float));

    for (i = 0; i < x->x_numFrames; i++)
        for (j = 0; j < x->x_featureLength; j++)
            x->x_instances[i].data[j] = 0.0;
}

static void binWrangler_numFrames (t_binWrangler* x, t_float num)
{
    t_attributeIdx i, j;

    num = (num < 1) ? 1 : num;

    x->x_listOut = (t_atom *)t_resizebytes (x->x_listOut, x->x_featureLength * x->x_numFrames * sizeof (t_atom), (x->x_featureLength*num) * sizeof (t_atom));

    // free the database memory
    for (i = 0; i < x->x_numFrames; i++)
        t_freebytes (x->x_instances[i].data, x->x_featureLength * sizeof (float));

    t_freebytes (x->x_instances, x->x_numFrames * sizeof (t_instance));

    x->x_currentFrame = 0;
    x->x_numFrames = num;

    for (i = 0; i < x->x_featureLength * x->x_numFrames; i++)
        SETFLOAT (x->x_listOut + i, 0.0);

    x->x_instances = (t_instance *)t_getbytes (x->x_numFrames * sizeof (t_instance));

    for (i = 0; i < x->x_numFrames; i++)
        x->x_instances[i].data = (float *)t_getbytes (x->x_featureLength * sizeof (float));

    for (i = 0; i < x->x_numFrames; i++)
        for (j = 0; j < x->x_featureLength; j++)
            x->x_instances[i].data[j] = 0.0;
}

static void binWrangler_length (t_binWrangler* x, t_float len)
{
    t_attributeIdx i, j;

    len = (len < 1) ? 1 : len;

    x->x_listOut = (t_atom *)t_resizebytes (x->x_listOut, x->x_featureLength * x->x_numFrames * sizeof (t_atom), (len*x->x_numFrames) * sizeof (t_atom));

    // free the database memory
    for (i = 0; i < x->x_numFrames; i++)
        t_freebytes (x->x_instances[i].data, x->x_featureLength * sizeof (float));

    t_freebytes (x->x_instances, x->x_numFrames * sizeof (t_instance));

    x->x_instances = (t_instance *)t_getbytes (x->x_numFrames * sizeof (t_instance));

    x->x_featureLength = len;
    x->x_currentFrame = 0;

    for (i = 0; i < x->x_featureLength * x->x_numFrames; i++)
        SETFLOAT (x->x_listOut + i, 0.0);

    for (i = 0; i < x->x_numFrames; i++)
        x->x_instances[i].data = (float *)t_getbytes (x->x_featureLength * sizeof (float));

    for (i = 0; i < x->x_numFrames; i++)
        for (j = 0; j < x->x_featureLength; j++)
            x->x_instances[i].data[j] = 0.0;
}

static void binWrangler_spew (t_binWrangler* x, t_floatarg s)
{
    s = (s <= 0) ? 0 : s;
    s = (s >= 1) ? 1 : s;
    x->x_spew = s;
}

static void* binWrangler_new (t_float numFrames, t_float length, t_float spew)
{
    t_binWrangler* x = (t_binWrangler *)pd_new (binWrangler_class);
    t_attributeIdx i, j;

    x->x_binList = outlet_new (&x->x_obj, gensym ("list"));

    x->x_objSymbol = gensym ("binWrangler");

    length = (length<1)?1:length;
    x->x_featureLength = length;
    numFrames = (numFrames < 1) ? 1 : numFrames;
    x->x_numFrames = numFrames;
    x->x_currentFrame = 0;
    spew = (spew > 1) ? 1 : spew;
    spew = (spew < 0) ? 0 : spew;
    x->x_spew = spew;


    x->x_listOut = (t_atom *)t_getbytes (x->x_featureLength * x->x_numFrames * sizeof (t_atom));
    x->x_instances = (t_instance *)t_getbytes (x->x_numFrames * sizeof (t_instance));

    for (i = 0; i < x->x_featureLength * x->x_numFrames; i++)
        SETFLOAT (x->x_listOut + i, 0.0);

    for (i = 0; i < x->x_numFrames; i++)
        x->x_instances[i].data = (float *)t_getbytes (x->x_featureLength * sizeof (float));

    for (i = 0; i < x->x_numFrames; i++)
        for (j = 0; j < x->x_featureLength; j++)
            x->x_instances[i].data[j] = 0.0;

    return (void *)x;
}

static void binWrangler_free (t_binWrangler* x)
{
    t_attributeIdx i;

    // free listOut memory
    t_freebytes (x->x_listOut, x->x_featureLength * x->x_numFrames * sizeof (t_atom));

    // free the database memory
    for (i = 0; i < x->x_currentFrame; i++)
        t_freebytes (x->x_instances[i].data, x->x_featureLength * sizeof (float));

    t_freebytes (x->x_instances, x->x_numFrames * sizeof (t_instance));

}

void binWrangler_setup (void) {

    binWrangler_class = class_new (gensym ("binWrangler"),
        (t_newmethod)binWrangler_new,
        (t_method)binWrangler_free,
        sizeof (t_binWrangler),
        CLASS_DEFAULT,
        A_DEFFLOAT,
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addcreator (
        (t_newmethod)binWrangler_new,
        gensym ("timbreIDLib/binWrangler"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addlist (binWrangler_class, (t_method)binWrangler_accum);

    class_addmethod (
        binWrangler_class,
        (t_method)binWrangler_clear,
        gensym ("clear"),
        0
    );

    class_addmethod (
        binWrangler_class,
        (t_method)binWrangler_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        binWrangler_class,
        (t_method)binWrangler_numFrames,
        gensym ("num_frames"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        binWrangler_class,
        (t_method)binWrangler_length,
        gensym ("length"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        binWrangler_class,
        (t_method)binWrangler_spew,
        gensym ("spew"),
        A_DEFFLOAT,
        0
    );
}
