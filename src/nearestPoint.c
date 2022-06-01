/*

nearestPoint - Compare an incoming coordinate with a database of coordinates and output the nearest point along with its distance from the input.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class* nearestPoint_class;

typedef struct _nearestPoint
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_instance* x_instances;
    t_attributeData* x_attributeData;
    t_attributeIdx x_dimensions;
    t_instanceIdx x_numInstances;
    t_instanceIdx x_numMatches;
    t_outlet* x_nearest;
    t_outlet* x_nearestDist;

} t_nearestPoint;


/* ------------------------ nearestPoint -------------------------------- */

static void nearestPoint_add (t_nearestPoint* x, t_symbol* s, int argc, t_atom* argv)
{
    t_attributeIdx i, dimensions;
    t_instanceIdx pointIdx;

    pointIdx = x->x_numInstances;
    dimensions = argc;

    if (dimensions == x->x_dimensions)
    {
        x->x_instances = (t_instance *)t_resizebytes (x->x_instances, x->x_numInstances * sizeof (t_instance), (x->x_numInstances + 1) * sizeof (t_instance));

        x->x_instances[pointIdx].data = (t_float *)t_getbytes (dimensions * sizeof (t_float));

        x->x_numInstances++;

        for (i = 0; i < dimensions; i++)
            x->x_instances[pointIdx].data[i] = atom_getfloat (argv + i);
    }
    else
        pd_error (x, "%s: dimensionality mismatch. input ignored", x->x_objSymbol->s_name);
}


static void nearestPoint_nearest (t_nearestPoint* x, t_symbol* s, int argc, t_atom* argv)
{
    t_float dist, *inputBuffer, *instanceBuffer, *weights;
    t_attributeIdx j, dimensions;
    t_instanceIdx i;

    if (x->x_numInstances)
    {
        dimensions = argc;

        if (dimensions == x->x_dimensions)
        {
            inputBuffer = (t_float *)t_getbytes (dimensions * sizeof (t_float));
            instanceBuffer = (t_float *)t_getbytes (dimensions * sizeof (t_float));
            weights = (t_float *)t_getbytes (dimensions * sizeof (t_float));

            for (i = 0; i < x->x_numInstances; i++)
            {
                for (j = 0; j < dimensions; j++)
                {
                    x->x_attributeData[j].inputData = atom_getfloat (argv + j);
                    inputBuffer[j] = x->x_attributeData[j].inputData;
                    instanceBuffer[j] = x->x_instances[i].data[j];
                    weights[j] = x->x_attributeData[j].weight;
                }

                dist = 0.0;
                dist = tIDLib_euclidDist (dimensions, inputBuffer, instanceBuffer, weights, true);

                x->x_instances[i].knnInfo.dist = x->x_instances[i].knnInfo.safeDist = dist; // store the distance
                x->x_instances[i].knnInfo.idx = i; // store the idx
            };

            tIDLib_sortKnnInfo (x->x_numMatches, x->x_numInstances, UINT_MAX, x->x_instances);

            for (i = 0; i < x->x_numMatches; i++)
            {
                outlet_float (x->x_nearestDist, x->x_instances[i].knnInfo.safeDist);
                outlet_float (x->x_nearest, x->x_instances[i].knnInfo.idx);
            }

            // free local memory
            t_freebytes (inputBuffer, dimensions * sizeof (t_float));
            t_freebytes (instanceBuffer, dimensions * sizeof (t_float));
            t_freebytes (weights, dimensions * sizeof (t_float));
        }
        else
            pd_error (x, "%s: dimensionality mismatch.", x->x_objSymbol->s_name);
    }
}


static void nearestPoint_dimensions (t_nearestPoint* x, t_floatarg dim)
{
    if (x->x_numInstances > 0)
    {
        pd_error (x, "%s: clear all coordinates before changing dimensionality.", x->x_objSymbol->s_name);
        return;
    }
    else if (dim <= 0)
        x->x_dimensions = 1;
    else
        x->x_dimensions = dim;

    post ("%s dimensionality: %i.", x->x_objSymbol->s_name, x->x_dimensions);
}


static void nearestPoint_num_matches (t_nearestPoint* x, t_floatarg n)
{
    n = (n > x->x_numInstances) ? x->x_numInstances : n;
    n = (n < 1) ? 1 : n;

    x->x_numMatches = n;

    post ("%s num_matches: %i.", x->x_objSymbol->s_name, x->x_numMatches);
}


static void nearestPoint_print (t_nearestPoint* x)
{
    post ("%s dimensions: %i", x->x_objSymbol->s_name, x->x_dimensions);
    post ("%s number of points: %i", x->x_objSymbol->s_name, x->x_numInstances);
    post ("%s num_matches: %i", x->x_objSymbol->s_name, x->x_numMatches);
}


static void nearestPoint_clear (t_nearestPoint* x)
{
    t_instanceIdx i;

    for (i = 0; i < x->x_numInstances; i++)
        t_freebytes (x->x_instances[i].data, x->x_dimensions * sizeof (t_float));

    x->x_instances = (t_instance *)t_resizebytes (x->x_instances, x->x_numInstances * sizeof (t_instance), 0);

    x->x_numInstances = 0;
}


static void* nearestPoint_new (t_float dim)
{
    t_nearestPoint* x = (t_nearestPoint *)pd_new (nearestPoint_class);
    t_attributeIdx i;

    x->x_nearest = outlet_new (&x->x_obj, &s_float);
    x->x_nearestDist = outlet_new (&x->x_obj, &s_float);
    inlet_new (&x->x_obj, &x->x_obj.ob_pd, gensym ("list"), gensym ("nearest"));

    x->x_objSymbol = gensym ("nearestPoint");

    if (dim)
        x->x_dimensions = dim;
    else
        x->x_dimensions = 2;

    x->x_numInstances = 0;
    x->x_numMatches = 1;

    x->x_instances = (t_instance *)t_getbytes (0);

    // resize feature input buffer to default dimensions
    x->x_attributeData = (t_attributeData *)t_getbytes (x->x_dimensions * sizeof (t_attributeData));

    // initialize feature input buffer
    for (i = 0; i < x->x_dimensions; i++)
    {
        x->x_attributeData[i].inputData = 0.0;
        x->x_attributeData[i].weight = 1.0;
    }

    post ("%s dimensionality: %i.", x->x_objSymbol->s_name, x->x_dimensions);

    return (x);
}


static void nearestPoint_free (t_nearestPoint* x)
{
    t_instanceIdx i;

    for (i = 0; i < x->x_numInstances; i++)
        t_freebytes (x->x_instances[i].data, x->x_dimensions * sizeof (t_float));

    t_freebytes (x->x_instances, x->x_numInstances * sizeof (t_instance));
    t_freebytes (x->x_attributeData, x->x_dimensions * sizeof (t_attributeData));
}


void nearestPoint_setup (void)
{
    nearestPoint_class =
    class_new (
        gensym ("nearestPoint"),
        (t_newmethod)nearestPoint_new,
        (t_method)nearestPoint_free,
        sizeof (t_nearestPoint),
        CLASS_DEFAULT,
        A_DEFFLOAT,
        0
    );

    class_addcreator (
        (t_newmethod)nearestPoint_new,
        gensym ("timbreIDLib/nearestPoint"),
        A_DEFFLOAT,
        0
    );

    class_addlist (
        nearestPoint_class,
        (t_method)nearestPoint_add
    );

    class_addmethod (
        nearestPoint_class,
        (t_method)nearestPoint_nearest,
        gensym ("nearest"),
        A_GIMME,
        0
    );

    class_addmethod (
        nearestPoint_class,
        (t_method)nearestPoint_dimensions,
        gensym ("dimensions"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        nearestPoint_class,
        (t_method)nearestPoint_num_matches,
        gensym ("num_matches"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        nearestPoint_class,
        (t_method)nearestPoint_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        nearestPoint_class,
        (t_method)nearestPoint_clear,
        gensym ("clear"),
        0
    );
}
