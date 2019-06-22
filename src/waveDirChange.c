/*

waveDirChange

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *waveDirChange_class;

typedef struct _waveDirChange
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
	t_sampIdx x_window;
	t_float *x_analysisBuffer;
	t_word *x_vec;
	t_symbol *x_arrayName;
	t_sampIdx x_arrayPoints;
    t_outlet *x_changes;
} t_waveDirChange;


/* ------------------------ waveDirChange -------------------------------- */

static void waveDirChange_analyze(t_waveDirChange *x, t_floatarg start, t_floatarg n)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j, oldWindow, window, startSamp, endSamp;
		t_float changes, *flagsBuf, minVal, maxVal;

		startSamp = (start<0)?0:start;

		if(n)
			endSamp = startSamp + n-1;
		else
			endSamp = startSamp + x->x_window-1;

		if(endSamp >= x->x_arrayPoints)
			endSamp = x->x_arrayPoints-1;

		window = endSamp-startSamp+1;

		if(endSamp <= startSamp)
		{
			post("%s: bad range of samples.", x->x_objSymbol->s_name);
			return;
		}

		if(x->x_window != window)
		{			
			oldWindow = x->x_window;
	
			// window must be at least 4 points long
			if(window<MINWINDOWSIZE)
			{
				window = WINDOWSIZEDEFAULT;
				post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
			}

			// hang on to these values for next time
			x->x_window = window;

			endSamp = startSamp + x->x_window-1;
			if(endSamp > x->x_arrayPoints)
				endSamp = x->x_arrayPoints-1;

			x->x_analysisBuffer = (t_float *)t_resizebytes(x->x_analysisBuffer, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

		}

		// construct analysis window
		for(i=0, j=startSamp; j<=endSamp; i++, j++)
			x->x_analysisBuffer[i] = x->x_vec[j].w_float;

		changes = minVal = maxVal = 0.0;
		flagsBuf = (t_float *)t_getbytes(window*sizeof(t_float));

		for(i=0; i<window; i++)
			flagsBuf[i] = 0.0;
			
		tIDLib_peaksValleys(window, x->x_analysisBuffer, flagsBuf, &minVal, &maxVal);

		for(i=0; i<window; i++)
			changes += fabsf(flagsBuf[i]);
		
		// free local memory
		t_freebytes(flagsBuf, window*sizeof(t_float));
		
		outlet_float(x->x_changes, changes);
	}
}


// analyze the whole damn array
static void waveDirChange_bang(t_waveDirChange *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx window, startSamp;
		startSamp = 0;
		window = x->x_arrayPoints;
		waveDirChange_analyze(x, startSamp, window);
	}
}


static void waveDirChange_set(t_waveDirChange *x, t_symbol *s)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
		pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
	else
	    x->x_arrayName = s;
}


static void waveDirChange_print(t_waveDirChange *x)
{
	post("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
}


static void waveDirChange_samplerate(t_waveDirChange *x, t_floatarg sr)
{
	if(sr<MINSAMPLERATE)
		x->x_sr = MINSAMPLERATE;
	else
		x->x_sr = sr;
}


static void *waveDirChange_new(t_symbol *s, int argc, t_atom *argv)
{
    t_waveDirChange *x = (t_waveDirChange *)pd_new(waveDirChange_class);
//	t_garray *a;

	x->x_changes = outlet_new(&x->x_obj, &s_float);

	// store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
	x->x_objSymbol = s;
	
	switch(argc)
	{
		case 1:
			x->x_arrayName = atom_getsymbol(argv);
			/*
			if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
				pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
			else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
				pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
			*/
			break;
	
		case 0:
			post("%s: no array specified.", x->x_objSymbol->s_name);
			// a bogus array name to trigger the safety check in _analyze()
			x->x_arrayName = gensym("NOARRAYSPECIFIED");
			break;
		
		default:
			x->x_arrayName = atom_getsymbol(argv);
			/*
			if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
				pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
			else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
				pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
			*/
			post("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
			break;
	}
	
	x->x_sr = SAMPLERATEDEFAULT;
	x->x_window = WINDOWSIZEDEFAULT;

	x->x_analysisBuffer = (t_sample *)t_getbytes(x->x_window*sizeof(t_sample));

    return (x);
}


static void waveDirChange_free(t_waveDirChange *x)
{
	// free the input buffer memory
    t_freebytes(x->x_analysisBuffer, x->x_window*sizeof(t_sample));
}


void waveDirChange_setup(void)
{
    waveDirChange_class =
    class_new(
    	gensym("waveDirChange"),
    	(t_newmethod)waveDirChange_new,
    	(t_method)waveDirChange_free,
        sizeof(t_waveDirChange),
        CLASS_DEFAULT,
        A_GIMME,
		0
    );

	class_addcreator(
		(t_newmethod)waveDirChange_new,
		gensym("timbreIDLib/waveDirChange"),
		A_GIMME,
		0
	);
	
	class_addbang(waveDirChange_class, waveDirChange_bang);

	class_addmethod(
		waveDirChange_class,
        (t_method)waveDirChange_analyze,
		gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		waveDirChange_class,
		(t_method)waveDirChange_set,
		gensym("set"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		waveDirChange_class,
		(t_method)waveDirChange_print,
		gensym("print"),
		0
	);

	class_addmethod(
		waveDirChange_class,
        (t_method)waveDirChange_samplerate,
		gensym("samplerate"),
		A_DEFFLOAT,
		0
	);
}
