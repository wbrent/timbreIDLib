/*

tempo~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

/*

TODO:

should address the fact that analysis hop isn't restricted to be a multiple of x->x_n, and yet the analysis routine will always be called on block boundaries because it occurs when x_dspTick*x_n is greater or equal to x_hop

Preprocesssing of onset peaks: find each peak and retain its value. then add a user-defined amount of padding on either side of it that ramps down linearly to the belowThreshDefault value. this is an entirely experimental idea that could wait for later

*/


#include "tIDLib.h"

static t_class* tempo_tilde_class;

typedef struct _tempo_tilde
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_bool x_debug;
    t_float x_sr;
    t_float x_n;
    t_uInt x_hop;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_bool x_powerSpectrum;
    t_bool x_squaredDiff;

    t_float x_growthThresh;
    t_float x_belowThreshDefault;
    t_float x_loFreq;
    t_float x_hiFreq;
    t_binIdx x_loBin;
    t_binIdx x_hiBin;
    t_float x_onsetsBufDur;
    t_uInt x_onsetsBufSize;
    t_uShortInt x_numHarm;
    t_float x_loTempo;
    t_float x_hiTempo;
    t_uInt x_loTempoIdx;
    t_uInt x_hiTempoIdx;
    t_float x_lastGoodTempo;
    t_float x_maxOnsetPeakVal;
    t_bool x_maxOnsetUpdate;
    t_float x_lastFirstOnsetPeakVal;
    double x_lastDspTime;
    t_uInt x_dspTicks;
    t_uShortInt x_tempoBufferIdx;
    t_uShortInt x_tempoBufferSize;
    t_float* x_tempoBuffer;
    t_sample* x_signalBuffer;
    t_sample* x_fftwInForwardWindow;
    t_sample* x_fftwInBackWindow;
    t_float* x_onsetsBuffer;
    fftwf_complex* x_fftwOutForwardWindow;
    fftwf_complex* x_fftwOutBackWindow;
    fftwf_plan x_fftwPlanForwardWindow;
    fftwf_plan x_fftwPlanBackWindow;
    t_float* x_blackman;
    t_float* x_cosine;
    t_float* x_hamming;
    t_float* x_hann;
    t_atom* x_listOut;
    t_outlet* x_onsetList;
    t_outlet* x_tempoRaw;
    t_outlet* x_conf;
    t_outlet* x_tempo;
    t_float x_f;

} t_tempo_tilde;


static void tempo_tilde_resizeOnsetsBuffer (t_tempo_tilde* x, t_floatarg newSize)
{
    t_uInt minSize, i;

    minSize = x->x_loTempoIdx*(x->x_numHarm+1);

    if (newSize<minSize)
    {
        newSize = minSize;
        post ("%s: requested onset buffer duration is too short based on tempo and harmonic settings. resizing to minimum required duration.", x->x_objSymbol->s_name);
    }

    x->x_onsetsBuffer = (t_float *)t_resizebytes (x->x_onsetsBuffer, x->x_onsetsBufSize * sizeof (t_float), newSize * sizeof (t_float));
    x->x_listOut = (t_atom *)t_resizebytes (x->x_listOut, x->x_onsetsBufSize * sizeof (t_atom), newSize * sizeof (t_atom));

    x->x_onsetsBufSize = newSize;
    x->x_onsetsBufDur = ((x->x_onsetsBufSize*x->x_hop)/x->x_sr)*1000.0;

    for (i = 0; i < x->x_onsetsBufSize; i++)
        x->x_onsetsBuffer[i] = x->x_belowThreshDefault;

    post ("%s onset buffer size updated: %0.2f ms, %i frames.", x->x_objSymbol->s_name, x->x_onsetsBufDur, x->x_onsetsBufSize);
}


/* ------------------------ tempo~ -------------------------------- */

static void tempo_tilde_analyze (t_tempo_tilde* x)
{
    t_sampIdx i, j, startIdx, window, windowHalf;
    t_sInt periodIdx;
    t_uInt onsetBufferShiftedSize;
    t_uLongInt modeCount;
    t_float growth, tempo, tempoMode, tempoConf, peakSlope, *windowFuncPtr, *onsetsBufferShifted, *yValues, maxYValue;
    t_bool hpsFlag;


// SPECTRUM CAPTURE

    window = x->x_window;
    windowHalf = x->x_windowHalf;

    // construct forward analysis window
    for (i = 0, j = 0; i < window; i++, j++)
        x->x_fftwInForwardWindow[i] = x->x_signalBuffer[window + j];

    // construct back analysis window x->x_hop samples earlier
    for (i = 0, j = 0; i < window; i++, j++)
        x->x_fftwInBackWindow[i] = x->x_signalBuffer[window - x->x_hop + j];

    windowFuncPtr = x->x_blackman;

    switch (x->x_windowFunction)
    {
        case rectangular:
            break;
        case blackman:
            windowFuncPtr = x->x_blackman;
            break;
        case cosine:
            windowFuncPtr = x->x_cosine;
            break;
        case hamming:
            windowFuncPtr = x->x_hamming;
            break;
        case hann:
            windowFuncPtr = x->x_hann;
            break;
        default:
            windowFuncPtr = x->x_blackman;
            break;
    };

    // if x_windowFunction == 0, skip the windowing (rectangular)
    if (x->x_windowFunction != rectangular)
        for (i = 0; i < window; i++, windowFuncPtr++)
            x->x_fftwInForwardWindow[i] *= *windowFuncPtr;

    fftwf_execute (x->x_fftwPlanForwardWindow);

    // put the result of power calc back in x_fftwIn
    tIDLib_power (windowHalf + 1, x->x_fftwOutForwardWindow, x->x_fftwInForwardWindow);

    if ( !x->x_powerSpectrum)
        tIDLib_mag (windowHalf + 1, x->x_fftwInForwardWindow);

    switch (x->x_windowFunction)
    {
        case rectangular:
            break;
        case blackman:
            windowFuncPtr = x->x_blackman;
            break;
        case cosine:
            windowFuncPtr = x->x_cosine;
            break;
        case hamming:
            windowFuncPtr = x->x_hamming;
            break;
        case hann:
            windowFuncPtr = x->x_hann;
            break;
        default:
            windowFuncPtr = x->x_blackman;
            break;
    };

    // if x_windowFunction == 0, skip the windowing (rectangular)
    if (x->x_windowFunction != rectangular)
        for (i = 0; i < window; i++, windowFuncPtr++)
            x->x_fftwInBackWindow[i] *= *windowFuncPtr;

    fftwf_execute (x->x_fftwPlanBackWindow);

    // put the result of power calc back in x_fftwIn
    tIDLib_power (windowHalf + 1, x->x_fftwOutBackWindow, x->x_fftwInBackWindow);

    if ( !x->x_powerSpectrum)
        tIDLib_mag (windowHalf + 1, x->x_fftwInBackWindow);

// END SPECTRUM CAPTURE


    hpsFlag = false;


// MEASURE GROWTH

    growth=0;

    for (i = x->x_loBin; i <= x->x_hiBin; i++)
    {
        t_float diff;

        diff = x->x_fftwInForwardWindow[i] - x->x_fftwInBackWindow[i];

        // only look at growth (positive flux)
        if (diff>0.0)
        {
            if (x->x_squaredDiff)
                diff = diff * diff;
            else
                diff = fabsf (diff);
        }
        else
            diff = 0.0;

        growth += diff;
    }

    // set the maxPeakVal to the most recent growth measurement if it's greater
    // but only do this if the x_maxOnsetUpdate flag is on
    if (x->x_maxOnsetUpdate && growth>x->x_maxOnsetPeakVal)
    {
        x->x_maxOnsetPeakVal = growth;
        if (x->x_debug)
            post ("%s: updated max onset peak: %0.2f", x->x_objSymbol->s_name, x->x_maxOnsetPeakVal);
    }

    // shift the contents backwards
    for (i = 0; i < x->x_onsetsBufSize - 1; i++)
        x->x_onsetsBuffer[i] = x->x_onsetsBuffer[i + 1];

    if (growth>=((x->x_growthThresh/100.0)*x->x_maxOnsetPeakVal))
        x->x_onsetsBuffer[x->x_onsetsBufSize - 1] = growth;
    else
        x->x_onsetsBuffer[x->x_onsetsBufSize - 1] = x->x_belowThreshDefault; //set below thresh growth to a default value that is NOT ZERO, which would ruin the running product in the HPS algo

    startIdx=0;

    // fast forward to the beginning of the first peak
    while (x->x_onsetsBuffer[startIdx]==x->x_belowThreshDefault)
    {
        startIdx++;

        if (startIdx>=x->x_onsetsBufSize - 1)
            break;
    }

    // now, move to the peak of that first peak
    if (startIdx>0)
        peakSlope = x->x_onsetsBuffer[startIdx] - x->x_onsetsBuffer[startIdx - 1];
    else
        peakSlope = 0.0;

    while (peakSlope>0.0)
    {
        startIdx++;

        if (startIdx>=x->x_onsetsBufSize - 1)
            break;

        peakSlope = x->x_onsetsBuffer[startIdx] - x->x_onsetsBuffer[startIdx - 1];
    }

    // at this point, we're one idx past the peak, so back up a step
    // BUT: shouldn't back up here in the special case of the peak was already perfectly aligned. if that's the case, peakSlope would be exactly zero from startIdx>0 check above
    if (peakSlope!=0.0)
        startIdx--;

    // TODO: should we make this onsetsBufSize - 10 or something, just so that we don't get a buffer size of 0, which could cause problems with the onsetsBufferShifted and yValues arrays below?
    // safety to make sure we stay in bounds
    startIdx = (startIdx>x->x_onsetsBufSize - 1)?x->x_onsetsBufSize - 1:startIdx;

    if (x->x_onsetsBuffer[startIdx] != x->x_lastFirstOnsetPeakVal)
    {
        hpsFlag = true;
        x->x_lastFirstOnsetPeakVal = x->x_onsetsBuffer[startIdx];
    }

    onsetBufferShiftedSize = x->x_onsetsBufSize-startIdx;

    onsetsBufferShifted = (t_float *)t_getbytes (onsetBufferShiftedSize * sizeof (t_float));

    for (i=startIdx, j = 0; i < x->x_onsetsBufSize; i++, j++)
    {
        SETFLOAT (x->x_listOut+j, x->x_onsetsBuffer[i]);
        onsetsBufferShifted[j] = x->x_onsetsBuffer[i];
    }

    for (; j < x->x_onsetsBufSize; j++)
        SETFLOAT (x->x_listOut+j, x->x_belowThreshDefault);

// END MEASURE GROWTH




// HPS STUFF

    if (hpsFlag)
    {
        maxYValue = -1;

        yValues = (t_float *)t_getbytes (onsetBufferShiftedSize * sizeof (t_float));

        // init yValues arrays to zero
        for (i = 0; i < onsetBufferShiftedSize; i++)
            yValues[i] = 0.0;

        periodIdx = tIDLib_hps (onsetsBufferShifted, onsetBufferShiftedSize, x->x_hiTempoIdx, x->x_loTempoIdx, x->x_numHarm, yValues, &maxYValue, x->x_debug);

        // OUTLET 4: shifted onset peak buffer
        outlet_list (x->x_onsetList, 0, x->x_onsetsBufSize, x->x_listOut);

        if (periodIdx == -1)
        {
            outlet_float (x->x_tempoRaw, x->x_lastGoodTempo);

            if (x->x_debug)
                post ("%s: bad HPS result", x->x_objSymbol->s_name);
        }
        else if (maxYValue<0.0)
        {
            outlet_float (x->x_tempoRaw, x->x_lastGoodTempo);

            if (x->x_debug)
                post ("%s: no Y peak. max Y value: %f", x->x_objSymbol->s_name, maxYValue);
        }
        else
        {
            tempo = periodIdx*x->x_hop; // frames to samples
            tempo /= x->x_sr; // samples to seconds
            tempo = 60.0f/tempo; // seconds to BPM

            // OUTLET 3: raw tempo per HPS call
            outlet_float (x->x_tempoRaw, tempo);

            x->x_tempoBuffer[x->x_tempoBufferIdx] = roundf (tempo);

            // modeCount must be t_uLongInt for tIDLib_mode()
            modeCount = 0;
            tempoMode = tIDLib_mode (x->x_tempoBuffer, x->x_tempoBufferSize, &modeCount);

            tempoConf = (modeCount/(t_float)x->x_tempoBufferSize)*100.0f;

            // OUTLET 2: tempo confidence
            outlet_float (x->x_conf, tempoConf);

            // OUTLET 1: main tempo output
            outlet_float (x->x_tempo, tempoMode);

            x->x_tempoBufferIdx = (x->x_tempoBufferIdx+1) % x->x_tempoBufferSize;
            x->x_lastGoodTempo = tempo;
        }

        // free local memory
        t_freebytes (yValues, onsetBufferShiftedSize * sizeof (t_float));
    }

    t_freebytes (onsetsBufferShifted, onsetBufferShiftedSize * sizeof (t_float));
}


static void tempo_tilde_debug (t_tempo_tilde* x, t_floatarg d)
{
    d = (d<0)?0:d;
    d = (d>1)?1:d;
    x->x_debug = d;

    post ("%s debug mode: %i", x->x_objSymbol->s_name, x->x_debug);
}


static void tempo_tilde_print (t_tempo_tilde* x)
{
    post ("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr / x->x_overlap));
    post ("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
    post ("%s window: %i", x->x_objSymbol->s_name, x->x_window);
    post ("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
    post ("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
    post ("%s hop: %i", x->x_objSymbol->s_name, x->x_hop);
    post ("%s squared difference: %i", x->x_objSymbol->s_name, x->x_squaredDiff);
    post ("%s onset buffer size: %0.2f ms, %i frames", x->x_objSymbol->s_name, x->x_onsetsBufDur, x->x_onsetsBufSize);
    post ("%s tempo range: %0.2f through %0.2f BPM, frame index %i through %i.", x->x_objSymbol->s_name, x->x_loTempo, x->x_hiTempo, x->x_loTempoIdx, x->x_hiTempoIdx);
    post ("%s harmonics: %i.", x->x_objSymbol->s_name, x->x_numHarm);
    post ("%s growth thresh: %0.2f%% of max peak value = %f.", x->x_objSymbol->s_name, x->x_growthThresh, (x->x_growthThresh/100.0)*x->x_maxOnsetPeakVal);
    post ("%s frequency range: %0.2f through %0.2f Hz, bin %i through %i.", x->x_loFreq, x->x_hiFreq, x->x_objSymbol->s_name, x->x_loBin, x->x_hiBin);
    post ("%s below threshold default value: %f.", x->x_objSymbol->s_name, x->x_belowThreshDefault);
    post ("%s max onset update: %i", x->x_objSymbol->s_name, x->x_maxOnsetUpdate);
    post ("%s current max onset peak value: %f", x->x_objSymbol->s_name, x->x_maxOnsetPeakVal);
    post ("%s tempo buffer size %i.", x->x_objSymbol->s_name, x->x_tempoBufferSize);
    post ("%s debug: %i.", x->x_objSymbol->s_name, x->x_debug);
}


static void tempo_tilde_window (t_tempo_tilde* x, t_floatarg w)
{
    t_sampIdx i, window, windowHalf;

    window = w;

    if (window < TID_MINWINDOWSIZE)
    {
        window = TID_WINDOWSIZEDEFAULT;
        post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
    }

    if (x->x_hop > window)
    {
        post ("%s change in window size caused frame hop to be less than %i samples apart. Setting frame hop to half of current window size instead.", x->x_objSymbol->s_name, window);
        x->x_hop = window * 0.5;
    }

    windowHalf = window * 0.5;

    x->x_signalBuffer = (t_sample *)t_resizebytes (x->x_signalBuffer, (x->x_window * 2) * sizeof (t_sample), (window * 2) * sizeof (t_sample));
    x->x_fftwInForwardWindow = (t_sample *)t_resizebytes (x->x_fftwInForwardWindow, x->x_window * sizeof (t_sample), window * sizeof (t_sample));
    x->x_fftwInBackWindow = (t_sample *)t_resizebytes (x->x_fftwInBackWindow, x->x_window * sizeof (t_sample), window * sizeof (t_sample));

    x->x_blackman = (t_float *)t_resizebytes (x->x_blackman, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_cosine = (t_float *)t_resizebytes (x->x_cosine, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_hamming = (t_float *)t_resizebytes (x->x_hamming, x->x_window * sizeof (t_float), window * sizeof (t_float));
    x->x_hann = (t_float *)t_resizebytes (x->x_hann, x->x_window * sizeof (t_float), window * sizeof (t_float));

    x->x_window = window;
    x->x_windowHalf = windowHalf;

    // free the FFTW output buffer, and re-malloc according to new window
    fftwf_free (x->x_fftwOutForwardWindow);
    fftwf_free (x->x_fftwOutBackWindow);

    // destroy old plan, which depended on x->x_window
    fftwf_destroy_plan (x->x_fftwPlanForwardWindow);
    fftwf_destroy_plan (x->x_fftwPlanBackWindow);

    // allocate new FFTW output buffer memory
    x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);
    x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);

    // create a new plan
    x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
    x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

     // we're supposed to initialize the input array after we create the plan
    for (i = 0; i < x->x_window; i++)
     {
        x->x_fftwInForwardWindow[i] = 0.0;
        x->x_fftwInBackWindow[i] = 0.0;
    }

    // initialize signal buffer
    for (i = 0; i < (x->x_window * 2); i++)
        x->x_signalBuffer[i] = 0.0;

    // re-init window functions
    tIDLib_blackmanWindow (x->x_blackman, x->x_window);
    tIDLib_cosineWindow (x->x_cosine, x->x_window);
    tIDLib_hammingWindow (x->x_hamming, x->x_window);
    tIDLib_hannWindow (x->x_hann, x->x_window);

    post ("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void tempo_tilde_hop (t_tempo_tilde* x, t_floatarg h)
{
    t_float loTempo, hiTempo;

    if (h > x->x_window)
    {
        post ("%s analysis hop cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
        x->x_hop = x->x_windowHalf;
    }
    else if (h <= 0)
    {
        post ("%s analysis hop must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
        x->x_hop = x->x_windowHalf;
    }
    else
        x->x_hop = h;

    // recalc tempo indices based on new hop
    loTempo = x->x_loTempo;
    loTempo = 60.0f/loTempo; // BPM to seconds
    loTempo *= x->x_sr; // seconds to samples
    loTempo /= x->x_hop;//samples to frames
    x->x_loTempoIdx = roundf (loTempo);

       hiTempo = x->x_hiTempo;
    hiTempo = 60.0f/hiTempo; // BPM to seconds
    hiTempo *= x->x_sr; // seconds to samples
    hiTempo /= x->x_hop;//samples to frames
    x->x_hiTempoIdx = roundf (hiTempo);

    // since x_loTempoIdx changed, need to resize the onset buffer
    tempo_tilde_resizeOnsetsBuffer (x, x->x_loTempoIdx*(x->x_numHarm+1));

    post ("%s analysis hop: %i", x->x_objSymbol->s_name, x->x_hop);
}


static void tempo_tilde_overlap (t_tempo_tilde* x, t_floatarg o)
{
    // this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr / x->x_overlap;
    x->x_overlap = (o < 1) ? 1 : o;

    post ("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void tempo_tilde_squaredDiff (t_tempo_tilde* x, t_floatarg sd)
{
    sd = (sd < 0) ? 0 : sd;
    sd = (sd > 1) ? 1 : sd;
    x->x_squaredDiff = sd;

    if (x->x_squaredDiff)
        post ("%s using squared spectral difference", x->x_objSymbol->s_name);
    else
        post ("%s using absolute value of spectral difference", x->x_objSymbol->s_name);
}


static void tempo_tilde_windowFunction (t_tempo_tilde* x, t_floatarg f)
{
    f = (f < 0) ? 0 : f;
    f = (f > 4) ? 4 : f;
    x->x_windowFunction = f;

    switch (x->x_windowFunction)
    {
        case rectangular:
            post ("%s window function: rectangular.", x->x_objSymbol->s_name);
            break;
        case blackman:
            post ("%s window function: blackman.", x->x_objSymbol->s_name);
            break;
        case cosine:
            post ("%s window function: cosine.", x->x_objSymbol->s_name);
            break;
        case hamming:
            post ("%s window function: hamming.", x->x_objSymbol->s_name);
            break;
        case hann:
            post ("%s window function: hann.", x->x_objSymbol->s_name);
            break;
        default:
            break;
    };
}


static void tempo_tilde_powerSpectrum (t_tempo_tilde* x, t_floatarg spec)
{
    spec = (spec < 0) ? 0 : spec;
    spec = (spec > 1) ? 1 : spec;
    x->x_powerSpectrum = spec;

    if (x->x_powerSpectrum)
        post ("%s using power spectrum", x->x_objSymbol->s_name);
    else
        post ("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void tempo_tilde_freqRange (t_tempo_tilde* x, t_floatarg loFreq, t_floatarg hiFreq)
{
    loFreq = (loFreq<0)?0:loFreq;
    loFreq = (loFreq>(x->x_sr * 0.5))?x->x_sr * 0.5:loFreq;

    x->x_loFreq = loFreq;
    x->x_hiFreq = hiFreq;

    x->x_loBin = tIDLib_freq2bin (loFreq, x->x_window, x->x_sr);

    hiFreq = (hiFreq<0)?0:hiFreq;
    hiFreq = (hiFreq>(x->x_sr * 0.5))?x->x_sr * 0.5:hiFreq;

    x->x_hiBin = tIDLib_freq2bin (hiFreq, x->x_window, x->x_sr);

    post ("%s frequency range: %0.2fHz through %0.2fHz (bin %i through %i).", x->x_objSymbol->s_name, loFreq, hiFreq, x->x_loBin, x->x_hiBin);
}


static void tempo_tilde_numHarm (t_tempo_tilde* x, t_floatarg h)
{
    x->x_numHarm = (h<2)?2:h;

    // since numHarm changed, we need to resize the onsets buffer
    tempo_tilde_resizeOnsetsBuffer (x, x->x_loTempoIdx*(x->x_numHarm+1));

    post ("%s HPS number of harmonics: %i", x->x_objSymbol->s_name, x->x_numHarm);
}


static void tempo_tilde_tempoRange (t_tempo_tilde* x, t_floatarg loTempo, t_floatarg hiTempo)
{
    loTempo = (loTempo<10)?10:loTempo;
    loTempo = (loTempo>500)?500:loTempo;
    x->x_loTempo = loTempo;

    loTempo = 60.0f/loTempo; // BPM to seconds
    loTempo *= x->x_sr; // seconds to samples
    loTempo /= x->x_hop;//samples to frames
    x->x_loTempoIdx = roundf (loTempo);

    hiTempo = (hiTempo<10)?10:hiTempo;
    hiTempo = (hiTempo>500)?500:hiTempo;
       x->x_hiTempo = hiTempo;

    hiTempo = 60.0f/hiTempo; // BPM to seconds
    hiTempo *= x->x_sr; // seconds to samples
    hiTempo /= x->x_hop;//samples to frames
    x->x_hiTempoIdx = roundf (hiTempo);

    // since loTempoIdx changed, we need to resize the onsets buffer
    tempo_tilde_resizeOnsetsBuffer (x, x->x_loTempoIdx*(x->x_numHarm+1));

    post ("%s tempo range: %0.2fBPM through %0.2fBPM (frameIdx %i through %i).", x->x_objSymbol->s_name, x->x_loTempo, x->x_hiTempo, x->x_loTempoIdx, x->x_hiTempoIdx);
}


static void tempo_tilde_onsetsBufDur (t_tempo_tilde* x, t_floatarg dur)
{
    t_float newSizeSamps;
    t_uInt newSize;

    dur = (dur<100)?100:dur;

    // ms to samples
    newSizeSamps = (dur*x->x_sr/1000.0f);

    // samples to frames
    newSize = newSizeSamps/(t_float)x->x_hop;

    tempo_tilde_resizeOnsetsBuffer (x, newSize);
}


static void tempo_tilde_thresh (t_tempo_tilde* x, t_floatarg thresh)
{
    thresh = (thresh<0)?0:thresh;
    x->x_growthThresh = thresh;

    post ("%s growth thresh: %0.2f%%.", x->x_objSymbol->s_name, x->x_growthThresh);
}


static void tempo_tilde_belowThreshDefault (t_tempo_tilde* x, t_floatarg val)
{
    x->x_belowThreshDefault = val;

    post ("%s below threshold default value: %f.", x->x_objSymbol->s_name, x->x_belowThreshDefault);
}


static void tempo_tilde_resetMaxOnsetPeakVal (t_tempo_tilde* x)
{
    x->x_maxOnsetPeakVal = -1;

    post ("%s maximum onset peak reset", x->x_objSymbol->s_name);
}


static void tempo_tilde_maxOnsetUpdate (t_tempo_tilde* x, t_floatarg update)
{
    update = (update<0)?0:update;
    update = (update>1)?1:update;
    x->x_maxOnsetUpdate = update;

    post ("%s maximum onset peak updating: %i", x->x_objSymbol->s_name, x->x_maxOnsetUpdate);
}


static void tempo_tilde_tempoBufferSize (t_tempo_tilde* x, t_floatarg n)
{
    t_uShortInt i;

    n = (n<5)?5:n;

    x->x_tempoBuffer = (t_float *)t_resizebytes (x->x_tempoBuffer, x->x_tempoBufferSize * sizeof (t_float), n * sizeof (t_float));

    x->x_tempoBufferSize = n;

    for (i = 0; i < x->x_tempoBufferSize; i++)
        x->x_tempoBuffer[i] = -1;

    post ("%s tempo buffer size: %i", x->x_objSymbol->s_name, x->x_tempoBufferSize);
}


static void tempo_tilde_tempoBufferConst (t_tempo_tilde* x, t_floatarg val)
{
    t_uShortInt i;

    for (i = 0; i < x->x_tempoBufferSize; i++)
        x->x_tempoBuffer[i] = val;

    post ("%s tempo buffer filled with a constant of %0.2f", x->x_objSymbol->s_name, val);
}


static void tempo_tilde_tempoBufferClear (t_tempo_tilde* x)
{
    t_uShortInt i;

    for (i = 0; i < x->x_tempoBufferSize; i++)
        x->x_tempoBuffer[i] = -1;

    post ("%s tempo buffer reset", x->x_objSymbol->s_name);
}


static void* tempo_tilde_new (t_symbol* s, int argc, t_atom* argv)
{
    t_tempo_tilde* x = (t_tempo_tilde *)pd_new (tempo_tilde_class);
    t_float hopFloat, loTempo, hiTempo;
    t_sampIdx i;

    x->x_tempo = outlet_new (&x->x_obj, &s_float);
    x->x_conf = outlet_new (&x->x_obj, &s_float);
    x->x_tempoRaw = outlet_new (&x->x_obj, &s_float);
    x->x_onsetList = outlet_new (&x->x_obj, gensym ("list"));

    // store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
    x->x_objSymbol = s;

    switch (argc)
    {
        case 2:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            hopFloat = atom_getfloat (argv + 1);
            if (hopFloat > x->x_window)
            {
                post ("%s WARNING: analysis hop cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
                x->x_hop = x->x_window * 0.5;
            }
            else if (hopFloat < 0)
            {
                post ("%s WARNING: analysis hop must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
                x->x_hop = x->x_window * 0.5;
            }
            else
                x->x_hop = hopFloat;
            break;

        case 1:
            x->x_window = atom_getfloat (argv);
            if (x->x_window < TID_MINWINDOWSIZE)
            {
                x->x_window = TID_WINDOWSIZEDEFAULT;
                post ("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, TID_MINWINDOWSIZE, TID_WINDOWSIZEDEFAULT);
            }

            x->x_hop = x->x_window * 0.5;
            break;

        case 0:
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_hop = x->x_window * 0.5;
            break;

        default:
            post ("%s WARNING: Too many arguments supplied. Using default window size of %i, and analysis hop of %i.", x->x_objSymbol->s_name, TID_WINDOWSIZEDEFAULT, (t_sampIdx)(TID_WINDOWSIZEDEFAULT * 0.5));
            x->x_window = TID_WINDOWSIZEDEFAULT;
            x->x_hop = x->x_window * 0.5;
            break;
    }

    x->x_windowHalf = x->x_window * 0.5;
    x->x_loBin = 0;
    x->x_hiBin = x->x_windowHalf;
    x->x_growthThresh = 20.0; // this is percent of max peak found
    x->x_maxOnsetUpdate = true;
    x->x_belowThreshDefault = 1.0;
    x->x_numHarm = 6;

    x->x_sr = TID_SAMPLERATEDEFAULT;
    x->x_n = TID_BLOCKSIZEDEFAULT;
    x->x_overlap = 1;

    x->x_loFreq = 0;
    x->x_hiFreq = x->x_sr * 0.5;

    x->x_loTempo = 40;
    x->x_hiTempo = 240;
    x->x_lastGoodTempo = -1;
    x->x_maxOnsetPeakVal = -1.0;
    x->x_lastFirstOnsetPeakVal = -1.0;

    loTempo = 60.0f/x->x_loTempo; // BPM to seconds
    loTempo *= x->x_sr; // seconds to samples
    loTempo /= x->x_hop;//samples to frames
    x->x_loTempoIdx = roundf (loTempo);

    hiTempo = 60.0f/x->x_hiTempo; // BPM to seconds
    hiTempo *= x->x_sr; // seconds to samples
    hiTempo /= x->x_hop;//samples to frames
    x->x_hiTempoIdx = roundf (hiTempo);

    x->x_windowFunction = blackman;
    x->x_powerSpectrum = false;
    x->x_lastDspTime = clock_getlogicaltime();
    x->x_squaredDiff = true;
    x->x_dspTicks = 0;
    x->x_tempoBufferIdx = 0;
    x->x_tempoBufferSize = 20;
    x->x_debug = false;

    x->x_onsetsBufSize = x->x_loTempoIdx*(x->x_numHarm+1);
    x->x_onsetsBufDur = ((x->x_onsetsBufSize*x->x_hop)/x->x_sr)*1000.0;

    x->x_signalBuffer = (t_sample *)t_getbytes ((x->x_window * 2) * sizeof (t_sample));
    x->x_fftwInForwardWindow = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));
    x->x_fftwInBackWindow = (t_sample *)t_getbytes (x->x_window * sizeof (t_sample));
    x->x_onsetsBuffer = (t_float *)t_getbytes (x->x_onsetsBufSize * sizeof (t_float));
    x->x_listOut = (t_atom *)t_getbytes (x->x_onsetsBufSize * sizeof (t_atom));

     for (i = 0; i < (x->x_window * 2); i++)
        x->x_signalBuffer[i] = 0.0;

     for (i = 0; i < x->x_onsetsBufSize; i++)
        x->x_onsetsBuffer[i] = x->x_belowThreshDefault;

      x->x_blackman = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
      x->x_cosine = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
      x->x_hamming = (t_float *)t_getbytes (x->x_window * sizeof (t_float));
      x->x_hann = (t_float *)t_getbytes (x->x_window * sizeof (t_float));

     // initialize signal windowing functions
    tIDLib_blackmanWindow (x->x_blackman, x->x_window);
    tIDLib_cosineWindow (x->x_cosine, x->x_window);
    tIDLib_hammingWindow (x->x_hamming, x->x_window);
    tIDLib_hannWindow (x->x_hann, x->x_window);

    // set up the FFTW output buffer. Is there no function to initialize it?
    x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);
    x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex (x->x_windowHalf + 1);

    // DFT plan
    x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
    x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d (x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

    // we're supposed to initialize the input array after we create the plan
     for (i = 0; i < x->x_window; i++)
     {
        x->x_fftwInForwardWindow[i] = 0.0;
        x->x_fftwInBackWindow[i] = 0.0;
    }

      x->x_tempoBuffer = (t_float *)t_getbytes (x->x_tempoBufferSize * sizeof (t_float));

     for (i = 0; i < x->x_tempoBufferSize; i++)
        x->x_tempoBuffer[i] = -1;

    return (x);
}


static t_int* tempo_tilde_perform (t_int* w)
{
    t_uShortInt n;
    t_sampIdx i;

    t_tempo_tilde* x = (t_tempo_tilde *)(w[1]);

    t_sample* in = (t_float *)(w[2]);
    n = w[3];

     // shift signal buffer contents back.
    for (i = 0; i < (x->x_window * 2-n); i++)
        x->x_signalBuffer[i] = x->x_signalBuffer[i + n];

    // write new block to end of signal buffer.
    for (i = 0; i < n; i++)
        x->x_signalBuffer[(x->x_window * 2-n)+i] = in[i];

    x->x_dspTicks++;

     if (x->x_dspTicks*n >= x->x_hop)
     {
         x->x_dspTicks = 0;
         tempo_tilde_analyze (x);
     }

    x->x_lastDspTime = clock_getlogicaltime();

    return (w + 4);
}


static void tempo_tilde_dsp (t_tempo_tilde* x, t_signal** sp)
{
    dsp_add (
        tempo_tilde_perform,
        3,
        x,
        sp[0]->s_vec,
        sp[0]->s_n
    );

// compare sr to stored sr and update if different
    if (sp[0]->s_sr != x->x_sr * x->x_overlap)
    {
        t_sampIdx i;
        t_float loTempo, hiTempo;

        x->x_sr = sp[0]->s_sr / x->x_overlap;

        // reset dspTicks to start fresh at the new samplerate
        x->x_dspTicks = 0;

        // since we're resetting dspTicks, we might as well clear the audio buffer and start fresh at the new samplerate
        for (i = 0; i < (x->x_window * 2); i++)
            x->x_signalBuffer[i] = 0.0;

        // if samplerate changes, we need to update loTempoIdx, hiTempoIdx, and therefore resize the onsets buffer
        loTempo = 60.0f/x->x_loTempo; // BPM to seconds
        loTempo *= x->x_sr; // seconds to samples
        loTempo /= x->x_hop;//samples to frames
        x->x_loTempoIdx = roundf (loTempo);

        hiTempo = 60.0f/x->x_hiTempo; // BPM to seconds
        hiTempo *= x->x_sr; // seconds to samples
        hiTempo /= x->x_hop;//samples to frames
        x->x_hiTempoIdx = roundf (hiTempo);

        // since x_loTempoIdx changed, need to resize the onset buffer
        tempo_tilde_resizeOnsetsBuffer (x, x->x_loTempoIdx*(x->x_numHarm+1));

        post ("%s samplerate: %0.2f", x->x_objSymbol->s_name, x->x_sr);
    };

// compare n to stored n and update if different
    if (sp[0]->s_n != x->x_n)
    {
        x->x_n = sp[0]->s_n;
        x->x_lastDspTime = clock_getlogicaltime();
    }
};

static void tempo_tilde_free (t_tempo_tilde* x)
{
    // free the input buffer memory
    t_freebytes (x->x_signalBuffer, (x->x_window * 2) * sizeof (t_sample));
    t_freebytes (x->x_onsetsBuffer, x->x_onsetsBufSize * sizeof (t_float));

    // free FFTW stuff
    t_freebytes (x->x_fftwInForwardWindow, x->x_window * sizeof (t_sample));
    t_freebytes (x->x_fftwInBackWindow, x->x_window * sizeof (t_sample));
    fftwf_free (x->x_fftwOutForwardWindow);
    fftwf_free (x->x_fftwOutBackWindow);
    fftwf_destroy_plan (x->x_fftwPlanForwardWindow);
    fftwf_destroy_plan (x->x_fftwPlanBackWindow);

    // free the window memory
    t_freebytes (x->x_blackman, x->x_window * sizeof (t_float));
    t_freebytes (x->x_cosine, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hamming, x->x_window * sizeof (t_float));
    t_freebytes (x->x_hann, x->x_window * sizeof (t_float));

    t_freebytes (x->x_tempoBuffer, x->x_tempoBufferSize * sizeof (t_float));
    t_freebytes (x->x_listOut, x->x_onsetsBufSize * sizeof (t_atom));
}

void tempo_tilde_setup (void)
{
    tempo_tilde_class =
    class_new (
        gensym ("tempo~"),
        (t_newmethod)tempo_tilde_new,
        (t_method)tempo_tilde_free,
        sizeof (t_tempo_tilde),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    class_addcreator (
        (t_newmethod)tempo_tilde_new,
        gensym ("timbreIDLib/tempo~"),
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN (tempo_tilde_class, t_tempo_tilde, x_f);

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_debug,
        gensym ("debug"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_window,
        gensym ("window"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_hop,
        gensym ("hop"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_overlap,
        gensym ("overlap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_squaredDiff,
        gensym ("squared_diff"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_windowFunction,
        gensym ("window_function"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_powerSpectrum,
        gensym ("power_spectrum"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_freqRange,
        gensym ("freq_range"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_numHarm,
        gensym ("harmonics"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_tempoRange,
        gensym ("tempo_range"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_onsetsBufDur,
        gensym ("onset_buf_duration"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_thresh,
        gensym ("thresh"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_belowThreshDefault,
        gensym ("below_thresh_default"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_resetMaxOnsetPeakVal,
        gensym ("reset_max_onset_peak"),
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_maxOnsetUpdate,
        gensym ("max_onset_peak_update"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_tempoBufferSize,
        gensym ("tempo_buffer_size"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_tempoBufferConst,
        gensym ("tempo_buffer_const"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_tempoBufferClear,
        gensym ("tempo_buffer_clear"),
        0
    );

    class_addmethod (
        tempo_tilde_class,
        (t_method)tempo_tilde_dsp,
        gensym ("dsp"),
        0
    );
}
