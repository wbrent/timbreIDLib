/*

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

typedef struct timbreIDLib
{
    t_object t_obj;
} t_timbreIDLib;

t_class* timbreIDLib_class;


static void* timbreIDLib_new (void)
{
    t_timbreIDLib* x = (t_timbreIDLib *)pd_new (timbreIDLib_class);

    return (x);
}

void timbreIDLib_setup (void)
{
    timbreIDLib_class = class_new (gensym ("timbreIDLib"), timbreIDLib_new, 0, sizeof (t_timbreIDLib), CLASS_PD, 0);

    post ("timbreIDLib version %s", TID_VERSION);

    attackTime_setup();
    attackTime_tilde_setup();
    bark_setup();
    bark2freq_setup();
    barkSpec_setup();
    barkSpecBrightness_setup();
    barkSpecBrightness_tilde_setup();
    barkSpecCentroid_setup();
    barkSpecCentroid_tilde_setup();
    barkSpecFlatness_setup();
    barkSpecFlatness_tilde_setup();
    barkSpecFlux_setup();
    barkSpecFlux_tilde_setup();
    barkSpecIrregularity_setup();
    barkSpecIrregularity_tilde_setup();
    barkSpecKurtosis_setup();
    barkSpecKurtosis_tilde_setup();
    barkSpecRolloff_setup();
    barkSpecRolloff_tilde_setup();
    barkSpecSkewness_setup();
    barkSpecSkewness_tilde_setup();
    barkSpecSlope_setup();
    barkSpecSlope_tilde_setup();
    barkSpecSpread_setup();
    barkSpecSpread_tilde_setup();
    barkSpec_tilde_setup();
    bark_tilde_setup();
    bfcc_setup();
    bfcc_tilde_setup();
    bin2freq_setup();
    binWrangler_setup();
    cepstrum_setup();
    cepstrumPitch_setup();
    cepstrumPitch_tilde_setup();
    cepstrum_tilde_setup();
    chroma_setup();
    chroma_tilde_setup();
    dct_setup();
    dct_tilde_setup();
    energy_setup();
    energy_tilde_setup();
    energyEntropy_setup();
    energyEntropy_tilde_setup();
    featureAccum_setup();
    featureDelta_setup();
    featureNorm_setup();
    freq2bark_setup();
    freq2bin_setup();
    freq2mel_setup();
    magSpec_setup();
    magSpec_tilde_setup();
    maxSample_setup();
    maxSampleDelta_setup();
    maxSampleDelta_tilde_setup();
    maxSample_tilde_setup();
    mel2freq_setup();
    melSpec_setup();
    melSpec_tilde_setup();
    mfcc_setup();
    mfcc_tilde_setup();
    minSample_setup();
    minSampleDelta_setup();
    minSampleDelta_tilde_setup();
    minSample_tilde_setup();
    nearestPoint_setup();
    peakSample_setup();
    peakSample_tilde_setup();
    phaseSpec_setup();
    phaseSpec_tilde_setup();
    sampleBuffer_tilde_setup();
    specBrightness_setup();
    specBrightness_tilde_setup();
    specCentroid_setup();
    specCentroid_tilde_setup();
    specFlatness_setup();
    specFlatness_tilde_setup();
    specFlux_setup();
    specFlux_tilde_setup();
    specHarmonicity_setup();
    specHarmonicity_tilde_setup();
    specIrregularity_setup();
    specIrregularity_tilde_setup();
    specKurtosis_setup();
    specKurtosis_tilde_setup();
    specRolloff_setup();
    specRolloff_tilde_setup();
    specSkewness_setup();
    specSkewness_tilde_setup();
    specSlope_setup();
    specSlope_tilde_setup();
    specSpread_setup();
    specSpread_tilde_setup();
    tempo_tilde_setup();
    tID_fft_setup();
    tID_fft_tilde_setup();
    tID_mean_setup();
    tID_std_setup();
    tabletool_setup();
    timbreID_setup();
    waveSlope_setup();
    waveSlope_tilde_setup();
    waveNoise_setup();
    waveNoise_tilde_setup();
    zeroCrossing_setup();
    zeroCrossing_tilde_setup();
}
