/*

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

version 0.9.0A, May 20, 2021

*/

#include <math.h>  // for fabs(), cos(), etc
#include <stdlib.h>  // for abs()
//#include <stdio.h>
//#include <stdbool.h>
#include <string.h>  // for memset() and strlen()
#include <limits.h>  // for INT_MAX, etc
#include <float.h>
#include <fftw3.h>
#include "m_pd.h"

#define TIDVERSION "0.9.0A"

// choose either FFTW_MEASURE or FFTW_ESTIMATE here.
#define FFTWPLANNERFLAG FFTW_ESTIMATE

#define MINBARKSPACING 0.1
#define MAXBARKSPACING 6.0
#define BARKSPACINGDEFAULT 0.5
#define MINMELSPACING 5.0
#define MAXMELSPACING 1000.0
#define MELSPACINGDEFAULT 100.0
#define MAXBARKS 26.0
#define MAXMELS 3962.0
#define MAXBARKFREQ 22855.4
#define MAXMELFREQ 22843.6
#define WINDOWSIZEDEFAULT 1024
#define MINWINDOWSIZE 4
#define SAMPLERATEDEFAULT 44100
#define MINSAMPLERATE 64
#define BLOCKSIZEDEFAULT 64
#define NUMWEIGHTPOINTS 29
#define MAXTIDTEXTSTRING 100000
#define PBINRANGEBUFSIZE 64
// this can just be 64 because the most octaves we'd have even with a pitch below MIDI 21 is about 10, so we only need around 20 pairs of bin indices at most.

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


typedef signed char t_sChar; // -128 to 127
typedef unsigned char t_uChar; // 0 to 255
typedef signed int t_sInt; // -2,147,483,648 to 2,147,483,647
typedef signed long int t_sLongInt; // -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807
typedef unsigned short int t_uShortInt; // 0 to 65,535
typedef unsigned int t_uInt; // 0 to 4,294,967,295
typedef unsigned long int t_uLongInt; // 0 to 18,446,744,073,709,551,615

// indices to sample arrays. useful for window sizes as well. non-negative only
typedef t_uLongInt t_sampIdx; // gets us ~27 hours at 44.1kHz. this is plenty and unsigned long long is not needed
typedef t_sampIdx t_binIdx; // should be the same as sampIdx if we're using the former for window size N. non-negative only
typedef t_uShortInt t_filterIdx; // can be small because min bark spacing allowed forces filterbanks to be small. non-negative only
// indices to attributes. non-negative only
typedef t_uInt t_attributeIdx; // gets us over 4 billion attributes. should be plenty
// indices to instances in database. members of a cluster. non-negative only
typedef t_uInt t_instanceIdx; // again, gets us over 4 billion instances


typedef enum
{
    false = 0,
    true
} t_bool;

typedef enum
{
    rectangular = 0,
    blackman,
    cosine,
    hamming,
    hann
} t_windowFunction;

typedef enum
{
    freq2barkFormula0 = 0,
    freq2barkFormula1,
    freq2barkFormula2
} t_freq2barkFormula;

typedef enum
{
    bark2freqFormula0 = 0,
    bark2freqFormula1,
    bark2freqFormula2
} t_bark2freqFormula;

typedef enum
{
    freq2melFormula0 = 0,
    freq2melFormula1
} t_freq2melFormula;

typedef enum
{
    mel2freqFormula0 = 0,
    mel2freqFormula1
} t_mel2freqFormula;

typedef enum
{
    mFlux = 0,
    mGrowth,
    mDecay
} t_fluxMode;

typedef struct filter
{
    t_float *filter;
    t_binIdx filterSize;
    t_binIdx indices[2];
    t_float filterFreqs[2];
} t_filter;

typedef struct knnInfo
{
    t_float dist;
    t_float safeDist;
    t_instanceIdx idx;
    t_instanceIdx cluster;
} t_knnInfo;

typedef struct instance
{
    t_float *data;
    t_attributeIdx length;
    t_instanceIdx clusterMembership;
    t_knnInfo knnInfo;
} t_instance;

typedef struct cluster
{
    t_instanceIdx *members;
    t_instanceIdx numMembers;
    t_uInt votes;
} t_cluster;

typedef struct normData
{
    t_float max;
    t_float min;
    t_float normScalar;
} t_normData;

typedef struct attributeData
{
    t_float inputData;
    t_normData normData;
    t_float weight;
    t_attributeIdx order;
    t_symbol *name;
} t_attributeData;


/* ---------------- conversion functions ---------------------- */
t_float tIDLib_freq2bin(t_float freq, t_float n, t_float sr);
t_float tIDLib_bin2freq(t_float bin, t_float n, t_float sr);
t_float tIDLib_freq2bark(t_float freq);
t_float tIDLib_bark2freq(t_float bark);
t_float tIDLib_freq2mel(t_float freq);
t_float tIDLib_mel2freq(t_float mel);
/* ---------------- END conversion functions ---------------------- */


/* ---------------- utility functions ---------------------- */
void tIDLib_linspace(t_float *ramp, t_float start, t_float finish, t_binIdx n);
t_sChar tIDLib_signum(t_float sample);
t_float tIDLib_fitLineSlope(t_sampIdx n, t_float *input);
t_float tIDLib_hps(t_float *data, t_uInt n, t_float loIdx, t_float hiIdx, t_uShortInt numHarm, t_float *yValues, t_float *maxYValue, t_bool debug);
t_float tIDLib_mode(t_float *data, t_uLongInt n, t_uLongInt *countOut);
void tIDLib_bubbleSort(t_sampIdx n, t_float *list);
void tIDLib_knnInfoBubbleSort(t_uShortInt n, t_instance *instances);
void tIDLib_sortKnnInfo(t_uShortInt k, t_instanceIdx numInstances, t_instanceIdx prevMatch, t_instance *instances);
t_float tIDLib_dotProd(t_attributeIdx n, t_float *v1, t_float *v2);
t_float tIDLib_euclidDist(t_attributeIdx n, t_float *v1, t_float *v2, t_float *weights, t_bool sqroot);
t_float tIDLib_taxiDist(t_attributeIdx n, t_float *v1, t_float *v2, t_float *weights);
t_float tIDLib_corr(t_attributeIdx n, t_float *v1, t_float *v2);
void tIDLib_peaksValleys(t_sampIdx n, t_float *data, t_float *flags, t_float *minVal, t_float *maxVal);
/* ---------------- END utility functions ---------------------- */


/* ---------------- filterbank functions ---------------------- */
t_binIdx tIDLib_nearestBinIndex(t_float target, t_float *binFreqs, t_binIdx n);
t_filterIdx tIDLib_getBarkBoundFreqs(t_float **filterFreqs, t_filterIdx oldSizeFilterFreqs, t_float spacing, t_float sr);
t_filterIdx tIDLib_getMelBoundFreqs(t_float **filterFreqs, t_filterIdx oldSizeFilterFreqs, t_float spacing, t_float sr);
void tIDLib_createFilterbank(t_float *filterFreqs, t_filter **filterbank, t_filterIdx oldNumFilters, t_filterIdx newNumFilters, t_float window, t_float sr);
void tIDLib_specFilterBands(t_binIdx n, t_filterIdx numFilters, t_float *spectrum, t_filter *filterbank, t_bool normalize);
void tIDLib_filterbankMultiply(t_float *spectrum, t_bool normalize, t_bool filterAvg, t_filter *filterbank, t_filterIdx numFilters);
void tIDLib_cosineTransform(t_float *output, t_sample *input, t_filterIdx numFilters);
/* ---------------- END filterbank functions ---------------------- */


/* ---------------- stat computation functions ---------------------- */
t_float tIDLib_computeCentroid(t_binIdx n, t_float *spectrum, t_float *freqList, t_float energySum);
t_float tIDLib_computeSpread(t_binIdx n, t_float *spectrum, t_float *freqList, t_float energySum, t_float centroid);
t_float tIDLib_computeSkewness(t_binIdx n, t_float *spectrum, t_float *freqList, t_float energySum, t_float centroid, t_float spread);
t_float tIDLib_computeKurtosis(t_binIdx n, t_float *spectrum, t_float *freqList, t_float energySum, t_float centroid, t_float spread);
/* ---------------- END stat computation functions ---------------------- */


/* ---------------- windowing buffer functions ---------------------- */
void tIDLib_blackmanWindow(t_float *wPtr, t_sampIdx n);
void tIDLib_cosineWindow(t_float *wPtr, t_sampIdx n);
void tIDLib_hammingWindow(t_float *wPtr, t_sampIdx n);
void tIDLib_hannWindow(t_float *wPtr, t_sampIdx n);
/* ---------------- END windowing buffer functions ---------------------- */


/* ---------------- dsp utility functions ---------------------- */
t_float tIDLib_ampDB(t_sampIdx n, t_sample *input);
void tIDLib_peakSample(t_sampIdx n, t_float *input, t_sampIdx *peakIdx, t_float *peakVal);
t_sampIdx tIDLib_findAttackStartSamp(t_sampIdx n, t_float *input, t_float sampDeltaThresh, t_uShortInt numSampsThresh);
t_float tIDLib_zeroCrossingRate(t_sampIdx n, t_sample *input, t_bool normalize);
t_uInt tIDLib_getPitchBinRanges(t_binIdx *binRanges, t_float thisPitch, t_float loFreq, t_float hiFreq, t_float pitchTolerance, t_sampIdx n, t_float sr);
void tIDLib_power(t_binIdx n, void *fftw_out, t_float *powBuf);
void tIDLib_mag(t_binIdx n, t_float *input);
void tIDLib_normal(t_binIdx n, t_float *input);
void tIDLib_normalPeak(t_binIdx n, t_float *input);
void tIDLib_log(t_binIdx n, t_float *input);
/* ---------------- END dsp utility functions ---------------------- */


/* ---------------- external setup function declarations ---------------------- */
void attackTime_setup(void);
void attackTime_tilde_setup(void);
void bark_setup(void);
void bark2freq_setup(void);
void barkSpec_setup(void);
void barkSpecBrightness_setup(void);
void barkSpecBrightness_tilde_setup(void);
void barkSpecCentroid_setup(void);
void barkSpecCentroid_tilde_setup(void);
void barkSpecFlatness_setup(void);
void barkSpecFlatness_tilde_setup(void);
void barkSpecFlux_setup(void);
void barkSpecFlux_tilde_setup(void);
void barkSpecIrregularity_setup(void);
void barkSpecIrregularity_tilde_setup(void);
void barkSpecKurtosis_setup(void);
void barkSpecKurtosis_tilde_setup(void);
void barkSpecRolloff_setup(void);
void barkSpecRolloff_tilde_setup(void);
void barkSpecSkewness_setup(void);
void barkSpecSkewness_tilde_setup(void);
void barkSpecSlope_setup(void);
void barkSpecSlope_tilde_setup(void);
void barkSpecSpread_setup(void);
void barkSpecSpread_tilde_setup(void);
void barkSpec_tilde_setup(void);
void bark_tilde_setup(void);
void bfcc_setup(void);
void bfcc_tilde_setup(void);
void bin2freq_setup(void);
void binWrangler_setup(void);
void cepstrum_setup(void);
void cepstrumPitch_setup(void);
void cepstrumPitch_tilde_setup(void);
void cepstrum_tilde_setup(void);
void chroma_setup(void);
void chroma_tilde_setup(void);
void dct_setup(void);
void dct_tilde_setup(void);
void featureAccum_setup(void);
void featureDelta_setup(void);
void featureNorm_setup(void);
void freq2bark_setup(void);
void freq2bin_setup(void);
void freq2mel_setup(void);
void magSpec_setup(void);
void magSpec_tilde_setup(void);
void maxSample_setup(void);
void maxSampleDelta_setup(void);
void maxSampleDelta_tilde_setup(void);
void maxSample_tilde_setup(void);
void mel2freq_setup(void);
void melSpec_setup(void);
void melSpec_tilde_setup(void);
void mfcc_setup(void);
void mfcc_tilde_setup(void);
void minSample_setup(void);
void minSampleDelta_setup(void);
void minSampleDelta_tilde_setup(void);
void minSample_tilde_setup(void);
void nearestPoint_setup(void);
void peakSample_setup(void);
void peakSample_tilde_setup(void);
void phaseSpec_setup(void);
void phaseSpec_tilde_setup(void);
void sampleBuffer_tilde_setup(void);
void specBrightness_setup(void);
void specBrightness_tilde_setup(void);
void specCentroid_setup(void);
void specCentroid_tilde_setup(void);
void specFlatness_setup(void);
void specFlatness_tilde_setup(void);
void specFlux_setup(void);
void specFlux_tilde_setup(void);
void specHarmonicity_setup(void);
void specHarmonicity_tilde_setup(void);
void specIrregularity_setup(void);
void specIrregularity_tilde_setup(void);
void specKurtosis_setup(void);
void specKurtosis_tilde_setup(void);
void specRolloff_setup(void);
void specRolloff_tilde_setup(void);
void specSkewness_setup(void);
void specSkewness_tilde_setup(void);
void specSlope_setup(void);
void specSlope_tilde_setup(void);
void specSpread_setup(void);
void specSpread_tilde_setup(void);
void tempo_tilde_setup(void);
void tID_fft_setup(void);
void tID_fft_tilde_setup(void);
void tID_mean_setup(void);
void tID_std_setup(void);
void tabletool_setup(void);
void timbreID_setup(void);
void waveSlope_setup(void);
void waveSlope_tilde_setup(void);
void waveDirChange_setup(void);
void waveDirChange_tilde_setup(void);
void zeroCrossing_setup(void);
void zeroCrossing_tilde_setup(void);
/* ---------------- END external setup function declarations ---------------------- */
