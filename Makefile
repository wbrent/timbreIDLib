# Makefile for timbreIDLib

# specify a location for Pd if desired
# PDDIR = /home/yourname/somedirectory/pd-0.49-0

lib.name = timbreIDLib

# specify the location and name of the FFTW library
ldlibs = -L/usr/local/lib -lfftw3f

# specify the location of FFTW header file
cflags = -Iinclude -I/usr/local/include

lib.setup.sources = src/$(lib.name).c

class.sources = src/attackTime.c src/attackTime~.c src/autoCorrPitch.c src/autoCorrPitch~.c src/bark.c src/bark2freq.c src/barkSpec.c src/barkSpecBrightness.c src/barkSpecBrightness~.c src/barkSpecCentroid.c src/barkSpecCentroid~.c src/barkSpecFlatness.c src/barkSpecFlatness~.c src/barkSpecFlux.c src/barkSpecFlux~.c src/barkSpecIrregularity.c src/barkSpecIrregularity~.c src/barkSpecKurtosis.c src/barkSpecKurtosis~.c src/barkSpecRolloff.c src/barkSpecRolloff~.c src/barkSpecSkewness.c src/barkSpecSkewness~.c src/barkSpecSlope.c src/barkSpecSlope~.c src/barkSpecSpread.c src/barkSpecSpread~.c src/barkSpec~.c src/bark~.c src/bfcc.c src/bfcc~.c src/bin2freq.c src/binWrangler.c src/cepstrum.c src/cepstrum~.c src/cepstrumPitch.c src/cepstrumPitch~.c src/chroma.c src/chroma~.c src/dct.c src/dct~.c src/energy.c src/energy~.c src/energyEntropy.c src/energyEntropy~.c src/featureAccum.c src/featureDelta.c src/featureNorm.c src/freq2bark.c src/freq2bin.c src/freq2mel.c src/harmonicRatio.c src/harmonicRatio~.c src/magSpec.c src/magSpec~.c src/maxSample.c src/maxSample~.c src/maxSampleDelta.c src/maxSampleDelta~.c src/mel2freq.c src/melSpec.c src/melSpec~.c src/mfcc.c src/mfcc~.c src/minSample.c src/minSample~.c src/minSampleDelta.c src/minSampleDelta~.c src/nearestPoint.c src/peakSample.c src/peakSample~.c src/phaseSpec.c src/phaseSpec~.c src/sampleBuffer~.c src/specBrightness.c src/specBrightness~.c src/specCentroid.c src/specCentroid~.c src/specFlatness.c src/specFlatness~.c src/specFlux.c src/specFlux~.c src/specHarmonicity.c src/specHarmonicity~.c src/specIrregularity.c src/specIrregularity~.c src/specKurtosis.c src/specKurtosis~.c src/specRolloff.c src/specRolloff~.c src/specSkewness.c src/specSkewness~.c src/specSlope.c src/specSlope~.c src/specSpread.c src/specSpread~.c src/tempo~.c src/tID_fft.c src/tID_fft~.c src/tID_mean.c src/tID_std.c src/timbreID.c src/tabletool.c src/waveSlope.c src/waveSlope~.c src/waveNoise.c src/waveNoise~.c src/zeroCrossing.c src/zeroCrossing~.c

common.sources = src/tIDLib.c

datafiles = help/attackTime-help.pd help/attackTime~-help.pd help/autoCorrPitch-help.pd help/autoCorrPitch~-help.pd help/bark-help.pd help/bark2freq-help.pd help/barkSpec-help.pd help/barkSpecBrightness-help.pd help/barkSpecBrightness~-help.pd help/barkSpecCentroid-help.pd help/barkSpecCentroid~-help.pd help/barkSpecFlatness-help.pd help/barkSpecFlatness~-help.pd help/barkSpecFlux-help.pd help/barkSpecFlux~-help.pd help/barkSpecIrregularity-help.pd help/barkSpecIrregularity~-help.pd help/barkSpecKurtosis-help.pd help/barkSpecKurtosis~-help.pd help/barkSpecRolloff-help.pd help/barkSpecRolloff~-help.pd help/barkSpecSkewness-help.pd help/barkSpecSkewness~-help.pd help/barkSpecSlope-help.pd help/barkSpecSlope~-help.pd help/barkSpecSpread-help.pd help/barkSpecSpread~-help.pd help/barkSpec~-help.pd help/bark~-help.pd help/bfcc-help.pd help/bfcc~-help.pd help/bin2freq-help.pd help/binWrangler-help.pd help/cepstrum-help.pd help/cepstrum~-help.pd help/cepstrumPitch-help.pd help/cepstrumPitch~-help.pd help/chroma-help.pd help/chroma~-help.pd help/dct-help.pd help/dct~-help.pd help/energy-help.pd help/energy~-help.pd help/energyEntropy-help.pd help/energyEntropy~-help.pd help/featureAccum-help.pd help/featureNorm-help.pd help/featureDelta-help.pd help/freq2bark-help.pd help/freq2bin-help.pd help/freq2mel-help.pd help/harmonicRatio-help.pd help/harmonicRatio~-help.pd help/magSpec-help.pd help/magSpec~-help.pd help/maxSample-help.pd help/maxSample~-help.pd help/maxSampleDelta-help.pd help/maxSampleDelta~-help.pd help/mel2freq-help.pd help/melSpec-help.pd help/melSpec~-help.pd help/mfcc-help.pd help/mfcc~-help.pd help/minSample-help.pd help/minSample~-help.pd help/minSampleDelta-help.pd help/minSampleDelta~-help.pd help/nearestPoint-help.pd help/peakSample-help.pd help/peakSample~-help.pd help/phaseSpec-help.pd help/phaseSpec~-help.pd help/sampleBuffer~-help.pd help/specBrightness-help.pd help/specBrightness~-help.pd help/specCentroid-help.pd help/specCentroid~-help.pd help/specFlatness-help.pd help/specFlatness~-help.pd help/specFlux-help.pd help/specFlux~-help.pd help/specHarmonicity-help.pd help/specHarmonicity~-help.pd help/specIrregularity-help.pd help/specIrregularity~-help.pd help/specKurtosis-help.pd help/specKurtosis~-help.pd help/specRolloff-help.pd help/specRolloff~-help.pd help/specSkewness-help.pd help/specSkewness~-help.pd help/specSlope-help.pd help/specSlope~-help.pd help/specSpread-help.pd help/specSpread~-help.pd help/tempo~-help.pd help/tID_fft-help.pd help/tID_fft~-help.pd help/tID_mean-help.pd help/tID_std-help.pd help/tIDLib-help.pd help/timbreID-help.pd help/tabletool-help.pd help/waveSlope-help.pd help/waveSlope~-help.pd help/waveNoise-help.pd help/waveNoise~-help.pd help/zeroCrossing-help.pd help/zeroCrossing~-help.pd README.md LICENSE INSTALL.txt

# build a multi-object library
make-lib-executable=yes

# provide the path to pd-lib-builder
PDLIBBUILDER_DIR=./pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
