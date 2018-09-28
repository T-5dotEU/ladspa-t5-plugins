/* t5_3band_parameq_with_shelves.c

   Free software by Juergen Herrmann, t-5@t-5.eu. Do with it, whatever you
   want. No warranty. None, whatsoever. Also see license.txt .

   This LADSPA plugin provides a three band parametric equalizer with
   shelving low- and highpass filters based on biquad coefficients

   This file has poor memory protection. Failures during malloc() will
   not recover nicely. */

/*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <ladspa.h>
#include "helpers.h"

/*****************************************************************************/

#define SF_INPUT       0
#define SF_OUTPUT      1
#define SF_LOW_F       2
#define SF_LOW_G       3
#define SF_LOW_Q       4
#define SF_P1_F        5
#define SF_P1_G        6
#define SF_P1_Q        7
#define SF_P2_F        8
#define SF_P2_G        9
#define SF_P2_Q       10
#define SF_P3_F       11
#define SF_P3_G       12
#define SF_P3_Q       13
#define SF_HIGH_F     14
#define SF_HIGH_G     15
#define SF_HIGH_Q     16
#define SF_GAIN       17
#define SF_MMAPFNAME  18
#define PORTCOUNT     19

/*****************************************************************************/

/* Instance data for the ThreeBandParametricEqWithShelves filter */
typedef struct {

    long m_created_ns;
    time_t m_created_s;

    LADSPA_Data * m_mmapArea;

    LADSPA_Data m_fSampleRate;
    // previous input samples of biquad filters
    LADSPA_Data m_fxnm1_LOW;
    LADSPA_Data m_fxnm2_LOW;
    LADSPA_Data m_fxnm1_P1;
    LADSPA_Data m_fxnm2_P1;
    LADSPA_Data m_fxnm1_P2;
    LADSPA_Data m_fxnm2_P2;
    LADSPA_Data m_fxnm1_P3;
    LADSPA_Data m_fxnm2_P3;
    LADSPA_Data m_fxnm1_HIGH;
    LADSPA_Data m_fxnm2_HIGH;
    // previous output samples of biquad filters
    LADSPA_Data m_fynm1_LOW;
    LADSPA_Data m_fynm2_LOW;
    LADSPA_Data m_fynm1_P1;
    LADSPA_Data m_fynm2_P1;
    LADSPA_Data m_fynm1_P2;
    LADSPA_Data m_fynm2_P2;
    LADSPA_Data m_fynm1_P3;
    LADSPA_Data m_fynm2_P3;
    LADSPA_Data m_fynm1_HIGH;
    LADSPA_Data m_fynm2_HIGH;
    // port pointers
    LADSPA_Data * m_pfInput;
    LADSPA_Data * m_pfOutput;
    LADSPA_Data * m_pfLowF;
    LADSPA_Data * m_pfLowG;
    LADSPA_Data * m_pfLowQ;
    LADSPA_Data * m_pfP1F;
    LADSPA_Data * m_pfP1G;
    LADSPA_Data * m_pfP1Q;
    LADSPA_Data * m_pfP2F;
    LADSPA_Data * m_pfP2G;
    LADSPA_Data * m_pfP2Q;
    LADSPA_Data * m_pfP3F;
    LADSPA_Data * m_pfP3G;
    LADSPA_Data * m_pfP3Q;
    LADSPA_Data * m_pfHighF;
    LADSPA_Data * m_pfHighG;
    LADSPA_Data * m_pfHighQ;
    LADSPA_Data * m_pfGain;
    LADSPA_Data * m_pfMmapFname;

} ThreeBandParametricEqWithShelves;

/* Helpers... ****************************************************************/

BiquadCoeffs calcCoeffsLowShelf(float f, float g, float q, float samplerate) {
    BiquadCoeffs coeffs;
    float w0 = 2.0 * M_PI * f / samplerate;
    float alpha = sin(w0) / (2.0 * q);
    float A = pow(10, g / 40.0);
    float cs = cos(w0);
    float norm = 1 / ((A+1.0) + (A-1.0)*cs + 2.0*sqrt(A)*alpha);
    coeffs.b0 = norm * (    A*( (A+1.0) - (A-1.0)*cs + 2.0*sqrt(A)*alpha ));
    coeffs.b1 = norm * (2.0*A*( (A-1.0) - (A+1.0)*cs                     ));
    coeffs.b2 = norm * (    A*( (A+1.0) - (A-1.0)*cs - 2.0*sqrt(A)*alpha ));
    coeffs.a1 = norm * ( -2.0*( (A-1.0) + (A+1.0)*cs                     ));
    coeffs.a2 = norm * (        (A+1.0) + (A-1.0)*cs - 2.0*sqrt(A)*alpha);
    return coeffs;
}

BiquadCoeffs calcCoeffsPeaking(float f, float g, float q, float samplerate) {
    BiquadCoeffs coeffs;
    float w0 = 2.0 * M_PI * f / samplerate;
    float alpha = sin(w0) / (2.0 * q);
    float A = pow(10, g / 40.0);
    float cs = cos(w0);
    float norm = 1 / (1.0 + alpha / A);
    coeffs.b0 = norm * (1.0 + alpha * A);
    coeffs.b1 = norm * (-2.0 * cs);
    coeffs.b2 = norm * (1.0 - alpha * A);
    coeffs.a1 = norm * (-2.0 * cs);
    coeffs.a2 = norm * (1.0 - alpha / A);
    return coeffs;
}

BiquadCoeffs calcCoeffsHighShelf(float f, float g, float q, float samplerate) {
    BiquadCoeffs coeffs;
    float w0 = 2.0 * M_PI * f / samplerate;
    float alpha = sin(w0) / (2.0 * q);
    float A = pow(10, g / 40.0);
    float cs = cos(w0);
    float norm = 1 / ((A+1.0) - (A-1.0)*cs + 2.0*sqrt(A)*alpha);
    coeffs.b0 = norm * (     A*( (A+1.0) + (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha ));
    coeffs.b1 = norm * (-2.0*A*( (A-1.0) + (A+1.0)*cos(w0)                     ));
    coeffs.b2 = norm * (     A*( (A+1.0) + (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha ));
    coeffs.a1 = norm * (   2.0*( (A-1.0) - (A+1.0)*cos(w0)                     ));
    coeffs.a2 = norm * (         (A+1.0) - (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha);
    return coeffs;
}

void setupMmapFileForThreeBandParametricEqWithShelves(ThreeBandParametricEqWithShelves * psInstance) {
    // setup shared memory area to enable parametrization at runtime from external processes
    TimeMmapStruct ret;
    ret = setupMmapFile("3BandParamEqWithShelves", *(psInstance->m_pfMmapFname), PORTCOUNT); 
    psInstance->m_mmapArea = ret.mmap;
    psInstance->m_created_s = ret.s;
    psInstance->m_created_ns = ret.ns;
}

/*****************************************************************************/

/* Construct a new plugin instance. */
LADSPA_Handle instantiateThreeBandParametricEqWithShelves(const LADSPA_Descriptor * Descriptor,
                                  unsigned long SampleRate) {
    ThreeBandParametricEqWithShelves * psInstance;
    psInstance = (ThreeBandParametricEqWithShelves *)malloc(sizeof(ThreeBandParametricEqWithShelves));
    if (psInstance) {
        psInstance->m_fSampleRate = (LADSPA_Data)SampleRate;
        psInstance->m_mmapArea = NULL;
    }  
    return psInstance;
}

/*****************************************************************************/

/* Initialise and activate a plugin instance. */
void activateThreeBandParametricEqWithShelves(LADSPA_Handle Instance) {
    ThreeBandParametricEqWithShelves * psInstance;
    psInstance = (ThreeBandParametricEqWithShelves *)Instance;
    psInstance->m_fxnm1_LOW = 0;
    psInstance->m_fxnm2_LOW = 0;
    psInstance->m_fynm1_LOW = 0;
    psInstance->m_fynm2_LOW = 0;
    psInstance->m_fxnm1_P1 = 0;
    psInstance->m_fxnm2_P1 = 0;
    psInstance->m_fynm1_P1 = 0;
    psInstance->m_fynm2_P1 = 0;
    psInstance->m_fxnm1_P2 = 0;
    psInstance->m_fxnm2_P2 = 0;
    psInstance->m_fynm1_P2 = 0;
    psInstance->m_fynm2_P2 = 0;
    psInstance->m_fxnm1_P3 = 0;
    psInstance->m_fxnm2_P3 = 0;
    psInstance->m_fynm1_P3 = 0;
    psInstance->m_fynm2_P3 = 0;
    psInstance->m_fxnm1_HIGH = 0;
    psInstance->m_fxnm2_HIGH = 0;
    psInstance->m_fynm1_HIGH = 0;
    psInstance->m_fynm2_HIGH = 0;
}

/*****************************************************************************/

/* Connect a port to a data location.  */
void connectPortToThreeBandParametricEqWithShelves(LADSPA_Handle Instance,
                                                   unsigned long Port,
                                                   LADSPA_Data * DataLocation) {
  
    ThreeBandParametricEqWithShelves * psInstance;
    psInstance = (ThreeBandParametricEqWithShelves *)Instance;

    switch (Port) {
    case SF_INPUT:
        psInstance->m_pfInput = DataLocation;
        break;
    case SF_OUTPUT:
        psInstance->m_pfOutput = DataLocation;
        break;
    case SF_LOW_F:
        psInstance->m_pfLowF = DataLocation;
        break;
    case SF_LOW_G:
        psInstance->m_pfLowG = DataLocation;
        break;
    case SF_LOW_Q:
        psInstance->m_pfLowQ = DataLocation;
        break;
    case SF_P1_F:
        psInstance->m_pfP1F = DataLocation;
        break;
    case SF_P1_G:
        psInstance->m_pfP1G = DataLocation;
        break;
    case SF_P1_Q:
        psInstance->m_pfP1Q = DataLocation;
        break;
    case SF_P2_F:
        psInstance->m_pfP2F = DataLocation;
        break;
    case SF_P2_G:
        psInstance->m_pfP2G = DataLocation;
        break;
    case SF_P2_Q:
        psInstance->m_pfP2Q = DataLocation;
        break;
    case SF_P3_F:
        psInstance->m_pfP3F = DataLocation;
        break;
    case SF_P3_G:
        psInstance->m_pfP3G = DataLocation;
        break;
    case SF_P3_Q:
        psInstance->m_pfP3Q = DataLocation;
        break;
    case SF_HIGH_F:
        psInstance->m_pfHighF = DataLocation;
        break;
    case SF_HIGH_G:
        psInstance->m_pfHighG = DataLocation;
        break;
    case SF_HIGH_Q:
        psInstance->m_pfHighQ = DataLocation;
        break;
    case SF_GAIN:
        psInstance->m_pfGain = DataLocation;
        break;
    case SF_MMAPFNAME:
        psInstance->m_pfMmapFname = DataLocation;
        break;
    }
}

/*****************************************************************************/

/* Run the filter algorithm for a block of SampleCount samples. */
void runThreeBandParametricEqWithShelves(LADSPA_Handle Instance,
                                         unsigned long SampleCount) {

    LADSPA_Data * pfInput;
    LADSPA_Data * pfOutput;
    ThreeBandParametricEqWithShelves * psInstance;
    BiquadCoeffs coeffs_LOW, coeffs_P1, coeffs_P2, coeffs_P3, coeffs_HIGH;
    unsigned long lSampleIndex;
    LADSPA_Data unchanged = 0.0;
    LADSPA_Data changed;
    LADSPA_Data * mmptr;
    float fGainFactor;
    float xn, yn; // xn/yn holds currently processed input/output samples.
    float xnm1, xnm2, ynm1, ynm2; // holds previously processed input/output samples.
    // get ThreeBandParametricEqWithShelves Instance
    psInstance = (ThreeBandParametricEqWithShelves *)Instance;
    // get input and output buffers
    pfInput = psInstance->m_pfInput;
    pfOutput = psInstance->m_pfOutput;
    // memcpy parameters over from mmapped area
    mmptr = psInstance->m_mmapArea;
    if (mmptr != NULL) {
        memcpy(&changed, mmptr, sizeof(LADSPA_Data));
        if (changed != 0.0) {
            mmptr += 1;
            memcpy(psInstance->m_pfLowF, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfLowG, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfLowQ, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP1F, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP1G, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP1Q, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP2F, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP2G, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP2Q, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP3F, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP3G, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfP3Q, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfHighF, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfHighG, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfHighQ, mmptr, sizeof(LADSPA_Data));
            mmptr += 1;
            memcpy(psInstance->m_pfGain, mmptr, sizeof(LADSPA_Data));
        }
        // reset changed flag to re-enable parametrization by control inputs
        memcpy(psInstance->m_mmapArea, &unchanged, sizeof(LADSPA_Data));
    } else if (*(psInstance->m_pfMmapFname) != 0.0) {
        setupMmapFileForThreeBandParametricEqWithShelves(psInstance);
    }
    // calculate coeffs and gain factor
    coeffs_LOW = calcCoeffsLowShelf(*(psInstance->m_pfLowF),
                                    *(psInstance->m_pfLowG),
                                    *(psInstance->m_pfLowQ),
                                    psInstance->m_fSampleRate);
    coeffs_P1 = calcCoeffsPeaking(*(psInstance->m_pfP1F),
                                  *(psInstance->m_pfP1G),
                                  *(psInstance->m_pfP1Q),
                                  psInstance->m_fSampleRate);
    coeffs_P2 = calcCoeffsPeaking(*(psInstance->m_pfP2F),
                                  *(psInstance->m_pfP2G),
                                  *(psInstance->m_pfP2Q),
                                  psInstance->m_fSampleRate);
    coeffs_P3 = calcCoeffsPeaking(*(psInstance->m_pfP3F),
                                  *(psInstance->m_pfP3G),
                                  *(psInstance->m_pfP3Q),
                                  psInstance->m_fSampleRate);
    coeffs_HIGH = calcCoeffsHighShelf(*(psInstance->m_pfHighF),
                                      *(psInstance->m_pfHighG),
                                      *(psInstance->m_pfHighQ),
                                      psInstance->m_fSampleRate);
    fGainFactor = dbToGainFactor(*(psInstance->m_pfGain));
    // FILTER PROCESSING FOR LOW SHELF EQ ///////////////////////////////////////
    // get preciosuly processed samples
    xnm1 = psInstance->m_fxnm1_LOW;
    xnm2 = psInstance->m_fxnm2_LOW;
    ynm1 = psInstance->m_fynm1_LOW;
    ynm2 = psInstance->m_fynm2_LOW;
    // apply biquad calculation to buffers
    for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
        xn = *(pfInput++);
        yn = (coeffs_LOW.b0 * xn   + coeffs_LOW.b1 * xnm1 + coeffs_LOW.b2 * xnm2 - 
              coeffs_LOW.a1 * ynm1 - coeffs_LOW.a2 * ynm2);
        xnm2 = xnm1;
        xnm1 = xn;
        ynm2 = ynm1;
        ynm1 = yn;
        *(pfOutput++) = yn;
    }
    // store preciously calculated samples in BiQuad instance for later
    psInstance->m_fxnm1_LOW = xnm1;
    psInstance->m_fxnm2_LOW = xnm2;
    psInstance->m_fynm1_LOW = ynm1;
    psInstance->m_fynm2_LOW = ynm2;
    // FILTER PROCESSING FOR PEAKING EQ 1 ///////////////////////////////////////
    // reset output buffer pointer
    pfOutput = psInstance->m_pfOutput;
    // get preciosuly processed samples
    xnm1 = psInstance->m_fxnm1_P1;
    xnm2 = psInstance->m_fxnm2_P1;
    ynm1 = psInstance->m_fynm1_P1;
    ynm2 = psInstance->m_fynm2_P1;
    // apply biquad calculation to buffers
    for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
        xn = *(pfOutput);
        yn = (coeffs_P1.b0 * xn   + coeffs_P1.b1 * xnm1 + coeffs_P1.b2 * xnm2 - 
              coeffs_P1.a1 * ynm1 - coeffs_P1.a2 * ynm2);
        xnm2 = xnm1;
        xnm1 = xn;
        ynm2 = ynm1;
        ynm1 = yn;
        *(pfOutput++) = yn;
    }
    // store preciously calculated samples in BiQuad instance for later
    psInstance->m_fxnm1_P1 = xnm1;
    psInstance->m_fxnm2_P1 = xnm2;
    psInstance->m_fynm1_P1 = ynm1;
    psInstance->m_fynm2_P1 = ynm2;
    // FILTER PROCESSING FOR PEAKING EQ 2 ///////////////////////////////////////
    // reset output buffer pointer
    pfOutput = psInstance->m_pfOutput;
    // get preciosuly processed samples
    xnm1 = psInstance->m_fxnm1_P2;
    xnm2 = psInstance->m_fxnm2_P2;
    ynm1 = psInstance->m_fynm1_P2;
    ynm2 = psInstance->m_fynm2_P2;
    // apply biquad calculation to buffers
    for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
        xn = *(pfOutput);
        yn = (coeffs_P2.b0 * xn   + coeffs_P2.b1 * xnm1 + coeffs_P2.b2 * xnm2 - 
              coeffs_P2.a1 * ynm1 - coeffs_P2.a2 * ynm2);
        xnm2 = xnm1;
        xnm1 = xn;
        ynm2 = ynm1;
        ynm1 = yn;
        *(pfOutput++) = yn;
    }
    // store preciously calculated samples in BiQuad instance for later
    psInstance->m_fxnm1_P2 = xnm1;
    psInstance->m_fxnm2_P2 = xnm2;
    psInstance->m_fynm1_P2 = ynm1;
    psInstance->m_fynm2_P2 = ynm2;
    // FILTER PROCESSING FOR PEAKING EQ 3 ///////////////////////////////////////
    // reset output buffer pointer
    pfOutput = psInstance->m_pfOutput;
    // get preciosuly processed samples
    xnm1 = psInstance->m_fxnm1_P3;
    xnm2 = psInstance->m_fxnm2_P3;
    ynm1 = psInstance->m_fynm1_P3;
    ynm2 = psInstance->m_fynm2_P3;
    // apply biquad calculation to buffers
    for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
        xn = *(pfOutput);
        yn = (coeffs_P3.b0 * xn   + coeffs_P3.b1 * xnm1 + coeffs_P3.b2 * xnm2 - 
              coeffs_P3.a1 * ynm1 - coeffs_P3.a2 * ynm2);
        xnm2 = xnm1;
        xnm1 = xn;
        ynm2 = ynm1;
        ynm1 = yn;
        *(pfOutput++) = yn;
    }
    // store preciously calculated samples in BiQuad instance for later
    psInstance->m_fxnm1_P3 = xnm1;
    psInstance->m_fxnm2_P3 = xnm2;
    psInstance->m_fynm1_P3 = ynm1;
    psInstance->m_fynm2_P3 = ynm2;
    // FILTER PROCESSING FOR HIGH SHELF EQ //////////////////////////////////////
    // reset output buffer pointer
    pfOutput = psInstance->m_pfOutput;
    // get preciosuly processed samples
    xnm1 = psInstance->m_fxnm1_HIGH;
    xnm2 = psInstance->m_fxnm2_HIGH;
    ynm1 = psInstance->m_fynm1_HIGH;
    ynm2 = psInstance->m_fynm2_HIGH;
    // apply biquad calculation to buffers
    for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
        xn = *(pfOutput);
        yn = (coeffs_HIGH.b0 * xn   + coeffs_HIGH.b1 * xnm1 + coeffs_HIGH.b2 * xnm2 - 
              coeffs_HIGH.a1 * ynm1 - coeffs_HIGH.a2 * ynm2);
        xnm2 = xnm1;
        xnm1 = xn;
        ynm2 = ynm1;
        ynm1 = yn;
        *(pfOutput++) = yn * fGainFactor;
    }
    // store preciously calculated samples in BiQuad instance for later
    psInstance->m_fxnm1_HIGH = xnm1;
    psInstance->m_fxnm2_HIGH = xnm2;
    psInstance->m_fynm1_HIGH = ynm1;
    psInstance->m_fynm2_HIGH = ynm2;
}

/*****************************************************************************/

/* Throw away a ThreeBandParametricEqWithShelves instance. */
void cleanupThreeBandParametricEqWithShelves(LADSPA_Handle Instance) {
    ThreeBandParametricEqWithShelves * psInstance;
    psInstance = (ThreeBandParametricEqWithShelves *)Instance;
    if (psInstance->m_mmapArea != NULL) {
        cleanupMmapFile("3BandParamEqWithShelves",
                    *(psInstance->m_pfMmapFname),
                    psInstance->m_created_s,
                    psInstance->m_created_ns);
        free(Instance);
    }
}

/*****************************************************************************/

LADSPA_Descriptor * g_psThreeBandParametricEqWithShelvesInstanceDescriptor = NULL;

/*****************************************************************************/

/* _init() is called automatically when the plugin library is first loaded. */
void _init() {

    char ** pcPortNames;
    LADSPA_PortDescriptor * piPortDescriptors;
    LADSPA_PortRangeHint * psPortRangeHints;
    
    g_psThreeBandParametricEqWithShelvesInstanceDescriptor
      = (LADSPA_Descriptor *)malloc(sizeof(LADSPA_Descriptor));

    if (g_psThreeBandParametricEqWithShelvesInstanceDescriptor != NULL) {
    
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->UniqueID
            = 5541;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->Label
            = strdup("3band_parameq_with_shelves");
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->Properties
            = LADSPA_PROPERTY_HARD_RT_CAPABLE;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->Name 
            = strdup("T5's 3-Band Parametric with Shelves");
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->Maker
            = strdup("Juergen Herrmann (t-5@t-5.eu)");
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->Copyright
            = strdup("3-clause BSD licence");
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->PortCount
            = PORTCOUNT;
        piPortDescriptors
            = (LADSPA_PortDescriptor *)calloc(PORTCOUNT, sizeof(LADSPA_PortDescriptor));
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->PortDescriptors
            = (const LADSPA_PortDescriptor *)piPortDescriptors;
        piPortDescriptors[SF_INPUT]
            = LADSPA_PORT_INPUT | LADSPA_PORT_AUDIO;
        piPortDescriptors[SF_OUTPUT]
            = LADSPA_PORT_OUTPUT | LADSPA_PORT_AUDIO;
        piPortDescriptors[SF_LOW_F]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_LOW_G]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_LOW_Q]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P1_F]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P1_G]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P1_Q]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P2_F]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P2_G]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P2_Q]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P3_F]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P3_G]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_P3_Q]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_HIGH_F]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_HIGH_G]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_HIGH_Q]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_GAIN]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        piPortDescriptors[SF_MMAPFNAME]
            = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
        pcPortNames
            = (char **)calloc(PORTCOUNT, sizeof(char *));
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->PortNames 
            = (const char **)pcPortNames;
        pcPortNames[SF_INPUT]
            = strdup("Input");
        pcPortNames[SF_OUTPUT]
            = strdup("Output");
        pcPortNames[SF_LOW_F]
            = strdup("Low Shelf Frequency [Hz]");
        pcPortNames[SF_LOW_G]
            = strdup("Low Shelf Gain [dB]");
        pcPortNames[SF_LOW_Q]
            = strdup("Low Shelf Q");
        pcPortNames[SF_P1_F]
            = strdup("Peaking EQ 1 Frequency [Hz]");
        pcPortNames[SF_P1_G]
            = strdup("Peaking EQ 1 Gain [dB]");
        pcPortNames[SF_P1_Q]
            = strdup("Peaking EQ 1 Q");
        pcPortNames[SF_P2_F]
            = strdup("Peaking EQ 2 Frequency [Hz]");
        pcPortNames[SF_P2_G]
            = strdup("Peaking EQ 2 Gain [dB]");
        pcPortNames[SF_P2_Q]
            = strdup("Peaking EQ 2 Q");
        pcPortNames[SF_P3_F]
            = strdup("Peaking EQ 3 Frequency [Hz]");
        pcPortNames[SF_P3_G]
            = strdup("Peaking EQ 3 Gain [dB]");
        pcPortNames[SF_P3_Q]
            = strdup("Peaking EQ 3 Q");
        pcPortNames[SF_HIGH_F]
            = strdup("High Shelf Frequency [Hz]");
        pcPortNames[SF_HIGH_G]
            = strdup("High Shelf Gain [dB]");
        pcPortNames[SF_HIGH_Q]
            = strdup("High Shelf Q");
        pcPortNames[SF_GAIN]
            = strdup("Overall Gain [dB]");
        pcPortNames[SF_MMAPFNAME]
            = strdup("MMAP-Filename-Part");
        psPortRangeHints = ((LADSPA_PortRangeHint *)
            calloc(PORTCOUNT, sizeof(LADSPA_PortRangeHint)));
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->PortRangeHints
            = (const LADSPA_PortRangeHint *)psPortRangeHints;
        // Low Shelf --------------------------------------------------------- */
        psPortRangeHints[SF_LOW_F].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_SAMPLE_RATE
            | LADSPA_HINT_LOGARITHMIC
            | LADSPA_HINT_DEFAULT_440);
        psPortRangeHints[SF_LOW_F].LowerBound 
            = 0;
        psPortRangeHints[SF_LOW_F].UpperBound
            = 0.5;
        psPortRangeHints[SF_LOW_G].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_0);
        psPortRangeHints[SF_LOW_G].LowerBound 
            = -12;
        psPortRangeHints[SF_LOW_G].UpperBound
            = 12;
        psPortRangeHints[SF_LOW_Q].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_1);
        psPortRangeHints[SF_LOW_Q].LowerBound 
            = 0.1;
        psPortRangeHints[SF_LOW_Q].UpperBound
            = 10;
        // Peaking Parametric EQ 1 -------------------------------------------- */
        psPortRangeHints[SF_P1_F].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_SAMPLE_RATE
            | LADSPA_HINT_LOGARITHMIC
            | LADSPA_HINT_DEFAULT_440);
        psPortRangeHints[SF_P1_F].LowerBound 
            = 0;
        psPortRangeHints[SF_P1_F].UpperBound
            = 0.5;
        psPortRangeHints[SF_P1_G].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_0);
        psPortRangeHints[SF_P1_G].LowerBound 
            = -12;
        psPortRangeHints[SF_P1_G].UpperBound
            = 12;
        psPortRangeHints[SF_P1_Q].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_1);
        psPortRangeHints[SF_P1_Q].LowerBound 
            = 0.1;
        psPortRangeHints[SF_P1_Q].UpperBound
            = 10;
        // Peaking Parametric EQ 2 -------------------------------------------- */
        psPortRangeHints[SF_P2_F].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_SAMPLE_RATE
            | LADSPA_HINT_LOGARITHMIC
            | LADSPA_HINT_DEFAULT_440);
        psPortRangeHints[SF_P2_F].LowerBound 
            = 0;
        psPortRangeHints[SF_P2_F].UpperBound
            = 0.5;
        psPortRangeHints[SF_P2_G].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_0);
        psPortRangeHints[SF_P2_G].LowerBound 
            = -12;
        psPortRangeHints[SF_P2_G].UpperBound
            = 12;
        psPortRangeHints[SF_P2_Q].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_1);
        psPortRangeHints[SF_P2_Q].LowerBound 
            = 0.1;
        psPortRangeHints[SF_P2_Q].UpperBound
            = 10;
        // Peaking Parametric EQ 3 -------------------------------------------- */
        psPortRangeHints[SF_P3_F].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_SAMPLE_RATE
            | LADSPA_HINT_LOGARITHMIC
            | LADSPA_HINT_DEFAULT_440);
        psPortRangeHints[SF_P3_F].LowerBound 
            = 0;
        psPortRangeHints[SF_P3_F].UpperBound
            = 0.5;
        psPortRangeHints[SF_P3_G].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_0);
        psPortRangeHints[SF_P3_G].LowerBound 
            = -12;
        psPortRangeHints[SF_P3_G].UpperBound
            = 12;
        psPortRangeHints[SF_P3_Q].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_1);
        psPortRangeHints[SF_P3_Q].LowerBound 
            = 0.1;
        psPortRangeHints[SF_P3_Q].UpperBound
            = 10;
        // High Shelf -------------------------------------------------------- */
        psPortRangeHints[SF_HIGH_F].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_SAMPLE_RATE
            | LADSPA_HINT_LOGARITHMIC
            | LADSPA_HINT_DEFAULT_440);
        psPortRangeHints[SF_HIGH_F].LowerBound 
            = 0;
        psPortRangeHints[SF_HIGH_F].UpperBound
            = 0.5;
        psPortRangeHints[SF_HIGH_G].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_0);
        psPortRangeHints[SF_HIGH_G].LowerBound 
            = -12;
        psPortRangeHints[SF_HIGH_G].UpperBound
            = 12;
        psPortRangeHints[SF_HIGH_Q].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_1);
        psPortRangeHints[SF_HIGH_Q].LowerBound 
            = 0.1;
        psPortRangeHints[SF_HIGH_Q].UpperBound
            = 10;
        // Gain ------------------------------------------------------------ */
        psPortRangeHints[SF_GAIN].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_0);
        psPortRangeHints[SF_GAIN].LowerBound 
            = -12;
        psPortRangeHints[SF_GAIN].UpperBound
            = 12;
        // MMAP Filename --------------------------------------------------- */
        psPortRangeHints[SF_MMAPFNAME].HintDescriptor
            = (LADSPA_HINT_BOUNDED_BELOW 
            | LADSPA_HINT_BOUNDED_ABOVE
            | LADSPA_HINT_DEFAULT_0);
        psPortRangeHints[SF_MMAPFNAME].LowerBound 
            = 0;
        psPortRangeHints[SF_MMAPFNAME].UpperBound
            = 10000000000;
        // In- and Outpus -------------------------------------------------- */
        psPortRangeHints[SF_INPUT].HintDescriptor
            = 0;
        psPortRangeHints[SF_OUTPUT].HintDescriptor
            = 0;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->instantiate 
            = instantiateThreeBandParametricEqWithShelves;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->connect_port 
            = connectPortToThreeBandParametricEqWithShelves;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->activate
            = activateThreeBandParametricEqWithShelves;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->run
            = runThreeBandParametricEqWithShelves;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->run_adding
            = NULL;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->set_run_adding_gain
            = NULL;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->deactivate
            = NULL;
        g_psThreeBandParametricEqWithShelvesInstanceDescriptor->cleanup
            = cleanupThreeBandParametricEqWithShelves;
    }
}
  
/*****************************************************************************/

void deleteDescriptor(LADSPA_Descriptor * psDescriptor) {
    unsigned long lIndex;
    if (psDescriptor) {
        free((char *)psDescriptor->Label);
        free((char *)psDescriptor->Name);
        free((char *)psDescriptor->Maker);
        free((char *)psDescriptor->Copyright);
        free((LADSPA_PortDescriptor *)psDescriptor->PortDescriptors);
        for (lIndex = 0; lIndex < psDescriptor->PortCount; lIndex++)
            free((char *)(psDescriptor->PortNames[lIndex]));
        free((char **)psDescriptor->PortNames);
        free((LADSPA_PortRangeHint *)psDescriptor->PortRangeHints);
        free(psDescriptor);
    }
}

/*****************************************************************************/

/* _fini() is called automatically when the library is unloaded. */
void _fini() {
    deleteDescriptor(g_psThreeBandParametricEqWithShelvesInstanceDescriptor);
}

/*****************************************************************************/

/* Return a descriptor of the requested plugin types. */
const LADSPA_Descriptor * ladspa_descriptor(unsigned long Index) {
    /* Return the requested descriptor or null if the index is out of range. */
    switch (Index) {
    case 0:
        return g_psThreeBandParametricEqWithShelvesInstanceDescriptor;
    default:
        return NULL;
    }
}

/*****************************************************************************/

/* EOF */
