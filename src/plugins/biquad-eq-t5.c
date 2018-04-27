/* biquad-eq-t5.c

   Free software by Juergen Herrmann. Do with it, whatever you want.
   No warranty. None, whatsoever :)

   This LADSPA plugin provides a three band parametric equalizer with
   shelving low- and highpass filters based on biquad coefficients

   This file has poor memory protection. Failures during malloc() will
   not recover nicely. */

/*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

/*****************************************************************************/

#include <ladspa.h>

/*****************************************************************************/

#define SF_INPUT   0
#define SF_OUTPUT  1
#define SF_LOW_F   2
#define SF_LOW_G   3
#define SF_LOW_Q   4
#define SF_P1_F    5
#define SF_P1_G    6
#define SF_P1_Q    7
#define SF_P2_F    8
#define SF_P2_G    9
#define SF_P2_Q    10
#define SF_P3_F    11
#define SF_P3_G    12
#define SF_P3_Q    13
#define SF_HIGH_F  14
#define SF_HIGH_G  15
#define SF_HIGH_Q  16
#define SF_GAIN    17
#define PORTCOUNT  18

/*****************************************************************************/

/* Instance data for the biquad eq filter */
typedef struct {

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

} BiquadEq;

/* biquad coefficients */
typedef struct {

  float a0;
  float a1;
  float a2;
  float b0;
  float b1;
  float b2;

} BiquadCoeffs;

/* Helpers... ****************************************************************/

float dbToGainFactor(float db) {
  return pow(10, db/20.0);
}

BiquadCoeffs calcCoeffsLowShelf(float f, float g, float q, float samplerate) {
  BiquadCoeffs coeffs;
  float w0 = 2.0 * M_PI * f / samplerate;
  float alpha = sin(w0) / (2.0 * q);
  float A = pow(10, g / 40.0);
  coeffs.b0 =     A*( (A+1.0) - (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha );
  coeffs.b1 = 2.0*A*( (A-1.0) - (A+1.0)*cos(w0)                     );
  coeffs.b2 =     A*( (A+1.0) - (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha );
  coeffs.a0 =         (A+1.0) + (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha;
  coeffs.a1 =  -2.0*( (A-1.0) + (A+1.0)*cos(w0)                     );
  coeffs.a2 =         (A+1.0) + (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha;
  return coeffs;
}

BiquadCoeffs calcCoeffsPeaking(float f, float g, float q, float samplerate) {
  BiquadCoeffs coeffs;
  float w0 = 2.0 * M_PI * f / samplerate;
  float alpha = sin(w0) / (2.0 * q);
  float A = pow(10, g / 40.0);
  coeffs.b0 = 1.0 + alpha * A;
  coeffs.b1 = -2.0 * cos(w0);
  coeffs.b2 = 1.0 - alpha * A;
  coeffs.a0 = 1.0 + alpha / A;
  coeffs.a1 = -2.0 * cos(w0);
  coeffs.a2 = 1.0 - alpha / A;
  return coeffs;
}

BiquadCoeffs calcCoeffsHighShelf(float f, float g, float q, float samplerate) {
  BiquadCoeffs coeffs;
  float w0 = 2.0 * M_PI * f / samplerate;
  float alpha = sin(w0) / (2.0 * q);
  float A = pow(10, g / 40.0);
  coeffs.b0 =      A*( (A+1.0) + (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha );
  coeffs.b1 = -2.0*A*( (A-1.0) + (A+1.0)*cos(w0)                     );
  coeffs.b2 =      A*( (A+1.0) + (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha );
  coeffs.a0 =          (A+1.0) - (A-1.0)*cos(w0) + 2.0*sqrt(A)*alpha;
  coeffs.a1 =    2.0*( (A-1.0) - (A+1.0)*cos(w0)                     );
  coeffs.a2 =          (A+1.0) - (A-1.0)*cos(w0) - 2.0*sqrt(A)*alpha;
  return coeffs;
}

/*****************************************************************************/

/* Construct a new plugin instance. */
LADSPA_Handle instantiateBiquadEq(const LADSPA_Descriptor * Descriptor,
                    unsigned long SampleRate) {
  BiquadEq * psBiquadEq;
  psBiquadEq = (BiquadEq *)malloc(sizeof(BiquadEq));
  if (psBiquadEq) {
    psBiquadEq->m_fSampleRate = (LADSPA_Data)SampleRate;
  }
  return psBiquadEq;
}

/*****************************************************************************/

/* Initialise and activate a plugin instance. */
void activateBiquadEq(LADSPA_Handle Instance) {
  BiquadEq * psBiquadEq;
  psBiquadEq = (BiquadEq *)Instance;
  psBiquadEq->m_fxnm1_LOW = 0;
  psBiquadEq->m_fxnm2_LOW = 0;
  psBiquadEq->m_fynm1_LOW = 0;
  psBiquadEq->m_fynm1_LOW = 0;
  psBiquadEq->m_fxnm1_P1 = 0;
  psBiquadEq->m_fxnm2_P1 = 0;
  psBiquadEq->m_fynm1_P1 = 0;
  psBiquadEq->m_fynm1_P1 = 0;
  psBiquadEq->m_fxnm1_P2 = 0;
  psBiquadEq->m_fxnm2_P2 = 0;
  psBiquadEq->m_fynm1_P2 = 0;
  psBiquadEq->m_fynm1_P2 = 0;
  psBiquadEq->m_fxnm1_P3 = 0;
  psBiquadEq->m_fxnm2_P3 = 0;
  psBiquadEq->m_fynm1_P3 = 0;
  psBiquadEq->m_fynm1_P3 = 0;
  psBiquadEq->m_fxnm1_HIGH = 0;
  psBiquadEq->m_fxnm2_HIGH = 0;
  psBiquadEq->m_fynm1_HIGH = 0;
  psBiquadEq->m_fynm1_HIGH = 0;
}

/*****************************************************************************/

/* Connect a port to a data location.  */
void connectPortToBiquadEq(LADSPA_Handle Instance,
                           unsigned long Port,
                           LADSPA_Data * DataLocation) {
  
  BiquadEq * psBiquadEq;
  psBiquadEq = (BiquadEq *)Instance;

  switch (Port) {
  case SF_INPUT:
    psBiquadEq->m_pfInput = DataLocation;
    break;
  case SF_OUTPUT:
    psBiquadEq->m_pfOutput = DataLocation;
    break;
  case SF_LOW_F:
    psBiquadEq->m_pfLowF = DataLocation;
    break;
  case SF_LOW_G:
    psBiquadEq->m_pfLowG = DataLocation;
    break;
  case SF_LOW_Q:
    psBiquadEq->m_pfLowQ = DataLocation;
    break;
  case SF_P1_F:
    psBiquadEq->m_pfP1F = DataLocation;
    break;
  case SF_P1_G:
    psBiquadEq->m_pfP1G = DataLocation;
    break;
  case SF_P1_Q:
    psBiquadEq->m_pfP1Q = DataLocation;
    break;
  case SF_P2_F:
    psBiquadEq->m_pfP2F = DataLocation;
    break;
  case SF_P2_G:
    psBiquadEq->m_pfP2G = DataLocation;
    break;
  case SF_P2_Q:
    psBiquadEq->m_pfP2Q = DataLocation;
    break;
  case SF_P3_F:
    psBiquadEq->m_pfP3F = DataLocation;
    break;
  case SF_P3_G:
    psBiquadEq->m_pfP3G = DataLocation;
    break;
  case SF_P3_Q:
    psBiquadEq->m_pfP3Q = DataLocation;
    break;
  case SF_HIGH_F:
    psBiquadEq->m_pfHighF = DataLocation;
    break;
  case SF_HIGH_G:
    psBiquadEq->m_pfHighG = DataLocation;
    break;
  case SF_HIGH_Q:
    psBiquadEq->m_pfHighQ = DataLocation;
    break;
  case SF_GAIN:
    psBiquadEq->m_pfGain = DataLocation;
    break;
  }
}

/*****************************************************************************/

/* Run the filter algorithm for a block of SampleCount samples. */
void runBiquadEq(LADSPA_Handle Instance, unsigned long SampleCount) {

  LADSPA_Data * pfInput;
  LADSPA_Data * pfOutput;
  BiquadEq * psBiquadEq;
  BiquadCoeffs coeffs_LOW, coeffs_P1, coeffs_P2, coeffs_P3, coeffs_HIGH;
  unsigned long lSampleIndex;
  float fGainFactor;
  float xn, yn; // xn/yn holds currently processed input/output sample.
  float xnm1, xnm2, ynm1, ynm2; // holds previously processed input/output samples.
  // get BiquadEq Instance
  psBiquadEq = (BiquadEq *)Instance;
  // get input and output buffers
  pfInput = psBiquadEq->m_pfInput;
  pfOutput = psBiquadEq->m_pfOutput;
  // calculate coeffs and gain factor
  coeffs_LOW = calcCoeffsLowShelf(*(psBiquadEq->m_pfLowF),
	  		          *(psBiquadEq->m_pfLowG),
			          *(psBiquadEq->m_pfLowQ),
			          psBiquadEq->m_fSampleRate);
  coeffs_P1 = calcCoeffsPeaking(*(psBiquadEq->m_pfP1F),
	  		        *(psBiquadEq->m_pfP1G),
			        *(psBiquadEq->m_pfP1Q),
			        psBiquadEq->m_fSampleRate);
  coeffs_P2 = calcCoeffsPeaking(*(psBiquadEq->m_pfP2F),
	  		        *(psBiquadEq->m_pfP2G),
			        *(psBiquadEq->m_pfP2Q),
			        psBiquadEq->m_fSampleRate);
  coeffs_P3 = calcCoeffsPeaking(*(psBiquadEq->m_pfP3F),
	  		        *(psBiquadEq->m_pfP3G),
			        *(psBiquadEq->m_pfP3Q),
			        psBiquadEq->m_fSampleRate);
  coeffs_HIGH = calcCoeffsHighShelf(*(psBiquadEq->m_pfHighF),
	  	  	            *(psBiquadEq->m_pfHighG),
			            *(psBiquadEq->m_pfHighQ),
			            psBiquadEq->m_fSampleRate);
  fGainFactor = dbToGainFactor(*(psBiquadEq->m_pfGain));
  // FILTER PROCESSING FOR LOW SHELF EQ ///////////////////////////////////////
  // get preciosuly processed samples
  xnm1 = psBiquadEq->m_fxnm1_LOW;
  xnm2 = psBiquadEq->m_fxnm2_LOW;
  ynm1 = psBiquadEq->m_fynm1_LOW;
  ynm2 = psBiquadEq->m_fynm1_LOW;
  // apply biquad calculation to buffers
  for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
    xn = *(pfInput++);
    yn = (coeffs_LOW.b0 * xn   + coeffs_LOW.b1 * xnm1 + coeffs_LOW.b2 * xnm2 - 
	  coeffs_LOW.a1 * ynm1 - coeffs_LOW.a2 * ynm2) / coeffs_LOW.a0;
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;
    *(pfOutput++) = yn;
  }
  // store preciously calculated samples in BiQuad instance for later
  psBiquadEq->m_fxnm1_P1 = xnm1;
  psBiquadEq->m_fxnm2_P1 = xnm2;
  psBiquadEq->m_fynm1_P1 = ynm1;
  psBiquadEq->m_fynm1_P1 = ynm2;
  // FILTER PROCESSING FOR PEAKING EQ 1 ///////////////////////////////////////
  // reset output buffer pointer
  pfOutput = psBiquadEq->m_pfOutput;
  // get preciosuly processed samples
  xnm1 = psBiquadEq->m_fxnm1_P1;
  xnm2 = psBiquadEq->m_fxnm2_P1;
  ynm1 = psBiquadEq->m_fynm1_P1;
  ynm2 = psBiquadEq->m_fynm1_P1;
  // apply biquad calculation to buffers
  for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
    xn = *(pfOutput);
    yn = (coeffs_P1.b0 * xn   + coeffs_P1.b1 * xnm1 + coeffs_P1.b2 * xnm2 - 
	  coeffs_P1.a1 * ynm1 - coeffs_P1.a2 * ynm2) / coeffs_P1.a0;
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;
    *(pfOutput++) = yn;
  }
  // store preciously calculated samples in BiQuad instance for later
  psBiquadEq->m_fxnm1_P1 = xnm1;
  psBiquadEq->m_fxnm2_P1 = xnm2;
  psBiquadEq->m_fynm1_P1 = ynm1;
  psBiquadEq->m_fynm1_P1 = ynm2;
  // FILTER PROCESSING FOR PEAKING EQ 2 ///////////////////////////////////////
  // reset output buffer pointer
  pfOutput = psBiquadEq->m_pfOutput;
  // get preciosuly processed samples
  xnm1 = psBiquadEq->m_fxnm1_P2;
  xnm2 = psBiquadEq->m_fxnm2_P2;
  ynm1 = psBiquadEq->m_fynm1_P2;
  ynm2 = psBiquadEq->m_fynm1_P2;
  // apply biquad calculation to buffers
  for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
    xn = *(pfOutput);
    yn = (coeffs_P2.b0 * xn   + coeffs_P2.b1 * xnm1 + coeffs_P2.b2 * xnm2 - 
	  coeffs_P2.a1 * ynm1 - coeffs_P2.a2 * ynm2) / coeffs_P2.a0;
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;
    *(pfOutput++) = yn;
  }
  // store preciously calculated samples in BiQuad instance for later
  psBiquadEq->m_fxnm1_P2 = xnm1;
  psBiquadEq->m_fxnm2_P2 = xnm2;
  psBiquadEq->m_fynm1_P2 = ynm1;
  psBiquadEq->m_fynm1_P2 = ynm2;
  // FILTER PROCESSING FOR PEAKING EQ 3 ///////////////////////////////////////
  // reset output buffer pointer
  pfOutput = psBiquadEq->m_pfOutput;
  // get preciosuly processed samples
  xnm1 = psBiquadEq->m_fxnm1_P3;
  xnm2 = psBiquadEq->m_fxnm2_P3;
  ynm1 = psBiquadEq->m_fynm1_P3;
  ynm2 = psBiquadEq->m_fynm1_P3;
  // apply biquad calculation to buffers
  for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
    xn = *(pfOutput);
    yn = (coeffs_P3.b0 * xn   + coeffs_P3.b1 * xnm1 + coeffs_P3.b2 * xnm2 - 
	  coeffs_P3.a1 * ynm1 - coeffs_P3.a2 * ynm2) / coeffs_P3.a0;
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;
    *(pfOutput++) = yn;
  }
  // store preciously calculated samples in BiQuad instance for later
  psBiquadEq->m_fxnm1_P3 = xnm1;
  psBiquadEq->m_fxnm2_P3 = xnm2;
  psBiquadEq->m_fynm1_P3 = ynm1;
  psBiquadEq->m_fynm1_P3 = ynm2;
  // FILTER PROCESSING FOR HIGH SHELF EQ //////////////////////////////////////
  // reset output buffer pointer
  pfOutput = psBiquadEq->m_pfOutput;
  // get preciosuly processed samples
  xnm1 = psBiquadEq->m_fxnm1_HIGH;
  xnm2 = psBiquadEq->m_fxnm2_HIGH;
  ynm1 = psBiquadEq->m_fynm1_HIGH;
  ynm2 = psBiquadEq->m_fynm1_HIGH;
  // apply biquad calculation to buffers
  for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
    xn = *(pfOutput);
    yn = (coeffs_HIGH.b0 * xn   + coeffs_HIGH.b1 * xnm1 + coeffs_HIGH.b2 * xnm2 - 
	  coeffs_HIGH.a1 * ynm1 - coeffs_HIGH.a2 * ynm2) / coeffs_HIGH.a0;
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;
    *(pfOutput++) = yn * fGainFactor;
  }
  // store preciously calculated samples in BiQuad instance for later
  psBiquadEq->m_fxnm1_HIGH = xnm1;
  psBiquadEq->m_fxnm2_HIGH = xnm2;
  psBiquadEq->m_fynm1_HIGH = ynm1;
  psBiquadEq->m_fynm1_HIGH = ynm2;
}

/*****************************************************************************/

/* Throw away a filter instance. */
void cleanupBiquadEq(LADSPA_Handle Instance) {
  free(Instance);
}

/*****************************************************************************/

LADSPA_Descriptor * g_psBiquadEqDescriptor = NULL;

/*****************************************************************************/

/* _init() is called automatically when the plugin library is first loaded. */
void 
_init() {

  char ** pcPortNames;
  LADSPA_PortDescriptor * piPortDescriptors;
  LADSPA_PortRangeHint * psPortRangeHints;
  
  g_psBiquadEqDescriptor
    = (LADSPA_Descriptor *)malloc(sizeof(LADSPA_Descriptor));

  if (g_psBiquadEqDescriptor != NULL) {
  
    g_psBiquadEqDescriptor->UniqueID
      = 5;
    g_psBiquadEqDescriptor->Label
      = strdup("biquad-eq-t5");
    g_psBiquadEqDescriptor->Properties
      = LADSPA_PROPERTY_HARD_RT_CAPABLE;
    g_psBiquadEqDescriptor->Name 
      = strdup("3-Band Biquad-EQ with Shelves");
    g_psBiquadEqDescriptor->Maker
      = strdup("Juergen Herrmann (t-5@gmx.de)");
    g_psBiquadEqDescriptor->Copyright
      = strdup("None");
    g_psBiquadEqDescriptor->PortCount
      = PORTCOUNT;
    piPortDescriptors
      = (LADSPA_PortDescriptor *)calloc(PORTCOUNT, sizeof(LADSPA_PortDescriptor));
    g_psBiquadEqDescriptor->PortDescriptors
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
    pcPortNames
      = (char **)calloc(PORTCOUNT, sizeof(char *));
    g_psBiquadEqDescriptor->PortNames 
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
    psPortRangeHints = ((LADSPA_PortRangeHint *)
			calloc(PORTCOUNT, sizeof(LADSPA_PortRangeHint)));
    g_psBiquadEqDescriptor->PortRangeHints
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
    // In- and Outpus -------------------------------------------------- */
    psPortRangeHints[SF_INPUT].HintDescriptor
      = 0;
    psPortRangeHints[SF_OUTPUT].HintDescriptor
      = 0;
    g_psBiquadEqDescriptor->instantiate 
      = instantiateBiquadEq;
    g_psBiquadEqDescriptor->connect_port 
      = connectPortToBiquadEq;
    g_psBiquadEqDescriptor->activate
      = activateBiquadEq;
    g_psBiquadEqDescriptor->run
      = runBiquadEq;
    g_psBiquadEqDescriptor->run_adding
      = NULL;
    g_psBiquadEqDescriptor->set_run_adding_gain
      = NULL;
    g_psBiquadEqDescriptor->deactivate
      = NULL;
    g_psBiquadEqDescriptor->cleanup
      = cleanupBiquadEq;
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
  deleteDescriptor(g_psBiquadEqDescriptor);
}

/*****************************************************************************/

/* Return a descriptor of the requested plugin type. There is one
   plugin types available in this library. */
const LADSPA_Descriptor * ladspa_descriptor(unsigned long Index) {
  /* Return the requested descriptor or null if the index is out of
     range. */
  switch (Index) {
  case 0:
    return g_psBiquadEqDescriptor;
  default:
    return NULL;
  }
}

/*****************************************************************************/

/* EOF */
