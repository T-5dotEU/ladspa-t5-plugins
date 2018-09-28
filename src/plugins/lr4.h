//#include "helpers.h"

/*****************************************************************************/

#define SF_INPUT       0
#define SF_OUTPUT      1
#define SF_F           2
#define SF_GAIN        3
#define SF_MMAPFNAME   4
#define PORTCOUNT      5

/*****************************************************************************/

/* Instance data for the Lr4(Low|High)Pass filter */
typedef struct {

    long m_created_ns;
    time_t m_created_s;

    LADSPA_Data * m_mmapArea;

    LADSPA_Data m_fSampleRate;
    // previous input samples of biquad filters
    LADSPA_Data m_fxnm1_pass1;
    LADSPA_Data m_fxnm2_pass1;
    LADSPA_Data m_fxnm1_pass2;
    LADSPA_Data m_fxnm2_pass2;
    // previous output samples of biquad filters
    LADSPA_Data m_fynm1_pass1;
    LADSPA_Data m_fynm2_pass1;
    LADSPA_Data m_fynm1_pass2;
    LADSPA_Data m_fynm2_pass2;
    // port pointers
    LADSPA_Data * m_pfInput;
    LADSPA_Data * m_pfOutput;
    LADSPA_Data * m_pfF;
    LADSPA_Data * m_pfGain;
    LADSPA_Data * m_pfMmapFname;

} Lr4LowHighPass;

/* Construct a new plugin instance. */
LADSPA_Handle instantiateLr4LowHighPass(const LADSPA_Descriptor * Descriptor,
                                        unsigned long SampleRate) {
    Lr4LowHighPass * psInstance;
    psInstance = (Lr4LowHighPass *)malloc(sizeof(Lr4LowHighPass));
    if (psInstance) {
        psInstance->m_fSampleRate = (LADSPA_Data)SampleRate;
        psInstance->m_mmapArea = NULL;
    }  
    return psInstance;
}

/* Initialise and activate a plugin instance. */
void activateLr4LowHighPass(LADSPA_Handle Instance);
void activateLr4LowHighPass(LADSPA_Handle Instance) {
    Lr4LowHighPass * psInstance;
    psInstance = (Lr4LowHighPass *)Instance;
    psInstance->m_fxnm1_pass1 = 0;
    psInstance->m_fxnm2_pass1 = 0;
    psInstance->m_fynm1_pass1 = 0;
    psInstance->m_fynm2_pass1 = 0;
    psInstance->m_fxnm1_pass2 = 0;
    psInstance->m_fxnm2_pass2 = 0;
    psInstance->m_fynm1_pass2 = 0;
    psInstance->m_fynm2_pass2 = 0;
}


/* Connect a port to a data location.  */
void connectPortToLr4LowHighPass(LADSPA_Handle Instance,
                                 unsigned long Port,
                                 LADSPA_Data * DataLocation);
void connectPortToLr4LowHighPass(LADSPA_Handle Instance,
                                 unsigned long Port,
                                 LADSPA_Data * DataLocation) {
  
    Lr4LowHighPass * psInstance;
    psInstance = (Lr4LowHighPass *)Instance;

    switch (Port) {
    case SF_INPUT:
        psInstance->m_pfInput = DataLocation;
        break;
    case SF_OUTPUT:
        psInstance->m_pfOutput = DataLocation;
        break;
    case SF_F:
        psInstance->m_pfF = DataLocation;
        break;
    case SF_GAIN:
        psInstance->m_pfGain = DataLocation;
        break;
    case SF_MMAPFNAME:
        psInstance->m_pfMmapFname = DataLocation;
        break;
    }
}

void deleteLr4LowHighPassDescriptor(LADSPA_Descriptor * psDescriptor);
void deleteLr4LowHighPassDescriptor(LADSPA_Descriptor * psDescriptor) {
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

/* Run the filter algorithm for a block of SampleCount samples. */
void runLr4LowHighPass(LADSPA_Handle Instance, unsigned long SampleCount, BiquadCoeffs coeffs);
void runLr4LowHighPass(LADSPA_Handle Instance, unsigned long SampleCount, BiquadCoeffs coeffs) {

  LADSPA_Data * pfInput;
  LADSPA_Data * pfOutput;
  Lr4LowHighPass * psInstance;
  unsigned long lSampleIndex;
  LADSPA_Data unchanged = 0.0;
  LADSPA_Data changed;
  LADSPA_Data * mmptr;
  float fGainFactor;
  float xn, yn; // xn/yn holds currently processed input/output samples.
  float xnm1, xnm2, ynm1, ynm2; // holds previously processed input/output samples.
  // get Lr4LowHighPass Instance
  psInstance = (Lr4LowHighPass *)Instance;
  // get input and output buffers
  pfInput = psInstance->m_pfInput;
  pfOutput = psInstance->m_pfOutput;
  // memcpy parameters over from mmapped area
  mmptr = psInstance->m_mmapArea;
  if (mmptr != NULL) {
    memcpy(&changed, mmptr, sizeof(LADSPA_Data));
    if (changed != 0.0) {
      mmptr += 1;
      memcpy(psInstance->m_pfF, mmptr, sizeof(LADSPA_Data));
      mmptr += 1;
      memcpy(psInstance->m_pfGain, mmptr, sizeof(LADSPA_Data));
    }
    // reset changed flag to re-enable parametrization by control inputs
    memcpy(psInstance->m_mmapArea, &unchanged, sizeof(LADSPA_Data));
  }
  fGainFactor = dbToGainFactor(*(psInstance->m_pfGain));
  // get preciosuly processed samples
  xnm1 = psInstance->m_fxnm1_pass1;
  xnm2 = psInstance->m_fxnm2_pass1;
  ynm1 = psInstance->m_fynm1_pass1;
  ynm2 = psInstance->m_fynm2_pass1;
  // FILTER PROCESSING First Pass /////////////////////////////////////////////
  for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
    xn = *(pfInput++);
    yn = (coeffs.b0 * xn   + coeffs.b1 * xnm1 + coeffs.b2 * xnm2 - 
          coeffs.a1 * ynm1 - coeffs.a2 * ynm2);
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;
    *(pfOutput++) = yn;
  }
  // store preciously calculated samples in BiQuad instance for later
  psInstance->m_fxnm1_pass1 = xnm1;
  psInstance->m_fxnm2_pass1 = xnm2;
  psInstance->m_fynm1_pass1 = ynm1;
  psInstance->m_fynm2_pass1 = ynm2;
  // reset output buffer pointer
  pfOutput = psInstance->m_pfOutput;
  // get preciosuly processed samples
  xnm1 = psInstance->m_fxnm1_pass2;
  xnm2 = psInstance->m_fxnm2_pass2;
  ynm1 = psInstance->m_fynm1_pass2;
  ynm2 = psInstance->m_fynm2_pass2;
  // FILTER PROCESSING Second Pass ////////////////////////////////////////////
  for (lSampleIndex = 0; lSampleIndex < SampleCount; lSampleIndex++) {
    xn = *(pfOutput);
    yn = (coeffs.b0 * xn   + coeffs.b1 * xnm1 + coeffs.b2 * xnm2 - 
          coeffs.a1 * ynm1 - coeffs.a2 * ynm2);
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;
    *(pfOutput++) = yn * fGainFactor;
  }
  // store preciously calculated samples in BiQuad instance for later
  psInstance->m_fxnm1_pass2 = xnm1;
  psInstance->m_fxnm2_pass2 = xnm2;
  psInstance->m_fynm1_pass2 = ynm1;
  psInstance->m_fynm2_pass2 = ynm2;
}

