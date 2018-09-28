/* t5_lr4_lowpass.c

   Free software by Juergen Herrmann, t-5@t-5.eu. Do with it, whatever you
   want. No warranty. None, whatsoever. Also see license.txt .

   This LADSPA plugin provides a Linkwitz-Riley 24dB/octave low pass filter.

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
#include "lr4.h"

/* Helpers... ****************************************************************/

BiquadCoeffs calcCoeffsLr4Lowpass(float f, float samplerate) {
    BiquadCoeffs coeffs;
    float w0 = 2 * M_PI * f / samplerate;
    float alpha = sin(w0) / 2 / 0.7071067811865476; // Butterworth characteristic, Q = 0.707...
    float cs = cos(w0);
    float norm = 1 / (1 + alpha);
    coeffs.b0 = (1 - cs) / 2 * norm;
    coeffs.b1 = (1 - cs) * norm;
    coeffs.b2 = coeffs.b0;
    coeffs.a1 = -2 * cs * norm;
    coeffs.a2 = (1 - alpha) * norm;
    return coeffs;
}

void setupMmapFileForLr4Lowpass(Lr4LowHighPass * psInstance) {
    // setup shared memory area to enable parametrization at runtime from external processes
    TimeMmapStruct ret;
    ret = setupMmapFile("Lr4Lowpass", *(psInstance->m_pfMmapFname), PORTCOUNT); 
    psInstance->m_mmapArea = ret.mmap;
    psInstance->m_created_s = ret.s;
    psInstance->m_created_ns = ret.ns;
}

/*****************************************************************************/

/* Run the filter algorithm for a block of SampleCount samples. */
void runLr4Lowpass(LADSPA_Handle Instance, unsigned long SampleCount) {
    BiquadCoeffs coeffs;
    Lr4LowHighPass * psInstance;
    psInstance = (Lr4LowHighPass *)Instance;
    if (psInstance->m_mmapArea == NULL && *(psInstance->m_pfMmapFname) != 0.0) {
        setupMmapFileForLr4Lowpass(psInstance);
    }
    coeffs = calcCoeffsLr4Lowpass(*(psInstance->m_pfF), psInstance->m_fSampleRate);
    runLr4LowHighPass(Instance, SampleCount, coeffs);}

/*****************************************************************************/

/* Throw away a Lr4LowHighPass instance. */
void cleanupLr4Lowpass(LADSPA_Handle Instance) {
  Lr4LowHighPass * psInstance;
  psInstance = (Lr4LowHighPass *)Instance;
  if (psInstance->m_mmapArea != NULL) {
    cleanupMmapFile("Lr4Lowpass",
                    *(psInstance->m_pfMmapFname),
                    psInstance->m_created_s,
                    psInstance->m_created_ns);
    free(Instance);
  }
}

/*****************************************************************************/

LADSPA_Descriptor * g_psLr4LowpassInstanceDescriptor = NULL;

/*****************************************************************************/

/* _init() is called automatically when the plugin library is first loaded. */
void _init() {

  char ** pcPortNames;
  LADSPA_PortDescriptor * piPortDescriptors;
  LADSPA_PortRangeHint * psPortRangeHints;
  
  g_psLr4LowpassInstanceDescriptor
    = (LADSPA_Descriptor *)malloc(sizeof(LADSPA_Descriptor));

  if (g_psLr4LowpassInstanceDescriptor != NULL) {
  
    g_psLr4LowpassInstanceDescriptor->UniqueID
      = 5542;
    g_psLr4LowpassInstanceDescriptor->Label
      = strdup("lr4_lowpass");
    g_psLr4LowpassInstanceDescriptor->Properties
      = LADSPA_PROPERTY_HARD_RT_CAPABLE;
    g_psLr4LowpassInstanceDescriptor->Name 
      = strdup("T5's LR-4 Low Pass");
    g_psLr4LowpassInstanceDescriptor->Maker
      = strdup("Juergen Herrmann (t-5@t-5.eu)");
    g_psLr4LowpassInstanceDescriptor->Copyright
      = strdup("3-clause BSD licence");
    g_psLr4LowpassInstanceDescriptor->PortCount
      = PORTCOUNT;
    piPortDescriptors
      = (LADSPA_PortDescriptor *)calloc(PORTCOUNT, sizeof(LADSPA_PortDescriptor));
    g_psLr4LowpassInstanceDescriptor->PortDescriptors
      = (const LADSPA_PortDescriptor *)piPortDescriptors;
    piPortDescriptors[SF_INPUT]
      = LADSPA_PORT_INPUT | LADSPA_PORT_AUDIO;
    piPortDescriptors[SF_OUTPUT]
      = LADSPA_PORT_OUTPUT | LADSPA_PORT_AUDIO;
    piPortDescriptors[SF_F]
      = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
    piPortDescriptors[SF_GAIN]
      = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
    piPortDescriptors[SF_MMAPFNAME]
      = LADSPA_PORT_INPUT | LADSPA_PORT_CONTROL;
    pcPortNames
      = (char **)calloc(PORTCOUNT, sizeof(char *));
    g_psLr4LowpassInstanceDescriptor->PortNames 
      = (const char **)pcPortNames;
    pcPortNames[SF_INPUT]
      = strdup("Input");
    pcPortNames[SF_OUTPUT]
      = strdup("Output");
    pcPortNames[SF_F]
      = strdup("Cutoff Frequency [Hz]");
    pcPortNames[SF_GAIN]
      = strdup("Overall Gain [dB]");
    pcPortNames[SF_MMAPFNAME]
      = strdup("MMAP-Filename-Part");
    psPortRangeHints = ((LADSPA_PortRangeHint *)
      calloc(PORTCOUNT, sizeof(LADSPA_PortRangeHint)));
    g_psLr4LowpassInstanceDescriptor->PortRangeHints
      = (const LADSPA_PortRangeHint *)psPortRangeHints;
    psPortRangeHints[SF_F].HintDescriptor
      = (LADSPA_HINT_BOUNDED_BELOW 
       | LADSPA_HINT_BOUNDED_ABOVE
       | LADSPA_HINT_SAMPLE_RATE
       | LADSPA_HINT_LOGARITHMIC
       | LADSPA_HINT_DEFAULT_440);
    psPortRangeHints[SF_F].LowerBound 
      = 0;
    psPortRangeHints[SF_F].UpperBound
      = 0.5;
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
    g_psLr4LowpassInstanceDescriptor->instantiate 
      = instantiateLr4LowHighPass;
    g_psLr4LowpassInstanceDescriptor->connect_port 
      = connectPortToLr4LowHighPass;
    g_psLr4LowpassInstanceDescriptor->activate
      = activateLr4LowHighPass;
    g_psLr4LowpassInstanceDescriptor->run
      = runLr4Lowpass;
    g_psLr4LowpassInstanceDescriptor->run_adding
      = NULL;
    g_psLr4LowpassInstanceDescriptor->set_run_adding_gain
      = NULL;
    g_psLr4LowpassInstanceDescriptor->deactivate
      = NULL;
    g_psLr4LowpassInstanceDescriptor->cleanup
      = cleanupLr4Lowpass;
  }
}
  
/*****************************************************************************/

/* _fini() is called automatically when the library is unloaded. */
void _fini() {
  deleteLr4LowHighPassDescriptor(g_psLr4LowpassInstanceDescriptor);
}

/*****************************************************************************/

/* Return a descriptor of the requested plugin types. */
const LADSPA_Descriptor * ladspa_descriptor(unsigned long Index) {
  /* Return the requested descriptor or null if the index is out of range. */
  switch (Index) {
  case 0:
    return g_psLr4LowpassInstanceDescriptor;
  default:
    return NULL;
  }
}

/*****************************************************************************/

/* EOF */
