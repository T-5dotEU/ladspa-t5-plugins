/* helpers.h

   Free software by Juergen Herrmann, t-5@t-5.eu. Do with it, whatever you
   want. No warranty. None, whatsoever. Also see license.txt .

*/

/*****************************************************************************/

#include <math.h>
#include <ladspa.h>

/* biquad coefficients */
typedef struct {

  float a1;
  float a2;
  float b0;
  float b1;
  float b2;

} BiquadCoeffs;

/* s/ns return value */
typedef struct {
    long s;
    long ns;
    LADSPA_Data * mmap;
} TimeMmapStruct;

/* Helpers... ****************************************************************/

float dbToGainFactor(float db) {
  return pow(10, db/20.0);
}

void cleanupMmapFile(char pluginname[], float mmapfname, long s, long ns);
void cleanupMmapFile(char pluginname[], float mmapfname, long s, long ns) {
    char name[255];
    sprintf(name,
            "/dev/shm/t5_%s_%u_%011lu.%09lu",
            pluginname,
            (int)round(mmapfname),
            s,
            ns);
    int ret = remove(name);
    if (ret != 0) {
      printf("Error removing mmap file %s\n", name);
    }
}

TimeMmapStruct setupMmapFile(char pluginname[], float mmapfname, int portcount);
TimeMmapStruct setupMmapFile(char pluginname[], float mmapfname, int portcount) {
    TimeMmapStruct ret;
    char name[255];
    long ns;
    time_t s;
    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);
    s = spec.tv_sec;
    ns = spec.tv_nsec;
    sprintf(name,
            "/dev/shm/t5_%s_%u_%011lu.%09lu",
            pluginname,
            (int)round(mmapfname),
            s,
            ns);
    int fd = open(name, O_RDWR | O_CREAT, 0600);
    int fret = ftruncate(fd, (portcount-1) * sizeof(LADSPA_Data));
    if (fret != 0) {
        printf("ERROR: could not truncate mmaped file %s", name);
    }
    ret.mmap = (LADSPA_Data *) mmap(NULL, portcount * sizeof(LADSPA_Data),
                PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    ret.s = s;
    ret.ns = ns;
    return ret;
}