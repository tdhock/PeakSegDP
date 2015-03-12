/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiSampleSeg.h"
#include <stdio.h>

int
multiSampleSeg(
  struct Profile **samples,
  int n_samples,
  int n_bins,
  int *optimal_start_end // array of length 2.
  ) {
  int sample_i, coverage_i;
  struct Profile *profile;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    printf("sample_i=%d profile_size=%d\n", sample_i, profile->n_entries);
    /* for(coverage_i=0; coverage_i < n_entries; coverage_i++){ */
    /*   profile->chromStart[0]; */
    /* } */
  }
  optimal_start_end[0] = 5;
  optimal_start_end[1] = 10;
  return 0;
}

