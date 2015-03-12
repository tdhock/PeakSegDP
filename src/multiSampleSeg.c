/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiSampleSeg.h"

int
multiSampleSeg(
  struct Profile **samples,
  int n_samples,
  int n_bins,
  int *optimal_start_end // array of length 2.
  ) {
  int sample_i;
  struct Profile *profile;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
  }
  optimal_start_end[0] = 5;
  optimal_start_end[1] = 10;
  return 0;
}

