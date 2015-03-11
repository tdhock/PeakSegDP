#include "multiSampleSeg.h"

int
multiSampleSeg(
  Profile **samples,
  int n_samples,
  int n_bins,
  int *optimal_start_end // array of length 2.
  ) {
  int sample_i;
  Profile *profile;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
  }
  return 0;
}

