/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiSampleSeg.h"
#include <stdio.h>
#include <stdlib.h>

int
multiSampleSeg(
  struct Profile **samples,
  int n_samples,
  int n_bins,
  int *optimal_start_end // array of length 2.
  ) {
  int sample_i, coverage_i, min_chromEnd, max_chromStart,
    chromStart, chromEnd;
  struct Profile *profile;
  profile = samples[0];
  min_chromEnd = get_max_chromEnd(profile);
  max_chromStart = get_min_chromStart(profile);
  for(sample_i=1; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    chromStart = get_min_chromStart(profile);
    if(max_chromStart < chromStart){
      max_chromStart = chromStart;
    }
    chromEnd = get_max_chromEnd(profile);
    if(chromEnd < min_chromEnd){
      min_chromEnd = chromEnd;
    }
    //printf("sample_i=%d profile_size=%d\n", sample_i, profile->n_entries);
    //for(coverage_i=0; coverage_i < profile->n_entries; coverage_i++){
    //profile->chromStart[0];
    //}
  }
  int bases = min_chromEnd - max_chromStart;
  if(bases <= n_bins){
    return ERROR_FEWER_BASES_THAN_BINS;
  }
  int bin_i, bin_chromStart, bin_chromEnd;
  int *bin_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    coverage_i=0;
    for(bin_i=0; bin_i < n_bins; bin_i++){
      bin_chromStart = bin_i * bases/n_bins + max_chromStart;
      bin_chromEnd = (bin_i+1)*bases/n_bins + max_chromStart;
      printf("[%3d] (%10d, %10d]\n", bin_i, bin_chromStart, bin_chromEnd);
    }
  }
  free(bin_mat);
  optimal_start_end[0] = max_chromStart;
  optimal_start_end[1] = min_chromEnd;
  return 0;
}

