/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiSampleSeg.h"
#include "binSum.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
  int bases_per_bin = bases/n_bins;
  if(bases <= n_bins){
    return ERROR_FEWER_BASES_THAN_BINS;
  }
  int *count_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int *cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  double *mean1_mat = (double*) malloc(n_bins * n_samples * sizeof(double));
  double *loss1_mat = (double*) malloc(n_bins * n_samples * sizeof(double));
  int status;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    count_mat + n_bins*sample_i,
		    bases_per_bin, n_bins, max_chromStart);
  }//for sample_i
  int bin_i, offset;
  int *count_vec, *cumsum_vec, cumsum_value;
  double *mean_vec, *loss_vec, mean_value, loss_value;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    cumsum_value = 0;
    offset = n_bins * sample_i;
    count_vec = count_mat + offset;
    cumsum_vec = cumsum_mat + offset;
    mean_vec = mean1_mat + offset;
    loss_vec = loss1_mat + offset;
    for(bin_i=0; bin_i < n_bins; bin_i++){
      cumsum_value += count_vec[bin_i];
      cumsum_vec[bin_i] = cumsum_value;
      mean_value = ((double) cumsum_value) / ((double)bin_i+1);
      mean_vec[bin_i] = mean_value;
      if(cumsum_value == 0){
	loss_vec[bin_i] = 0.0;
      }else{
	loss_vec[bin_i] = cumsum_value * (1-log(mean_value));
      }
      printf("[%3d,%3d]=%d %f %f\n", sample_i, bin_i, 
	     cumsum_value, mean_value, 
	     ((double)count_vec[bin_i])/((double)bases_per_bin));
    }
  }
  free(count_mat);
  free(cumsum_mat);
  free(mean1_mat);
  free(loss1_mat);
  optimal_start_end[0] = max_chromStart;
  optimal_start_end[1] = min_chromEnd;
  return 0;
}

