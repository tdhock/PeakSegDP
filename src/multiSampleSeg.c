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
  // sample_*_mat variables are matrices n_bins x n_samples (in
  // contrast to model_*_mat which are n_bins x n_segments=3).
  int *sample_count_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int *sample_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  double *sample_mean1_mat = (double*) malloc(
    n_bins * n_samples * sizeof(double));
  double *sample_loss1_mat = (double*) malloc(
    n_bins * n_samples * sizeof(double));
  int status;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    sample_count_mat + n_bins*sample_i,
		    bases_per_bin, n_bins, max_chromStart);
  }//for sample_i
  int bin_i, offset;
  int *count_vec, *cumsum_vec, cumsum_value;
  double *mean_vec, *loss_vec, mean_value, loss_value;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    cumsum_value = 0;
    offset = n_bins * sample_i;
    count_vec = sample_count_mat + offset;
    cumsum_vec = sample_cumsum_mat + offset;
    mean_vec = sample_mean1_mat + offset;
    loss_vec = sample_loss1_mat + offset;
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
      /* printf("[%3d,%3d]=%d %f %f\n", sample_i, bin_i,  */
      /* 	     cumsum_value, mean_value,  */
      /* 	     ((double)count_vec[bin_i])/((double)bases_per_bin)); */
    }
  }
  int maxSegments = 3;
  double *model_loss_mat = (double*) malloc(
    n_bins * maxSegments * sizeof(double));
  for(bin_i=0; bin_i < n_bins; bin_i++){
    model_loss_mat[bin_i] = 0;
    for(sample_i=0; sample_i < n_samples; sample_i++){
      model_loss_mat[bin_i] += sample_loss1_mat[bin_i + n_bins*sample_i];
    }
  }
  int segment_i, LastSeg_FirstIndex, LastSeg_LastIndex, best_FirstIndex;
  double prev_loss, min_loss, LastSeg_loss, candidate_loss;
  for(segment_i=1; segment_i < maxSegments; segment_i++){
    for(LastSeg_LastIndex=segment_i; 
	LastSeg_LastIndex < n_bins; 
	LastSeg_LastIndex++){
      min_loss = INFINITY;
      for(LastSeg_FirstIndex=segment_i; 
	  LastSeg_FirstIndex <= LastSeg_LastIndex;
	  LastSeg_FirstIndex++){
	prev_loss = model_loss_mat[LastSeg_FirstIndex-1 + n_bins*(segment_i-1)];
	LastSeg_loss = 0.0;
	for(sample_i=0; sample_i < n_samples; sample_i++){
	  cumsum_vec = sample_cumsum_mat + n_bins * sample_i;
	  cumsum_value = cumsum_vec[LastSeg_LastIndex] - 
	    cumsum_vec[LastSeg_FirstIndex-1];
	  mean_value = ((double)cumsum_value)/
	    ((double)LastSeg_LastIndex-LastSeg_FirstIndex+1);
	  LastSeg_loss += cumsum_value *(1-log(mean_value));
	}
	candidate_loss = LastSeg_loss + prev_loss;
	if(candidate_loss < min_loss){
	  min_loss = candidate_loss;
	  best_FirstIndex = LastSeg_FirstIndex;
	}
      }
      printf("segment_i=%d First=%d loss=%f\n",
	     segment_i, best_FirstIndex, min_loss);
      model_loss_mat[LastSeg_LastIndex + n_bins*segment_i] = min_loss;
    }
  }
  free(model_loss_mat);
  free(sample_count_mat);
  free(sample_cumsum_mat);
  free(sample_mean1_mat);
  free(sample_loss1_mat);
  optimal_start_end[0] = max_chromStart;
  optimal_start_end[1] = min_chromEnd;
  return 0;
}

