/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiSampleSeg.h"
#include "binSum.h"
#include <stdio.h>
#include <stdlib.h>

int get_min_chromStart(struct Profile *profile){
  return profile->chromStart[0];
}

int get_max_chromEnd(struct Profile *profile){
  return profile->chromEnd[profile->n_entries-1];
}

double OptimalPoissonLoss(int cumsum_value, double mean_value){
  if(cumsum_value == 0){
    return 0.0;
  }
  return cumsum_value * (1-log(mean_value));
}

/*
  solver for the DPA recursion: 
  L_{s,t} = min_{t' < t} L_{s-1, t'} + c_(t', t]
  where L is the optimal loss in s segments up to t
  and c is the optimal loss of segment (t', t]

  - PrevSegs_loss_vec is the L_{s-1} vector,
  - sample_cumsum_mat is the matrix of cumsums
    which is used to compute the optimal loss c.
  - first_possible_index = s.
  - last_possible_index = t.
  - LastSeg_FirstIndex = t'.
 */
void get_best_FirstIndex(
  double *PrevSegs_loss_vec, // n_data
  int n_data,
  int *sample_cumsum_mat, // n_data x n_samples
  int n_samples,
  int first_possible_index,
  int last_possible_index,
  int *best_FirstIndex, //output
  double *best_loss //output
  ){
  int sample_i;
  int* cumsum_vec;
  int LastSeg_FirstIndex;
  double candidate_loss, cumsum_value, mean_value;
  *best_loss = INFINITY;
  for(LastSeg_FirstIndex = first_possible_index; 
      LastSeg_FirstIndex <= last_possible_index;
      LastSeg_FirstIndex++){
    // start with previous segment optimal loss.
    candidate_loss = PrevSegs_loss_vec[LastSeg_FirstIndex-1];
    for(sample_i=0; sample_i < n_samples; sample_i++){
      cumsum_vec = sample_cumsum_mat + n_data * sample_i;
      /*  
	  the optimal loss of segment (t', t]
	  can be computed in O(1) time given the vector of cumsums.
      */
      cumsum_value = cumsum_vec[last_possible_index] - 
	cumsum_vec[LastSeg_FirstIndex-1];
      mean_value = ((double)cumsum_value)/
	((double)last_possible_index-LastSeg_FirstIndex+1);
      candidate_loss += OptimalPoissonLoss(cumsum_value, mean_value);
    }
    if(candidate_loss < *best_loss){
      //printf("LastSeg_FirstIndex=%d\n", LastSeg_FirstIndex);
      *best_loss = candidate_loss;
      *best_FirstIndex = LastSeg_FirstIndex;
    }
  }
}

/*
  Fast heuristic solver for the multiple samples with 1 shared peak
  problem. The sub-optimality parameter bin_factor is used to
  down-sample the profiles, which results in a faster segmentation
  problem.
 */
int
multiSampleSegHeuristic(
  struct Profile **samples,
  int n_samples,
  int bin_factor,
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
  if(bases/bin_factor < 4){
    /*
      4 is smallest the number of data points for which the 3-segment
      optimization problem is not trivial.

      If we don't have at least this many data points for the first
      bin step, than we stop with an error.
    */
    return ERROR_BIN_FACTOR_TOO_LARGE;
  }
  int bases_per_bin = 1;
  while(bases/bases_per_bin/bin_factor >= 4){
    bases_per_bin *= bin_factor;
  }
  int n_bins = bases / bases_per_bin;
  printf("n_bins=%d bases_per_bin=%d\n", n_bins, bases_per_bin);
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
      loss_vec[bin_i] = OptimalPoissonLoss(cumsum_value, mean_value);
      /* printf("[%3d,%3d]=%d %f %f\n", sample_i, bin_i,  */
      /* 	     cumsum_value, mean_value,  */
      /* 	     ((double)count_vec[bin_i])/((double)bases_per_bin)); */
    }
  }
  /*
    First step of DPA: compute optimal loss for 1 segment up to data
    point t, for all data points = bins.
  */
  double *seg1_loss_vec = (double*) malloc(n_bins * sizeof(double));
  for(bin_i=0; bin_i < n_bins; bin_i++){
    seg1_loss_vec[bin_i] = 0;
    for(sample_i=0; sample_i < n_samples; sample_i++){
      seg1_loss_vec[bin_i] += sample_loss1_mat[bin_i + n_bins*sample_i];
    }
  }
  /* 
     Second step of DPA: compute optimal loss in 2 segments up to data
     point t, for all data points = bins.
   */
  double *seg12_loss_vec = (double*) malloc(n_bins * sizeof(double));
  int *seg2_first_vec = (int*) malloc(n_bins * sizeof(int));
  int seg2_FirstIndex, seg2_LastIndex, best_FirstIndex;
  double seg1_loss, min_loss, seg2_loss, candidate_loss;
  for(seg2_LastIndex=1; 
      seg2_LastIndex < n_bins; 
      seg2_LastIndex++){
    get_best_FirstIndex(
      seg1_loss_vec,
      n_bins,
      sample_cumsum_mat,
      n_samples,
      1, // first_possible_index
      seg2_LastIndex, // last_possible_index
      seg2_first_vec + seg2_LastIndex,
      seg12_loss_vec + seg2_LastIndex);
  }
  /*
    For the best segmentation in 3 segments up to n_bins-1, the first
    index of the 3rd segment is computed using the following code. No
    need for a for loop on seg3_LastIndex, since we only want to know
    the optimal model which ends at the last data point = bin (we are
    not continuing the DPA past 3 segments).
  */
  int seg3_FirstIndex;
  double best_loss;
  get_best_FirstIndex(
    seg12_loss_vec,
    n_bins,
    sample_cumsum_mat,
    n_samples,
    2, // first_possible_index
    n_bins-1, // last_possible_index
    &seg3_FirstIndex,
    &best_loss);
  seg2_LastIndex = seg3_FirstIndex-1;
  seg2_FirstIndex = seg2_first_vec[seg2_LastIndex];
  printf("[0,%d] [%d,%d] [%d,%d]\n",
	 seg2_FirstIndex-1,
	 seg2_FirstIndex, seg2_LastIndex,
	 seg3_FirstIndex, n_bins-1);
  int peakStart = max_chromStart + bases_per_bin * seg2_FirstIndex;
  int peakEnd = max_chromStart + bases_per_bin * seg3_FirstIndex;

  //Now we zoom in, and search on the left and right bins.
  int *left_count_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int *right_count_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  printf("next bases/bin=%d\n", bases_per_bin / bin_factor);

  //cleanup!
  free(left_count_mat);
  free(right_count_mat);
  free(seg12_loss_vec);
  free(seg1_loss_vec);
  free(seg2_first_vec);
  free(sample_count_mat);
  free(sample_cumsum_mat);
  free(sample_mean1_mat);
  free(sample_loss1_mat);
  optimal_start_end[0] = peakStart;
  optimal_start_end[1] = peakEnd;
  return 0;
}


/*
  Implements base-pair level DPA = Dynamic Programming Algorithm (not
  constrained), in order to recover the most likely Poisson model with
  the same 3 segments (but different mean values) on each of n_samples
  profiles.
 */
int
multiSampleSegOptimal(
  struct Profile **samples,
  int n_samples,
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
  }
  int n_bases = min_chromEnd - max_chromStart;

  // sample_*_mat variables are matrices n_bases x n_samples (in
  // contrast to model_*_mat which are n_bases x n_segments=3).
  int *sample_count_mat = (int*) malloc(n_bases * n_samples * sizeof(int));
  int *sample_cumsum_mat = (int*) malloc(n_bases * n_samples * sizeof(int));
  double *sample_mean1_mat = (double*) malloc(
    n_bases * n_samples * sizeof(double));
  double *sample_loss1_mat = (double*) malloc(
    n_bases * n_samples * sizeof(double));
  int status;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    sample_count_mat + n_bases*sample_i,
		    1, n_bases, max_chromStart);
  }//for sample_i
  int base_i, offset;
  int *count_vec, *cumsum_vec, cumsum_value;
  double *mean_vec, *loss_vec, mean_value, loss_value;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    cumsum_value = 0;
    offset = n_bases * sample_i;
    count_vec = sample_count_mat + offset;
    cumsum_vec = sample_cumsum_mat + offset;
    mean_vec = sample_mean1_mat + offset;
    loss_vec = sample_loss1_mat + offset;
    for(base_i=0; base_i < n_bases; base_i++){
      cumsum_value += count_vec[base_i];
      cumsum_vec[base_i] = cumsum_value;
      mean_value = ((double) cumsum_value) / ((double)base_i+1);
      mean_vec[base_i] = mean_value;
      loss_vec[base_i] = OptimalPoissonLoss(cumsum_value, mean_value);
    }
  }
  int maxSegments = 3;
  double *model_loss_mat = (double*) malloc(
    n_bases * maxSegments * sizeof(double));
  int *model_first_mat = (int*) malloc(
    n_bases * maxSegments * sizeof(int));
  for(base_i=0; base_i < n_bases; base_i++){
    model_loss_mat[base_i] = 0;
    for(sample_i=0; sample_i < n_samples; sample_i++){
      model_loss_mat[base_i] += sample_loss1_mat[base_i + n_bases*sample_i];
    }
  }
  int segment_i, LastSeg_FirstIndex, LastSeg_LastIndex, best_FirstIndex;
  double prev_loss, min_loss, LastSeg_loss, candidate_loss;
  for(segment_i=1; segment_i < maxSegments; segment_i++){
    for(LastSeg_LastIndex=segment_i; 
	LastSeg_LastIndex < n_bases; 
	LastSeg_LastIndex++){
      min_loss = INFINITY;
      for(LastSeg_FirstIndex=segment_i; 
	  LastSeg_FirstIndex <= LastSeg_LastIndex;
	  LastSeg_FirstIndex++){
	prev_loss = model_loss_mat[
	  LastSeg_FirstIndex-1 + n_bases*(segment_i-1)];
	LastSeg_loss = 0.0;
	for(sample_i=0; sample_i < n_samples; sample_i++){
	  cumsum_vec = sample_cumsum_mat + n_bases * sample_i;
	  cumsum_value = cumsum_vec[LastSeg_LastIndex] - 
	    cumsum_vec[LastSeg_FirstIndex-1];
	  mean_value = ((double)cumsum_value)/
	    ((double)LastSeg_LastIndex-LastSeg_FirstIndex+1);
	  LastSeg_loss += OptimalPoissonLoss(cumsum_value, mean_value);
	}
	candidate_loss = LastSeg_loss + prev_loss;
	if(candidate_loss < min_loss){
	  min_loss = candidate_loss;
	  best_FirstIndex = LastSeg_FirstIndex;
	}
      }
      model_first_mat[LastSeg_LastIndex + n_bases*segment_i] = best_FirstIndex;
      model_loss_mat[LastSeg_LastIndex + n_bases*segment_i] = min_loss;
    }
  }
  // For the best segmentation in 3 segments up to n_bases-1, the first
  // index of the 3rd segment is
  int seg3_FirstIndex = model_first_mat[n_bases*3-1];
  int seg2_LastIndex = seg3_FirstIndex-1;
  int seg2_FirstIndex = model_first_mat[n_bases*1+seg2_LastIndex];
  int peakStart = max_chromStart + seg2_FirstIndex;
  int peakEnd = max_chromStart + seg3_FirstIndex;

  free(model_loss_mat);
  free(model_first_mat);
  free(sample_count_mat);
  free(sample_cumsum_mat);
  free(sample_mean1_mat);
  free(sample_loss1_mat);
  optimal_start_end[0] = peakStart;
  optimal_start_end[1] = peakEnd;
  return 0;
}

