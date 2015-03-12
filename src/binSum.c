#include "binSum.h"

int binSum
(int *profile_chromStart, 
 int *profile_chromEnd, 
 int *profile_coverage, 
 int n_profiles,
 int *bin_total, 
 int bin_size,
 int n_bins, 
 int bin_chromStart){
  int profile_i = 0, bin_i = 0;
  for(bin_i = 0; bin_i < n_bins; bin_i++){
    bin_total[bin_i] = 0;
  }
  // bin_chromStart gives the base before the first position that we
  // want to count, for example 1000 means we want to start counting at
  // 1001. so we should ignore profile entries of (0, 10], (0, 1000], 
  // (999, 1000], but start counting (0, 1001], (1000, 1001], (1000, 1002].
  bin_i = 0;
  while(profile_chromEnd[profile_i] <= bin_chromStart){
    profile_i ++;
  }
  int count_until, bases, bin_add, profile_add;
  int begin_count_after = bin_chromStart;
  int bin_end = bin_chromStart + bin_size;
  while(bin_i < n_bins && profile_i < n_profiles){
    // at this point there are two cases.
    if(bin_end <= profile_chromEnd[profile_i]){
      // 1. the profile segment continues to the end of this bin,
      //    so add profile_coverage * bin_size to this bin total.
      // -profile----]
      //             (-----------
      //          bin]
      //         bin]
      //  bin]
      count_until = bin_end;
      if(bin_end == profile_chromEnd[profile_i]){
	profile_add = 1; // done adding from this profile segment.
      }else{
	profile_add = 0; // not done adding from this profile segment.
      }
      bin_add = 1;
    }else{      
      // 2. the profile segment ends before this bin ends.
      // -profile----]
      //             (-----------
      //           bin]
      //             bin]
      count_until = profile_chromEnd[profile_i];
      profile_add = 1; // done adding from this profile segment.
      bin_add = 0; // not done adding to this bin total.
    }
    bases = count_until - begin_count_after;
    bin_total[bin_i] += profile_coverage[profile_i] * bases;
    // setup next iteration.
    begin_count_after = count_until;
    profile_i += profile_add;
    if(bin_add){
      bin_i++;
      bin_end += bin_size;
    }
  }
  while(bin_i < n_bins){
    //printf("bin_i=%d bin_total[bin_i]=%d\n", bin_i, bin_total[bin_i]);
    if(bin_total[bin_i] == 0){
      bin_total[bin_i] = -1;
    }
    bin_i++;
  }
  return 0;
}
