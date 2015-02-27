#include "binSum.h"
#include <R.h>

void binSum_interface(
  int *profile_chromStart,
  int *profile_chromEnd,
  int *profile_coverage,
  int *n_profiles,
  int *bin_total,
  int *bin_size,
  int *n_bins,
  int *bin_chromStart){
  int status;
  status = binSum(profile_chromStart, 
		  profile_chromEnd,
		  profile_coverage,
		  *n_profiles,
		  bin_total,
		  *bin_size,
		  *n_bins,
		  *bin_chromStart);
  if(status != 0){
    error("error code %d", status);
  }
}
