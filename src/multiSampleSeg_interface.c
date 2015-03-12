/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiSampleSeg.h"
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

SEXP
multiSampleSeg_interface(
  SEXP profile_list,
  SEXP n_bins
  ){
  int n_profiles = length(profile_list);
  int profile_i;
  SEXP df;
  struct Profile *profile;
  struct Profile **samples = malloc(n_profiles * sizeof(struct Profile *));
  for(profile_i=0; profile_i < n_profiles; profile_i++){
    df = VECTOR_ELT(profile_list, profile_i);
    profile = malloc(sizeof(struct Profile));
    profile->chromStart = INTEGER(VECTOR_ELT(df, 0));
    profile->chromEnd = INTEGER(VECTOR_ELT(df, 1));
    profile->coverage = INTEGER(VECTOR_ELT(df, 2));
    profile->n_entries = length(df);
    samples[profile_i] = profile;
  }
  int status, n_bins_int;
  n_bins_int = INTEGER(n_bins)[0];
  SEXP optimal;
  PROTECT(optimal = allocVector(INTSXP, 2));
  int *optimal_ptr = INTEGER(optimal);
  status = multiSampleSeg(samples, n_profiles, n_bins_int, optimal_ptr);
  // TODO: if known status codes...
  if(status != 0){
    error("unrecognized error code %d", status);
  }
  for(profile_i=0; profile_i < n_profiles; profile_i++){
    free(samples[profile_i]);
  }
  free(samples);
  UNPROTECT(1);
  return optimal;
}
