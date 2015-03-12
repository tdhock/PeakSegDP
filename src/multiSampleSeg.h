/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_FEWER_BASES_THAN_BINS 1

struct Profile {
  int *chromStart;
  int *chromEnd;
  int *coverage;
  int n_entries;
};

int multiSampleSeg(
  struct Profile **samples,
  int n_samples,
  int n_bins,
  int *optimal_start_end // array of length 2.
  );

int get_min_chromStart(struct Profile *profile){
  return profile->chromStart[0];
}

int get_max_chromEnd(struct Profile *profile){
  return profile->chromEnd[profile->n_entries-1];
}
