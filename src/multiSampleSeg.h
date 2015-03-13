/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <math.h>

#define ERROR_BIN_FACTOR_TOO_LARGE 1

struct Profile {
  int *chromStart;
  int *chromEnd;
  int *coverage;
  int n_entries;
};

int multiSampleSegOptimal(
  struct Profile **samples,
  int n_samples,
  int *optimal_start_end // array of length 2.
  );

int multiSampleSegHeuristic(
  struct Profile **samples,
  int n_samples,
  int bin_factor,
  int *optimal_start_end // array of length 2.
  );

int get_min_chromStart(struct Profile *profile);
int get_max_chromEnd(struct Profile *profile);
double OptimalPoissonLoss(int cumsum_value, double mean_value);


