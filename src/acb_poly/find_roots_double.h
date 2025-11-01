#ifndef ACB_POLY_FIND_ROOTS_DOUBLE_H
#define ACB_POLY_FIND_ROOTS_DOUBLE_H

#include "ulong_extras.h"

#ifdef __cplusplus
extern "C" {
#endif

slong double_cpoly_partition_pivot(double* z_r, double* z_i, slong n);

void double_cpoly_horner(double* results_r, double* results_i,
                         const double* values_r, const double* values_i, slong n_start, slong n_end,
                         const double* coefficients_r, const double* coefficients_i, slong n);

void double_cpoly_weierstrass(double* results_r, double* results_i,
                              double lc_r, double lc_i,
                              const double* twice_values_r, const double* twice_values_i,
                              slong n_start, slong n_end, slong n);

void double_cpoly_find_roots(double * z, const double * p, slong n, slong max_iter, int verbose);

#ifdef __cplusplus
}
#endif

#endif /* ACB_POLY_FIND_ROOTS_DOUBLE_H */
