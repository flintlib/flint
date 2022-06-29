
#ifndef ACB_MODULAR_H
#define ACB_MODULAR_H

#include <stdio.h>
#include "acb.h"

/* #ifdef __cplusplus
extern "C" {
#endif */


/* General comments:
   - Throughout, the matrix tau is assumed to be square of size g
   - In case a computation fails, output values are set to the full complex plane, according to conventions in acb_modular.h
   - A suffix sqr indicates that we compute squares of theta values
   - A suffix proj indicates that we compute theta values up to common scaling, and derivatives of those
   - A suffix half indicates that theta values are taken at tau/2
   - A suffix all indicates that theta values are computed for all characteristics (a,b), not only for a=0; in any case theta values are ordered w.r.t the characteristic
   - A suffix ind indicates that we compute a single theta value
   - A suffix const indicates that we restrict to theta constants (z=0). If not present, th vector has double length, and contains theta constants, then regular theta values; "proj" is understood for each half independently.
   - A suffix jet indicates that we compute successive derivatives with respect to z. We return a vector of matrices as follows: one matrix per derivation order; in each of these, a row of the matrix contains partial derivatives of a fixed theta value
   - Order of suffixes: const, half, all/ind, proj, sqr, jet, then algorithm type
   - Characteristics (a,b) are encoded as ulongs; first half is a, second half is b
   - Mixed algorithm is guaranteed to be uniform quasi-linear time on the Siegel fundamental domain if g <= 2
   - Vector lengths for theta values are 2^g in general, 2^(2g) in case of _all
   - Elements of the modular group Sp_2g(ZZ) are encoded as fmpz_mat's
   - Output values (not marked 'const') are always assumed to be initialized, with correct lengths in case of vectors
   - Naive algorithms are based on sums over ellipsoids. They accept an 'enum' argument, for optimization in case many values are needed for a given tau (may be NULL)
   - 
*/

/* Argument reduction */

void acb_theta_reduce_tau(fmpz_mat_t m, acb_mat_t tau, slong prec);


/* Naive algorithms */

void acb_theta_tail_ubound(arb_t B, const arb_t R, const arb_mat_t Y, slong p, slong g, slong prec);

void acb_theta_naive_radius(arb_t R, const arb_mat_t Y, const acb_t tau, slong p, const arb_t epsilon, slong prec);

/* Todo: define a structure for easy ellipsoid enumeration */

typedef struct
{
  arb_mat_struct Y;
}
acb_theta_enum_struct;

typedef acb_theta_enum_t acb_theta_enum_struct[1];

void acb_theta_enum_init(acb_theta_enum_t en);

/* Init, clear, set, get lattice vectors in some ordering, addition chains? */

void acb_theta_naive(acb_ptr th, acb_srcptr z, const acb_mat_t tau, const acb_theta_enum_t en, slong prec);

void acb_theta_const_naive(acb_ptr th, const acb_mat_t tau, const acb_theta_enum_t en, slong prec);

void acb_theta_all_naive(acb_ptr th, acb_srcptr z, const acb_mat_t tau, const acb_theta_enum_t en, slong prec);

void acb_theta_const_all_naive(acb_ptr th, const acb_mat_t tau, const acb_theta_enum_t en, slong prec);

void acb_theta_ind_naive(acb_t th, ulong ab, acb_srcptr z, const acb_mat_t tau, const acb_theta_enum_t en, slong prec);

void acb_theta_const_ind_naive(acb_t th, ulong ab, const acb_mat_t tau, const acb_theta_enum_t en, slong prec);


slong acb_theta_nb_partials(slong ord, slong nvars);

void acb_theta_partial(slong* tup, slong k, slong ord, slong nvars);

slong acb_theta_partial_index(slong* tup, slong ord, slong nvars);


void acb_theta_jet_naive(acb_mat_struct* th, acb_srcptr z, const acb_mat_t tau, const acb_theta_enum_t en, slong ord, slong prec);

void acb_theta_const_jet_naive(acb_mat_struct* dth, const acb_mat_t tau, const acb_theta_enum_t en, slong ord, slong prec);


/* AGM algorithms */

void acb_theta_borchardt(acb_t r, acb_srcptr a, acb_srcptr sqrt_bad, slong nb_bad, slong g, slong prec);

void acb_theta_ext_borchardt(acb_t r, acb_srcptr a, acb_srcptr b, acb_srcptr sqrt_a_bad, acb_srcptr sqrt_b_bad, slong nb_bad, slong g, slong prec);

void acb_theta_half_proj_agm(acb_ptr th, acb_ptr dth, fmpz_mat_struct* gamma, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_const_half_proj_agm(acb_ptr th, acb_ptr dth, fmpz_mat_struct* gamma, const acb_mat_t tau, slong prec);

void acb_theta_all_sqr_agm(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_const_all_sqr_agm(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);


/* Mixed naive-AGM algorithms */

void acb_theta_all_sqr(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_const_all_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Conversions */

void acb_theta_const_all_from_sqr(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_all_from_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Finite difference algorithms */

void acb_theta_derivatives(acb_ptr dth, const acb_mat_t tau, const acb_mat_t dtau, acb_ptr z, acb_ptr dz, slong order, slong prec);

  
/* #ifdef __cplusplus
}
#endif */

#endif
