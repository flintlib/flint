
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
   - A suffix const indicates that we restrict to theta constants (z=0)
   - Order of suffixes: const, half, all, proj, sqr, then algorithm type
   - Characteristics (a,b) are encoded as ulongs; first half is a, second half is b
   - Mixed algorithm is guaranteed to be uniform quasi-linear time on the Siegel fundamental domain if g <= 2
   - Vector lengths for theta values are 2^g in general, 2^(2g) in case of _all
   - Elements of the modular group Sp_2g(ZZ) are encoded as fmpz_mat's
   - Output values (not marked 'const') are always assumed to be initialized, with correct lengths in case of vectors
*/


/* Naive algorithms */

void acb_theta_naive(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_const_naive(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_all_naive(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_const_all_naive(acb_ptr th, const acb_mat_t tau, slong prec);


/* AGM algorithms */

void acb_theta_borchardt(acb_t r, acb_srcptr a, acb_srcptr sqrt_bad, slong nb_bad, slong g, slong prec);

void acb_theta_ext_borchardt(acb_t r, acb_srcptr a, acb_srcptr b, acb_srcptr sqrt_a_bad, acb_srcptr sqrt_b_bad, slong nb_bad, slong g, slong prec);

void _acb_theta_half_proj_agm(acb_ptr th, acb_ptr dth, fmpz_mat_struct* gamma, const acb_mat_t tau, acb_srcptr z, slong prec);

void _acb_theta_const_half_proj_agm(acb_ptr th, acb_ptr dth, fmpz_mat_struct* gamma, const acb_mat_t tau, slong prec);

void acb_theta_all_sqr_agm(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_const_all_sqr_agm(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);


/* Mixed naive-AGM algorithms */

void acb_theta_all_sqr(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_const_all_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Finite difference algorithms */

void acb_theta_derivatives(acb_ptr dth, const acb_mat_t tau, const acb_mat_t dtau, acb_ptr z, acb_ptr dz, slong order, slong prec);

  
/* #ifdef __cplusplus
}
#endif */

#endif
