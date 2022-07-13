
#ifndef ACB_MODULAR_H
#define ACB_MODULAR_H

#include <stdio.h>
#include "flint/fmpz_mat.h"
#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"

/* #ifdef __cplusplus
extern "C" {
#endif */


/* General comments:
   - In case a computation fails, output values are set to NaNs if possible, otherwise abort
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
*/


/* Argument reduction */

void acb_theta_reduce_tau(fmpz_mat_t m, acb_mat_t tau, slong prec);


/* Ellipsoids for naive algorithms */

struct arb_eld_struct
{
  slong dim;
  slong ambient_dim;
  slong* last_coords;
  arb_struct* offset;
  arb_struct normsqr;
  
  arb_struct ctr;
  arb_struct rad;
  slong min, mid, max, step;
  struct arb_eld_struct* rchildren;
  slong nr;
  struct arb_eld_struct* lchildren;
  slong nl;
  slong nb_pts;
};

typedef struct arb_eld_struct arb_eld_t[1];

#define arb_eld_dim(E) ((E)->dim)
#define arb_eld_ambient_dim(E) ((E)->ambient_dim)
#define arb_eld_coord(E, k) ((E)->last_coords[(k) - arb_eld_dim(E)])
#define arb_eld_offset(E) ((E)->offset)
#define arb_eld_normsqr(E) (&(E)->normsqr)
#define arb_eld_ctr(E) (&(E)->ctr)
#define arb_eld_rad(E) (&(E)->rad)
#define arb_eld_min(E) ((E)->min)
#define arb_eld_mid(E) ((E)->mid)
#define arb_eld_max(E) ((E)->max)
#define arb_eld_step(E) ((E)->step)
#define arb_eld_rchild(E, k) (&(E)->rchildren[(k)])
#define arb_eld_lchild(E, k) (&(E)->lchildren[(k)])
#define arb_eld_nr(E) ((E)->nr)
#define arb_eld_nl(E) ((E)->nl)
#define arb_eld_nb_pts(E) ((E)->nb_pts)

void arb_eld_init(arb_eld_t E, slong d, slong g);

void arb_eld_clear(arb_eld_t E);

void arb_eld_init_children(arb_eld_t E, slong nr, slong nl);

void arb_eld_interval(slong* nmin, slong* nmid, slong* nmax,
		      const arb_t ctr, const arb_t rad, int a, slong prec);

void arb_eld_next_normsqr(arb_t next_normsqr, const arb_t normsqr, const arb_t gamma,
			  const arb_t ctr, slong k, slong prec);

void arb_eld_fill(arb_eld_t E, const arb_mat_t Y, const arb_t normsqr,
		  arb_srcptr offset, slong* last_coords, ulong a, slong prec);

void arb_eld_points(slong* pts, const arb_eld_t E);


/* Naive algorithms */

void acb_theta_naive_tail(arf_t B, const arf_t R, const arb_mat_t Y, slong p, slong prec);

void acb_theta_naive_radius(arf_t R, const arb_mat_t Y, slong p, const arf_t epsilon, slong prec);


void acb_theta_naive(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_const_naive(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_all_naive(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_const_all_naive(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_ind_naive(acb_t th, ulong ab, acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_const_ind_naive(acb_t th, ulong ab, const acb_mat_t tau, slong prec);


slong acb_theta_nb_partials(slong ord, slong nvars);

void acb_theta_partial(slong* tup, slong k, slong ord, slong nvars);

slong acb_theta_partial_index(slong* tup, slong ord, slong nvars);


void acb_theta_jet_naive(acb_mat_struct* th, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec);

void acb_theta_const_jet_naive(acb_mat_struct* dth, const acb_mat_t tau, slong ord, slong prec);


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
