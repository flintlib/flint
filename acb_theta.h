
#ifndef ACB_THETA_H
#define ACB_THETA_H

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
   - A suffix sqr indicates squares of theta values
   - A suffix proj indicates theta values up to common scaling, and derivatives of those
   - A suffix half indicates theta values taken at tau/2
   - A suffix all indicates theta values for all characteristics (a,b), not only a=0
   - A suffix ind indicates a single theta value
   - A suffix const indicates theta constants (z=0). If not present, we compute both theta constants and regular theta values; "proj" is understood for each half independently.
   - A suffix jet indicates successive derivatives with respect to z. Return a vector of matrices as follows: one matrix per derivation order; in each of these, a row of the matrix contains partial derivatives of a fixed theta value
   - Order: naive/agm, all/ind, const, half,  proj, sqr, jet
   - Characteristics (a,b) are encoded as ulongs; first half is a, second half is b
*/

/* Extras for arb_mat's and acb_mat's */

void arb_randtest_pos(arb_t x, flint_rand_t state, slong prec, slong mag_bits);

void acb_mat_get_real(arb_mat_t re, const acb_mat_t tau);

void acb_mat_get_imag(arb_mat_t im, const acb_mat_t tau);

void acb_mat_set_arb_arb(acb_mat_t z, const arb_mat_t re, const arb_mat_t im);

void arb_mat_randtest_cho(arb_mat_t r, flint_rand_t state, slong prec, slong mag_bits);

void arb_mat_randtest_sym_pos(arb_mat_t r, flint_rand_t state, slong prec, slong mag_bits);

int arb_mat_is_nonsymmetric(const arb_mat_t m);

void arb_mat_pos_lambda(arb_t lambda, const arb_mat_t m, slong prec);

void arb_mat_reduce(arb_mat_t r, fmpz_mat_t u, const arb_mat_t m, slong prec);


/* Extras for fmpz_mat's */

void fmpz_mat_get_a(fmpz_mat_t a, const fmpz_mat_t m);

void fmpz_mat_get_b(fmpz_mat_t a, const fmpz_mat_t m);

void fmpz_mat_get_c(fmpz_mat_t a, const fmpz_mat_t m);

void fmpz_mat_get_d(fmpz_mat_t a, const fmpz_mat_t m);

void fmpz_mat_set_abcd(fmpz_mat_t m,
		       const fmpz_mat_t a, const fmpz_mat_t b,
		       const fmpz_mat_t c, const fmpz_mat_t d);

void fmpz_mat_J(fmpz_mat_t m);

int fmpz_mat_is_scalar(const fmpz_mat_t m);

int fmpz_mat_is_sp(const fmpz_mat_t m);

int fmpz_mat_is_gsp(const fmpz_mat_t m);

void fmpz_mat_diag_sp(fmpz_mat_t m, const fmpz_mat_t u);

void fmpz_mat_trig_sp(fmpz_mat_t m, const fmpz_mat_t s);

void fmpz_mat_randtest_sp(fmpz_mat_t m, flint_rand_t state, slong bits);

void fmpz_mat_siegel_fd(fmpz_mat_t m, slong j);


/* Siegel space */

void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits);

void acb_siegel_cocycle(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t tau, slong prec);

void acb_siegel_transform(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t tau, slong prec);

int acb_siegel_is_real_reduced(const acb_mat_t tau, const arb_t tol, slong prec);

int acb_siegel_not_real_reduced(const acb_mat_t tau, slong prec);

void acb_siegel_reduce_real(acb_mat_t w, fmpz_mat_t u, const acb_mat_t tau, slong prec);

void acb_siegel_reduce(acb_mat_t w, fmpz_mat_t m, const acb_mat_t tau, slong prec);

int acb_siegel_is_in_fundamental_domain(const acb_mat_t tau, slong prec);

void acb_theta_duplication(acb_ptr th2, acb_srcptr th, slong g, slong prec);

void acb_theta_duplication_all(acb_ptr th2, acb_srcptr th, slong g, slong prec);


/* Ellipsoids for naive algorithms */

struct acb_theta_eld_struct
{
  slong dim;
  slong ambient_dim;
  slong* last_coords;
  arb_struct* offset;
  arb_struct normsqr;
  
  arb_struct ctr;
  arb_struct rad;
  slong min, mid, max, step;
  struct acb_theta_eld_struct* rchildren;
  slong nr;
  struct acb_theta_eld_struct* lchildren;
  slong nl;
  slong nb_pts;
  slong* box;
};

typedef struct acb_theta_eld_struct acb_theta_eld_t[1];

#define acb_theta_eld_dim(E) ((E)->dim)
#define acb_theta_eld_ambient_dim(E) ((E)->ambient_dim)
#define acb_theta_eld_coord(E, k) ((E)->last_coords[(k) - acb_theta_eld_dim(E)])
#define acb_theta_eld_offset(E) ((E)->offset)
#define acb_theta_eld_normsqr(E) (&(E)->normsqr)
#define acb_theta_eld_ctr(E) (&(E)->ctr)
#define acb_theta_eld_rad(E) (&(E)->rad)
#define acb_theta_eld_min(E) ((E)->min)
#define acb_theta_eld_mid(E) ((E)->mid)
#define acb_theta_eld_max(E) ((E)->max)
#define acb_theta_eld_step(E) ((E)->step)
#define acb_theta_eld_rchild(E, k) (&(E)->rchildren[(k)])
#define acb_theta_eld_lchild(E, k) (&(E)->lchildren[(k)])
#define acb_theta_eld_nr(E) ((E)->nr)
#define acb_theta_eld_nl(E) ((E)->nl)
#define acb_theta_eld_nb_pts(E) ((E)->nb_pts)
#define acb_theta_eld_box(E, k) ((E)->box[(k)])

void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g);

void acb_theta_eld_clear(acb_theta_eld_t E);

void acb_theta_eld_init_children(acb_theta_eld_t E, slong nr, slong nl);

void acb_theta_eld_interval(slong* min, slong* mid, slong* max,
			    const arb_t ctr, const arb_t rad, int a, slong prec);

void acb_theta_eld_next_normsqr(arb_t next_normsqr, const arb_t normsqr, const arb_t gamma,
				const arb_t ctr, slong k, slong prec);

void acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t Y, const arb_t normsqr,
			arb_srcptr offset, slong* last_coords, ulong a, slong prec);

void acb_theta_eld_points(slong* pts, const acb_theta_eld_t E);

int acb_theta_eld_contains(const acb_theta_eld_t E, slong* pt);

void acb_theta_eld_print(const acb_theta_eld_t E);

/* Choice of radii and precisions in naive algorithms */

#define ACB_THETA_ELD_DEFAULT_PREC 50
#define ACB_THETA_NAIVE_EPS_2EXP 0
#define ACB_THETA_NAIVE_FULLPREC_ADDLOG 1.5
#define ACB_THETA_NAIVE_NEWPREC_MARGIN 1.0

void acb_theta_naive_tail(arf_t B, const arf_t R, const arb_mat_t Y, slong p, slong prec);

void acb_theta_naive_radius(arf_t R, const arb_mat_t Y, slong p, const arf_t epsilon, slong prec);

slong acb_theta_naive_newprec(slong prec, slong coord, slong dist, slong max_dist,
			      slong step, slong ord);

slong acb_theta_naive_fullprec(const acb_theta_eld_t E, slong prec);

/* Precomputations for naive algorithms */
/* For this to work, we assume that step is 1 or 2 and constant among ellipsoid layers */

typedef struct
{
  slong g;
  acb_mat_struct exp_mat;
  slong* box;
  slong step;
  slong* indices;
  acb_ptr sqr_powers;
  slong nb;
} acb_theta_precomp_struct;

typedef acb_theta_precomp_struct acb_theta_precomp_t[1];

#define acb_theta_precomp_g(D) ((D)->g)
#define acb_theta_precomp_exp_mat(D) (&(D)->exp_mat)
#define acb_theta_precomp_box(D, k) ((D)->box[(k)])
#define acb_theta_precomp_sqr_pow(D, k, i) (&(D)->sqr_powers[(i) + (D)->indices[(k)]])

void acb_theta_precomp_init(acb_theta_precomp_t D, slong g);

void acb_theta_precomp_clear(acb_theta_precomp_t D);

void acb_theta_precomp_set(acb_theta_precomp_t D, const acb_mat_t tau,
			   const acb_theta_eld_t E, slong prec);


/* Generic code for naive algorithms */

/* All naive algorithms enumerate points in ellipsoids and arrive at a
   point where an exponential term is computed. What to do with it is
   encoded in a function of the form:

   void acb_theta_worker(acb_ptr th, const acb_t term, slong* coords, slong g,
                         ulong ab, slong ord, slong prec, slong fullprec)
*/

typedef void (*acb_theta_naive_worker_t)(acb_ptr, const acb_t, slong*, slong,
					 ulong, slong, slong, slong);

/* Comments on input:
   - E is an ellipsoid of dim 1
   - each exponential term is of the form cofactor * lin^k * (some k^2-power);
     last factor is precomputed in D
   - prec is the current relative precision for this slice
   - rest of the data is passed on to the worker. */

void acb_theta_naive_worker_dim1(acb_ptr th,
				 const acb_theta_eld_t E, const acb_theta_precomp_t D,
				 const acb_t lin, const acb_t cofactor,
				 ulong ab, slong ord, slong prec, slong fullprec,
				 acb_theta_naive_worker_t worker_dim0);

/* Comments on input:
   - (k,j)-th entry of lin_powers is \exp(\pi i n_j tau_{k,j}/4) for k <= d, j > d
   - entries of lin_powers with j > d are const, others can be modified for recursion
   - k-th entry of exp_z is \exp(\pi i z_k)
   - cofactor is the common part for all exponential terms in current slice
   - rest of the data is passed on recursively as given, except that prec is adjusted 
     depending on chosen slice
   In case dim=1, fall back to worker_dim1
*/

void acb_theta_naive_worker_rec(acb_ptr th, acb_mat_t lin_powers,
				const acb_theta_eld_t E, const acb_theta_precomp_t D,
				acb_srcptr exp_z, const acb_t cofactor,
				ulong ab, slong ord, slong prec, slong fullprec,
				acb_theta_naive_worker_t worker_dim0);


/* Naive algorithms */

void acb_theta_naive_term(acb_t exp, const acb_mat_t tau, acb_srcptr z,
			  ulong ab, slong* coords, slong prec);

void acb_theta_naive_ellipsoid(acb_theta_eld_t E, arf_t epsilon,
			       ulong ab, int all, int unif, slong ord,
			       acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_naive(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_naive_const(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_naive_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_naive_all_const(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_naive_ind(acb_t th, ulong ab, acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_naive_ind_const(acb_t th, ulong ab, const acb_mat_t tau, slong prec);


slong acb_theta_nb_partials(slong ord, slong nvars);

void acb_theta_partial(slong* tup, slong k, slong ord, slong nvars);

slong acb_theta_partial_index(slong* tup, slong ord, slong nvars);


void acb_theta_jet_naive(acb_mat_struct* th, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec);

void acb_theta_const_jet_naive(acb_mat_struct* dth, const acb_mat_t tau, slong ord, slong prec);


/* Upper bounds on theta constants and their derivatives */

void acb_theta_neighborhood(arf_t rad, arf_t bound, acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_neighborhood_const(arf_t rad, arf_t bound, acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_theta_cauchy(arf_t bound_der, const arf_t rad, const arf_t bound, slong ord, slong prec);


/* AGM sequences */

#define ACB_THETA_AGM_LOWPREC 20

void acb_theta_agm_hadamard(acb_ptr r, acb_ptr s, slong g, slong prec);

void acb_theta_agm_sqrt_lowprec(acb_t r, const acb_t x, const acb_t r0, slong prec);

void acb_theta_agm_step_sqrt(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr r0, slong g, slong prec);

void acb_theta_agm_step_good(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_ext_step_sqrt(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_ext_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr r0, slong g, slong prec);

void acb_theta_agm_ext_step_good(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm(acb_t r, acb_srcptr a, acb_srcptr all_r0, slong nb_bad,
		   slong nb_total, slong g, slong prec);

void acb_theta_agm_ext(acb_t r, acb_srcptr a, acb_srcptr all_r0, slong nb_bad,
		       slong nb_total, slong g, slong prec);


/* Newton and convergence data */

slong acb_theta_agm_nb_bad_steps(const acb_mat_t tau, const fmpz_mat_t N, slong prec);

void acb_theta_agm_collect(acb_ptr all_r0, arb_t M0, arb_t minf, arb_ptr mi,
			   slong nb_bad, const acb_mat_t tau, const fmpz_mat_t N, slong prec);

void acb_theta_agm_ext_collect(acb_ptr all_r0, slong* nb_bad, slong* nb_total,
			       arb_t rho, arb_t M0, arb_t minf,
			       arb_ptr mi, acb_srcptr z, const acb_mat_t tau,
			       const fmpz_mat_t N, slong prec);

void acb_theta_agm_setup(fmpz_mat_struct* Ni, arb_t rho, arb_t M, arb_t Binv,
			 const acb_mat_t tau, slong prec);

void acb_theta_agm_ext_setup(fmpz_mat_struct* Ni, arb_t rho, arb_t M, arb_t Binv,
			     acb_srcptr z, const acb_mat_t tau, slong prec);

/* Ideas:
   - attempt setup at default prec, if not, double (up to prec/4?)
   - necessary to write generic newton?
*/

void acb_theta_agm_half_proj(acb_ptr th, acb_ptr dth, fmpz_mat_struct* gamma, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_agm_const_half_proj(acb_ptr th, acb_ptr dth, fmpz_mat_struct* gamma, const acb_mat_t tau, slong prec);

void acb_theta_agm_all_sqr(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_agm_all_const_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Mixed naive-AGM algorithms */

void acb_theta_all_sqr(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_all_const_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Conversions */

void acb_theta_all_const_from_sqr(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_all_from_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Finite difference algorithms */

void acb_theta_derivatives(acb_ptr dth, const acb_mat_t tau, const acb_mat_t dtau, acb_ptr z, acb_ptr dz, slong order, slong prec);

  
/* #ifdef __cplusplus
}
#endif */

#endif
