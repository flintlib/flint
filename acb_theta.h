
#ifndef ACB_THETA_H
#define ACB_THETA_H

#include <stdio.h>
#include "flint/fmpz_mat.h"
#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Extras for arb_mat's and acb_mat's */

void arb_randtest_pos(arb_t x, flint_rand_t state, slong prec, slong mag_bits);

void acb_randtest_disk(acb_t x, const acb_t ctr, const arf_t rad,
	flint_rand_t state, slong prec);

void acb_mat_get_real(arb_mat_t re, const acb_mat_t mat);

void acb_mat_get_imag(arb_mat_t im, const acb_mat_t mat);

void acb_mat_set_arb_arb(acb_mat_t mat, const arb_mat_t re,
	const arb_mat_t im);

void arb_mat_add_error_arf(arb_mat_t mat, const arf_t err);

void arb_mat_randtest_cho(arb_mat_t mat, flint_rand_t state, slong prec,
	slong mag_bits);

void arb_mat_randtest_sym_pos(arb_mat_t mat, flint_rand_t state, slong prec,
	slong mag_bits);

int arb_mat_is_nonsymmetric(const arb_mat_t mat);

void arb_mat_pos_lambda(arb_t lambda, const arb_mat_t mat, slong prec);

void arb_mat_pos_radius(arf_t rad, const arb_mat_t mat, slong prec);

void arb_mat_reduce(arb_mat_t R, fmpz_mat_t U, const arb_mat_t M, slong prec);

void acb_mat_ninf(arb_t norm, const acb_mat_t mat, slong prec);


/* Extras for fmpz_mat's */

void fmpz_mat_get_a(fmpz_mat_t res, const fmpz_mat_t mat);

void fmpz_mat_get_b(fmpz_mat_t res, const fmpz_mat_t mat);

void fmpz_mat_get_c(fmpz_mat_t res, const fmpz_mat_t mat);

void fmpz_mat_get_d(fmpz_mat_t res, const fmpz_mat_t mat);

void fmpz_mat_set_abcd(fmpz_mat_t mat, const fmpz_mat_t a, const fmpz_mat_t b,
	const fmpz_mat_t c, const fmpz_mat_t d);

void fmpz_mat_J(fmpz_mat_t mat);

int fmpz_mat_is_scalar(const fmpz_mat_t mat);

int fmpz_mat_is_sp(const fmpz_mat_t mat);

int fmpz_mat_is_gsp(const fmpz_mat_t mat);

void fmpz_mat_diag_sp(fmpz_mat_t mat, const fmpz_mat_t U);

void fmpz_mat_trig_sp(fmpz_mat_t mat, const fmpz_mat_t S);

void fmpz_mat_randtest_sp(fmpz_mat_t mat, flint_rand_t state, slong bits);

void fmpz_mat_siegel_fund(fmpz_mat_t mat, slong j);


/* Siegel space */

void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec,
	slong mag_bits);

void acb_siegel_randtest_fund(acb_mat_t tau, flint_rand_t state,
	slong prec);

void acb_siegel_cocycle(acb_mat_t res, const fmpz_mat_t mat,
	const acb_mat_t tau, slong prec);

void acb_siegel_transform(acb_mat_t res, const fmpz_mat_t mat,
	const acb_mat_t tau, slong prec);

int acb_siegel_is_real_reduced(const acb_mat_t tau, const arf_t eps,
	slong prec);

int acb_siegel_not_real_reduced(const acb_mat_t tau, slong prec);

void acb_siegel_reduce_real(acb_mat_t res, fmpz_mat_t mat, const acb_mat_t tau,
	slong prec);

void acb_siegel_reduce(acb_mat_t res, fmpz_mat_t mat, const acb_mat_t tau,
	slong prec);

int acb_siegel_is_reduced(const acb_mat_t tau, const arf_t eps, slong prec);


/* AGM sequences */

#define ACB_THETA_AGM_LOWPREC 20

void acb_theta_agm_hadamard(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_sqrt_lowprec(acb_t r, const acb_t a, const acb_t root,
	slong prec);

void acb_theta_agm_step_sqrt(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr roots, slong g,
	slong prec);

void acb_theta_agm_step_good(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_ext_step_sqrt(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_ext_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr roots,
        slong g, slong prec);

void acb_theta_agm_ext_step_good(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm(acb_t r, acb_srcptr a, acb_srcptr all_roots,
	const arf_t rel_err, slong nb_bad, slong nb_good, slong g, slong prec);

void acb_theta_agm_ext(acb_t r, acb_srcptr a, acb_srcptr all_roots,
	const arf_t rel_err, slong nb_bad, slong nb_good, slong g, slong prec);

slong acb_theta_agm_nb_bad_steps(const acb_mat_t tau, slong prec);

slong acb_theta_agm_nb_good_steps(arf_t rel_err, slong g, slong prec);


/* Transformation formulas */

slong acb_theta_char_dot(ulong a, ulong b, slong g);

slong acb_theta_dot(ulong a, slong* n, slong g);

void acb_theta_duplication(acb_ptr th2, acb_srcptr th, slong g, slong prec);

void acb_theta_duplication_all(acb_ptr th2, acb_srcptr th, slong g,
	slong prec);

void acb_theta_duplication_ext(acb_ptr th2, acb_srcptr th, slong g,
        slong prec);

void acb_theta_duplication_all_ext(acb_ptr tr2, acb_srcptr th, slong g,
        slong prec);

ulong acb_theta_transform_image_char(fmpz_t eps, ulong ab,
	const fmpz_mat_t mat);

void acb_theta_transform_sqr_proj(acb_ptr res, acb_srcptr th2,
        const fmpz_mat_t mat, slong prec);


/* Ellipsoids for naive algorithms */

struct acb_theta_eld_struct
{
  slong dim;
  slong ambient_dim;
  slong* last_coords;
  slong min, mid, max;
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
#define acb_theta_eld_min(E) ((E)->min)
#define acb_theta_eld_mid(E) ((E)->mid)
#define acb_theta_eld_max(E) ((E)->max)
#define acb_theta_eld_nr(E) ((E)->nr)
#define acb_theta_eld_nl(E) ((E)->nl)
#define acb_theta_eld_rchild(E, k) (&(E)->rchildren[(k)])
#define acb_theta_eld_lchild(E, k) (&(E)->lchildren[(k)])
#define acb_theta_eld_nb_pts(E) ((E)->nb_pts)
#define acb_theta_eld_box(E, k) ((E)->box[(k)])

void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g);

void acb_theta_eld_clear(acb_theta_eld_t E);

void acb_theta_eld_interval(slong* min, slong* mid, slong* max,
        const arb_t ctr, const arf_t rad, int a, slong prec);

void acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t Y, const arf_t R2,
        arb_srcptr offset, slong* last_coords, ulong a, slong prec);

void acb_theta_eld_points(slong* pts, const acb_theta_eld_t E);

int acb_theta_eld_contains(const acb_theta_eld_t E, slong* pt);

void acb_theta_eld_print(const acb_theta_eld_t E);


/* Choice of radii and precisions in naive algorithms */

#define ACB_THETA_ELD_DEFAULT_PREC 50
#define ACB_THETA_NAIVE_EPS_2EXP 0
#define ACB_THETA_NAIVE_FULLPREC_ADDLOG 1.1
#define ACB_THETA_NAIVE_NEWPREC_MARGIN 1.0

void acb_theta_naive_tail(arf_t bound, const arf_t R2, const arb_mat_t Y,
        slong ord, slong prec);

void acb_theta_naive_radius(arf_t R2, const arb_mat_t Y, slong ord,
        const arf_t eps, slong prec);

void acb_theta_naive_ellipsoid(acb_theta_eld_t E, arf_t eps, ulong ab,
        int all, int unif, slong ord, acb_srcptr z, const acb_mat_t tau,
        slong prec);

slong acb_theta_naive_newprec(slong prec, slong coord, slong dist,
        slong max_dist, slong ord);

slong acb_theta_naive_fullprec(const acb_theta_eld_t E, slong prec);


/* Precomputations for naive algorithms */

typedef struct
{
    acb_mat_struct exp_mat;
    acb_ptr sqr_powers;
    slong* indices;
    acb_ptr exp_z;
    slong nb_z;
} acb_theta_precomp_struct;

typedef acb_theta_precomp_struct acb_theta_precomp_t[1];

#define acb_theta_precomp_exp_mat(D) (&(D)->exp_mat)
#define acb_theta_precomp_sqr_pow(D, k, j) \
    (&(D)->sqr_powers[(j) + (D)->indices[(k)]])
#define acb_theta_precomp_exp_z(D, k, j) \
    (&(D)->exp_z[(k) * acb_mat_nrows(acb_theta_precomp_exp_mat(D)) + (j)])
#define acb_theta_precomp_nb_z(D) ((D)->nb_z)

void acb_theta_precomp_init(acb_theta_precomp_t D, slong nb_z, slong g);

void acb_theta_precomp_clear(acb_theta_precomp_t D);

void acb_theta_precomp_set(acb_theta_precomp_t D, acb_srcptr z,
        const acb_mat_t tau, const acb_theta_eld_t E, slong prec);


/* Naive algorithms */

typedef void (*acb_theta_naive_worker_t)(acb_ptr, const acb_t, slong*, slong,
        ulong, slong, slong, slong);

void acb_theta_naive_worker(acb_ptr th, slong nb, const arf_t epsilon,
        const acb_theta_eld_t E, const acb_theta_precomp_t D, slong k,
        ulong ab, slong ord, slong prec, acb_theta_naive_worker_t worker_dim0);

ulong acb_theta_naive_a(slong* coords, slong g);

void acb_theta_naive(acb_ptr th, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_naive_const(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_naive_const_proj(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_naive_ext(acb_ptr th, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_naive_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_naive_all_const(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_naive_all_ext(acb_ptr th, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_naive_ind(acb_t th, ulong ab, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_naive_ind_const(acb_t th, ulong ab, const acb_mat_t tau,
        slong prec);

slong acb_theta_nb_partials(slong ord, slong nvars);

void acb_theta_partial(slong* tup, slong k, slong ord, slong nvars);

slong acb_theta_partial_index(slong* tup, slong ord, slong nvars);

void acb_theta_jet_naive(acb_mat_struct* th, acb_srcptr z, const acb_mat_t tau,
        slong ord, slong prec);

void acb_theta_const_jet_naive(acb_mat_struct* dth, const acb_mat_t tau,
        slong ord, slong prec);


/* Upper bounds on theta constants and their derivatives */

void acb_theta_bound(arf_t rad, arf_t bound, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_bound_const(arf_t rad, arf_t bound, const acb_mat_t tau,
        slong prec);

void acb_theta_cauchy(arf_t bound_der, const arf_t rad, const arf_t bound,
        slong ord, slong dim, slong prec);

/* Context for Newton iterations */

#define ACB_THETA_AGM_NB_MATRIX_SETUPS 10
#define ACB_THETA_AGM_BASEPREC 2000
#define ACB_THETA_AGM_BASEPREC_MAXQ 4
#define ACB_THETA_AGM_GUARD 5

typedef struct
{
  slong g, nb;
  fmpz_mat_struct* matrices;
  slong* nb_bad_steps;
  acb_ptr* roots;
  arf_struct** mi;
  arf_struct* M0;
  arf_struct* minf;
  arf_struct rho, max, inv_der;
} acb_theta_agm_ctx_struct;

typedef acb_theta_agm_ctx_struct acb_theta_agm_ctx_t[1];

#define acb_theta_agm_ctx_g(ctx) ((ctx)->g)
#define acb_theta_agm_ctx_nb(ctx) ((ctx)->nb)
#define acb_theta_agm_ctx_matrix(ctx, k) (&(ctx)->matrices[(k)])
#define acb_theta_agm_ctx_nb_bad_steps(ctx, k) ((ctx)->nb_bad_steps[(k)])
#define acb_theta_agm_ctx_roots(ctx, k) ((ctx)->roots[(k)])
#define acb_theta_agm_ctx_mi(ctx, k) ((ctx)->mi[(k)])
#define acb_theta_agm_ctx_M0(ctx, k) (&(ctx)->M0[(k)])
#define acb_theta_agm_ctx_minf(ctx, k) (&(ctx)->minf[(k)])
#define acb_theta_agm_ctx_rho(ctx) (&(ctx)->rho)
#define acb_theta_agm_ctx_max(ctx) (&(ctx)->max)
#define acb_theta_agm_ctx_inv_der(ctx) (&(ctx)->inv_der)

void acb_theta_agm_ctx_init(acb_theta_agm_ctx_t ctx, slong g, slong n);

void acb_theta_agm_ctx_clear(acb_theta_agm_ctx_t ctx);

void acb_theta_agm_ctx_reset_steps(acb_theta_agm_ctx_t ctx, slong k, slong m);

void acb_theta_agm_ctx_set_matrix(acb_theta_agm_ctx_t ctx, slong k, const acb_mat_t tau,
				  const fmpz_mat_t N, slong prec);

void acb_theta_agm_ctx_matrices(fmpz_mat_struct* Ni, slong k, slong g);

void acb_theta_agm_ctx_set_all(acb_theta_agm_ctx_t ctx, const acb_mat_t tau, slong prec);

int acb_theta_agm_ctx_is_valid(const acb_theta_agm_ctx_t ctx);


/* Newton iterations */

void acb_theta_newton_eval(acb_ptr r, acb_srcptr th, const acb_theta_agm_ctx_t ctx, slong prec);

void acb_theta_newton_fd(acb_ptr r, acb_mat_t fd, acb_srcptr th, const arb_t eta,
			 const acb_theta_agm_ctx_t ctx, slong prec);

void acb_theta_newton_logs(slong* log_max, slong* log_rho, slong* log_B1, slong* log_B2,
			   slong* log_B3, const acb_theta_agm_ctx_t ctx);

slong acb_theta_newton_start(acb_ptr start, acb_ptr im, arf_t err, const acb_mat_t tau,
			     const acb_theta_agm_ctx_t ctx, slong prec);

slong acb_theta_newton_step(acb_ptr next, acb_srcptr current, acb_srcptr im,
			    const acb_theta_agm_ctx_t, slong prec);

void acb_theta_newton_run(acb_ptr r, const acb_mat_t tau, const acb_theta_agm_ctx_t ctx,
			  slong prec);


/* AGM/Newton algorithms for theta functions */

void acb_theta_newton_const_half_proj(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_newton_all_sqr(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_newton_all_const_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Mixed naive-AGM algorithms */

void acb_theta_all_sqr(acb_ptr th, const acb_mat_t tau, acb_srcptr z, slong prec);

void acb_theta_all_const_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Conversions */

void acb_theta_all_const_from_sqr(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_all_from_sqr(acb_ptr th, const acb_mat_t tau, slong prec);


/* Finite difference algorithms */

  
#ifdef __cplusplus
}
#endif

#endif
