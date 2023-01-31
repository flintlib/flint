/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ACB_THETA_H
#define ACB_THETA_H

#include <stdio.h>
#include "flint/fmpz_mat.h"
#include "flint/fmpz_lll.h"
#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

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

slong fmpz_mat_nb_siegel_fund(slong g);

void fmpz_mat_siegel_fund(fmpz_mat_t mat, slong j);

/* Siegel space */

void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec,
    slong mag_bits);

void acb_siegel_randtest_fund(acb_mat_t tau, flint_rand_t state,
    slong prec);

void acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec,
    slong mag_bits);

void acb_siegel_cocycle(acb_mat_t res, const fmpz_mat_t mat,
    const acb_mat_t tau, slong prec);

void acb_siegel_transform(acb_mat_t res, const fmpz_mat_t mat,
    const acb_mat_t tau, slong prec);

void acb_siegel_transform_ext(acb_ptr r, acb_mat_t w, const fmpz_mat_t mat,
    acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_siegel_reduce_imag(fmpz_mat_t mat, const acb_mat_t tau, slong prec);

void acb_siegel_reduce_real(fmpz_mat_t mat, const acb_mat_t tau, slong prec);

void acb_siegel_reduce(acb_mat_t res, fmpz_mat_t mat, const acb_mat_t tau,
    slong prec);

/* AGM sequences */

#define ACB_THETA_AGM_LOWPREC 50

void acb_theta_agm_hadamard(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_sqrt_lowprec(acb_t r, const acb_t a, const acb_t root,
    slong prec);

void acb_theta_agm_step_sqrt(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr roots, slong g,
    slong prec);

void acb_theta_agm_step_good(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_step_last(acb_t r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_ext_step_sqrt(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_ext_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr roots,
    slong g, slong prec);

void acb_theta_agm_ext_step_good(acb_ptr r, acb_srcptr a, slong g, slong prec);

void acb_theta_agm_ext_step_last(acb_t r, const acb_t s, acb_srcptr a, slong g,
    slong prec);

void acb_theta_agm_max_abs(arb_t max, acb_srcptr a, slong nb, slong prec);

void acb_theta_agm_min_abs(arb_t min, acb_srcptr a, slong nb, slong prec);

void acb_theta_agm_abs_dist(arb_t eps, acb_srcptr a, slong nb, slong lowprec,
    slong prec);

void acb_theta_agm_rel_dist(arb_t eps, acb_srcptr a, slong nb, slong lowprec,
    slong prec);

void acb_theta_agm_radius(arf_t rad, const arf_struct* mi, const arf_struct* Mi,
    const arf_t abs_dist, slong nb, slong prec);

void acb_theta_agm_conv_rate(arf_t c, arf_t r, const arf_t eps, slong prec);

slong acb_theta_agm_nb_good_steps(const arf_t c, const arf_t r, slong prec);

void acb_theta_agm(acb_t r, acb_srcptr a, acb_srcptr roots, slong nb_bad,
    slong g, slong prec);

void acb_theta_agm_ext_conv_rate(arf_t c1, arf_t c2, arf_t r, const arf_t eps,
    const arf_t m, const arf_t M, slong prec);

void acb_theta_agm_ext_rel_err(arf_t err, const arf_t c2, const arf_t r,
    slong nb_good, slong prec);

void acb_theta_agm_ext(acb_t r, acb_t s, acb_srcptr a, acb_srcptr roots,
    slong nb_bad, slong g, slong prec);

slong acb_theta_agm_nb_bad_steps(const acb_mat_t tau, slong prec);

slong acb_theta_agm_ext_nb_bad_steps(acb_srcptr z, const acb_mat_t tau,
    slong prec);

void acb_theta_agm_roots(acb_ptr roots, const acb_mat_t tau, slong nb_bad,
    slong prec);

void acb_theta_agm_ext_roots(acb_ptr roots, acb_srcptr z, const acb_mat_t tau,
    slong nb_bad, slong prec);

/* Transformation formulas */

slong acb_theta_char_dot(ulong a, ulong b, slong g);

slong acb_theta_dot(ulong a, slong* n, slong g);

void acb_theta_vecsqr(acb_ptr th2, acb_srcptr th, slong n, slong prec);

void acb_theta_dupl_const(acb_ptr th2, acb_srcptr th, slong g, slong prec);

void acb_theta_dupl_all_const(acb_ptr th2, acb_srcptr th, slong g, slong prec);

void acb_theta_dupl(acb_ptr th2, acb_srcptr th, slong g, slong prec);

void acb_theta_dupl_all(acb_ptr th2, acb_srcptr th, slong g, slong prec);

void acb_theta_dupl_z(acb_ptr r, acb_srcptr th, slong g, slong prec);

ulong acb_theta_transform_image_char(fmpz_t eps, ulong ab,
    const fmpz_mat_t mat);

void acb_theta_transform_proj(acb_ptr res, acb_srcptr th,
        const fmpz_mat_t mat, slong prec);

void acb_theta_transform_sqr_proj(acb_ptr res, acb_srcptr th2,
        const fmpz_mat_t mat, slong prec);

void acb_theta_transform_all_sqr_proj(acb_ptr res, acb_srcptr th2,
    const fmpz_mat_t mat, slong prec);
    
void acb_theta_transform_scal_const(acb_t scal, const acb_mat_t tau,
        const fmpz_mat_t mat, slong k2, slong prec);

void acb_theta_transform_scal(acb_t scal_z, acb_t scal_0, acb_srcptr z,
        const acb_mat_t tau, const fmpz_mat_t mat, slong k2, slong prec);

void acb_theta_dupl_radius(arf_t rho, const arf_t r, acb_srcptr th, slong nb,
        slong prec);

void acb_theta_transform_sqr_radius(arf_t rho, const arf_t r, acb_srcptr th2,
        const fmpz_mat_t mat, slong prec);

void acb_theta_dupl_transform_radius_const(arf_t rho, const arf_t r,
        acb_srcptr th, const fmpz_mat_t mat, slong prec);

void acb_theta_dupl_transform_radius(arf_t rho, const arf_t r,
        acb_srcptr th, const fmpz_mat_t mat, slong prec);

/* Naive algorithms */

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

void acb_theta_eld_round(slong* r, const arb_mat_t v);

void acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t Y, const arf_t R2,
        arb_srcptr offset, slong* last_coords, ulong a, slong prec);

void acb_theta_eld_points(slong* pts, const acb_theta_eld_t E);

int acb_theta_eld_contains(const acb_theta_eld_t E, slong* pt);

void acb_theta_eld_print(const acb_theta_eld_t E);

#define ACB_THETA_ELD_DEFAULT_PREC 50
#define ACB_THETA_NAIVE_EPS_2EXP 0
#define ACB_THETA_NAIVE_FULLPREC_ADDLOG 1.1
#define ACB_THETA_NAIVE_NEWPREC_MARGIN 1.0

void acb_theta_naive_tail(arf_t bound, const arf_t R2, const arb_mat_t Y,
        slong ord, slong prec);

void acb_theta_naive_radius(arf_t R2, const arb_mat_t Y, slong ord,
        const arf_t eps, slong prec);

void acb_theta_naive_ellipsoid(acb_theta_eld_t E, arf_struct* eps, acb_ptr c,
        acb_ptr new_z, ulong ab, int all, slong ord, acb_srcptr z, slong nb_z,
        const acb_mat_t tau, slong prec);

slong acb_theta_naive_newprec(slong prec, slong coord, slong dist,
        slong max_dist, slong ord);

slong acb_theta_naive_fullprec(const acb_theta_eld_t E, slong prec);

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

typedef void (*acb_theta_naive_worker_t)(acb_ptr, const acb_t, slong*, slong,
        ulong, slong, slong, slong);

void acb_theta_naive_worker(acb_ptr th, slong nb, const acb_t c,
        const arf_t eps, const acb_theta_eld_t E,
        const acb_theta_precomp_t D, slong k, ulong ab, slong ord, slong prec,
        acb_theta_naive_worker_t worker_dim0);

ulong acb_theta_naive_a(slong* coords, slong g);

void acb_theta_naive(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau,
        slong prec);

void acb_theta_naive_proj(acb_ptr th, acb_srcptr z, slong nb_z,
        const acb_mat_t tau, slong prec);

void acb_theta_naive_const(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_naive_const_proj(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_naive_all(acb_ptr th, acb_srcptr z, slong nb_z,
        const acb_mat_t tau, slong prec);

void acb_theta_naive_all_const(acb_ptr th, const acb_mat_t tau, slong prec);

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


/* Conversions */

void acb_theta_all_const_from_sqr(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_all_from_sqr(acb_ptr th, const acb_mat_t tau, slong prec);

void acb_theta_renormalize_const_sqr(acb_t scal, acb_srcptr th2,
        const acb_mat_t tau, slong prec);

void acb_theta_renormalize_sqr(acb_t scal_z, acb_t scal_0, acb_srcptr th2_z,
        acb_srcptr th2_0, acb_srcptr z, const acb_mat_t tau, slong prec);

slong acb_theta_k2(const fmpz_mat_t mat);

/* Newton iterations */

void acb_theta_bound(arf_t rad, arf_t bound, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_bound_const(arf_t rad, arf_t bound, const acb_mat_t tau,
        slong prec);

void acb_theta_cauchy(arf_t bound_der, const arf_t rad, const arf_t bound,
        slong ord, slong dim, slong prec);

#define ACB_THETA_AGM_NB_MATRIX_SETUPS 10
#define ACB_THETA_AGM_BASEPREC 2000
#define ACB_THETA_AGM_BASEPREC_MAXQ 4
#define ACB_THETA_AGM_GUARD 5

typedef struct
{
    int is_ext;
    slong dim;
    acb_mat_struct tau;
    acb_struct* z;
    acb_struct* th;
    slong nb;
    
    fmpz_mat_struct* mat;
    slong* k2;
    ulong* ab;
    fmpz* eps;
    slong* nb_bad;
    acb_ptr* roots;
    
    slong log_th;
    slong log_rho;
    slong log_M;
    slong log_B1;
    slong log_B2;
    slong log_B3;
    
} acb_theta_agm_ctx_struct;

typedef acb_theta_agm_ctx_struct acb_theta_agm_ctx_t[1];

#define acb_theta_agm_ctx_is_ext(ctx) ((ctx)->is_ext)
#define acb_theta_agm_ctx_dim(ctx) ((ctx)->dim)
#define acb_theta_agm_ctx_tau(ctx) (&(ctx)->tau)
#define acb_theta_agm_ctx_z(ctx) ((ctx)->z)
#define acb_theta_agm_ctx_th(ctx) ((ctx)->th)
#define acb_theta_agm_ctx_g(ctx) (acb_mat_nrows(acb_theta_agm_ctx_tau(ctx)))
#define acb_theta_agm_ctx_nb(ctx) ((ctx)->nb)

#define acb_theta_agm_ctx_mat(ctx, k) (&(ctx)->mat[(k)])
#define acb_theta_agm_ctx_k2(ctx, k) ((ctx)->k2[(k)])
#define acb_theta_agm_ctx_ab(ctx, k) ((ctx)->ab[(k)])
#define acb_theta_agm_ctx_eps(ctx, k) (&(ctx)->eps[(k)])
#define acb_theta_agm_ctx_nb_bad(ctx, k) ((ctx)->nb_bad[(k)])
#define acb_theta_agm_ctx_roots(ctx, k) ((ctx)->roots[(k)])

#define acb_theta_agm_ctx_log_th(ctx) ((ctx)->log_th)
#define acb_theta_agm_ctx_log_rho(ctx) ((ctx)->log_rho)
#define acb_theta_agm_ctx_log_M(ctx) ((ctx)->log_M)
#define acb_theta_agm_ctx_log_B1(ctx) ((ctx)->log_B1)
#define acb_theta_agm_ctx_log_B2(ctx) ((ctx)->log_B2)
#define acb_theta_agm_ctx_log_B3(ctx) ((ctx)->log_B3)

void acb_theta_agm_ctx_init_internal(acb_theta_agm_ctx_t ctx, slong nb,
        slong g);

void acb_theta_agm_ctx_init(acb_theta_agm_ctx_t ctx, const acb_mat_t tau);

void acb_theta_agm_ctx_init_ext(acb_theta_agm_ctx_t ctx, acb_srcptr z,
        const acb_mat_t tau);

void acb_theta_agm_ctx_clear(acb_theta_agm_ctx_t ctx);

void acb_theta_agm_ctx_reset_steps(acb_theta_agm_ctx_t ctx, slong k, slong m);

int acb_theta_agm_ctx_set(acb_theta_agm_ctx_t ctx, slong prec);

void acb_theta_newton_eval(acb_ptr r, acb_srcptr th,
        const acb_theta_agm_ctx_t ctx, slong prec);

void acb_theta_newton_fd(acb_ptr r, acb_mat_t fd, acb_srcptr th,
        const arb_t eta, const acb_theta_agm_ctx_t ctx, slong prec);

void acb_theta_newton_run(acb_ptr r, const acb_theta_agm_ctx_t ctx, slong prec);

void acb_theta_newton_const_half_proj(acb_ptr th, const acb_mat_t tau,
        slong prec);

void acb_theta_newton_const_sqr(acb_ptr th2, const acb_mat_t tau, slong prec);

void acb_theta_newton_all_const_sqr(acb_ptr th2, const acb_mat_t tau,
        slong prec);

void acb_theta_newton_half_proj(acb_ptr th, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_newton_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau,
        slong prec);

void acb_theta_newton_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau,
        slong prec);

/* Mixed Newton/naive algorithms */

#define ACB_THETA_NAIVE_CONST_THRESHOLD 500
#define ACB_THETA_NAIVE_THRESHOLD 500
#define ACB_THETA_BALANCE_LOWPREC_MUL 10
#define ACB_THETA_BALANCE_THRESHOLD 4
#define ACB_THETA_REDUCE_Z 4

void acb_theta_balance(acb_ptr z2, acb_mat_t tau2, fmpz_mat_t mat,
    acb_srcptr z, const acb_mat_t tau, slong j);

int acb_theta_is_balanced(slong* j0, const acb_mat_t tau, slong prec);

slong acb_theta_balance_lowprec(slong g, slong prec);

void acb_theta_all_const_sqr(acb_ptr th2, const acb_mat_t tau, slong prec);

void acb_theta_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau,
    slong prec);

#ifdef __cplusplus
}
#endif

#endif
