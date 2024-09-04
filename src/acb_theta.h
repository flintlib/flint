/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_THETA_H
#define ACB_THETA_H

#include "fmpz_mat.h"
#include "acb_types.h"
#include "acb_theta_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ACB_THETA_LOW_PREC 32

/* The Siegel modular group */

FLINT_FORCE_INLINE slong
sp2gz_dim(const fmpz_mat_t mat)
{
    return fmpz_mat_nrows(mat) / 2;
}

void sp2gz_set_blocks(fmpz_mat_t mat, const fmpz_mat_t alpha, const fmpz_mat_t beta,
    const fmpz_mat_t gamma, const fmpz_mat_t delta);
void sp2gz_j(fmpz_mat_t mat);
void sp2gz_block_diag(fmpz_mat_t mat, const fmpz_mat_t U);
void sp2gz_trig(fmpz_mat_t mat, const fmpz_mat_t S);
void sp2gz_embed(fmpz_mat_t res, const fmpz_mat_t mat);
void sp2gz_restrict(fmpz_mat_t res, const fmpz_mat_t mat);

slong sp2gz_nb_fundamental(slong g);
void sp2gz_fundamental(fmpz_mat_t mat, slong j);

int sp2gz_is_correct(const fmpz_mat_t mat);
int sp2gz_is_j(const fmpz_mat_t mat);
int sp2gz_is_block_diag(const fmpz_mat_t mat);
int sp2gz_is_trig(const fmpz_mat_t mat);
int sp2gz_is_embedded(fmpz_mat_t res, const fmpz_mat_t mat);

void sp2gz_inv(fmpz_mat_t inv, const fmpz_mat_t mat);
fmpz_mat_struct * sp2gz_decompose(slong * nb, const fmpz_mat_t mat);

void sp2gz_randtest(fmpz_mat_t mat, flint_rand_t state, slong bits);

/* The Siegel half space */

void acb_siegel_cocycle(acb_mat_t c, const fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_transform_cocycle_inv(acb_mat_t w, acb_mat_t c, acb_mat_t cinv,
    const fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_transform(acb_mat_t w, const fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_transform_z(acb_ptr r, acb_mat_t w, const fmpz_mat_t mat,
    acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_siegel_cho(arb_mat_t C, const acb_mat_t tau, slong prec);
void acb_siegel_yinv(arb_mat_t Yinv, const acb_mat_t tau, slong prec);

void acb_siegel_reduce(fmpz_mat_t mat, const acb_mat_t tau, slong prec);
int acb_siegel_is_reduced(const acb_mat_t tau, slong tol_exp, slong prec);

void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits);
void acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits);
void acb_siegel_randtest_vec(acb_ptr z, flint_rand_t state, slong g, slong prec);
void acb_siegel_randtest_vec_reduced(acb_ptr z, flint_rand_t state,
    const acb_mat_t tau, int exact, slong prec);

/* Theta characteristics */

void acb_theta_char_get_slong(slong * n, ulong a, slong g);
ulong acb_theta_char_get_a(const slong * n, slong g);
void acb_theta_char_get_arb(arb_ptr v, ulong a, slong g);
void acb_theta_char_get_acb(acb_ptr v, ulong a, slong g);

slong acb_theta_char_dot(ulong a, ulong b, slong g);
slong acb_theta_char_dot_slong(ulong a, const slong * n, slong g);
void acb_theta_char_dot_acb(acb_t x, ulong a, acb_srcptr z, slong g, slong prec);

int acb_theta_char_is_even(ulong ab, slong g);
int acb_theta_char_is_goepel(ulong ch1, ulong ch2, ulong ch3, ulong ch4, slong g);
int acb_theta_char_is_syzygous(ulong ch1, ulong ch2, ulong ch3, slong g);

/* Ellipsoids in summation algorithms */

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
#define acb_theta_eld_nb_border(E) ((E)->nb_border)
#define acb_theta_eld_box(E, k) ((E)->box[(k)])

void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g);
void acb_theta_eld_clear(acb_theta_eld_t E);

int acb_theta_eld_set(acb_theta_eld_t E, const arb_mat_t C, const arf_t R2, arb_srcptr v);
void acb_theta_eld_points(slong * pts, const acb_theta_eld_t E);
void acb_theta_eld_border(slong * pts, const acb_theta_eld_t E);
int acb_theta_eld_contains(const acb_theta_eld_t E, const slong * pt);
void acb_theta_eld_print(const acb_theta_eld_t E);

/* Distances */

void acb_theta_dist_pt(arb_t d, arb_srcptr v, const arb_mat_t C, const slong * n, slong prec);
void acb_theta_dist_lat(arb_t d, arb_srcptr v, const arb_mat_t C, slong prec);
void acb_theta_dist_a0(arb_ptr d, acb_srcptr z, const acb_mat_t tau, slong prec);
slong acb_theta_dist_addprec(const arb_t d);

/* AGM steps */

void acb_theta_agm_hadamard(acb_ptr res, acb_srcptr a, slong g, slong prec);
void acb_theta_agm_sqrt(acb_ptr res, acb_srcptr a, acb_srcptr roots, slong nb, slong prec);
void acb_theta_agm_mul(acb_ptr res, acb_srcptr a1, acb_srcptr a2, slong g, slong prec);
void acb_theta_agm_mul_tight(acb_ptr res, acb_srcptr a0, acb_srcptr a,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec);

/* The transformation formula */

ulong acb_theta_transform_char(slong * e, const fmpz_mat_t mat, ulong ab);
void acb_theta_transform_sqrtdet(acb_t res, const acb_mat_t tau, slong prec);
slong acb_theta_transform_kappa(acb_t sqrtdet, const fmpz_mat_t mat,
    const acb_mat_t tau, slong prec);
slong acb_theta_transform_kappa2(const fmpz_mat_t mat);
void acb_theta_transform_proj(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th,
    int sqr, slong prec);

/* Toolbox for derivatives */

slong acb_theta_jet_nb(slong ord, slong g);
slong acb_theta_jet_total_order(const slong * tup, slong g);
void acb_theta_jet_tuples(slong * tups, slong ord, slong g);
slong acb_theta_jet_index(const slong * tup, slong g);

void acb_theta_jet_mul(acb_ptr res, acb_srcptr v1, acb_srcptr v2, slong ord,
    slong g, slong prec);
void acb_theta_jet_compose(acb_ptr res, acb_srcptr v, const acb_mat_t N,
    slong ord, slong prec);
void acb_theta_jet_exp_pi_i(acb_ptr res, arb_srcptr a, slong ord, slong g, slong prec);
void acb_theta_jet_exp_qf(acb_ptr res, acb_srcptr z, const acb_mat_t N, slong ord, slong prec);

void acb_theta_jet_error_bounds(arb_ptr err, acb_srcptr z, const acb_mat_t tau,
    acb_srcptr dth, slong ord, slong prec);

/* Error bounds in summation algorithms */

void acb_theta_naive_radius(arf_t R2, arf_t eps, const arb_mat_t C, slong ord, slong prec);
void acb_theta_jet_naive_radius(arf_t R2, arf_t eps, const arb_mat_t C, arb_srcptr v,
    slong ord, slong prec);
void acb_theta_jet_ql_bounds(arb_t c, arb_t rho, acb_srcptr z, const acb_mat_t tau, slong ord);
void acb_theta_jet_ql_radius(arf_t eps, arf_t err, const arb_t c, const arb_t rho,
    slong ord, slong g, slong prec);

void acb_theta_naive_term(acb_t res, acb_srcptr z, const acb_mat_t tau, const slong * tup,
    const slong * n, slong prec);

/* Context structures for theta function evaluations */

void acb_theta_ctx_exp_inv(acb_t exp_inv, const acb_t exp, const acb_t x, int is_real, slong prec);
void acb_theta_ctx_sqr_inv(acb_t sqr_inv, const acb_t inv, const acb_t sqr, int is_real, slong prec);

#define acb_theta_ctx_g(ctx) (acb_mat_nrows(&(ctx)->tau))
#define acb_theta_ctx_tau(ctx) (&(ctx)->tau)
#define acb_theta_ctx_y(ctx) (&(ctx)->Y)
#define acb_theta_ctx_yinv(ctx) (&(ctx)->Yinv)
#define acb_theta_ctx_exp_tau_div_4(ctx) ((ctx)->exp_tau_div_4)
#define acb_theta_ctx_exp_tau_div_2(ctx) ((ctx)->exp_tau_div_2)
#define acb_theta_ctx_exp_tau(ctx) ((ctx)->exp_tau)
#define acb_theta_ctx_cho(ctx) (&(ctx)->C)
#define acb_theta_ctx_choinv(ctx) (&(ctx)->Cinv)
#define acb_theta_ctx_exp_tau_div_4_inv(ctx) ((ctx)->exp_tau_div_4_inv)
#define acb_theta_ctx_exp_tau_div_2_inv(ctx) ((ctx)->exp_tau_div_2_inv)
#define acb_theta_ctx_exp_tau_inv(ctx) ((ctx)->exp_tau_inv)
#define acb_theta_ctx_exp_tau_a_div_2(ctx, a) ((ctx)->exp_tau_a_div_2 + (a) * g)
#define acb_theta_ctx_exp_tau_a(ctx, a) ((ctx)->exp_tau_a + (a) * g)
#define acb_theta_ctx_exp_tau_a_div_2_inv(ctx, a) ((ctx)->exp_tau_a_div_2_inv + (a) * g)
#define acb_theta_ctx_exp_tau_a_inv(ctx, a) ((ctx)->exp_tau_a_inv + (a) * g)
#define acb_theta_ctx_exp_a_tau_a_div_4(ctx, a) ((ctx)->exp_a_tau_a_div_4 + (a))

void acb_theta_ctx_tau_init(acb_theta_ctx_tau_t ctx, slong g);
void acb_theta_ctx_tau_clear(acb_theta_ctx_tau_t ctx);
void acb_theta_ctx_tau_set(acb_theta_ctx_tau_t ctx, const acb_mat_t tau, slong prec);
void acb_theta_ctx_tau_copy(acb_theta_ctx_tau_t res, const acb_theta_ctx_tau_t ctx);
int acb_theta_ctx_tau_overlaps(const acb_theta_ctx_tau_t ctx1, const acb_theta_ctx_tau_t ctx2);
void acb_theta_ctx_tau_dupl(acb_theta_ctx_tau_t ctx, slong prec);

#define acb_theta_ctx_z(ctx) ((ctx)->z)
#define acb_theta_ctx_exp_z(ctx) ((ctx)->exp_z)
#define acb_theta_ctx_exp_z_inv(ctx) ((ctx)->exp_z_inv)
#define acb_theta_ctx_exp_2z(ctx) ((ctx)->exp_2z)
#define acb_theta_ctx_exp_2z_inv(ctx) ((ctx)->exp_2z_inv)
#define acb_theta_ctx_c(ctx) (&(ctx)->c)
#define acb_theta_ctx_u(ctx) (&(ctx)->u)
#define acb_theta_ctx_uinv(ctx) (&(ctx)->uinv)
#define acb_theta_ctx_r(ctx) ((ctx)->r)
#define acb_theta_ctx_v(ctx) ((ctx)->v)
#define acb_theta_ctx_is_real(ctx) ((ctx)->is_real)

void acb_theta_ctx_z_init(acb_theta_ctx_z_t ctx, slong g);
void acb_theta_ctx_z_clear(acb_theta_ctx_z_t ctx);
acb_theta_ctx_z_struct * acb_theta_ctx_z_vec_init(slong nb, slong g);
void acb_theta_ctx_z_vec_clear(acb_theta_ctx_z_struct * vec, slong nb);
void acb_theta_ctx_z_set(acb_theta_ctx_z_t ctx, acb_srcptr z, const acb_theta_ctx_tau_t ctx_tau, slong prec);
void acb_theta_ctx_z_copy(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx);
int acb_theta_ctx_z_overlaps(const acb_theta_ctx_z_t ctx1, const acb_theta_ctx_z_t ctx2);

void acb_theta_ctx_z_add_real(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx,
    const acb_theta_ctx_z_t ctx_real, slong prec);
void acb_theta_ctx_z_dupl(acb_theta_ctx_z_t ctx, slong prec);
void acb_theta_ctx_z_shift_a0(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx,
    const acb_theta_ctx_tau_t ctx_tau, ulong a, slong prec);
void acb_theta_ctx_z_common_v(arb_ptr v, const acb_theta_ctx_z_struct * vec, slong nb, slong prec);

/* Summation algorithms */

void acb_theta_sum_work(acb_ptr th, slong len, acb_srcptr exp_zs, acb_srcptr exp_zs_inv,
    slong nb, const acb_mat_t exp_tau, const acb_mat_t exp_tau_inv, const acb_theta_eld_t E,
    slong ord, slong prec, acb_theta_sum_worker_t worker);

void acb_theta_sum_00_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec);
void acb_theta_sum_0b_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec);
void acb_theta_sum_jet_00_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec);
void acb_theta_sum_jet_all_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec);

void acb_theta_sum_00(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong prec);
void acb_theta_sum_0b(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong prec);
void acb_theta_sum_a0_tilde(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, slong prec);
void acb_theta_sum_all_tilde(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, slong prec);
void acb_theta_sum_jet_00(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong ord, slong prec);
void acb_theta_sum_jet_all(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong ord, slong prec);

/* Quasilinear algorithms on exact input */

int acb_theta_ql_nb_steps(slong * pattern, const arb_mat_t cho, slong prec);

int acb_theta_ql_lower_dim(acb_ptr * new_zs, acb_ptr * cofactors, slong ** pts,
    slong * nb, arf_t err, slong * fullprec, acb_srcptr z, const acb_mat_t tau,
    arb_srcptr distances, slong s, ulong a, slong prec);
void acb_theta_ql_recombine(acb_ptr th, acb_srcptr th0, acb_srcptr cofactors,
    const slong * pts, slong nb, const arf_t err, slong fullprec,
    slong s, ulong a, int all, slong g, slong prec);
int acb_theta_ql_setup(acb_ptr rts, acb_ptr rts0, acb_ptr t, slong * guard, slong * easy_steps,
    acb_srcptr zs, slong nb, const acb_mat_t tau, arb_srcptr distances,
    slong nb_steps, int all, slong prec);
void acb_theta_ql_steps(acb_ptr th, acb_ptr th_init, acb_srcptr rts,
    acb_srcptr rts_all, slong nb, slong nb_steps, arb_srcptr distances,
    const slong * easy_steps, int all, slong g, slong prec);
int acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    const slong * pattern, int all, int shifted_prec, slong prec);

void acb_theta_jet_ql_finite_diff(acb_ptr dth, const arf_t eps, const arf_t err,
    const arb_t rho, acb_srcptr val, slong ord, slong g, slong prec);

/* Main functions */

int acb_theta_reduce_tau(acb_ptr new_zs, acb_mat_t new_tau, fmpz_mat_t mat, acb_mat_t N,
    acb_mat_t ct, acb_ptr exps, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec);
int acb_theta_reduce_z(acb_ptr new_zs, arb_ptr rs, acb_ptr cs, acb_srcptr zs,
    slong nb, const acb_mat_t tau, slong prec);

void acb_theta_00_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong prec);
void acb_theta_one_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, ulong ab, slong prec);
void acb_theta_all_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, int sqr, slong prec);
void acb_theta_jet_00_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, slong prec);
void acb_theta_jet_one_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, ulong ab, slong prec);
void acb_theta_jet_all_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, slong prec);

void acb_theta_00(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong prec);
void acb_theta_all(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, int sqr, slong prec);
void acb_theta_jet_00(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, slong prec);
void acb_theta_jet_all(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, slong prec);

/* Genus 2 specifics */

#define ACB_THETA_G2_COV_NB 26

void acb_theta_g2_detk_symj(acb_poly_t res, const acb_mat_t m, const acb_poly_t f,
    slong k, slong j, slong prec);
void acb_theta_g2_transvectant(acb_poly_t res, const acb_poly_t g, const acb_poly_t h,
    slong m, slong n, slong k, slong prec);
void acb_theta_g2_transvectant_lead(acb_t r, const acb_poly_t g, const acb_poly_t h,
    slong m, slong n, slong k, slong prec);

slong acb_theta_g2_character(const fmpz_mat_t mat);

void acb_theta_g2_psi4(acb_t res, acb_srcptr th2, slong prec);
void acb_theta_g2_psi6(acb_t res, acb_srcptr th2, slong prec);
void acb_theta_g2_chi10(acb_t res, acb_srcptr th2, slong prec);
void acb_theta_g2_chi12(acb_t res, acb_srcptr th2, slong prec);
void acb_theta_g2_chi5(acb_t res, acb_srcptr th, slong prec);
void acb_theta_g2_chi35(acb_t res, acb_srcptr th, slong prec);
void acb_theta_g2_chi3_6(acb_poly_t res, acb_srcptr dth, slong prec);

void acb_theta_g2_sextic(acb_poly_t res, const acb_mat_t tau, slong prec);
void acb_theta_g2_sextic_chi5(acb_poly_t res, acb_t chi5, const acb_mat_t tau, slong prec);
void acb_theta_g2_covariants(acb_poly_struct * res, const acb_poly_t f, slong prec);
void acb_theta_g2_covariants_lead(acb_ptr res, const acb_poly_t f, slong prec);

#ifdef __cplusplus
}
#endif

#endif
