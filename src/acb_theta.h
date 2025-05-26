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

#include "fmpz_types.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ACB_THETA_LOW_PREC 32

/* The Siegel modular group */

FLINT_FORCE_INLINE slong
sp2gz_dim(const fmpz_mat_t mat)
{
    return (mat->r) / 2;
}

void sp2gz_set_blocks(fmpz_mat_t mat, const fmpz_mat_t alpha, const fmpz_mat_t beta,
    const fmpz_mat_t gamma, const fmpz_mat_t delta);
void sp2gz_j(fmpz_mat_t mat);
void sp2gz_block_diag(fmpz_mat_t mat, const fmpz_mat_t U);
void sp2gz_trig(fmpz_mat_t mat, const fmpz_mat_t S);
void sp2gz_embed(fmpz_mat_t res, const fmpz_mat_t mat);
void sp2gz_restrict(fmpz_mat_t res, const fmpz_mat_t mat);

FLINT_FORCE_INLINE slong
sp2gz_nb_fundamental(slong g)
{
    if (g == 1)
        return 1;
    if (g == 2)
        return 19;
    else
        return 19 * ((g * (g - 1)) / 2) + (1 << g);
}

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
void acb_siegel_cho_yinv(arb_mat_t cho, arb_mat_t yinv, const acb_mat_t tau, slong prec);

void acb_siegel_reduce(fmpz_mat_t mat, const acb_mat_t tau, slong prec);
int acb_siegel_is_reduced(const acb_mat_t tau, slong tol_exp, slong prec);

slong acb_siegel_kappa(acb_t sqrtdet, const fmpz_mat_t mat, const acb_mat_t tau, int sqr, slong prec);
slong acb_siegel_kappa2(const fmpz_mat_t mat);

void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits);
void acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits);
void acb_siegel_randtest_compact(acb_mat_t tau, flint_rand_t state, int exact, slong prec);
void acb_siegel_randtest_vec(acb_ptr z, flint_rand_t state, slong g, slong prec);
void acb_siegel_randtest_vec_reduced(acb_ptr zs, flint_rand_t state,
    slong nb, const acb_mat_t tau, int exact, slong prec);

/* Theta characteristics */

FLINT_FORCE_INLINE int
acb_theta_char_bit(ulong ch, slong j, slong n)
{
    return (ch >> (n - 1 - j)) & 1;
}

void acb_theta_char_get_arb(arb_ptr v, ulong a, slong g);
void acb_theta_char_get_acb(acb_ptr v, ulong a, slong g);
ulong acb_theta_char_set_slong_vec(const slong * vec, slong len);

slong acb_theta_char_dot(ulong a, ulong b, slong g);
slong acb_theta_char_dot_slong(ulong a, const slong * n, slong g);

FLINT_FORCE_INLINE int
acb_theta_char_is_even(ulong ab, slong g)
{
    return acb_theta_char_dot(ab >> g, ab, g) % 2 == 0;
}

void acb_theta_char_table(ulong * ch, slong * e, const fmpz_mat_t mat, ulong ab, int all);
void acb_theta_char_shuffle(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th,
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

/* Ellipsoids */

struct acb_theta_eld_struct
{
    slong dim, ambient_dim;
    slong * last_coords;
    slong min, mid, max, nr, nl;
    struct acb_theta_eld_struct * rchildren;
    struct acb_theta_eld_struct * lchildren;
    slong nb_pts, nb_border;
    slong * box;
};

typedef struct acb_theta_eld_struct acb_theta_eld_t[1];

void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g);
void acb_theta_eld_clear(acb_theta_eld_t E);

int acb_theta_eld_set(acb_theta_eld_t E, const arb_mat_t C, const arf_t R2, arb_srcptr v);

FLINT_FORCE_INLINE slong
acb_theta_eld_nb_pts(const acb_theta_eld_t E)
{
    return E->nb_pts;
}

void acb_theta_eld_points(slong * pts, const acb_theta_eld_t E);

FLINT_FORCE_INLINE slong
acb_theta_eld_box(const acb_theta_eld_t E, slong j)
{
    return E->box[j];
}

FLINT_FORCE_INLINE slong
acb_theta_eld_nb_border(const acb_theta_eld_t E)
{
    return E->nb_border;
}

void acb_theta_eld_border(slong * pts, const acb_theta_eld_t E);
int acb_theta_eld_contains(const acb_theta_eld_t E, const slong * pt);
void acb_theta_eld_print(const acb_theta_eld_t E);

void acb_theta_eld_distances(arb_ptr ds, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong prec);

/* Error bounds in summation algorithms */

void acb_theta_sum_radius(arf_t R2, arf_t eps, const arb_mat_t cho, slong ord, slong prec);
void acb_theta_sum_jet_radius(arf_t R2, arf_t eps, const arb_mat_t cho, arb_srcptr v,
    slong ord, slong prec);
void acb_theta_sum_term(acb_t res, acb_srcptr z, const acb_mat_t tau, const slong * tup,
    const slong * n, slong prec);
slong acb_theta_sum_addprec(const arb_t d);

/* Context structures in summation algorithms */

struct acb_theta_ctx_tau_struct
{
    slong g;
    int allow_shift;
    arb_mat_struct yinv;
    arb_mat_struct cho;

    acb_mat_t exp_tau_div_4;
    acb_mat_t exp_tau_div_2;
    acb_mat_t exp_tau;
    acb_mat_t exp_tau_div_4_inv;
    acb_mat_t exp_tau_div_2_inv;
    acb_mat_t exp_tau_inv;

    /* allow_shift only */
    acb_ptr exp_tau_a;
    acb_ptr exp_tau_a_inv;
    acb_ptr exp_a_tau_a_div_4;
};

typedef struct acb_theta_ctx_tau_struct acb_theta_ctx_tau_t[1];

typedef struct
{
    slong g;
    acb_ptr exp_z;
    acb_ptr exp_2z;
    acb_ptr exp_z_inv;
    acb_ptr exp_2z_inv;
    arb_ptr v;
    arb_struct u;
    arb_struct uinv;
    int is_real;
}
acb_theta_ctx_z_struct;

typedef acb_theta_ctx_z_struct acb_theta_ctx_z_t[1];

void acb_theta_ctx_tau_init(acb_theta_ctx_tau_t ctx, int allow_shift, slong g);
void acb_theta_ctx_tau_clear(acb_theta_ctx_tau_t ctx);
void acb_theta_ctx_z_init(acb_theta_ctx_z_t ctx, slong g);
void acb_theta_ctx_z_clear(acb_theta_ctx_z_t ctx);
acb_theta_ctx_z_struct * acb_theta_ctx_z_vec_init(slong nb, slong g);
void acb_theta_ctx_z_vec_clear(acb_theta_ctx_z_struct * vec, slong nb);

void acb_theta_ctx_exp_inv(acb_t exp_inv, const acb_t exp, const acb_t x, int is_real, slong prec);
void acb_theta_ctx_sqr_inv(acb_t sqr_inv, const acb_t inv, const acb_t sqr, int is_real, slong prec);

void acb_theta_ctx_tau_set(acb_theta_ctx_tau_t ctx, const acb_mat_t tau, slong prec);
void acb_theta_ctx_tau_dupl(acb_theta_ctx_tau_t ctx, slong prec);
int acb_theta_ctx_tau_overlaps(const acb_theta_ctx_tau_t ctx1, const acb_theta_ctx_tau_t ctx2);

void acb_theta_ctx_z_set(acb_theta_ctx_z_t ctx, acb_srcptr z, const acb_theta_ctx_tau_t ctx_tau, slong prec);
void acb_theta_ctx_z_dupl(acb_theta_ctx_z_t ctx, slong prec);
void acb_theta_ctx_z_add_real(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx,
    const acb_theta_ctx_z_t ctx_real, slong prec);
void acb_theta_ctx_z_common_v(arb_ptr v, const acb_theta_ctx_z_struct * vec, slong nb, slong prec);
int acb_theta_ctx_z_overlaps(const acb_theta_ctx_z_t ctx1, const acb_theta_ctx_z_t ctx2);

/* Summation algorithms */

typedef void (*acb_theta_sum_worker_t)(acb_ptr, acb_srcptr, acb_srcptr, const slong *,
    slong, const acb_t, const slong *, slong, slong, slong, slong);

void acb_theta_sum_sqr_pow(acb_ptr * sqr_pow, const acb_mat_t exp_tau, const acb_theta_eld_t E, slong prec);
void acb_theta_sum_work(acb_ptr th, slong len, acb_srcptr exp_z, acb_srcptr exp_z_inv,
    const acb_mat_t exp_tau, const acb_mat_t exp_tau_inv, const acb_ptr * sqr_pow,
    const acb_theta_eld_t E, slong ord, slong prec, acb_theta_sum_worker_t worker);
void acb_theta_sum(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, int all_a,
    int all_b, int tilde, slong prec);
void acb_theta_sum_jet(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong ord, int all_a, int all_b, slong prec);

/* AGM steps */

void acb_theta_agm_sqrt(acb_ptr res, acb_srcptr a, acb_srcptr roots, slong nb, slong prec);
void acb_theta_agm_mul(acb_ptr res, acb_srcptr a1, acb_srcptr a2, slong g, int all, slong prec);
void acb_theta_agm_mul_tight(acb_ptr res, acb_srcptr a0, acb_srcptr a,
    arb_srcptr d0, arb_srcptr d, slong g, int all, slong prec);

/* Quasilinear algorithms on reduced input */

int acb_theta_ql_nb_steps(slong * pattern, const acb_mat_t tau, int cst, slong prec);

int acb_theta_ql_lower_dim(acb_ptr * new_zs, acb_ptr * cofactors, slong ** pts,
    slong * nb, arf_t err, slong * fullprec, acb_srcptr z, const acb_mat_t tau,
    arb_srcptr distances, slong s, ulong a, slong prec);
void acb_theta_ql_recombine(acb_ptr th, acb_srcptr th0, acb_srcptr cofactors,
    const slong * pts, slong nb, const arf_t err, slong fullprec,
    slong s, ulong a, int all, slong g, slong prec);
int acb_theta_ql_setup(acb_ptr rts, acb_ptr rts_all, acb_ptr t, slong * guard, slong * easy_steps,
    acb_srcptr zs, slong nb, const acb_mat_t tau, arb_srcptr distances,
    slong nb_steps, int all, slong prec);
void acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    const slong * pattern, int all, int shifted_prec, slong prec);

void acb_theta_ql_local_bound(arb_t c, arb_t rho, acb_srcptr z, const acb_mat_t tau, slong ord);
void acb_theta_ql_jet_error(arb_ptr err, acb_srcptr z, const acb_mat_t tau,
    acb_srcptr dth, slong ord, slong prec);

void acb_theta_ql_jet_fd(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, int all, slong prec);
void acb_theta_ql_jet(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, int all, slong prec);

/* Reduction and main functions */

void acb_theta_jet_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, ulong ab, int all, int sqr, slong prec);

int acb_theta_reduce_tau(acb_ptr new_zs, acb_mat_t new_tau, fmpz_mat_t mat, acb_mat_t N,
    acb_mat_t ct, acb_ptr exps, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec);
int acb_theta_reduce_z(acb_ptr new_zs, arb_ptr rs, acb_ptr cs, acb_srcptr zs,
    slong nb, const acb_mat_t tau, slong prec);

void acb_theta_jet(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, ulong ab, int all, int sqr, slong prec);

FLINT_FORCE_INLINE void
acb_theta_one(acb_ptr th, acb_srcptr z, const acb_mat_t tau, ulong ab, slong prec)
{
    acb_theta_jet(th, z, 1, tau, 0, ab, 0, 0, prec);
}

FLINT_FORCE_INLINE void
acb_theta_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, int sqr, slong prec)
{
    acb_theta_jet(th, z, 1, tau, 0, 0, 1, sqr, prec);
}

/* Dimension 2 specifics */

void acb_theta_g2_detk_symj(acb_poly_t res, const acb_mat_t m, const acb_poly_t f,
    slong k, slong j, slong prec);
void acb_theta_g2_transvectant(acb_poly_t res, const acb_poly_t g, const acb_poly_t h,
    slong m, slong n, slong k, int lead, slong prec);
slong acb_theta_g2_character(const fmpz_mat_t mat);

void acb_theta_g2_even_weight(acb_t psi4, acb_t psi6, acb_t chi10, acb_t chi12,
    acb_srcptr th2, slong prec);
void acb_theta_g2_chi5(acb_t res, acb_srcptr th, slong prec);
void acb_theta_g2_chi35(acb_t res, acb_srcptr th, slong prec);
void acb_theta_g2_chi3_6(acb_poly_t res, acb_srcptr dth, slong prec);

void acb_theta_g2_sextic_chi5(acb_poly_t f, acb_t chi5, const acb_mat_t tau, slong prec);
void acb_theta_g2_covariants(acb_poly_struct * res, const acb_poly_t f, int lead, slong prec);

#ifdef __cplusplus
}
#endif

#endif
