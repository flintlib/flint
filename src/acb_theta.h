/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_THETA_H
#define ACB_THETA_H

#include "fmpz_mat.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ACB_THETA_LOW_PREC 32

/* The Siegel modular group */

static inline slong
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

/* Ellipsoids in naive algorithms */

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

/* Naive algorithms */

void acb_theta_naive_radius(arf_t R2, arf_t eps, const arb_mat_t C, slong ord, slong prec);
void acb_theta_naive_reduce(arb_ptr v, acb_ptr new_zs, arb_ptr as, acb_ptr cs, arb_ptr us,
    acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec);
void acb_theta_naive_term(acb_t res, acb_srcptr z, const acb_mat_t tau, const slong * tup,
    const slong * n, slong prec);

typedef void (*acb_theta_naive_worker_t)(acb_ptr, acb_srcptr, acb_srcptr, const slong *,
    slong, const acb_t, const slong *, slong, slong, slong, slong);

void acb_theta_naive_worker(acb_ptr th, slong len, acb_srcptr zs, slong nb,
    const acb_mat_t tau, const acb_theta_eld_t E, slong ord, slong prec,
    acb_theta_naive_worker_t worker);

void acb_theta_naive_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec);
void acb_theta_naive_0b(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec);

void acb_theta_naive_fixed_ab(acb_ptr th, ulong ab, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong prec);
void acb_theta_naive_fixed_a(acb_ptr th, ulong a, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong prec);
void acb_theta_naive_all(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec);

/* Naive algorithms for derivatives */

slong acb_theta_jet_nb(slong ord, slong g);
slong acb_theta_jet_total_order(const slong * tup, slong g);
void acb_theta_jet_tuples(slong * tups, slong ord, slong g);
slong acb_theta_jet_index(const slong * tup, slong g);

void acb_theta_jet_mul(acb_ptr res, acb_srcptr v1, acb_srcptr v2, slong ord,
    slong g, slong prec);
void acb_theta_jet_compose(acb_ptr res, acb_srcptr v, const acb_mat_t N,
    slong ord, slong prec);
void acb_theta_jet_exp_pi_i(acb_ptr res, arb_srcptr a, slong ord, slong g, slong prec);

void acb_theta_jet_naive_radius(arf_t R2, arf_t eps, const arb_mat_t C, arb_srcptr v,
    slong ord, slong prec);

void acb_theta_jet_naive_00(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec);
void acb_theta_jet_naive_fixed_ab(acb_ptr dth, ulong ab, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec);
void acb_theta_jet_naive_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec);

void acb_theta_jet_error_bounds(arb_ptr err, acb_srcptr z, const acb_mat_t tau,
    acb_srcptr dth, slong ord, slong prec);

/* Quasi-linear algorithms on the reduced domain */

void acb_theta_dist_pt(arb_t d, arb_srcptr v, const arb_mat_t C, const slong * n, slong prec);
void acb_theta_dist_lat(arb_t d, arb_srcptr v, const arb_mat_t C, slong prec);
void acb_theta_dist_a0(arb_ptr d, acb_srcptr z, const acb_mat_t tau, slong prec);
slong acb_theta_dist_addprec(const arb_t d);

void acb_theta_agm_hadamard(acb_ptr res, acb_srcptr a, slong g, slong prec);
void acb_theta_agm_sqrt(acb_ptr res, acb_srcptr a, acb_srcptr roots, slong nb, slong prec);
void acb_theta_agm_mul(acb_ptr res, acb_srcptr a1, acb_srcptr a2, slong g, slong prec);
void acb_theta_agm_mul_tight(acb_ptr res, acb_srcptr a0, acb_srcptr a,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec);

typedef int (*acb_theta_ql_worker_t)(acb_ptr, acb_srcptr, acb_srcptr,
    arb_srcptr, arb_srcptr, const acb_mat_t, slong, slong);

int acb_theta_ql_a0_naive(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong guard, slong prec);
int acb_theta_ql_a0_split(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d,
    const acb_mat_t tau, slong s, slong guard, slong prec, acb_theta_ql_worker_t worker);
int acb_theta_ql_a0_steps(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong nb_steps, slong s, slong guard,
    slong prec, acb_theta_ql_worker_t worker);
slong acb_theta_ql_a0_nb_steps(const arb_mat_t C, slong s, slong prec);
int acb_theta_ql_a0(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong guard, slong prec);

slong acb_theta_ql_reduce(acb_ptr x, acb_t c, arb_t u, slong * n1, acb_srcptr z,
    const acb_mat_t tau, slong prec);

void acb_theta_ql_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, int sqr, slong prec);

/* Quasi-linear algorithms for derivatives */

void acb_theta_jet_ql_bounds(arb_t c, arb_t rho, acb_srcptr z, const acb_mat_t tau, slong ord);
void acb_theta_jet_ql_radius(arf_t eps, arf_t err, const arb_t c, const arb_t rho,
    slong ord, slong g, slong prec);
void acb_theta_jet_ql_finite_diff(acb_ptr dth, const arf_t eps, const arf_t err,
    const arb_t rho, acb_srcptr val, slong ord, slong g, slong prec);

void acb_theta_jet_ql_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec);

/* Transformation formulas */

ulong acb_theta_transform_char(slong * e, const fmpz_mat_t mat, ulong ab);
void acb_theta_transform_sqrtdet(acb_t res, const acb_mat_t tau, slong prec);
slong acb_theta_transform_kappa(acb_t sqrtdet, const fmpz_mat_t mat,
    const acb_mat_t tau, slong prec);
slong acb_theta_transform_kappa2(const fmpz_mat_t mat);
void acb_theta_transform_proj(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th,
    int sqr, slong prec);

void acb_theta_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, int sqr, slong prec);
void acb_theta_jet_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec);

/* Genus 2 specifics */

#define ACB_THETA_G2_COV_NB 26

void acb_theta_g2_jet_naive_1(acb_ptr dth, const acb_mat_t tau, slong prec);
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
