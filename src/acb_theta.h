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
#include <math.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "arb.h"
#include "acb.h"
#include "acb_poly.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_modular.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ACB_THETA_LOW_PREC 32

/* The Siegel modular group */

static __inline__ slong
sp2gz_dim(const fmpz_mat_t mat)
{
    return fmpz_mat_nrows(mat) / 2;
}

void sp2gz_get_a(fmpz_mat_t res, const fmpz_mat_t mat);
void sp2gz_get_b(fmpz_mat_t res, const fmpz_mat_t mat);
void sp2gz_get_c(fmpz_mat_t res, const fmpz_mat_t mat);
void sp2gz_get_d(fmpz_mat_t res, const fmpz_mat_t mat);
void sp2gz_set_abcd(fmpz_mat_t mat, const fmpz_mat_t a, const fmpz_mat_t b,
    const fmpz_mat_t c, const fmpz_mat_t d);

void sp2gz_j(fmpz_mat_t mat);
void sp2gz_block_diag(fmpz_mat_t mat, const fmpz_mat_t U);
void sp2gz_trig(fmpz_mat_t mat, const fmpz_mat_t S);
slong sp2gz_nb_fundamental(slong g);
void sp2gz_fundamental(fmpz_mat_t mat, slong j);

void sp2gz_inv(fmpz_mat_t inv, const fmpz_mat_t mat);
int sp2gz_is_correct(const fmpz_mat_t mat);
void sp2gz_randtest(fmpz_mat_t mat, flint_rand_t state, slong bits);

/* The Siegel half space */

void acb_siegel_cocycle(acb_mat_t res, const fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_cocycle_det(acb_t det, const fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_cocycle_sqrtdet(acb_t r, const fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_transform(acb_mat_t res, const fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_transform_ext(acb_ptr r, acb_mat_t w, const fmpz_mat_t mat,
    acb_srcptr z, const acb_mat_t tau, slong prec);

void acb_siegel_reduce_imag(fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_reduce_real(fmpz_mat_t mat, const acb_mat_t tau, slong prec);
void acb_siegel_reduce(acb_mat_t res, fmpz_mat_t mat, const acb_mat_t tau, slong prec);

void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits);
void acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits);
void acb_siegel_randtest_nice(acb_mat_t tau, flint_rand_t state, slong prec);

/* Theta characteristics */

void acb_theta_char_get_slong(slong* n, ulong a, slong g);
ulong acb_theta_char_get_a(slong* n, slong g);
void acb_theta_char_get_acb(acb_ptr v, ulong a, slong g);
void acb_theta_char_get_arb(arb_ptr v, ulong a, slong g);

slong acb_theta_char_dot(ulong a, ulong b, slong g);
slong acb_theta_char_dot_slong(ulong a, slong* n, slong g);
void acb_theta_char_dot_acb(acb_t x, ulong a, acb_srcptr z, slong g, slong prec);

/* Ellipsoids in naive algorithms */

struct acb_theta_eld_struct
{
    slong dim;
    slong ambient_dim;
    slong* last_coords;
    slong min, mid, max, nr, nl;
    struct acb_theta_eld_struct* rchildren;
    struct acb_theta_eld_struct* lchildren;
    slong nb_pts, nb_border;
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
#define acb_theta_eld_nb_border(E) ((E)->nb_border)
#define acb_theta_eld_box(E, k) ((E)->box[(k)])

void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g);
void acb_theta_eld_clear(acb_theta_eld_t E);

void acb_theta_eld_interval(slong* min, slong* mid, slong* max,
    const arb_t ctr, const arf_t rad, slong prec);
void acb_theta_eld_cho(arb_mat_t cho, const acb_mat_t tau, slong prec);
void acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t cho, const arf_t R2,
    arb_srcptr offset, slong prec);
void acb_theta_eld_points(slong* pts, const acb_theta_eld_t E);
void acb_theta_eld_border(slong* pts, const acb_theta_eld_t E);
int acb_theta_eld_contains(const acb_theta_eld_t E, slong* pt);
void acb_theta_eld_print(const acb_theta_eld_t E);

/* Precomputations in naive algorithms */

typedef struct
{
    slong dim;
    acb_mat_struct exp_mat;
    acb_ptr sqr_powers;
    slong* indices;
    acb_ptr exp_z;
    slong nb_z;
} acb_theta_precomp_struct;

typedef acb_theta_precomp_struct acb_theta_precomp_t[1];

#define acb_theta_precomp_dim(D) ((D)->dim)
#define acb_theta_precomp_exp_mat(D) (&(D)->exp_mat)
#define acb_theta_precomp_sqr_pow(D, k, j) (&(D)->sqr_powers[(j) + (D)->indices[(k)]])
#define acb_theta_precomp_exp_z(D, k, j) (&(D)->exp_z[(k) * (D)->dim + (j)])
#define acb_theta_precomp_nb_z(D) ((D)->nb_z)

void acb_theta_precomp_init(acb_theta_precomp_t D, slong nb_z, slong g);
void acb_theta_precomp_clear(acb_theta_precomp_t D);
void acb_theta_precomp_set(acb_theta_precomp_t D, acb_srcptr z,
    const acb_mat_t tau, const acb_theta_eld_t E, slong prec);

/* Naive algorithms */

void acb_theta_naive_term(acb_t res, acb_srcptr z, const acb_mat_t tau,
    slong* n, slong prec);
void acb_theta_naive_tail(arf_t bound, const arf_t R2, const arb_mat_t cho,
    slong ord, slong prec);
void acb_theta_naive_radius(arf_t R2, const arb_mat_t cho, slong ord,
    const arf_t eps, slong prec);
void acb_theta_naive_reduce(arb_ptr offset, acb_ptr new_z, acb_ptr c, arb_ptr u,
    acb_srcptr z, slong nb_z, const acb_mat_t tau, const arb_mat_t cho, slong prec);
void acb_theta_naive_ellipsoid(acb_theta_eld_t E, acb_ptr new_z, acb_ptr c,
    arb_ptr u, slong ord, acb_srcptr z, slong nb_z, const acb_mat_t tau,
    const arf_t eps, slong prec);
slong acb_theta_naive_fullprec(const acb_theta_eld_t E, slong prec);

typedef void (*acb_theta_naive_worker_t)(acb_ptr, const acb_t, slong*, slong,
    slong, slong, slong);

void acb_theta_naive_worker(acb_ptr th, slong nb, const acb_t c, const arb_t u,
    const acb_theta_eld_t E, const acb_theta_precomp_t D, slong k, slong ord,
    slong prec, acb_theta_naive_worker_t worker_dim0);

void acb_theta_naive_00(acb_ptr th, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec);
void acb_theta_naive_0b(acb_ptr th, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec);

void acb_theta_naive_ind(acb_ptr th, ulong ab, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec);
void acb_theta_naive_fixed_a(acb_ptr th, ulong a, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec);
void acb_theta_naive_all(acb_ptr th, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec);

/* Quasi-linear algorithms on the reduced domain */

void acb_theta_dist_pt(arb_t d2, arb_srcptr v, const arb_mat_t cho, slong* pt, slong prec);
void acb_theta_dist_lat(arb_t d2, arb_srcptr v, const arb_mat_t cho, slong prec);
void acb_theta_dist_a0(arb_ptr dist, acb_srcptr z, const acb_mat_t tau, slong prec);
slong acb_theta_dist_addprec(const arb_t d2);

void acb_theta_agm_hadamard(acb_ptr r, acb_srcptr a, slong g, slong prec);
void acb_theta_agm_sqrt(acb_ptr r, acb_srcptr a, acb_srcptr roots, slong nb, slong prec);
void acb_theta_agm_mul(acb_ptr r, acb_srcptr a1, acb_srcptr a2, slong g, slong prec);
void acb_theta_agm_rel_mag_err(arf_t m, arf_t eps, acb_srcptr a,
    arb_srcptr dist, slong n, slong prec);
void acb_theta_agm_mul_tight(acb_ptr r, acb_srcptr a0, acb_srcptr a,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec);

#define ACB_THETA_QL_SPLIT 2
#define ACB_THETA_QL_TRY 100
/* See also acb_theta_ql_nb_steps for more tuning */

slong acb_theta_ql_nb_steps(const arb_mat_t cho, slong d, slong prec);
void acb_theta_ql_log_rescale(acb_t f, acb_srcptr z, const acb_mat_t tau, slong prec);
int acb_theta_ql_roots(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong nb_steps, slong guard, slong prec);
void acb_theta_ql_step_1(acb_ptr r, acb_srcptr th0, acb_srcptr th,
    acb_srcptr roots, arb_srcptr dist0, arb_srcptr dist, slong g, slong prec);
void acb_theta_ql_step_3(acb_ptr r, acb_srcptr th0, acb_srcptr th,
    acb_srcptr roots, arb_srcptr dist0, arb_srcptr dist, slong g, slong prec);

/* Use as worker(r, t, z, dist0, dist, tau, guard, prec). Computes theta_{a,0}
   at 0, t, 2t, z, z + t, z + 2t (less values if z = 0 or t = 0 or both) */

typedef int (*acb_theta_ql_worker_t)(acb_ptr, acb_srcptr, acb_srcptr,
    arb_srcptr, arb_srcptr, const acb_mat_t, slong, slong);

int acb_theta_ql_a0_naive(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec);
int acb_theta_ql_a0_split(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong d, slong guard, slong prec, acb_theta_ql_worker_t worker);
int acb_theta_ql_a0_steps(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec,
    acb_theta_ql_worker_t worker);
int acb_theta_ql_a0(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec);

void acb_theta_duplication(acb_ptr th2, acb_srcptr th0, acb_srcptr th,
    arb_srcptr dist0, arb_srcptr dist, slong g, slong prec);
void acb_theta_ql_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec);
void acb_theta_ql_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec);

/* Transformation formulas for theta functions */

ulong acb_theta_transform_char(fmpz_t eps, ulong ab, const fmpz_mat_t mat);
void acb_theta_transform_proj(acb_ptr res, acb_srcptr th, const fmpz_mat_t mat, slong prec);
void acb_theta_transform_sqr_proj(acb_ptr res, acb_srcptr th2, const fmpz_mat_t mat, slong prec);

slong acb_theta_transform_kappa(const fmpz_mat_t mat);
slong acb_theta_transform_k2(const fmpz_mat_t mat);
void acb_theta_transform(acb_ptr res, acb_srcptr th, acb_srcptr z,
    const acb_mat_t tau, const fmpz_mat_t mat, slong kappa, slong prec);
void acb_theta_transform_sqr(acb_ptr res, acb_srcptr th2, acb_srcptr z,
    const acb_mat_t tau, const fmpz_mat_t mat, slong k2, slong prec);

/* Main functions */

void acb_theta_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec);
void acb_theta_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec);

/* Derivatives/jets */

slong acb_theta_jet_nb(slong ord, slong g);
void acb_theta_jet_orders(slong* orders, slong ord, slong g);
slong acb_theta_jet_index(const slong* orders, slong g);

void acb_theta_jet_bounds(arb_t eps, arb_t c, arb_t rho, acb_srcptr z,
    const acb_mat_t tau, slong ord, slong hprec, slong prec);
void acb_theta_jet_fourier(acb_ptr res, acb_srcptr val, slong ord, slong g, slong prec);
void acb_theta_jet_fd(acb_ptr dth, const arb_t eps, const arb_t c,
    const arb_t rho, acb_srcptr val, slong ord, slong g, slong prec);

void acb_theta_jet_naive_all(acb_ptr dth, slong ord, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec);
void acb_theta_jet_all(acb_ptr dth, slong ord, acb_srcptr z, const acb_mat_t tau, slong prec);

#ifdef __cplusplus
}
#endif

#endif
