/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_DIRICHLET_H
#define ACB_DIRICHLET_H

#ifdef ACB_DIRICHLET_INLINES_C
#define ACB_DIRICHLET_INLINE
#else
#define ACB_DIRICHLET_INLINE static inline
#endif

#include "acb.h"
#include "dirichlet.h"

#ifdef __cplusplus
extern "C" {
#endif

void acb_dirichlet_powsum_term(acb_ptr res, arb_t log_prev, ulong * prev,
    const acb_t s, ulong k, int integer, int critical_line, slong len, slong prec);

void acb_dirichlet_powsum_sieved(acb_ptr z, const acb_t s, ulong n, slong len, slong prec);
void acb_dirichlet_powsum_smooth(acb_ptr z, const acb_t s, ulong n, slong len, slong prec);

void acb_dirichlet_zeta_bound(mag_t res, const acb_t s);
void acb_dirichlet_zeta_deriv_bound(mag_t der1, mag_t der2, const acb_t s);
void acb_dirichlet_zeta_rs_f_coeffs(acb_ptr c, const arb_t p, slong N, slong prec);
void acb_dirichlet_zeta_rs_d_coeffs(arb_ptr d, const arb_t sigma, slong k, slong prec);
void acb_dirichlet_zeta_rs_bound(mag_t err, const acb_t s, slong K);
void acb_dirichlet_zeta_rs_r(acb_t res, const acb_t s, slong K, slong prec);
void acb_dirichlet_zeta_rs(acb_t res, const acb_t s, slong K, slong prec);
void acb_dirichlet_zeta(acb_t res, const acb_t s, slong prec);

void acb_dirichlet_zeta_jet_rs(acb_ptr res, const acb_t s, slong len, slong prec);
void acb_dirichlet_zeta_jet(acb_t res, const acb_t s, int deflate, slong len, slong prec);

void acb_dirichlet_hurwitz(acb_t res, const acb_t s, const acb_t a, slong prec);

void acb_dirichlet_lerch_phi_integral(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec);
void acb_dirichlet_lerch_phi_direct(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec);
void acb_dirichlet_lerch_phi(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec);

void acb_dirichlet_stieltjes(acb_t res, const fmpz_t n, const acb_t a, slong prec);

typedef struct
{
    acb_struct s;
    mag_struct err;
    acb_ptr coeffs;
    int deflate;
    slong A;
    slong N;
    slong K;
}
acb_dirichlet_hurwitz_precomp_struct;

typedef acb_dirichlet_hurwitz_precomp_struct acb_dirichlet_hurwitz_precomp_t[1];

void acb_dirichlet_hurwitz_precomp_init(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, int deflate, slong A, slong K, slong N, slong prec);
void acb_dirichlet_hurwitz_precomp_init_num(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, int deflate, double num_eval, slong prec);
void acb_dirichlet_hurwitz_precomp_clear(acb_dirichlet_hurwitz_precomp_t pre);
void acb_dirichlet_hurwitz_precomp_bound(mag_t res, const acb_t s, slong A, slong K, slong N);
void acb_dirichlet_hurwitz_precomp_eval(acb_t res, const acb_dirichlet_hurwitz_precomp_t pre, ulong p, ulong q, slong prec);
void acb_dirichlet_hurwitz_precomp_choose_param(ulong * A, ulong * K, ulong * N, const acb_t s, double num_eval, slong prec);

void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s,
    const signed char * chi, int mod, int reciprocal, slong prec);

void acb_dirichlet_eta(acb_t res, const acb_t s, slong prec);

void acb_dirichlet_xi(acb_t res, const acb_t s, slong prec);

void acb_dirichlet_pairing(acb_t res, const dirichlet_group_t G, ulong m, ulong n, slong prec);
void acb_dirichlet_pairing_char(acb_t res, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b, slong prec);

/* precompute roots of unity */
typedef struct
{
    ulong order;
    ulong reduced_order;
    acb_t z;
    slong size;
    slong depth;
    acb_ptr * Z;
    int use_pow;
}
acb_dirichlet_roots_struct;

typedef acb_dirichlet_roots_struct acb_dirichlet_roots_t[1];

void acb_dirichlet_roots_init(acb_dirichlet_roots_t t, ulong order, slong num, slong prec);
void acb_dirichlet_roots_clear(acb_dirichlet_roots_t t);
void acb_dirichlet_root(acb_t z, const acb_dirichlet_roots_t t, ulong n, slong prec);

void acb_dirichlet_chi(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, ulong n, slong prec);
void acb_dirichlet_chi_vec(acb_ptr v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv, slong prec);

void acb_dirichlet_arb_quadratic_powers(arb_ptr v, slong nv, const arb_t x, slong prec);
void acb_dirichlet_qseries_arb(acb_t res, acb_srcptr a, const arb_t x, slong len, slong prec);
void acb_dirichlet_qseries_arb_powers_naive(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_roots_t z, slong len, slong prec);
void acb_dirichlet_qseries_arb_powers_smallorder(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_roots_t z, slong len, slong prec);

ulong acb_dirichlet_theta_length_d(ulong q, double x, slong prec);
ulong acb_dirichlet_theta_length(ulong q, const arb_t x, slong prec);
void mag_tail_kexpk2_arb(mag_t res, const arb_t a, ulong n);

void _acb_dirichlet_theta_argument_at_arb(arb_t xt, ulong q, const arb_t t, slong prec);
void acb_dirichlet_theta_arb(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t t, slong prec);
void acb_dirichlet_ui_theta_arb(acb_t res, const dirichlet_group_t G, ulong a, const arb_t t, slong prec);

void acb_dirichlet_gauss_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_factor(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_order2(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_root_number_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_root_number(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_si_poly_evaluate(acb_t res, slong * v, slong len, const acb_t z, slong prec);

void acb_dirichlet_jacobi_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);
ulong jacobi_one_prime(ulong p, ulong e, ulong pe, ulong cond);
void acb_dirichlet_jacobi_sum_factor(acb_t res,  const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum_gauss(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum_ui(acb_t res, const dirichlet_group_t G, ulong a, ulong b, slong prec);

void acb_dirichlet_l_euler_product(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_l_hurwitz(acb_t res, const acb_t s, const acb_dirichlet_hurwitz_precomp_t precomp,
    const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_l(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s, const acb_dirichlet_hurwitz_precomp_t precomp, const dirichlet_group_t G, slong prec);

void acb_dirichlet_l_jet(acb_ptr res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec);

void acb_dirichlet_l_fmpq_afe(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_l_fmpq(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void _acb_dirichlet_l_series(acb_ptr res, acb_srcptr s, slong slen,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    int deflate, slong len, slong prec);

void acb_dirichlet_l_series(acb_poly_t res, const acb_poly_t s,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    int deflate, slong len, slong prec);

void acb_dirichlet_hardy_theta(acb_ptr res, const acb_t t,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    slong len, slong prec);

void acb_dirichlet_hardy_z(acb_ptr res, const acb_t t,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    slong len, slong prec);

void _acb_dirichlet_hardy_theta_series(acb_ptr res, acb_srcptr s, slong slen, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec);
void acb_dirichlet_hardy_theta_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec);
void _acb_dirichlet_hardy_z_series(acb_ptr res, acb_srcptr s, slong slen, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec);
void acb_dirichlet_hardy_z_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec);

void acb_dirichlet_gram_point(arb_t res, const fmpz_t n, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
ulong acb_dirichlet_turing_method_bound(const fmpz_t p);
int _acb_dirichlet_definite_hardy_z(arb_t res, const arf_t t, slong *pprec);
void _acb_dirichlet_isolate_gram_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n);
void _acb_dirichlet_isolate_rosser_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n);
void _acb_dirichlet_isolate_turing_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n);
void acb_dirichlet_isolate_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n);
void _acb_dirichlet_refine_hardy_z_zero(arb_t res, const arf_t a, const arf_t b, slong prec);
void acb_dirichlet_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, slong prec);
void acb_dirichlet_zeta_zeros(acb_ptr res, const fmpz_t n, slong len, slong prec);
slong acb_dirichlet_platt_zeta_zeros(acb_ptr res, const fmpz_t n, slong len, slong prec);
void _acb_dirichlet_exact_zeta_nzeros(fmpz_t res, const arf_t t);
void acb_dirichlet_zeta_nzeros(arb_t res, const arb_t t, slong prec);
void acb_dirichlet_backlund_s(arb_t res, const arb_t t, slong prec);
void acb_dirichlet_backlund_s_bound(mag_t res, const arb_t t);
void acb_dirichlet_zeta_nzeros_gram(fmpz_t res, const fmpz_t n);
slong acb_dirichlet_backlund_s_gram(const fmpz_t n);

ACB_DIRICHLET_INLINE void
acb_dirichlet_hardy_z_zero(arb_t res, const fmpz_t n, slong prec)
{
    acb_dirichlet_hardy_z_zeros(res, n, 1, prec);
}

ACB_DIRICHLET_INLINE void
acb_dirichlet_zeta_zero(acb_t res, const fmpz_t n, slong prec)
{
    acb_dirichlet_zeta_zeros(res, n, 1, prec);
}

/* Platt zeta zeros */

typedef struct
{
    slong len;
    arb_ptr p;
    arb_struct Xa;
    arb_struct Xb;
}
acb_dirichlet_platt_c_precomp_struct;
typedef acb_dirichlet_platt_c_precomp_struct acb_dirichlet_platt_c_precomp_t[1];

typedef struct
{
    arb_struct c1;
    arb_struct c2;
}
acb_dirichlet_platt_i_precomp_struct;
typedef acb_dirichlet_platt_i_precomp_struct acb_dirichlet_platt_i_precomp_t[1];

typedef struct
{
    acb_dirichlet_platt_c_precomp_struct pre_c;
    acb_dirichlet_platt_i_precomp_struct pre_i;
}
acb_dirichlet_platt_ws_precomp_struct;
typedef acb_dirichlet_platt_ws_precomp_struct acb_dirichlet_platt_ws_precomp_t[1];

/* Platt C bound */

void acb_dirichlet_platt_c_precomp_init(acb_dirichlet_platt_c_precomp_t pre,
    slong sigma, const arb_t h, ulong k, slong prec);
void acb_dirichlet_platt_c_precomp_clear(acb_dirichlet_platt_c_precomp_t pre);
void acb_dirichlet_platt_c_bound_precomp(arb_t res,
    const acb_dirichlet_platt_c_precomp_t pre, slong sigma, const arb_t t0,
    const arb_t h, slong k, slong prec);
void acb_dirichlet_platt_c_bound(arb_t res,
    slong sigma, const arb_t t0, const arb_t h, slong k, slong prec);

/* Platt I bound */

void acb_dirichlet_platt_i_precomp_init(acb_dirichlet_platt_i_precomp_t pre,
        slong A, const arb_t H, slong sigma, slong prec);
void acb_dirichlet_platt_i_precomp_clear(acb_dirichlet_platt_i_precomp_t pre);
void acb_dirichlet_platt_i_bound_precomp(arb_t res,
    const acb_dirichlet_platt_i_precomp_t pre_i,
    const acb_dirichlet_platt_c_precomp_t pre_c,
    const arb_t t0, slong A, const arb_t H, slong sigma, slong prec);
void acb_dirichlet_platt_i_bound(arb_t res,
    const arb_t t0, slong A, const arb_t H, slong sigma, slong prec);

/* Platt Gaussian-windowed Whittaker-Shannon interpolation */

void acb_dirichlet_platt_ws_precomp_init(acb_dirichlet_platt_ws_precomp_t pre,
    slong A, const arb_t H, slong sigma, slong prec);
void acb_dirichlet_platt_ws_precomp_clear(acb_dirichlet_platt_ws_precomp_t pre);
void acb_dirichlet_platt_ws_interpolation_precomp(arb_t res, arf_t deriv,
    const acb_dirichlet_platt_ws_precomp_t pre, const arb_t t0,
    arb_srcptr p, const fmpz_t T, slong A, slong B, slong Ns_max,
    const arb_t H, slong sigma, slong prec);
void acb_dirichlet_platt_ws_interpolation(arb_t res, arf_t deriv,
    const arb_t t0, arb_srcptr p, const fmpz_t T, slong A, slong B,
    slong Ns_max, const arb_t H, slong sigma, slong prec);
void acb_dirichlet_platt_bound_C3(arb_t res, const arb_t t0, slong A,
    const arb_t H, slong Ns, slong prec);

void acb_dirichlet_platt_scaled_lambda(arb_t res, const arb_t t, slong prec);
void acb_dirichlet_platt_scaled_lambda_vec(arb_ptr res, const fmpz_t T,
    slong A, slong B, slong prec);

/* Platt lemma bounds of errors in the DFT grid evaluation of scaled Lambda */

void acb_dirichlet_platt_beta(arb_t res, const arb_t t, slong prec);
void acb_dirichlet_platt_lemma_32(arb_t out, const arb_t h, const arb_t t0,
    const arb_t x, slong prec);
void acb_dirichlet_platt_lemma_A5(arb_t out, slong B, const arb_t h, slong k,
    slong prec);
void acb_dirichlet_platt_lemma_A7(arb_t out, slong sigma, const arb_t t0,
    const arb_t h, slong k, slong A, slong prec);
void acb_dirichlet_platt_lemma_A9(arb_t out, slong sigma, const arb_t t0,
    const arb_t h, slong A, slong prec);
void acb_dirichlet_platt_lemma_A11(arb_t out, const arb_t t0, const arb_t h,
    slong B, slong prec);
void acb_dirichlet_platt_lemma_B1(arb_t out, slong sigma, const arb_t t0,
    const arb_t h, const fmpz_t J, slong prec);
void acb_dirichlet_platt_lemma_B2(arb_t out, slong K, const arb_t h,
    const arb_t xi, slong prec);

/* Platt DFT grid evaluation of scaled Lambda */

void acb_dirichlet_platt_multieval(arb_ptr out, const fmpz_t T,
    slong A, slong B, const arb_t h, const fmpz_t J, slong K,
    slong sigma, slong prec);
void acb_dirichlet_platt_multieval_threaded(arb_ptr out, const fmpz_t T,
    slong A, slong B, const arb_t h, const fmpz_t J, slong K,
    slong sigma, slong prec);

/* Platt Hardy Z zeros */

slong _acb_dirichlet_platt_local_hardy_z_zeros(
    arb_ptr res, const fmpz_t n, slong len,
    const fmpz_t T, slong A, slong B,
    const arb_t h, const fmpz_t J, slong K, slong sigma_grid,
    slong Ns_max, const arb_t H, slong sigma_interp, slong prec);
slong acb_dirichlet_platt_local_hardy_z_zeros(
    arb_ptr res, const fmpz_t n, slong len, slong prec);
slong acb_dirichlet_platt_hardy_z_zeros(
    arb_ptr res, const fmpz_t n, slong len, slong prec);

/* Discrete Fourier Transform */

void acb_dirichlet_dft_index(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec);
void acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec);

/* utils */

ACB_DIRICHLET_INLINE void
acb_vec_printd(acb_srcptr vec, slong len, slong digits)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_printd(vec + i, digits), flint_printf("\n");
}

#ifdef __cplusplus
}
#endif

#endif
