/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ACB_HYPGEOM_H
#define ACB_HYPGEOM_H

#include "acb.h"
#include "acb_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

void acb_hypgeom_rising_ui_forward(acb_t res, const acb_t x, ulong n, slong prec);
void acb_hypgeom_rising_ui_rs(acb_t res, const acb_t x, ulong n, ulong m, slong prec);
void acb_hypgeom_rising_ui_bs(acb_t res, const acb_t x, ulong n, slong prec);
void acb_hypgeom_rising_ui_rec(acb_t res, const acb_t x, ulong n, slong prec);
void acb_hypgeom_rising_ui(acb_t y, const acb_t x, ulong n, slong prec);
void acb_hypgeom_rising(acb_t y, const acb_t x, const acb_t n, slong prec);

void acb_hypgeom_rising_ui_jet_powsum(acb_ptr res, const acb_t x, ulong n, slong len, slong prec);
void acb_hypgeom_rising_ui_jet_rs(acb_ptr res, const acb_t x, ulong n, ulong m, slong len, slong prec);
void acb_hypgeom_rising_ui_jet_bs(acb_ptr res, const acb_t x, ulong n, slong len, slong prec);
void acb_hypgeom_rising_ui_jet(acb_ptr res, const acb_t x, ulong n, slong len, slong prec);

void acb_hypgeom_log_rising_ui_jet(acb_ptr res, const acb_t z, ulong r, slong len, slong prec);
void acb_hypgeom_log_rising_ui(acb_ptr res, const acb_t z, ulong r, slong prec);

void acb_hypgeom_gamma_stirling_sum_horner(acb_t s, const acb_t z, slong N, slong prec);
void acb_hypgeom_gamma_stirling_sum_improved(acb_t s, const acb_t z, slong N, slong K, slong prec);
void acb_hypgeom_gamma_stirling(acb_t res, const acb_t x, int reciprocal, slong prec);
int acb_hypgeom_gamma_taylor(acb_t res, const acb_t x, int reciprocal, slong prec);
void acb_hypgeom_gamma(acb_t y, const acb_t x, slong prec);
void acb_hypgeom_rgamma(acb_t y, const acb_t x, slong prec);

void acb_hypgeom_lgamma(acb_t y, const acb_t x, slong prec);

void acb_hypgeom_pfq_bound_factor(mag_t C,
    acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, ulong n);

slong acb_hypgeom_pfq_choose_n(acb_srcptr a, slong p,
                         acb_srcptr b, slong q, const acb_t z, slong prec);

void acb_hypgeom_pfq_sum_forward(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec);

void acb_hypgeom_pfq_sum_rs(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec);

void acb_hypgeom_pfq_sum_bs(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec);

void acb_hypgeom_pfq_sum_fme(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec);

void acb_hypgeom_pfq_sum(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec);

void acb_hypgeom_pfq_sum_bs_invz(acb_t s, acb_t t,
    acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec);

void acb_hypgeom_pfq_sum_invz(acb_t s, acb_t t, acb_srcptr a, slong p,
    acb_srcptr b, slong q, const acb_t z, const acb_t zinv, slong n, slong prec);

void acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec);

slong acb_hypgeom_pfq_series_choose_n(const acb_poly_struct * a, slong p,
                                const acb_poly_struct * b, slong q,
                                const acb_poly_t z, slong len, slong prec);

void
acb_hypgeom_pfq_series_sum_forward(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec);

void
acb_hypgeom_pfq_series_sum_bs(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec);

void
acb_hypgeom_pfq_series_sum_rs(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec);

void
acb_hypgeom_pfq_series_sum(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec);

void acb_hypgeom_pfq_series_direct(acb_poly_t res,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec);

void acb_hypgeom_pfq(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, int regularized, slong prec);

void acb_hypgeom_u_asymp(acb_t res, const acb_t a, const acb_t b,
    const acb_t z, slong n, slong prec);

void acb_hypgeom_u_1f1_series(acb_poly_t res,
    const acb_poly_t a, const acb_poly_t b, const acb_poly_t z,
    slong len, slong prec);

void acb_hypgeom_u_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec);
void acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec);

int acb_hypgeom_u_use_asymp(const acb_t z, slong prec);

void acb_hypgeom_m_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec);
void acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec);
void acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec);
void acb_hypgeom_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec);

void acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, slong prec);
void acb_hypgeom_bessel_j_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec);
void acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, slong prec);

void acb_hypgeom_bessel_i_0f1(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec);
void acb_hypgeom_bessel_i_asymp(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec);
void acb_hypgeom_bessel_i(acb_t res, const acb_t nu, const acb_t z, slong prec);
void acb_hypgeom_bessel_i_scaled(acb_t res, const acb_t nu, const acb_t z, slong prec);

void acb_hypgeom_bessel_k_0f1(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec);
void acb_hypgeom_bessel_k_0f1_series(acb_poly_t res, const acb_poly_t n, const acb_poly_t z, int scaled, slong len, slong prec);
void acb_hypgeom_bessel_k_asymp(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec);
void acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, slong prec);
void acb_hypgeom_bessel_k_scaled(acb_t res, const acb_t nu, const acb_t z, slong prec);

void acb_hypgeom_bessel_y(acb_t res, const acb_t nu, const acb_t z, slong prec);
void acb_hypgeom_bessel_jy(acb_t res1, acb_t res2, const acb_t nu, const acb_t z, slong prec);

void acb_hypgeom_0f1_asymp(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec);
void acb_hypgeom_0f1_direct(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec);
void acb_hypgeom_0f1(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec);

void acb_hypgeom_airy_bound(mag_t ai, mag_t aip, mag_t bi, mag_t bip, const acb_t z);
void acb_hypgeom_airy_asymp(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, slong n, slong prec);
void acb_hypgeom_airy_direct(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, slong n, slong prec);
void acb_hypgeom_airy(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, slong prec);
void acb_hypgeom_airy_jet(acb_ptr ai, acb_ptr bi, const acb_t z, slong len, slong prec);
void _acb_hypgeom_airy_series(acb_ptr ai, acb_ptr ai_prime, acb_ptr bi, acb_ptr bi_prime, acb_srcptr z, slong zlen, slong len, slong prec);
void acb_hypgeom_airy_series(acb_poly_t ai, acb_poly_t ai_prime, acb_poly_t bi, acb_poly_t bi_prime, const acb_poly_t z, slong len, slong prec);

void acb_hypgeom_coulomb(acb_t F, acb_t G, acb_t Hpos, acb_t Hneg, const acb_t l, const acb_t eta, const acb_t z, slong prec);
void acb_hypgeom_coulomb_jet(acb_ptr F, acb_ptr G, acb_ptr Hpos, acb_ptr Hneg, const acb_t l, const acb_t eta, const acb_t z, slong len, slong prec);
void _acb_hypgeom_coulomb_series(acb_ptr F, acb_ptr G, acb_ptr Hpos, acb_ptr Hneg, const acb_t l, const acb_t eta, acb_srcptr z, slong zlen, slong len, slong prec);
void acb_hypgeom_coulomb_series(acb_poly_t F, acb_poly_t G, acb_poly_t Hpos, acb_poly_t Hneg, const acb_t l, const acb_t eta, const acb_poly_t z, slong len, slong prec);

void acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s, const acb_t z, int modified, slong prec);
void acb_hypgeom_gamma_upper_1f1a(acb_t res, const acb_t s, const acb_t z, int modified, slong prec);
void acb_hypgeom_gamma_upper_1f1b(acb_t res, const acb_t s, const acb_t z, int modified, slong prec);
void acb_hypgeom_gamma_upper_singular(acb_t res, slong s, const acb_t z, int modified, slong prec);
void acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int modified, slong prec);

void _acb_hypgeom_gamma_upper_series(acb_ptr g, const acb_t s, acb_srcptr h, slong hlen, int regularized, slong n, slong prec);
void acb_hypgeom_gamma_upper_series(acb_poly_t g, const acb_t s, const acb_poly_t h, int regularized, slong n, slong prec);

void acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z, int modified, slong prec);

void _acb_hypgeom_gamma_lower_series(acb_ptr g, const acb_t s, acb_srcptr h, slong hlen, int regularized, slong n, slong prec);
void acb_hypgeom_gamma_lower_series(acb_poly_t g, const acb_t s, const acb_poly_t h, int regularized, slong n, slong prec);

void acb_hypgeom_beta_lower(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec);
void _acb_hypgeom_beta_lower_series(acb_ptr res, const acb_t a, const acb_t b, acb_srcptr z, slong zlen, int regularized, slong len, slong prec);
void acb_hypgeom_beta_lower_series(acb_poly_t res, const acb_t a, const acb_t b, const acb_poly_t z, int regularized, slong len, slong prec);

void acb_hypgeom_expint(acb_t res, const acb_t s, const acb_t z, slong prec);

void acb_hypgeom_erf_propagated_error(mag_t re, mag_t im, const acb_t z);
void acb_hypgeom_erf_1f1a(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_erf_1f1b(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_erf_asymp(acb_t res, const acb_t z, int complementary, slong prec, slong prec2);
void acb_hypgeom_erf(acb_t res, const acb_t z, slong prec);
void _acb_hypgeom_erf_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_hypgeom_erf_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec);

void acb_hypgeom_erfc(acb_t res, const acb_t z, slong prec);
void _acb_hypgeom_erfc_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_hypgeom_erfc_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec);

void acb_hypgeom_erfi(acb_t res, const acb_t z, slong prec);
void _acb_hypgeom_erfi_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_hypgeom_erfi_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec);

void acb_hypgeom_fresnel(acb_t res1, acb_t res2, const acb_t z, int normalized, slong prec);
void _acb_hypgeom_fresnel_series(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, int normalized, slong len, slong prec);
void acb_hypgeom_fresnel_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, int normalized, slong len, slong prec);

void acb_hypgeom_ei_asymp(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_ei_2f2(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_ei(acb_t res, const acb_t z, slong prec);
void _acb_hypgeom_ei_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_hypgeom_ei_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec);

void acb_hypgeom_si_asymp(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_si_1f2(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_si(acb_t res, const acb_t z, slong prec);
void _acb_hypgeom_si_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_hypgeom_si_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec);

void acb_hypgeom_ci_asymp(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_ci_2f3(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_ci(acb_t res, const acb_t z, slong prec);
void _acb_hypgeom_ci_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_hypgeom_ci_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec);

void acb_hypgeom_shi(acb_t res, const acb_t z, slong prec);
void _acb_hypgeom_shi_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_hypgeom_shi_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec);

void acb_hypgeom_chi_asymp(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_chi_2f3(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_chi(acb_t res, const acb_t z, slong prec);
void _acb_hypgeom_chi_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_hypgeom_chi_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec);

void acb_hypgeom_li(acb_t res, const acb_t z, int offset, slong prec);
void _acb_hypgeom_li_series(acb_ptr g, acb_srcptr h, slong hlen, int offset, slong len, slong prec);
void acb_hypgeom_li_series(acb_poly_t g, const acb_poly_t h, int offset, slong len, slong prec);

void acb_hypgeom_2f1_continuation(acb_t res0, acb_t res1,
    const acb_t a, const acb_t b, const acb_t c, const acb_t z0,
    const acb_t z1, const acb_t f0, const acb_t f1, slong prec);

void acb_hypgeom_2f1_series_direct(acb_poly_t res, const acb_poly_t a, const acb_poly_t b,
    const acb_poly_t c, const acb_poly_t z, int regularized, slong len, slong prec);
void acb_hypgeom_2f1_direct(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec);

void acb_hypgeom_2f1_transform(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, int which, slong prec);

void acb_hypgeom_2f1_transform_limit(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, int which, slong prec);

void acb_hypgeom_2f1_corner(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec);

int acb_hypgeom_2f1_choose(const acb_t z);

void acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec);

#define ACB_HYPGEOM_2F1_REGULARIZED 1
#define ACB_HYPGEOM_2F1_AB 2   /* a-b integer */
#define ACB_HYPGEOM_2F1_AC 4   /* a-c integer */
#define ACB_HYPGEOM_2F1_BC 8   /* b-c integer */
#define ACB_HYPGEOM_2F1_ABC 16  /* a+b-c integer */

void acb_hypgeom_legendre_p_uiui_rec(acb_t res, ulong n, ulong m, const acb_t z, slong prec);
void acb_hypgeom_legendre_p(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, slong prec);
void acb_hypgeom_legendre_q(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, slong prec);
void acb_hypgeom_jacobi_p(acb_t res, const acb_t n, const acb_t a, const acb_t b, const acb_t z, slong prec);
void acb_hypgeom_gegenbauer_c(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec);
void acb_hypgeom_laguerre_l(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec);
void acb_hypgeom_hermite_h(acb_t res, const acb_t n, const acb_t z, slong prec);
void acb_hypgeom_chebyshev_t(acb_t res, const acb_t n, const acb_t z, slong prec);
void acb_hypgeom_chebyshev_u(acb_t res, const acb_t n, const acb_t z, slong prec);
void acb_hypgeom_spherical_y(acb_t res, slong n, slong m, const acb_t theta, const acb_t phi, slong prec);

void acb_hypgeom_dilog_bernoulli(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_dilog_continuation(acb_t res, const acb_t a, const acb_t z, slong prec);
void acb_hypgeom_dilog_bitburst(acb_t res, acb_t z0, const acb_t z, slong prec);
void acb_hypgeom_dilog_transform(acb_t res, const acb_t z, int algorithm, slong prec);
void acb_hypgeom_dilog_zero_taylor(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_dilog_zero(acb_t res, const acb_t z, slong prec);
void acb_hypgeom_dilog(acb_t res, const acb_t z, slong prec);

#ifdef __cplusplus
}
#endif

#endif

