/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_SERIES_H
#define GR_SERIES_H

#include "flint.h"
#include "gr.h"
#include "dirichlet.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    gr_ctx_struct * base_ring;
    slong n;
    char * var;
}
series_mod_ctx_t;

typedef struct
{
    gr_ctx_struct * base_ring;
    slong prec;     /* default approximate truncation */
    char * var;
}
series_ctx_t;

#define GR_SERIES_CTX(ring_ctx) ((series_ctx_t *)((ring_ctx)))
#define GR_SERIES_ELEM_CTX(ring_ctx) (GR_SERIES_CTX(ring_ctx)->base_ring)
#define GR_SERIES_PREC(ring_ctx) (GR_SERIES_CTX(ring_ctx)->prec)

#define GR_SERIES_MOD_CTX(ring_ctx) ((series_mod_ctx_t *)((ring_ctx)))
#define GR_SERIES_MOD_ELEM_CTX(ring_ctx) (GR_SERIES_MOD_CTX(ring_ctx)->base_ring)
#define GR_SERIES_MOD_N(ring_ctx) (GR_SERIES_MOD_CTX(ring_ctx)->n)


#define GR_SERIES_ERR_EXACT WORD_MAX
#define GR_SERIES_ERR_MAX WORD_MAX / 4

typedef struct
{
    gr_poly_struct poly;
    slong error;
}
gr_series_struct;

typedef gr_series_struct gr_series_t[1];

#define GR_SERIES_POLY(x) (&((x)->poly))
#define GR_SERIES_ERROR(x) ((x)->error)

typedef struct
{
    gr_series_struct * entries;
    slong alloc;
    slong length;
}
gr_series_vec_struct;

typedef gr_series_vec_struct gr_series_vec_t[1];


void gr_series_ctx_init(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec);
void gr_series_ctx_clear(gr_ctx_t ctx);

int gr_series_ctx_write(gr_stream_t out, gr_ctx_t ctx);
truth_t gr_series_ctx_is_ring(gr_ctx_t ctx);
truth_t gr_series_ctx_is_commutative_ring(gr_ctx_t ctx);
truth_t gr_series_ctx_is_integral_domain(gr_ctx_t ctx);
truth_t gr_series_ctx_is_rational_vector_space(gr_ctx_t ctx);
truth_t gr_series_ctx_is_real_vector_space(gr_ctx_t ctx);
truth_t gr_series_ctx_is_complex_vector_space(gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_ctx_set_gen_name(gr_ctx_t ctx, const char * s);
WARN_UNUSED_RESULT int gr_series_ctx_set_gen_names(gr_ctx_t ctx, const char ** s);

void gr_series_init(gr_series_t res, gr_ctx_t ctx);
void gr_series_clear(gr_series_t res, gr_ctx_t ctx);
void gr_series_swap(gr_series_t x, gr_series_t y, gr_ctx_t ctx);
void gr_series_set_shallow(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);

slong _gr_series_get_error(const gr_series_t x, gr_ctx_t ctx);
truth_t _gr_series_is_exact(const gr_series_t x, gr_ctx_t ctx);
void _gr_series_set_error(gr_series_t x, slong err, gr_ctx_t ctx);
void _gr_series_make_exact(gr_series_t x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_series_randtest(gr_series_t res, flint_rand_t state, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_write(gr_stream_t out, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_zero(gr_series_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_one(gr_series_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_set(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_gen(gr_series_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_gens_recursive(gr_vec_t vec, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_neg(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_set_gr_poly(gr_series_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_set_scalar(gr_series_t res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_set_si(gr_series_t res, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_set_ui(gr_series_t res, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_set_fmpz(gr_series_t res, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_set_fmpq(gr_series_t res, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_set_other(gr_series_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
truth_t gr_series_is_zero(const gr_series_t x, gr_ctx_t ctx);
truth_t gr_series_is_one(const gr_series_t x, gr_ctx_t ctx);
truth_t gr_series_coeff_is_zero(const gr_series_t x, slong i, gr_ctx_t ctx);
truth_t gr_series_equal(const gr_series_t x, const gr_series_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_add(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_sub(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mul(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_inv(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_div(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_divexact(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_sqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_rsqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_series_exp(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_log(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_tan(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_asin(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_acos(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_atan(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_asinh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_acosh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_atanh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);

/* currently arb/acb wrappers only */
WARN_UNUSED_RESULT int gr_series_gamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_rgamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_lgamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_digamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_erf(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_erfc(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_erfi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_exp_integral_ei(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_cos_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_cosh_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_sin_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_sinh_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_series_fresnel(gr_series_t res1, gr_series_t res2, const gr_series_t x, int normalized, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_fresnel_s(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_fresnel_c(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_airy(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_airy_ai(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_airy_ai_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_airy_bi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_airy_bi_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_log_integral(gr_series_t res, const gr_series_t x, int offset, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_gamma_upper(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_gamma_lower(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_beta_lower(gr_series_t res, const gr_series_t a, const gr_series_t b, const gr_series_t x, int regularized, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_polylog(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_hurwitz_zeta(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_dirichlet_l(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_dirichlet_hardy_theta(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_dirichlet_hardy_z(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_jacobi_theta(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_jacobi_theta_1(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_jacobi_theta_2(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_jacobi_theta_3(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_jacobi_theta_4(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_agm1(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_elliptic_k(gr_series_t res, const gr_series_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_weierstrass_p(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_hypgeom_pfq(gr_series_t res, const gr_series_vec_t a, const gr_series_vec_t b, const gr_series_t x, int regularized, gr_ctx_t ctx);

void gr_series_mod_ctx_init(gr_ctx_t ctx, gr_ctx_t base_ring, slong n);
void gr_series_mod_ctx_clear(gr_ctx_t ctx);
int gr_series_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx);
truth_t gr_series_mod_ctx_is_ring(gr_ctx_t ctx);
truth_t gr_series_mod_ctx_is_commutative_ring(gr_ctx_t ctx);
truth_t gr_series_mod_ctx_is_integral_domain(gr_ctx_t ctx);
truth_t gr_series_mod_ctx_is_rational_vector_space(gr_ctx_t ctx);
truth_t gr_series_mod_ctx_is_real_vector_space(gr_ctx_t ctx);
truth_t gr_series_mod_ctx_is_complex_vector_space(gr_ctx_t ctx);
truth_t gr_series_mod_ctx_is_field(gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_ctx_set_gen_name(gr_ctx_t ctx, const char * s);
WARN_UNUSED_RESULT int gr_series_mod_ctx_set_gen_names(gr_ctx_t ctx, const char ** s);
WARN_UNUSED_RESULT int gr_series_mod_gens_recursive(gr_vec_t vec, gr_ctx_t ctx);
void gr_series_mod_init(gr_poly_t res, gr_ctx_t ctx);
void gr_series_mod_clear(gr_poly_t res, gr_ctx_t ctx);
void gr_series_mod_swap(gr_poly_t x, gr_poly_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_randtest(gr_poly_t res, flint_rand_t state, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_write(gr_stream_t out, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_zero(gr_poly_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_one(gr_poly_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_gen(gr_poly_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_set(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_set_other(gr_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
void gr_series_mod_set_shallow(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_set_si(gr_poly_t res, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_set_ui(gr_poly_t res, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_set_fmpz(gr_poly_t res, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_set_fmpq(gr_poly_t res, const fmpq_t c, gr_ctx_t ctx);
truth_t gr_series_mod_is_zero(const gr_poly_t x, gr_ctx_t ctx);
truth_t gr_series_mod_is_one(const gr_poly_t x, gr_ctx_t ctx);
truth_t gr_series_mod_equal(const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_neg(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_add(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_sub(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_mul(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_inv(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_div(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_exp(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_log(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_sqrt(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_rsqrt(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_tan(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_asin(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_acos(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_atan(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_asinh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_acosh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_series_mod_atanh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif /* GR_SERIES_H */
