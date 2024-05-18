/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_MOD_H
#define MPN_MOD_H

#ifdef MPN_MOD_INLINES_C
#define MPN_MOD_INLINE
#else
#define MPN_MOD_INLINE static inline
#endif

#include "mpn_extras.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Single limbs are already dealt with well by nmod, and excluding them
   allows avoiding various special cases. */
#define MPN_MOD_MIN_LIMBS 2

/* This should be small enough that we can stack-allocate mpn_mods
   and temporary buffers a small multiple of the size.
   For bigger sizes, use fmpz_mod.
   Before resizing this, make sure that any lookup tables that go up
   to MPN_MOD_MAX_LIMBS (such as tuning tables) are large enough. */
#define MPN_MOD_MAX_LIMBS 16

typedef struct
{
    slong nlimbs;
    ulong d[MPN_MOD_MAX_LIMBS];
    ulong dinv[MPN_MOD_MAX_LIMBS];
    ulong dnormed[MPN_MOD_MAX_LIMBS];
    flint_bitcnt_t norm;
    truth_t is_prime;
}
_mpn_mod_ctx_struct;

#define MPN_MOD_CTX(ctx) ((_mpn_mod_ctx_struct *)(GR_CTX_DATA_AS_PTR(ctx)))
#define MPN_MOD_CTX_NLIMBS(ctx) (MPN_MOD_CTX(ctx)->nlimbs)
#define MPN_MOD_CTX_MODULUS(ctx) (MPN_MOD_CTX(ctx)->d)
#define MPN_MOD_CTX_MODULUS_NORMED(ctx) (MPN_MOD_CTX(ctx)->dnormed)
#define MPN_MOD_CTX_MODULUS_PREINV(ctx) (MPN_MOD_CTX(ctx)->dinv)
#define MPN_MOD_CTX_NORM(ctx) (MPN_MOD_CTX(ctx)->norm)
#define MPN_MOD_CTX_IS_PRIME(ctx) (MPN_MOD_CTX(ctx)->is_prime)
#define MPN_MOD_CTX_MODULUS_BITS(ctx) ((MPN_MOD_CTX_NLIMBS(ctx) - 1) * FLINT_BITS + (FLINT_BITS - MPN_MOD_CTX_NORM(ctx)))

MPN_MOD_INLINE int
mpn_mod_ctx_set_is_field(gr_ctx_t ctx, truth_t is_field)
{
    MPN_MOD_CTX_IS_PRIME(ctx) = is_field;
    return GR_SUCCESS;
}

/* Basic operations and arithmetic */

int gr_ctx_init_mpn_mod(gr_ctx_t ctx, const fmpz_t n);
int _gr_ctx_init_mpn_mod(gr_ctx_t ctx, nn_srcptr n, slong nlimbs);
void gr_ctx_init_mpn_mod_randtest(gr_ctx_t ctx, flint_rand_t state);

int mpn_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx);
void mpn_mod_ctx_clear(gr_ctx_t ctx);

MPN_MOD_INLINE truth_t
mpn_mod_ctx_is_field(gr_ctx_t ctx)
{
    return MPN_MOD_CTX_IS_PRIME(ctx);
}

MPN_MOD_INLINE void
mpn_mod_init(nn_ptr x, gr_ctx_t ctx)
{
    flint_mpn_zero(x, MPN_MOD_CTX_NLIMBS(ctx));
}

MPN_MOD_INLINE void
mpn_mod_clear(nn_ptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
}

MPN_MOD_INLINE void
mpn_mod_swap(nn_ptr x, nn_ptr y, gr_ctx_t ctx)
{
    slong i = 0, n = MPN_MOD_CTX_NLIMBS(ctx);
    for (i = 0; i < n; i++)
        FLINT_SWAP(ulong, x[i], y[i]);
}

MPN_MOD_INLINE int
mpn_mod_set(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
{
    flint_mpn_copyi(res, x, MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

MPN_MOD_INLINE int
mpn_mod_zero(nn_ptr res, gr_ctx_t ctx)
{
    flint_mpn_zero(res, MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

MPN_MOD_INLINE
int mpn_mod_one(nn_ptr res, gr_ctx_t ctx)
{
    res[0] = 1;
    flint_mpn_zero(res + 1, MPN_MOD_CTX_NLIMBS(ctx) - 1);
    return GR_SUCCESS;
}

int mpn_mod_set_ui(nn_ptr res, ulong x, gr_ctx_t ctx);
int mpn_mod_set_si(nn_ptr res, slong x, gr_ctx_t ctx);
int mpn_mod_neg_one(nn_ptr res, gr_ctx_t ctx);

int mpn_mod_set_mpn(nn_ptr res, nn_srcptr x, slong xn, gr_ctx_t ctx);
int mpn_mod_set_fmpz(nn_ptr res, const fmpz_t x, gr_ctx_t ctx);
int mpn_mod_set_other(nn_ptr res, gr_ptr v, gr_ctx_t v_ctx, gr_ctx_t ctx);
int mpn_mod_randtest(nn_ptr res, flint_rand_t state, gr_ctx_t ctx);
int mpn_mod_write(gr_stream_t out, nn_srcptr x, gr_ctx_t ctx);

int mpn_mod_get_fmpz(fmpz_t res, nn_srcptr x, gr_ctx_t ctx);

MPN_MOD_INLINE truth_t
mpn_mod_is_zero(nn_srcptr x, gr_ctx_t ctx)
{
    return flint_mpn_zero_p(x, MPN_MOD_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

MPN_MOD_INLINE truth_t
mpn_mod_is_one(nn_srcptr x, gr_ctx_t ctx)
{
    return (x[0] == 1 && flint_mpn_zero_p(x + 1, MPN_MOD_CTX_NLIMBS(ctx) - 1)) ? T_TRUE : T_FALSE;
}

truth_t mpn_mod_is_neg_one(gr_srcptr x, gr_ctx_t ctx);

MPN_MOD_INLINE truth_t
mpn_mod_equal(nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
{
    return flint_mpn_equal_p(x, y, MPN_MOD_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

int mpn_mod_neg(nn_ptr res, nn_srcptr x, gr_ctx_t ctx);
int mpn_mod_add(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx);
int mpn_mod_sub(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx);
int mpn_mod_add_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_sub_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_add_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_sub_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_add_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx);
int mpn_mod_sub_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx);

int mpn_mod_mul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx);

int mpn_mod_mul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_mul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_mul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx);
int mpn_mod_addmul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx);
int mpn_mod_addmul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_addmul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_addmul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx);
int mpn_mod_submul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx);
int mpn_mod_submul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_submul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_submul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx);

MPN_MOD_INLINE int
mpn_mod_sqr(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
{
    return mpn_mod_mul(res, x, x, ctx);
}

int mpn_mod_inv(nn_ptr res, nn_srcptr x, gr_ctx_t ctx);
int mpn_mod_div(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx);

/* Vector functions */

int _mpn_mod_vec_zero(nn_ptr res, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_clear(nn_ptr FLINT_UNUSED(res), slong FLINT_UNUSED(len), gr_ctx_t FLINT_UNUSED(ctx));
int _mpn_mod_vec_set(nn_ptr res, nn_srcptr x, slong len, gr_ctx_t ctx);
void _mpn_mod_vec_swap(nn_ptr vec1, nn_ptr vec2, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_neg(nn_ptr res, nn_srcptr x, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_add(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_sub(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_mul(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_mul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx);
int _mpn_mod_scalar_mul_vec(nn_ptr res, nn_srcptr y, nn_srcptr x, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_addmul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx);
int _mpn_mod_vec_dot(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, nn_srcptr vec2, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_dot_rev(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, nn_srcptr vec2, slong len, gr_ctx_t ctx);

/* Matrix algorithms */

int mpn_mod_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int mpn_mod_mat_mul_multi_mod(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int mpn_mod_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

int mpn_mod_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);
int mpn_mod_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx);
int mpn_mod_mat_lu_classical_delayed(slong * res_rank, slong * P, gr_mat_t A, const gr_mat_t A_in, int rank_check, gr_ctx_t ctx);
int mpn_mod_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);
int mpn_mod_mat_det(nn_ptr res, const gr_mat_t A, gr_ctx_t ctx);

/* Polynomial algorithms */

int _mpn_mod_poly_mullow_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx);
int _mpn_mod_poly_mullow_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, slong cutoff, gr_ctx_t ctx);
int _mpn_mod_poly_mullow_KS(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx);
int _mpn_mod_poly_mullow_fft_small(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx);
int _mpn_mod_poly_mullow(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx);

int _mpn_mod_poly_inv_series(nn_ptr Q, nn_srcptr B, slong lenB, slong len, gr_ctx_t ctx);
int _mpn_mod_poly_div_series(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, slong len, gr_ctx_t ctx);

int _mpn_mod_poly_divrem_basecase_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
int _mpn_mod_poly_divrem_basecase(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx);
int _mpn_mod_poly_divrem(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx);
int _mpn_mod_poly_div(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx);

int _mpn_mod_poly_gcd(nn_ptr G, slong * lenG, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx);
int _mpn_mod_poly_xgcd(slong * lenG, nn_ptr G, nn_ptr S, nn_ptr T, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
