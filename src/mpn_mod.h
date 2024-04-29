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

#include "flint.h"
#include "mpn_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_mat.h"

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
    mp_size_t nlimbs;
    mp_limb_t d[MPN_MOD_MAX_LIMBS];
    mp_limb_t dinv[MPN_MOD_MAX_LIMBS];
    mp_limb_t dnormed[MPN_MOD_MAX_LIMBS];
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
int _gr_ctx_init_mpn_mod(gr_ctx_t ctx, mp_srcptr n, mp_size_t nlimbs);
void gr_ctx_init_mpn_mod_randtest(gr_ctx_t ctx, flint_rand_t state);

int mpn_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx);
void mpn_mod_ctx_clear(gr_ctx_t ctx);

MPN_MOD_INLINE truth_t
mpn_mod_ctx_is_field(gr_ctx_t ctx)
{
    return MPN_MOD_CTX_IS_PRIME(ctx);
}

MPN_MOD_INLINE void
mpn_mod_init(mp_ptr x, gr_ctx_t ctx)
{
    flint_mpn_zero(x, MPN_MOD_CTX_NLIMBS(ctx));
}

MPN_MOD_INLINE void
mpn_mod_clear(mp_ptr FLINT_UNUSED(x), gr_ctx_t FLINT_UNUSED(ctx))
{
}

MPN_MOD_INLINE void
mpn_mod_swap(mp_ptr x, mp_ptr y, gr_ctx_t ctx)
{
    slong i = 0, n = MPN_MOD_CTX_NLIMBS(ctx);
    for (i = 0; i < n; i++)
        FLINT_SWAP(mp_limb_t, x[i], y[i]);
}

MPN_MOD_INLINE int
mpn_mod_set(mp_ptr res, mp_srcptr x, gr_ctx_t ctx)
{
    flint_mpn_copyi(res, x, MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

MPN_MOD_INLINE int
mpn_mod_zero(mp_ptr res, gr_ctx_t ctx)
{
    flint_mpn_zero(res, MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

MPN_MOD_INLINE
int mpn_mod_one(mp_ptr res, gr_ctx_t ctx)
{
    res[0] = 1;
    flint_mpn_zero(res + 1, MPN_MOD_CTX_NLIMBS(ctx) - 1);
    return GR_SUCCESS;
}

int mpn_mod_set_ui(mp_ptr res, ulong x, gr_ctx_t ctx);
int mpn_mod_set_si(mp_ptr res, slong x, gr_ctx_t ctx);
int mpn_mod_neg_one(mp_ptr res, gr_ctx_t ctx);

int mpn_mod_set_mpn(mp_ptr res, mp_srcptr x, mp_size_t xn, gr_ctx_t ctx);
int mpn_mod_set_fmpz(mp_ptr res, const fmpz_t x, gr_ctx_t ctx);
int mpn_mod_set_other(mp_ptr res, gr_ptr v, gr_ctx_t v_ctx, gr_ctx_t ctx);
int mpn_mod_randtest(mp_ptr res, flint_rand_t state, gr_ctx_t ctx);
int mpn_mod_write(gr_stream_t out, mp_srcptr x, gr_ctx_t ctx);

int mpn_mod_get_fmpz(fmpz_t res, mp_srcptr x, gr_ctx_t ctx);

MPN_MOD_INLINE truth_t
mpn_mod_is_zero(mp_srcptr x, gr_ctx_t ctx)
{
    return flint_mpn_zero_p(x, MPN_MOD_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

MPN_MOD_INLINE truth_t
mpn_mod_is_one(mp_srcptr x, gr_ctx_t ctx)
{
    return (x[0] == 1 && flint_mpn_zero_p(x + 1, MPN_MOD_CTX_NLIMBS(ctx) - 1)) ? T_TRUE : T_FALSE;
}

truth_t mpn_mod_is_neg_one(gr_srcptr x, gr_ctx_t ctx);

MPN_MOD_INLINE truth_t
mpn_mod_equal(mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    return flint_mpn_equal_p(x, y, MPN_MOD_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

int mpn_mod_neg(mp_ptr res, mp_srcptr x, gr_ctx_t ctx);
int mpn_mod_add(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx);
int mpn_mod_sub(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx);
int mpn_mod_add_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_sub_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_add_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_sub_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_add_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx);
int mpn_mod_sub_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx);

int mpn_mod_mul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx);

int mpn_mod_mul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_mul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_mul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx);
int mpn_mod_addmul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx);
int mpn_mod_addmul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_addmul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_addmul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx);
int mpn_mod_submul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx);
int mpn_mod_submul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx);
int mpn_mod_submul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx);
int mpn_mod_submul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx);

MPN_MOD_INLINE int
mpn_mod_sqr(mp_ptr res, mp_srcptr x, gr_ctx_t ctx)
{
    return mpn_mod_mul(res, x, x, ctx);
}

int mpn_mod_inv(mp_ptr res, mp_srcptr x, gr_ctx_t ctx);
int mpn_mod_div(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx);

/* Vector functions */

int _mpn_mod_vec_zero(mp_ptr res, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_clear(mp_ptr FLINT_UNUSED(res), slong FLINT_UNUSED(len), gr_ctx_t FLINT_UNUSED(ctx));
int _mpn_mod_vec_set(mp_ptr res, mp_srcptr x, slong len, gr_ctx_t ctx);
void _mpn_mod_vec_swap(mp_ptr vec1, mp_ptr vec2, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_neg(mp_ptr res, mp_srcptr x, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_add(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_sub(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_mul(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_mul_scalar(mp_ptr res, mp_srcptr x, slong len, mp_srcptr y, gr_ctx_t ctx);
int _mpn_mod_scalar_mul_vec(mp_ptr res, mp_srcptr y, mp_srcptr x, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_addmul_scalar(mp_ptr res, mp_srcptr x, slong len, mp_srcptr y, gr_ctx_t ctx);
int _mpn_mod_vec_dot(mp_ptr res, mp_srcptr initial, int subtract, mp_srcptr vec1, mp_srcptr vec2, slong len, gr_ctx_t ctx);
int _mpn_mod_vec_dot_rev(mp_ptr res, mp_srcptr initial, int subtract, mp_srcptr vec1, mp_srcptr vec2, slong len, gr_ctx_t ctx);

/* Matrix algorithms */

int mpn_mod_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int mpn_mod_mat_mul_multi_mod(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int mpn_mod_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

int mpn_mod_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);
int mpn_mod_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx);
int mpn_mod_mat_lu_classical_delayed(slong * res_rank, slong * P, gr_mat_t A, const gr_mat_t A_in, int rank_check, gr_ctx_t ctx);
int mpn_mod_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);
int mpn_mod_mat_det(mp_ptr res, const gr_mat_t A, gr_ctx_t ctx);

/* Polynomial algorithms */

int _mpn_mod_poly_mullow_classical(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, gr_ctx_t ctx);
int _mpn_mod_poly_mullow_karatsuba(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, slong cutoff, gr_ctx_t ctx);
int _mpn_mod_poly_mullow_KS(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, gr_ctx_t ctx);
int _mpn_mod_poly_mullow_fft_small(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, gr_ctx_t ctx);
int _mpn_mod_poly_mullow(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, gr_ctx_t ctx);

int _mpn_mod_poly_inv_series(mp_ptr Q, mp_srcptr B, slong lenB, slong len, gr_ctx_t ctx);
int _mpn_mod_poly_div_series(mp_ptr Q, mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, slong len, gr_ctx_t ctx);

int _mpn_mod_poly_divrem_basecase_preinv1(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, mp_srcptr invL, gr_ctx_t ctx);
int _mpn_mod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, gr_ctx_t ctx);
int _mpn_mod_poly_divrem(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, gr_ctx_t ctx);
int _mpn_mod_poly_div(mp_ptr Q, mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, gr_ctx_t ctx);

int _mpn_mod_poly_gcd(mp_ptr G, slong * lenG, mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, gr_ctx_t ctx);
int _mpn_mod_poly_xgcd(slong * lenG, mp_ptr G, mp_ptr S, mp_ptr T, mp_srcptr A, slong lenA, mp_srcptr B, slong lenB, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
