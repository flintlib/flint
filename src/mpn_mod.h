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

/* Helpers which actually belong in mpn_extras.h */

FLINT_FORCE_INLINE
int flint_mpn_equal_p(mp_srcptr x, mp_srcptr y, mp_size_t xsize)
{
    slong i;
    for (i = 0; i < xsize; i++)
    {
        if (x[i] != y[i])
            return 0;
    }
    return 1;
}

FLINT_FORCE_INLINE void
flint_mpn_negmod_n(mp_ptr res, mp_srcptr x, mp_srcptr m, mp_size_t n)
{
    if (flint_mpn_zero_p(x, n))
        flint_mpn_zero(res, n);
    else
        mpn_sub_n(res, m, x, n);
}

FLINT_FORCE_INLINE void
flint_mpn_addmod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
{
    mp_limb_t cy;
    cy = mpn_add_n(res, x, y, n);
    if (cy || mpn_cmp(res, m, n) >= 0)
        mpn_sub_n(res, res, m, n);
}

FLINT_FORCE_INLINE void
flint_mpn_submod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
{
    int cmp = (mpn_cmp(x, y, n) < 0);
    mpn_sub_n(res, x, y, n);
    if (cmp)
        mpn_add_n(res, res, m, n);
}

/* assumes m <= n and y < m */
FLINT_FORCE_INLINE void
flint_mpn_addmod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)
{
    mp_limb_t cy;
    cy = mpn_add(res, x, n, y, yn);
    if (cy || mpn_cmp(res, m, n) >= 0)
        mpn_sub_n(res, res, m, n);
}

/* assumes m <= n and y < m */
FLINT_FORCE_INLINE void
flint_mpn_submod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)
{
    int cmp = (flint_mpn_zero_p(x + yn, n - yn) && mpn_cmp(x, y, yn) < 0);
    mpn_sub(res, x, n, y, yn);
    if (cmp)
        mpn_add_n(res, res, m, n);
}

FLINT_FORCE_INLINE void
flint_mpn_negmod_2(mp_ptr res, mp_srcptr x, mp_srcptr m)
{
    if (x[0] == 0 && x[1] == 0)
        res[1] = res[0] = 0;
    else
        sub_ddmmss(res[1], res[0], m[1], m[0], x[1], x[0]);
}

FLINT_FORCE_INLINE void
flint_mpn_addmod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
{
    mp_limb_t cy;
    mp_limb_t m1 = m[1], m0 = m[0];
    add_sssaaaaaa(cy, res[1], res[0], 0, x[1], x[0], 0, y[1], y[0]);
    if (cy || (res[1] > m1 || (res[1] == m1 && res[0] >= m0)))
        sub_ddmmss(res[1], res[0], res[1], res[0], m1, m0);
}

FLINT_FORCE_INLINE void
flint_mpn_submod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
{
    int cmp;
    mp_limb_t m1 = m[1], m0 = m[0];
    mp_limb_t x1 = x[1], x0 = x[0];
    mp_limb_t y1 = y[1], y0 = y[0];
    cmp = (x1 < y1) || (x1 == y1 && x0 < y0);
    sub_ddmmss(res[1], res[0], x1, x0, y1, y0);
    if (cmp)
        add_ssaaaa(res[1], res[0], res[1], res[0], m1, m0);
}

#ifdef __cplusplus
}
#endif

#endif
