/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_MOD_H
#define N_MOD_H

#ifdef N_MOD_INLINES_C
# define N_MOD_INLINE
#else
# define N_MOD_INLINE static inline
#endif

#include "ulong_extras.h"
#include "n_mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* context management ********************************************************/

N_MOD_INLINE void n_mod_ctx_init(n_mod_ctx_ptr ctx, ulong n)
{
    unsigned int norm = flint_clz(n);

    ctx->nu = n;
    ctx->nn = n << norm;
    ctx->ninv = n_preinvert_limb_prenorm(n << norm);
    ctx->norm = norm;
    /* FIXME: Add more things here such as choice of algorithm */
}

N_MOD_INLINE void n_mod_ctx_clear(n_mod_ctx_ptr ctx) { }

void n_mod_ctx_init_rand(n_mod_ctx_ptr, flint_rand_t);

/* memory management *********************************************************/

N_MOD_INLINE int n_mod_is_canonical(ulong a, n_mod_ctx_srcptr ctx)
{
    return a < ctx->nu;
}

/* assignments and conversions ***********************************************/

/* FIXME */
N_MOD_INLINE ulong n_mod_set_uiui(ulong a1, ulong a0, n_mod_ctx_srcptr ctx)
{
    unsigned int norm = ctx->norm;
    ulong nn = ctx->nn, ninv = ctx->ninv;
    ulong q0, q1, r1, u0, u1;

    u1 = (a1 << norm) + ((norm == 0) ? UWORD(0) : a0 >> (FLINT_BITS - norm));
    u0 = a0 << norm;

    umul_ppmm(q1, q0, ninv, u1);
    add_ssaaaa(q1, q0, q1, q0, u1, u0);
    r1 = (u0 - (q1 + 1) * nn);
    if (r1 > q0)
        r1 += nn;

    return (r1 < nn) ? r1 >> norm : (r1 - nn) >> norm;
}

N_MOD_INLINE ulong n_mod_set_ui(ulong a, n_mod_ctx_srcptr ctx)
{
    if (a < ctx->nu)
        return a;
    else
        return n_mod_set_uiui(0, a, ctx);
}

N_MOD_INLINE ulong n_mod_set_uiuiui(ulong a2, ulong a1, ulong a0, n_mod_ctx_srcptr ctx)
{
}

/* randomisation *************************************************************/

N_MOD_INLINE ulong n_mod_rand(flint_rand_t state, n_mod_ctx_srcptr ctx)
{
    return n_mod_set_ui(n_randlimb(state), ctx);
}

/* addition ******************************************************************/

N_MOD_INLINE ulong n_mod_add(ulong a, ulong b, n_mod_ctx_srcptr ctx)
{
    ulong sum = a + b;
    return sum - ctx->nu + ((((slong) (sum - ctx->nu)) >> (FLINT_BITS - 1)) & ctx->nu);
}

N_MOD_INLINE ulong n_mod_sub(ulong a, ulong b, n_mod_ctx_srcptr ctx)
{
    ulong diff = a - b;
    return diff + ((((slong) diff) >> (FLINT_BITS - 1)) & ctx->nu);
}

N_MOD_INLINE ulong n_mod_neg(ulong a, n_mod_ctx_srcptr ctx)
{
    return (a != UWORD(0)) ? ctx->nu - a : 0;
}

/* multiplication ************************************************************/

N_MOD_INLINE ulong n_mod_mul(ulong a, ulong b, n_mod_ctx_srcptr ctx)
{
    ulong nn = ctx->nn, ninv = ctx->ninv;
    unsigned int norm = ctx->norm;
    ulong q1, q0, r, p1, p0;
    umul_ppmm(p1, p0, a, b << norm);
    umul_ppmm(q1, q0, ninv, p1);
    add_ssaaaa(q1, q0, q1, q0, p1, p0);
    r = p0 - (q1 + 1) * nn;
    if (r > q0)
        r += nn;
    return (r < nn ? r : r - nn) >> norm;
}

/* NOTE: Unchecked inverse. */
N_MOD_INLINE ulong n_mod_inv(ulong a, n_mod_ctx_srcptr ctx)
{
    nn_pair_t ret = _n_gcdinv(a, ctx->nu);
    return ret.m0;
}

N_MOD_INLINE ulong n_mod_div(ulong a, ulong b, n_mod_ctx_srcptr ctx)
{
    return n_mod_mul(a, n_mod_inv(b, ctx), ctx);
}

/* exponentiation ************************************************************/

/* FIXME: si, ui and fmpz methods makes sense */

#ifdef __cplusplus
}
#endif

#endif
