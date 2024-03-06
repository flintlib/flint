/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "fmpz.h"
#include "fmpz_mod.h"

void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t n)
{
    flint_bitcnt_t bits;

    if  (fmpz_sgn(n) <= 0)
    {
        flint_throw(FLINT_ERROR, "Exception in fmpz_mod_ctx_init: "
                                                    "Modulus is nonpositive.");
    }

    /* prepare for general case */
    fmpz_init_set(ctx->n, n);
    ctx->n_limbs[0] = 0;
    ctx->n_limbs[1] = 0;
    ctx->n_limbs[2] = 0;
    ctx->add_fxn = _fmpz_mod_addN;
    ctx->sub_fxn = _fmpz_mod_subN;
    ctx->mul_fxn = _fmpz_mod_mulN;
    ctx->ninv_huge = NULL;

    bits = fmpz_bits(n);
    if (bits <= FLINT_BITS)
    {
        ctx->add_fxn = _fmpz_mod_add1;
        ctx->sub_fxn = _fmpz_mod_sub1;
        ctx->mul_fxn = _fmpz_mod_mul1;
        nmod_init(&ctx->mod, fmpz_get_ui(n));
    }
    else if (bits <= 2*FLINT_BITS)
    {
        fmpz_get_ui_array(ctx->n_limbs, 3, n); /* n_limbs[2] will be 0 */

        /* n = 2^FLINT_BITS must be special case */
        if (ctx->n_limbs[1] == 1 && ctx->n_limbs[0] == 0)
        {
            ctx->add_fxn = _fmpz_mod_add2s;
            ctx->sub_fxn = _fmpz_mod_sub2s;
            ctx->mul_fxn = _fmpz_mod_mul2s;
        }
        else
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_one(t);
            fmpz_mul_2exp(t, t, 4*FLINT_BITS);
            fmpz_tdiv_q(t, t, n);
            fmpz_get_ui_array(ctx->ninv_limbs, 3, t);
            fmpz_clear(t);
            FLINT_ASSERT(ctx->ninv_limbs[2] != 0);
            ctx->add_fxn = _fmpz_mod_add2;
            ctx->sub_fxn = _fmpz_mod_sub2;
            ctx->mul_fxn = _fmpz_mod_mul2;
        }
    }
#if FLINT_HAVE_FFT_SMALL
    else if (bits >= 19000)
    {
        ctx->ninv_huge = flint_malloc(sizeof(fmpz_preinvn_struct));
        fmpz_preinvn_init(ctx->ninv_huge, n);
    }
#endif
}

void fmpz_mod_ctx_init_ui(fmpz_mod_ctx_t ctx, ulong n)
{
    fmpz_t n_;
    fmpz_init_set_ui(n_, n);
    fmpz_mod_ctx_init(ctx, n_);
    fmpz_clear(n_);
}

void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx)
{
    if (ctx->ninv_huge != NULL)
    {
        fmpz_preinvn_clear(ctx->ninv_huge);
        flint_free(ctx->ninv_huge);
    }

    fmpz_clear(ctx->n);
}

void fmpz_mod_ctx_init_rand_bits(fmpz_mod_ctx_t ctx,
                                   flint_rand_t state, flint_bitcnt_t max_bits)
{
    fmpz_t m;
    fmpz_init(m);
    fmpz_randtest_unsigned(m, state, max_bits);
    fmpz_add_ui(m, m, 2);
    fmpz_mod_ctx_init(ctx, m);
    fmpz_clear(m);
}

void fmpz_mod_ctx_init_rand_bits_prime(fmpz_mod_ctx_t ctx,
                                   flint_rand_t state, flint_bitcnt_t max_bits)
{
    fmpz_t m;
    fmpz_init(m);
    fmpz_randtest_unsigned(m, state, max_bits);
    fmpz_nextprime(m, m, 0);
    fmpz_mod_ctx_init(ctx, m);
    fmpz_clear(m);
}
