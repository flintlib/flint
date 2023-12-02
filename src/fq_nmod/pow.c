/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fq_nmod.h"

void _fq_nmod_pow(mp_limb_t *rop, const mp_limb_t *op, slong len, const fmpz_t e,
                  const fq_nmod_ctx_t ctx)
{
    const slong d = fq_nmod_ctx_degree(ctx);

    if (fmpz_is_zero(e))
    {
        rop[0] = WORD(1);
        _nmod_vec_zero(rop + 1, 2 * d - 1 - 1);
    }
    else if (fmpz_is_one(e))
    {
        _nmod_vec_set(rop, op, len);
        _nmod_vec_zero(rop + len, 2 * d - 1 - len);
    }
    else
    {
        ulong bit;
        mp_limb_t *v = _nmod_vec_init(2 * d - 1);
        mp_limb_t *R, *S, *T;

        _nmod_vec_zero(v, 2 * d - 1);
        _nmod_vec_zero(rop, 2 * d - 1);

        /*
           Set bits to the bitmask with a 1 one place lower than the msb of e
         */

        bit = fmpz_bits(e) - 2;

        /*
           Trial run without any polynomial arithmetic to determine the parity
           of the number of swaps;  then set R and S accordingly
         */

        {
            unsigned int swaps = 0U;
            ulong bit2 = bit;
            if (fmpz_tstbit(e, bit2))
                swaps = ~swaps;
            while (bit2--)
                if (!fmpz_tstbit(e, bit2))
                    swaps = ~swaps;

            if (swaps == 0U)
            {
                R = rop;
                S = v;
            }
            else
            {
                R = v;
                S = rop;
            }
        }

        /*
           We unroll the first step of the loop, referring to {op, len}
         */

        _nmod_poly_mul(R, op, len, op, len, ctx->mod);
        _fq_nmod_reduce(R, 2 * len - 1, ctx);

        if (fmpz_tstbit(e, bit))
        {
            _nmod_poly_mul(S, R, d, op, len, ctx->mod);
            _fq_nmod_reduce(S, d + len - 1, ctx);
            T = R;
            R = S;
            S = T;
        }

        while (bit--)
        {
            if (fmpz_tstbit(e, bit))
            {
                _nmod_poly_mul(S, R, d, R, d, ctx->mod);
                _fq_nmod_reduce(S, 2 * d - 1, ctx);
                _nmod_poly_mul(R, S, d, op, len, ctx->mod);
                _fq_nmod_reduce(R, d + len - 1, ctx);
            }
            else
            {
                _nmod_poly_mul(S, R, d, R, d, ctx->mod);
                _fq_nmod_reduce(S, 2 * d - 1, ctx);
                T = R;
                R = S;
                S = T;
            }
        }

        _nmod_vec_clear(v);
    }
}

void fq_nmod_pow(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t e, const fq_nmod_ctx_t ctx)
{
    if (fmpz_sgn(e) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fq_nmod_pow).  e < 0.\n");
    }

    if (fmpz_is_zero(e))
    {
        fq_nmod_one(rop, ctx);
    }
    else if (fq_nmod_is_zero(op, ctx))
    {
        fq_nmod_zero(rop, ctx);
    }
    else if (fmpz_is_one(e))
    {
        fq_nmod_set(rop, op, ctx);
    }
    else
    {
        const slong d = fq_nmod_ctx_degree(ctx);
        mp_limb_t *t;

        if (rop == op)
        {
            t = _nmod_vec_init(2 * d - 1);
        }
        else
        {
            nmod_poly_fit_length(rop, 2 * d - 1);
            t = rop->coeffs;
        }

        if (fmpz_cmpabs(e, fq_nmod_ctx_prime(ctx)) < 0)
            _fq_nmod_pow(t, op->coeffs, op->length, e, ctx);
        else
        {
            fmpz_t order, e_mod;
            fmpz_init(e_mod);
            fmpz_init(order);
            fq_nmod_ctx_order(order, ctx);
            fmpz_sub_ui(order, order, 1);
            fmpz_mod(e_mod, e, order);
            _fq_nmod_pow(t, op->coeffs, op->length, e_mod, ctx);
            fmpz_clear(order);
            fmpz_clear(e_mod);
        }

        if (rop == op)
        {
            _nmod_vec_clear(rop->coeffs);
            rop->coeffs = t;
            rop->alloc  = 2 * d - 1;
            rop->length = d;
        }
        else
        {
            _nmod_poly_set_length(rop, d);
        }
        _nmod_poly_normalise(rop);
    }
}


/* TODO: Move into separate function / optimize */
void fq_nmod_pow_ui(fq_nmod_t rop, const fq_nmod_t op, const ulong e, const fq_nmod_ctx_t ctx)
{
    fmpz_t exp;
    fmpz_init_set_ui(exp, e);
    fq_nmod_pow(rop, op, exp, ctx);
    fmpz_clear(exp);
}
