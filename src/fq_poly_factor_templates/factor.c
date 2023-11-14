/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#define ZASSENHAUS 0
#define BERLEKAMP 1
#define KALTOFEN 2

static inline void
__TEMPLATE(T, poly_factor1) (TEMPLATE(T, poly_factor_t) res,
                             const TEMPLATE(T, poly_t) f, int algorithm,
                             const TEMPLATE(T, ctx_t) ctx)
{
    if (algorithm == KALTOFEN)
        TEMPLATE(T, poly_factor_kaltofen_shoup) (res, f, ctx);
    else if (algorithm == ZASSENHAUS)
        TEMPLATE(T, poly_factor_cantor_zassenhaus) (res, f, ctx);
    else
        TEMPLATE(T, poly_factor_berlekamp) (res, f, ctx);
}

void
__TEMPLATE(T, poly_factor) (TEMPLATE(T, poly_factor_t) result,
                            TEMPLATE(T, t) leading_coeff,
                            const TEMPLATE(T, poly_t) input,
                            int algorithm, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) monic_input;
    TEMPLATE(T, poly_factor_t) sqfree_factors, factors;
    slong i, len;

    len = input->length;

    if (len <= 1)
    {
        if (len == 0)
        {
            TEMPLATE(T, zero) (leading_coeff, ctx);
            return;
        }
        else
        {
            TEMPLATE(T, set) (leading_coeff, input->coeffs, ctx);
        }
    }

    TEMPLATE(T, poly_get_coeff) (leading_coeff, input,
                                 TEMPLATE(T, poly_degree) (input, ctx), ctx);

    TEMPLATE(T, poly_init) (monic_input, ctx);
    TEMPLATE(T, poly_make_monic) (monic_input, input, ctx);

    if (len == 2)
    {
        TEMPLATE(T, poly_factor_insert) (result, monic_input, 1, ctx);
        TEMPLATE(T, poly_clear) (monic_input, ctx);
        TEMPLATE(T, set) (leading_coeff, input->coeffs + 1, ctx);
        return;
    }

    TEMPLATE(T, poly_factor_init) (sqfree_factors, ctx);
    TEMPLATE(T, poly_factor_squarefree) (sqfree_factors, monic_input, ctx);
    TEMPLATE(T, poly_clear) (monic_input, ctx);

    /* Run CZ on each of the square-free factors */
    for (i = 0; i < sqfree_factors->num; i++)
    {
        TEMPLATE(T, poly_factor_init) (factors, ctx);
        __TEMPLATE(T, poly_factor1) (factors, sqfree_factors->poly + i,
                                     algorithm, ctx);
        TEMPLATE(T, poly_factor_pow) (factors, sqfree_factors->exp[i], ctx);
        TEMPLATE(T, poly_factor_concat) (result, factors, ctx);
        TEMPLATE(T, poly_factor_clear) (factors, ctx);
    }

    TEMPLATE(T, poly_factor_clear) (sqfree_factors, ctx);

    return;
}

void
__TEMPLATE(T, poly_factor_deflation) (TEMPLATE(T, poly_factor_t) result,
                                      TEMPLATE(T, t) leading_coeff,
                                      const TEMPLATE(T, poly_t) input,
                                      int algorithm,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    ulong deflation;

    if (input->length <= 1)
    {
        if (input->length == 0)
        {
            TEMPLATE(T, zero) (leading_coeff, ctx);
            return;
        }
        else
        {
            TEMPLATE(T, set) (leading_coeff, input->coeffs, ctx);
            return;
        }
    }

    deflation = TEMPLATE(T, poly_deflation) (input, ctx);
    if (deflation == 1)
    {
        __TEMPLATE(T, poly_factor) (result, leading_coeff, input, algorithm,
                                    ctx);
        return;
    }
    else
    {
        TEMPLATE(T, poly_factor_t) def_res;
        TEMPLATE(T, poly_t) def;
        TEMPLATE(T, t) lc_dummy;

        TEMPLATE(T, init) (lc_dummy, ctx);
        TEMPLATE(T, poly_init) (def, ctx);
        TEMPLATE(T, poly_deflate) (def, input, deflation, ctx);
        TEMPLATE(T, poly_factor_init) (def_res, ctx);
        __TEMPLATE(T, poly_factor) (def_res, leading_coeff, def, algorithm,
                                    ctx);
        TEMPLATE(T, poly_clear) (def, ctx);

        for (i = 0; i < def_res->num; i++)
        {
            /* Inflate */
            TEMPLATE(T, poly_t) pol;
            TEMPLATE(T, poly_init) (pol, ctx);
            TEMPLATE(T, poly_inflate) (pol, def_res->poly + i, deflation, ctx);

            /* Factor inflation */
            if (def_res->exp[i] == 1)
                __TEMPLATE(T, poly_factor) (result, lc_dummy, pol,
                                            algorithm, ctx);
            else
            {
                TEMPLATE(T, poly_factor_t) t;
                TEMPLATE(T, poly_factor_init) (t, ctx);
                __TEMPLATE(T, poly_factor) (t, lc_dummy, pol, algorithm,
                                            ctx);
                TEMPLATE(T, poly_factor_pow) (t, def_res->exp[i], ctx);
                TEMPLATE(T, poly_factor_concat) (result, t, ctx);
                TEMPLATE(T, poly_factor_clear) (t, ctx);
            }
            TEMPLATE(T, poly_clear) (pol, ctx);
        }

        TEMPLATE(T, clear) (lc_dummy, ctx);
        TEMPLATE(T, poly_factor_clear) (def_res, ctx);
    }
}

void
TEMPLATE(T, poly_factor_with_berlekamp) (TEMPLATE(T, poly_factor_t) result,
                                         TEMPLATE(T, t) leading_coeff,
                                         const TEMPLATE(T, poly_t) input,
                                         const TEMPLATE(T, ctx_t) ctx)
{
    __TEMPLATE(T, poly_factor_deflation) (result, leading_coeff, input,
                                          BERLEKAMP, ctx);
}

void
TEMPLATE(T, poly_factor_with_cantor_zassenhaus) (
    TEMPLATE(T, poly_factor_t) result,
    TEMPLATE(T, t) leading_coeff,
    const TEMPLATE(T, poly_t) input,
    const TEMPLATE(T, ctx_t) ctx)
{
    __TEMPLATE(T, poly_factor_deflation) (result, leading_coeff, input,
                                          ZASSENHAUS, ctx);
}

void
TEMPLATE(T, poly_factor_with_kaltofen_shoup) (
    TEMPLATE(T, poly_factor_t) result,
    TEMPLATE(T, t) leading_coeff,
    const TEMPLATE(T, poly_t) input,
    const TEMPLATE(T, ctx_t) ctx)
{
    __TEMPLATE(T, poly_factor_deflation) (result, leading_coeff, input,
                                          KALTOFEN, ctx);
}

void
TEMPLATE(T, poly_factor) (TEMPLATE(T, poly_factor_t) result,
                          TEMPLATE(T, t) leading_coeff,
                          const TEMPLATE(T, poly_t) input,
                          const TEMPLATE(T, ctx_t) ctx)
{
    flint_bitcnt_t bits = fmpz_bits(TEMPLATE(T, ctx_prime) (ctx));
    slong n = TEMPLATE(T, poly_degree) (input, ctx);

    result->num = 0;

    if (n < 10 + 50 / bits)
        __TEMPLATE(T, poly_factor_deflation) (result, leading_coeff, input,
                                              ZASSENHAUS, ctx);
    else
        __TEMPLATE(T, poly_factor_deflation) (result, leading_coeff, input,
                                              KALTOFEN, ctx);
}


#endif
