/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

#define ZASSENHAUS 0
#define BERLEKAMP 1
#define KALTOFEN 2

static __inline__ void
__fq_poly_factor1(fq_poly_factor_t res, const fq_poly_t f, int algorithm,
                  const fq_ctx_t ctx)
{
    if (algorithm == KALTOFEN)
        fq_poly_factor_kaltofen_shoup(res, f, ctx);
    else if (algorithm == ZASSENHAUS)
        fq_poly_factor_cantor_zassenhaus(res, f, ctx);
    else
        fq_poly_factor_berlekamp(res, f, ctx);
}

void
__fq_poly_factor(fq_poly_factor_t result, fq_t leading_coeff,
                 const fq_poly_t input, int algorithm, const fq_ctx_t ctx)
{
    fq_poly_t monic_input;
    fq_poly_factor_t sqfree_factors, factors;
    slong i, len;

    len = input->length;

    if (len <= 1)
    {
        if (len == 0)
        {
            fq_zero(leading_coeff, ctx);
            return;
        }
        else
        {
            fq_set(leading_coeff, input->coeffs, ctx);
        }
    }

    fq_poly_get_coeff(leading_coeff, input, fq_poly_degree(input, ctx), ctx);

    fq_poly_init(monic_input, ctx);
    fq_poly_make_monic(monic_input, input, ctx);

    if (len == 2)
    {
        fq_poly_factor_insert(result, monic_input, 1, ctx);
        fq_poly_clear(monic_input, ctx);
        fq_set(leading_coeff, input->coeffs + 1, ctx);
        return;
    }

    fq_poly_factor_init(sqfree_factors, ctx);
    fq_poly_factor_squarefree(sqfree_factors, monic_input, ctx);
    fq_poly_clear(monic_input, ctx);

    /* Run CZ on each of the square-free factors */
    for (i = 0; i < sqfree_factors->num; i++)
    {
        fq_poly_factor_init(factors, ctx);
        __fq_poly_factor1(factors, sqfree_factors->poly + i, algorithm, ctx);
        fq_poly_factor_pow(factors, sqfree_factors->exp[i], ctx);
        fq_poly_factor_concat(result, factors, ctx);
        fq_poly_factor_clear(factors, ctx);
    }

    fq_poly_factor_clear(sqfree_factors, ctx);
    
    return;
}

void
__fq_poly_factor_deflation(fq_poly_factor_t result, fq_t leading_coeff,
                           const fq_poly_t input, int algorithm,
                           const fq_ctx_t ctx)
{
    slong i;
    ulong deflation;

    if (input->length <= 1)
    {
        if (input->length == 0)
        {
            fq_zero(leading_coeff, ctx);
            return;
        }
        else
        {
            fq_set(leading_coeff, input->coeffs, ctx);
            return;
        }
    }

    deflation = fq_poly_deflation(input, ctx);
    if (deflation == 1)
    {
        __fq_poly_factor(result, leading_coeff, input, algorithm, ctx);
        return;
    }
    else
    {
        fq_poly_factor_t def_res;
        fq_poly_t def;

        fq_poly_init(def, ctx);
        fq_poly_deflate(def, input, deflation, ctx);
        fq_poly_factor_init(def_res, ctx);
        __fq_poly_factor(def_res, leading_coeff, def, algorithm, ctx);
        fq_poly_clear(def, ctx);

        for (i = 0; i < def_res->num; i++)
        {
            /* Inflate */
            fq_poly_t pol;
            fq_poly_init(pol, ctx);
            fq_poly_inflate(pol, def_res->poly + i, deflation, ctx);

            /* Factor inflation */
            if (def_res->exp[i] == 1)
                __fq_poly_factor(result, leading_coeff, pol, algorithm, ctx);
            else
            {
                fq_poly_factor_t t;
                fq_poly_factor_init(t, ctx);
                __fq_poly_factor(t, leading_coeff, pol, algorithm, ctx);
                fq_poly_factor_pow(t, def_res->exp[i], ctx);
                fq_poly_factor_concat(result, t, ctx);
                fq_poly_factor_clear(t, ctx);
            }
            fq_poly_clear(pol, ctx);
        }

        fq_poly_factor_clear(def_res, ctx);
    }
}

void
fq_poly_factor_with_berlekamp(fq_poly_factor_t result, fq_t leading_coeff,
                              const fq_poly_t input, const fq_ctx_t ctx)
{
    __fq_poly_factor_deflation(result, leading_coeff, input, BERLEKAMP, ctx);
}

void
fq_poly_factor_with_cantor_zassenhaus(fq_poly_factor_t result, fq_t leading_coeff,
                                      const fq_poly_t input, const fq_ctx_t ctx)
{
    __fq_poly_factor_deflation(result, leading_coeff, input, ZASSENHAUS, ctx);
}

void
fq_poly_factor_with_kaltofen_shoup(fq_poly_factor_t result, fq_t leading_coeff,
                                   const fq_poly_t input, const fq_ctx_t ctx)
{
    __fq_poly_factor_deflation(result, leading_coeff, input, KALTOFEN, ctx);
}

void
fq_poly_factor(fq_poly_factor_t result, fq_t leading_coeff,
               const fq_poly_t input, const fq_ctx_t ctx)
{
    mp_bitcnt_t bits = fmpz_bits(fq_ctx_prime(ctx));
    slong n = fq_poly_degree(input, ctx);

    if (n < 10 + 50 / bits)
        __fq_poly_factor_deflation(result, leading_coeff, input,
                                   ZASSENHAUS, ctx);
    else
        __fq_poly_factor_deflation(result, leading_coeff, input,
                                   KALTOFEN, ctx);
}
