/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010, 2022 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

#define ZASSENHAUS 0
#define BERLEKAMP 1
#define KALTOFEN 2

static inline void
__nmod_poly_factor1(nmod_poly_factor_t res, const nmod_poly_t f, int algorithm)
{
    if (algorithm == KALTOFEN)
        nmod_poly_factor_kaltofen_shoup(res, f);
    else if (algorithm == ZASSENHAUS)
        nmod_poly_factor_cantor_zassenhaus(res, f);
    else
        nmod_poly_factor_berlekamp(res, f);
}

mp_limb_t
__nmod_poly_factor(nmod_poly_factor_t result,
                                    const nmod_poly_t input, int algorithm)
{
    nmod_poly_t monic_input;
    nmod_poly_factor_t sqfree_factors, factors;
    mp_limb_t leading_coeff;
    slong i, len;

    len = input->length;

    if (len <= 1)
    {
        if (len == 0)
            return 0;
        else
            return input->coeffs[0];
    }

    leading_coeff = *nmod_poly_lead(input);

    nmod_poly_init_mod(monic_input, input->mod);
    nmod_poly_make_monic(monic_input, input);

    if (len == 2)
    {
        nmod_poly_factor_insert(result, monic_input, 1);
        nmod_poly_clear(monic_input);
        return input->coeffs[1];
    }

    nmod_poly_factor_init(sqfree_factors);
    nmod_poly_factor_squarefree(sqfree_factors, monic_input);
    nmod_poly_clear(monic_input);

    /* Run CZ on each of the square-free factors */
    for (i = 0; i < sqfree_factors->num; i++)
    {
        nmod_poly_factor_init(factors);

        __nmod_poly_factor1(factors, sqfree_factors->p + i, algorithm);
        nmod_poly_factor_pow(factors, sqfree_factors->exp[i]);
        nmod_poly_factor_concat(result, factors);

        nmod_poly_factor_clear(factors);
    }

    nmod_poly_factor_clear(sqfree_factors);
    return leading_coeff;
}

mp_limb_t
__nmod_poly_factor_deflation(nmod_poly_factor_t result,
    const nmod_poly_t input, int algorithm)
{
    slong i;
    ulong deflation;

    if (input->length <= 1)
    {
        if (input->length == 0)
            return 0;
        else
            return input->coeffs[0];
    }

    deflation = nmod_poly_deflation(input);
    if (deflation == 1)
    {
        return __nmod_poly_factor(result, input, algorithm);
    }
    else
    {
        nmod_poly_factor_t def_res;
        nmod_poly_t def;
        mp_limb_t leading_coeff;

        nmod_poly_init_mod(def, input->mod);

        nmod_poly_deflate(def, input, deflation);
        nmod_poly_factor_init(def_res);
        leading_coeff = __nmod_poly_factor(def_res, def, algorithm);

        nmod_poly_clear(def);

        for (i = 0; i < def_res->num; i++)
        {
            /* Inflate */
            nmod_poly_t pol;

            nmod_poly_init_mod(pol, input->mod);

            nmod_poly_inflate(pol, def_res->p + i, deflation);

            /* Factor inflation */
            if (def_res->exp[i] == 1)
                __nmod_poly_factor(result, pol, algorithm);
            else
            {
                nmod_poly_factor_t t;

                nmod_poly_factor_init(t);

                __nmod_poly_factor(t, pol, algorithm);
                nmod_poly_factor_pow(t, def_res->exp[i]);
                nmod_poly_factor_concat(result, t);

                nmod_poly_factor_clear(t);
            }

            nmod_poly_clear(pol);
        }

        nmod_poly_factor_clear(def_res);

        return leading_coeff;
    }
}

mp_limb_t
nmod_poly_factor_with_berlekamp(nmod_poly_factor_t result,
    const nmod_poly_t input)
{
    return __nmod_poly_factor_deflation(result, input, BERLEKAMP);
}

mp_limb_t
nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t result,
    const nmod_poly_t input)
{
    return __nmod_poly_factor_deflation(result, input, ZASSENHAUS);
}

mp_limb_t
nmod_poly_factor_with_kaltofen_shoup(nmod_poly_factor_t result,
    const nmod_poly_t input)
{
    return __nmod_poly_factor_deflation(result, input, KALTOFEN);
}

mp_limb_t
nmod_poly_factor(nmod_poly_factor_t result, const nmod_poly_t input)
{
    mp_limb_t p = input->mod.n;
    unsigned int bits = FLINT_BIT_COUNT (p);
    slong n = nmod_poly_degree(input);

    result->num = 0;

    if (n < 10 + 50 / bits)
        return __nmod_poly_factor_deflation(result, input, ZASSENHAUS);
    else
        return __nmod_poly_factor_deflation(result, input, KALTOFEN);
}
