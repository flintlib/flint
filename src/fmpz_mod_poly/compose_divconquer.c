/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz_mod_poly.h"

/*
    Let i be such that 2^{i} < len1 <= 2^{i+1}.

    Note that the jth step of the recursion requires temporary space 
    of size no more than (len2 - 1)(2^j - 1) + 1.  Note the smallest 
    step j=0 doesn't require any temporary space and the largest step 
    has j = i, and hence the sum is 

        sum_{j=1}^i [(len2 - 1) (2^j - 1) + 1]

      = (len2 - 1)(2^{i+1} - 2) - (len2 - 2) i
 */

void _fmpz_mod_poly_compose_divconquer_recursive(fmpz *res, 
    const fmpz *poly1, slong len1, fmpz **pow2, slong len2, fmpz *v, 
    const fmpz_t p)
{
    if (len1 == 1)
    {
        fmpz_set(res, poly1);
    }
    else if (len1 == 2)
    {
        _fmpz_mod_poly_scalar_mul_fmpz(res, pow2[0], len2, poly1 + 1, p);
        fmpz_add(res, res, poly1);
        if (fmpz_cmpabs(res, p) >= 0)
            fmpz_sub(res, res, p);
    }
    else
    {
        const slong i = FLINT_BIT_COUNT(len1 - 1) - 1;
        fmpz *w = v + ((WORD(1) << i) - 1) * (len2 - 1) + 1;

        _fmpz_mod_poly_compose_divconquer_recursive(v, 
            poly1 + (WORD(1) << i), len1 - (WORD(1) << i), pow2, len2, w, p);

        _fmpz_mod_poly_mul(res, pow2[i], (len2 - 1) * (WORD(1) << i) + 1, 
                                v, (len2 - 1) * (len1 - (WORD(1) << i) - 1) + 1, p);

        _fmpz_mod_poly_compose_divconquer_recursive(v, poly1, WORD(1) << i, 
                                                       pow2, len2, w, p);

        _fmpz_mod_poly_add(res, res, (len2 - 1) * ((WORD(1) << i) - 1) + 1, 
                                  v, (len2 - 1) * ((WORD(1) << i) - 1) + 1, p);
    }
}

void _fmpz_mod_poly_compose_divconquer(fmpz *res, 
                                       const fmpz *poly1, slong len1, 
                                       const fmpz *poly2, slong len2, 
                                       const fmpz_t p)
{
    if (len1 == 1 || len2 == 0)
    {
        fmpz_set(res, poly1);
    }
    else
    {
        const slong k = FLINT_BIT_COUNT(len1 - 1);
        const slong lenV = len2 * ((WORD(1) << k) - 1) + k;
        const slong lenW = (len2 - 1) * ((WORD(1) << k) - 2) - (len2 - 2) * (k-1);
        slong i;
        fmpz *v, *w, **pow2;

        v    = _fmpz_vec_init(lenV + lenW);
        w    = v + lenV;
        pow2 = flint_malloc(k * sizeof(fmpz *));

        for (i = 0; i < k; i++)
        {
            pow2[i] = v + (len2 * ((WORD(1) << i) - 1) + i);
        }

        _fmpz_vec_set(pow2[0], poly2, len2);
        for (i = 1; i < k; i++)
        {
            _fmpz_mod_poly_sqr(pow2[i], 
                               pow2[i-1], (len2 - 1) * (WORD(1) << (i - 1)) + 1, p);
        }

        _fmpz_mod_poly_compose_divconquer_recursive(res, poly1, len1, 
                                                         pow2, len2, w, p);

        _fmpz_vec_clear(v, lenV + lenW);
        flint_free(pow2);
    }
}

void fmpz_mod_poly_compose_divconquer(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0)
    {
        fmpz_mod_poly_zero(res, ctx);
    }
    else if (len1 == 1 || len2 == 0)
    {
        fmpz_mod_poly_set_fmpz(res, poly1->coeffs, ctx);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;

        if ((res != poly1) && (res != poly2))
        {
            fmpz_mod_poly_fit_length(res, lenr, ctx);
            _fmpz_mod_poly_compose_divconquer(res->coeffs, poly1->coeffs, len1, 
                               poly2->coeffs, len2, fmpz_mod_ctx_modulus(ctx));
        }
        else
        {
            fmpz *t = _fmpz_vec_init(lenr);

            _fmpz_mod_poly_compose_divconquer(t, poly1->coeffs, len1,
                               poly2->coeffs, len2, fmpz_mod_ctx_modulus(ctx));
            _fmpz_vec_clear(res->coeffs, res->alloc);
            res->coeffs = t;
            res->alloc  = lenr;
            res->length = lenr;
        }

        _fmpz_mod_poly_set_length(res, lenr);
        _fmpz_mod_poly_normalise(res);
    }
}

