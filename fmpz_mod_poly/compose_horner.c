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

void _fmpz_mod_poly_compose_horner(fmpz *res, const fmpz *poly1, slong len1, 
                                              const fmpz *poly2, slong len2, 
                                              const fmpz_t p)
{
    if (len1 == 1 || len2 == 0)
    {
        fmpz_set(res, poly1);
    }
    else
    {
        const slong alloc = (len1 - 1) * (len2 - 1) + 1;
        slong i = len1 - 1, lenr = len2;
        fmpz * t = _fmpz_vec_init(alloc);
        
        /*
           Perform the first two steps as one, 
             "res = a(m) * poly2 + a(m-1)".
         */
        {
            _fmpz_mod_poly_scalar_mul_fmpz(res, poly2, len2, poly1 + i, p);
            i--;
            fmpz_add(res, res, poly1 + i);
            if (fmpz_cmpabs(res, p) >= 0)
                fmpz_sub(res, res, p);
        }
        while (i > 0)
        {
            i--;
            _fmpz_mod_poly_mul(t, res, lenr, poly2, len2, p);
            lenr += len2 - 1;
            _fmpz_mod_poly_add(res, t, lenr, poly1 + i, 1, p);
        }

        _fmpz_vec_clear(t, alloc);
    }
}

void fmpz_mod_poly_compose_horner(fmpz_mod_poly_t res,
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
            _fmpz_mod_poly_compose_horner(res->coeffs, poly1->coeffs, len1, 
                               poly2->coeffs, len2, fmpz_mod_ctx_modulus(ctx));
        }
        else
        {
            fmpz *t = _fmpz_vec_init(lenr);

            _fmpz_mod_poly_compose_horner(t, poly1->coeffs, len1,
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
