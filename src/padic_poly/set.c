/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

void padic_poly_set(padic_poly_t poly1, 
                    const padic_poly_t poly2, const padic_ctx_t ctx)
{
    if (poly1 != poly2)         /* Aliasing is trivial */
    {
        slong len2 = poly2->length, N1 = padic_poly_prec(poly1);

        if (len2 == 0 || poly2->val >= N1)
        {
            padic_poly_zero(poly1);
        }
        else 
        {
            padic_poly_fit_length(poly1, len2);
            _padic_poly_set_length(poly1, len2);
            poly1->val = poly2->val;

            if (N1 >= padic_poly_prec(poly2))  /* No reduction */
            {
                _fmpz_vec_set(poly1->coeffs, poly2->coeffs, len2);
            }
            else  /* Reduction necessary */
            {
                fmpz_t pow;
                int alloc;

                alloc = _padic_ctx_pow_ui(pow, N1 - poly1->val, ctx);

                _fmpz_vec_scalar_mod_fmpz(poly1->coeffs, poly2->coeffs, len2, pow);

                if (alloc)
                    fmpz_clear(pow);

                _padic_poly_normalise(poly1);
                /* Length cannot be zero, so no need to check */
                /* if (poly->length == 0)                     */
                /*     poly->val = 0;                         */
            }
        }
    }
}

