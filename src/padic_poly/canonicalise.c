/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

void _padic_poly_canonicalise(fmpz *poly, slong *v, slong len, const fmpz_t p)
{
    const slong min = _fmpz_vec_ord_p(poly, len, p);

    if (min == 0)
    {
        if (_fmpz_vec_is_zero(poly, len))
            *v = 0;
    }
    else  /* min > 0 */
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, p, min);
        _fmpz_vec_scalar_divexact_fmpz(poly, poly, len, pow);
        fmpz_clear(pow);

        *v += min;
    }
}

void padic_poly_canonicalise(padic_poly_t poly, const fmpz_t p)
{
    _padic_poly_canonicalise(poly->coeffs, &(poly->val), poly->length, p);
}

