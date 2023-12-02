/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_preinvert(fmpz * Binv, const fmpz * B, slong len)
{
    if (len == 1)  /* B is +-1 */
    {
        fmpz_set(Binv, B);
    }
    else
    {
        const slong alloc = len + FLINT_MAX(len, 3 * FMPZ_POLY_INV_NEWTON_CUTOFF);
        slong *a, i, m, n = len;
        fmpz *T, *W;

        T = _fmpz_vec_init(alloc);
        W = T + len;

        for (i = 1; (WORD(1) << i) < n; i++) ;

        a = (slong *) flint_malloc(i * sizeof(slong));
        a[i = 0] = n;
        while (n >= FMPZ_POLY_INV_NEWTON_CUTOFF)
            a[++i] = (n = (n + 1) / 2);

        if (len != n)
           _fmpz_poly_reverse(T, B, len, len); /* only reverse if it ... */

        /* Base case */
        {
            fmpz *Brev = W + 2 * FMPZ_POLY_INV_NEWTON_CUTOFF;
            if (len != n)
               _fmpz_poly_reverse(Brev, T, n, n); /* ... won't be undone */
            else Brev = (fmpz *) B;

            _fmpz_vec_zero(W, 2*n - 2);
            fmpz_one(W + (2*n - 2));
            _fmpz_poly_div_basecase(Binv, W, W, 2*n - 1, Brev, n, 0);
            _fmpz_poly_reverse(Binv, Binv, n, n);
        }

        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            _fmpz_poly_mullow(W, T, n, Binv, m, n);
            _fmpz_poly_mullow(Binv + m, Binv, m, W + m, n - m, n - m);
            _fmpz_vec_neg(Binv + m, Binv + m, n - m);
        }

        _fmpz_vec_clear(T, alloc);
        flint_free(a);
    }
}

void
fmpz_poly_preinvert(fmpz_poly_t B_inv, const fmpz_poly_t B)
{
    slong n = B->length;
    fmpz_poly_t temp;
    fmpz * Binv_coeffs;

    if (n == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_preinvert). Division by zero.\n");
    }

    if (B == B_inv)
    {
       fmpz_poly_init2(temp, n);
       Binv_coeffs = temp->coeffs;
    } else
    {
       fmpz_poly_fit_length(B_inv, n);
       Binv_coeffs = B_inv->coeffs;
    }

    _fmpz_poly_preinvert(Binv_coeffs, B->coeffs, n);


    if (B == B_inv)
    {
       _fmpz_poly_set_length(temp, n);
       fmpz_poly_swap(B_inv, temp);
       fmpz_poly_clear(temp);
    } else
       _fmpz_poly_set_length(B_inv, n);

    /* no need to normalise */
}
