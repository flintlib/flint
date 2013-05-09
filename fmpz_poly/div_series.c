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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void 
_fmpz_poly_div_series(fmpz * Q, const fmpz * A, const fmpz * B, len_t n)
{
    if (n == 1)
    {
        _fmpz_vec_set(Q, A, n);
    }
    else
    {
        fmpz * Binv = _fmpz_vec_init(n);

        _fmpz_poly_inv_series(Binv, B, n);
        _fmpz_poly_mullow(Q, A, n, Binv, n, n);

        _fmpz_vec_clear(Binv, n);
    }
}

void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A, 
                                         const fmpz_poly_t B, len_t n)
{
    fmpz *a, *b;
    ulong flags = 0UL;  /* 2^0 for a, 2^1 for b */

    if (Q == A)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_div_series(t, A, B, n);
        fmpz_poly_swap(Q, t);
        fmpz_poly_clear(t);
        return;
    }

    fmpz_poly_fit_length(Q, n);

    if (A->length >= n)
    {
        a = A->coeffs;
    }
    else
    {
        len_t i;
        a = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < A->length; i++)
            a[i] = A->coeffs[i];
        flint_mpn_zero((mp_ptr) a + A->length, n - A->length);
        flags |= 1UL;
    }

    if (B->length >= n)
    {
        b = B->coeffs;
    }
    else
    {
        len_t i;
        b = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < B->length; i++)
            b[i] = B->coeffs[i];
        flint_mpn_zero((mp_ptr) b + B->length, n - B->length);
        flags |= 2UL;
    }

    _fmpz_poly_div_series(Q->coeffs, a, b, n);

    _fmpz_poly_set_length(Q, n);
    _fmpz_poly_normalise(Q);

    if ((flags & 1UL))
        flint_free(a);
    if ((flags & 2UL))
        flint_free(b);
}

