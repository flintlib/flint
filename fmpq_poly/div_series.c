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
#include "fmpq_poly.h"

void 
_fmpq_poly_div_series(fmpz * Q, fmpz_t denQ, 
                      const fmpz * A, const fmpz_t denA, 
                      const fmpz * B, const fmpz_t denB, long n)
{
    fmpz * C = _fmpz_vec_init(n + 1);
    fmpz * denC = C + n;

    _fmpq_poly_inv_series(C, denC, B, denB, n);
    _fmpq_poly_mullow(Q, denQ, A, denA, n, C, denC, n, n);

    _fmpz_vec_clear(C, n + 1);
}

void fmpq_poly_div_series(fmpq_poly_t Q, const fmpq_poly_t A, 
                                         const fmpq_poly_t B, long n)
{
    fmpz *a, *b;
    ulong flags = 0UL;  /* 2^0 for a, 2^1 for b */

    if (Q == A)
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        fmpq_poly_div_series(t, A, B, n);
        fmpq_poly_swap(Q, t);
        fmpq_poly_clear(t);
        return;
    }

    fmpq_poly_fit_length(Q, n);

    if (A->length >= n)
    {
        a = A->coeffs;
    }
    else
    {
        long i;
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
        long i;
        b = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < B->length; i++)
            b[i] = B->coeffs[i];
        flint_mpn_zero((mp_ptr) b + B->length, n - B->length);
        flags |= 2UL;
    }

    _fmpq_poly_div_series(Q->coeffs, Q->den, a, A->den, b, B->den, n);

    _fmpq_poly_set_length(Q, n);
    fmpq_poly_canonicalise(Q);

    if ((flags & 1UL))
        flint_free(a);
    if ((flags & 2UL))
        flint_free(b);
}

