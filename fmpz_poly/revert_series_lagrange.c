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
   Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"


void
_fmpz_poly_revert_series_lagrange(fmpz * Qinv, const fmpz * Q, long n)
{
    long i;
    fmpz *R, *S, *T, *tmp;

    if (n <= 2)
    {
        _fmpz_vec_set(Qinv, Q, n);
        return;
    }

    R = _fmpz_vec_init(n - 1);
    S = _fmpz_vec_init(n - 1);
    T = _fmpz_vec_init(n - 1);

    fmpz_zero(Qinv);
    fmpz_set(Qinv + 1, Q + 1);

    _fmpz_poly_inv_series(R, Q + 1, n - 1);
    _fmpz_vec_set(S, R, n - 1);

    for (i = 2; i < n; i++)
    {
        _fmpz_poly_mullow(T, S, n - 1, R, n - 1, n - 1);
        fmpz_divexact_ui(Qinv + i, T + i - 1, i);
        tmp = S; S = T; T = tmp;
    }

    _fmpz_vec_clear(R, n - 1);
    _fmpz_vec_clear(S, n - 1);
    _fmpz_vec_clear(T, n - 1);
}

void
fmpz_poly_revert_series_lagrange(fmpz_poly_t Qinv,
                                        const fmpz_poly_t Q, long n)
{
    fmpz *Qcopy;
    int Qalloc;
    long Qlen = Q->length;

    if (Qlen < 2 || !fmpz_is_zero(Q->coeffs) || !fmpz_is_pm1(Q->coeffs + 1))
    {
        printf("Exception (fmpz_poly_revert_series_lagrange). Input must have \n"
               "zero constant term and +1 or -1 as coefficient of x^1.\n");
        abort();
    }

    if (Qlen >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        long i;
        Qcopy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < Qlen; i++)
            Qcopy[i] = Q->coeffs[i];
        for ( ; i < n; i++)
            Qcopy[i] = 0;
        Qalloc = 1;
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_revert_series_lagrange(Qinv->coeffs, Qcopy, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_revert_series_lagrange(t->coeffs, Qcopy, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }
    
    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);

    if (Qalloc)
        flint_free(Qcopy);
}
