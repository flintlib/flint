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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#define FLINT_NEWTON_INV_CUTOFF  32

void 
_fmpz_poly_newton_inv(fmpz * Qinv, fmpz * temp, const fmpz * Q, long n)
{
    long *a, i;

    if (n == 1)
    {
        fmpz_set_ui(Qinv, 1);
        return;
    }

    _fmpz_vec_zero(Qinv, n);

    for (i = 1; (1L << i) < n; i++) ;

    a = (long *) malloc(i * sizeof(long));
    a[i = 0] = n;
    while (n >= FLINT_NEWTON_INV_CUTOFF)
        a[++i] = (n = (n + 1) / 2);

    /* Base case */
    {
        fmpz *Qrev = temp + 2 * FLINT_NEWTON_INV_CUTOFF;

        _fmpz_poly_reverse(Qrev, Q, n, n);
        _fmpz_vec_zero(temp, 2*n - 2);
        fmpz_set_ui(temp + (2*n - 2), 1);
        _fmpz_poly_div_basecase(Qinv, temp, temp, 2*n - 1, Qrev, n);
        _fmpz_poly_reverse(Qinv, Qinv, n, n);
    }
    
    for (i--; i >= 0; i--)
    {
        n = a[i];

        _fmpz_poly_mullow_n(temp, Q, Qinv, n);
        fmpz_sub_ui(temp, temp, 1);
        _fmpz_poly_mullow_n(temp + n, temp, Qinv, n);
        _fmpz_vec_sub(Qinv, Qinv, temp + n, n);
    }

    free(a);
}

void fmpz_poly_newton_inv(fmpz_poly_t Qinv, const fmpz_poly_t Q, long n)
{
    fmpz *Qcopy;
    fmpz *T = _fmpz_vec_init(FLINT_MAX(2 * n, 3 * FLINT_NEWTON_INV_CUTOFF));
    int Qalloc;

    if (Q->length >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        long i;
        Qcopy = (fmpz *) malloc(n * sizeof(fmpz));
        for (i = 0; i < Q->length; i++)
            Qcopy[i] = Q->coeffs[i];
        for ( ; i < n; i++)
            Qcopy[i] = 0;
        Qalloc = 1;
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_newton_inv(Qinv->coeffs, T, Qcopy, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_newton_inv(t->coeffs, T, Qcopy, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }
    
    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
    _fmpz_vec_clear(T, FLINT_MAX(2 * n, 3 * FLINT_NEWTON_INV_CUTOFF));

    if (Qalloc)
        free(Qcopy);
}

