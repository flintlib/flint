/*============================================================================

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

    Copyright (C) 2011 Sebastian Pancratz
 
******************************************************************************/

#include "fmpz_mod_poly.h"
#include "padic_poly.h"

void padic_poly_inv_series(padic_poly_t Qinv, const padic_poly_t Q, long n, 
                           const padic_ctx_t ctx)
{
    fmpz_t cinv;

    fmpz_t pow;
    int palloc;

    fmpz *Qcopy;
    int Qalloc;

    if (Q->length == 0 || fmpz_is_zero(Q->coeffs + 0))
    {
        printf("Exception (padic_poly_inv_series):  Constant term is zero.\n");
        abort();
    }
    if (fmpz_divisible(Q->coeffs, ctx->p))
    {
        printf("Exception (padic_poly_inv_series):\n");
        printf("Valuation of constant term is not minimal.\n");
        abort();
    }

    if (- Q->val >= ctx->N)
    {
        padic_poly_zero(Qinv);
        return;
    }

    if (Q->length >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        long i;

        Qcopy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < Q->length; i++)
            Qcopy[i] = Q->coeffs[i];
        mpn_zero((mp_ptr) Qcopy + i, n - i);
        Qalloc = 1;
    }

    fmpz_init(cinv);
    fmpz_init(pow);

    _padic_inv(cinv, Q->coeffs, ctx->p, ctx->N + Q->val);
    palloc = _padic_ctx_pow_ui(pow, ctx->N + Q->val, ctx);

    if (Qinv != Q)
    {
        padic_poly_fit_length(Qinv, n);
        _fmpz_mod_poly_inv_series_newton(Qinv->coeffs, Qcopy, n, cinv, pow);
    }
    else
    {
        fmpz *t = _fmpz_vec_init(n);

        _fmpz_mod_poly_inv_series_newton(t, Qcopy, n, cinv, pow);

        _fmpz_vec_clear(Qinv->coeffs, Qinv->alloc);
        Qinv->coeffs = t;
        Qinv->alloc  = n;
        Qinv->length = n;
    }

    Qinv->val = - Q->val;

    _padic_poly_set_length(Qinv, n);
    _padic_poly_normalise(Qinv);

    fmpz_clear(cinv);
    if (palloc)
        fmpz_clear(pow);
    if (Qalloc)
        flint_free(Qcopy);
}

