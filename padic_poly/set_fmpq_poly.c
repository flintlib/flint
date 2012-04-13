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

#include "fmpz_vec.h"
#include "fmpq_poly.h"
#include "padic_poly.h"

void padic_poly_set_fmpq_poly(padic_poly_t f, 
                              const fmpq_poly_t g, const padic_ctx_t ctx)
{
    const long len = g->length;

    if (len == 0)
    {
        padic_poly_zero(f);
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);

        f->val = - fmpz_remove(t, g->den, ctx->p);

        if (f->val < ctx->N)
        {
            padic_poly_fit_length(f, len);
            _padic_poly_set_length(f, len);

            _padic_inv(t, t, ctx->p, ctx->N - f->val);
            _fmpz_vec_scalar_mul_fmpz(f->coeffs, g->coeffs, len, t);
            if (f->val == 0)
                padic_poly_canonicalise(f, ctx->p);

            padic_poly_reduce(f, ctx);
        }
        else
        {
            padic_poly_zero(f);
        }

        fmpz_clear(t);
    }
}

