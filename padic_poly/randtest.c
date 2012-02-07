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

    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "flint.h"
#include "padic_poly.h"

void padic_poly_randtest_val(padic_poly_t f, flint_rand_t state, 
                             long val, long len, const padic_ctx_t ctx)
{
    if (len == 0)
        return;

    if (val >= ctx->N)
    {
        padic_poly_zero(f);
    }
    else
    {
        long i;
        fmpz_t pow;
        int alloc;

        f->val = val;

        padic_poly_fit_length(f, len);

        alloc = _padic_ctx_pow_ui(pow, ctx->N - f->val, ctx);

        for (i = 0; i < len; i++)
            fmpz_randm(f->coeffs + i, state, pow);

        if (alloc)
            fmpz_clear(pow);

        _padic_poly_set_length(f, len);
        _padic_poly_normalise(f);

        padic_poly_canonicalise(f, ctx->p);
        padic_poly_reduce(f, ctx);
    }
}

void padic_poly_randtest(padic_poly_t f, flint_rand_t state, 
                         long len, const padic_ctx_t ctx)
{
    long min, max, val;

    if (ctx->N > 0)
    {
        min = - ((ctx->N + 9) / 10);
        max = ctx->N;
    }
    else if (ctx->N < 0)
    {
        min = ctx->N - ((-ctx->N + 9) / 10);
        max = ctx->N;
    }
    else  /* ctx->N == 0 */
    {
        min = -10;
        max = 0;
    }

    val = n_randint(state, max - min) + min;

    padic_poly_randtest_val(f, state, val, len, ctx);
}

void padic_poly_randtest_not_zero(padic_poly_t f, flint_rand_t state, 
                                  long len, const padic_ctx_t ctx)
{
    long i;

    if (len == 0)
    {
        printf("Exception (padic_poly_randtest_not_zero).  len == 0.\n");
        abort();
    }

    padic_poly_randtest(f, state, len, ctx);
    for (i = 0; !padic_poly_is_zero(f) && (i < 10); i++)
        padic_poly_randtest(f, state, len, ctx);

    if (padic_poly_is_zero(f))
    {
        padic_poly_fit_length(f, 1);
        _padic_poly_set_length(f, 1);
        fmpz_set_ui(f->coeffs + 0, 1);
        f->val = ctx->N - 1;
    }
}

