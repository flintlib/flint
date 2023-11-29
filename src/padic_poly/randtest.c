/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/


#include "flint.h"
#include "padic_poly.h"

void padic_poly_randtest_val(padic_poly_t f, flint_rand_t state,
                             slong val, slong len, const padic_ctx_t ctx)
{
    const slong N = padic_poly_prec(f);

    if (len == 0)
        return;

    if (val >= N)
    {
        padic_poly_zero(f);
    }
    else
    {
        slong i;
        fmpz_t pow;
        int alloc;

        f->val = val;

        padic_poly_fit_length(f, len);

        alloc = _padic_ctx_pow_ui(pow, N - f->val, ctx);

        for (i = 0; i < len; i++)
            fmpz_randm(f->coeffs + i, state, pow);

        if (alloc)
            fmpz_clear(pow);

        for (i = 0; i < len; i++)
            if (!fmpz_divisible(f->coeffs + i, ctx->p))
                break;
        if (i == len)
            fmpz_one(f->coeffs + n_randint(state, len));

        _padic_poly_set_length(f, len);
        _padic_poly_normalise(f);

        padic_poly_reduce(f, ctx);
    }
}

void padic_poly_randtest(padic_poly_t f, flint_rand_t state,
                         slong len, const padic_ctx_t ctx)
{
    const slong N = padic_poly_prec(f);

    slong min, max, val;

    if (N > 0)
    {
        min = - ((N + 9) / 10);
        max = N;
    }
    else if (N < 0)
    {
        min = N - ((-N + 9) / 10);
        max = N;
    }
    else  /* N == 0 */
    {
        min = -10;
        max = 0;
    }

    val = n_randint(state, max - min) + min;

    padic_poly_randtest_val(f, state, val, len, ctx);
}

void padic_poly_randtest_not_zero(padic_poly_t f, flint_rand_t state,
                                  slong len, const padic_ctx_t ctx)
{
    slong i;

    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (padic_poly_randtest_not_zero).  len == 0.\n");
    }

    padic_poly_randtest(f, state, len, ctx);
    for (i = 0; !padic_poly_is_zero(f) && (i < 10); i++)
        padic_poly_randtest(f, state, len, ctx);

    if (padic_poly_is_zero(f))
    {
        padic_poly_fit_length(f, 1);
        _padic_poly_set_length(f, 1);
        fmpz_one(f->coeffs + 0);
        f->val = f->N - 1;
    }
}

