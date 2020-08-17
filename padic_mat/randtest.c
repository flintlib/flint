/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_mat.h"

void padic_mat_randtest(padic_mat_t mat, flint_rand_t state, const padic_ctx_t ctx)
{
    if (!padic_mat_is_empty(mat))
    {
        const slong N = padic_mat_prec(mat);

        slong i, j, min, max;
        fmpz_t pow;

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
        else  /* ctx->N == 0 */
        {
            min = -10;
            max = 0;
        }

        padic_mat_val(mat) = n_randint(state, max - min) + min;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, N - padic_mat_val(mat));

        for (i = 0; i < padic_mat(mat)->r; i++)
            for (j = 0; j < padic_mat(mat)->c; j++)
                fmpz_randm(padic_mat_entry(mat, i, j), state, pow);

        fmpz_clear(pow);

        _padic_mat_canonicalise(mat, ctx);
    }
}

