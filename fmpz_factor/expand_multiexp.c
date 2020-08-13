/*
    Copyright (C) 2009, 2010 Paul Zimmermann
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_factor_eval_multiexp(fmpz_t res, const fmpz * p, const ulong * e, slong len)
{
    slong i, j;
    ulong mask, emax;
    fmpz * q;
    fmpz_t tmp;

    if (len <= 1)
    {
        if (len < 1)
            fmpz_one(res);
        else
            fmpz_pow_ui(res, p, e[0]);
        return;
    }

    q = flint_malloc(sizeof(fmpz) * len);

    emax = e[0];
    for (i = 1; i < len; i++)
        emax = FLINT_MAX(emax, e[i]);

    for (mask = 1; emax >= (mask << 1); mask <<= 1);

    fmpz_init(tmp);
    fmpz_one(res);

    while (mask)
    {
        for (i = 0, j = 0; i < len; i++)
            if (e[i] & mask)
                q[j++] = p[i];

        _fmpz_vec_prod(tmp, q, j);
        fmpz_mul(res, res, res);
        fmpz_mul(res, res, tmp);
        mask >>= 1;
    }

    fmpz_clear(tmp);
    flint_free(q);
}

void
fmpz_factor_expand_multiexp(fmpz_t n, const fmpz_factor_t factor)
{
    if (factor->num != 0 && factor->p[0] == 2)
    {
        _fmpz_factor_eval_multiexp(n, factor->p + 1, factor->exp + 1, factor->num - 1);
        fmpz_mul_2exp(n, n, factor->exp[0]);
    }
    else
        _fmpz_factor_eval_multiexp(n, factor->p, factor->exp, factor->num);

    fmpz_mul_si(n, n, factor->sign);
}
