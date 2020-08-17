/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "padic.h"
#include "ulong_extras.h"

void _padic_log_satoh(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)
{
    if (N < 16)
    {
        _padic_log_rectangular(z, y, v, p, N);
    }
    else
    {
        const slong k = n_sqrt(N);

        fmpz_t t, pk, pNk;

        fmpz_init(t);
        fmpz_init(pk);
        fmpz_init(pNk);
        fmpz_pow_ui(pk, p, k);
        fmpz_pow_ui(pNk, p, N + k);

        fmpz_sub_ui(t, y, 1);
        fmpz_neg(t, t);
        fmpz_powm(t, t, pk, pNk);
        fmpz_sub_ui(t, t, 1);
        fmpz_neg(t, t);

        /* TODO: Improve suggested valuation */
        _padic_log_rectangular(z, t, k + 1, p, N + k);

        fmpz_divexact(z, z, pk);

        fmpz_clear(t);
        fmpz_clear(pk);
        fmpz_clear(pNk);
    }
}

int padic_log_satoh(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    const fmpz *p = ctx->p;
    const slong N  = padic_prec(rop);

    if (padic_val(op) < 0)
    {
        return 0;
    }
    else
    {
        fmpz_t y;
        int ans;

        fmpz_init(y);

        padic_get_fmpz(y, op, ctx);
        fmpz_sub_ui(y, y, 1);
        fmpz_neg(y, y);

        if (fmpz_is_zero(y))
        {
            padic_zero(rop);
            ans = 1;
        }
        else
        {
            fmpz_t t;
            slong v;

            fmpz_init(t);
            v = fmpz_remove(t, y, ctx->p);
            fmpz_clear(t);

            if (v >= 2 || (!fmpz_equal_ui(p, 2) && v >= 1))
            {
                if (v >= N)
                {
                    padic_zero(rop);
                }
                else
                {
                    _padic_log_satoh(padic_unit(rop), y, v, p, N);
                    padic_val(rop) = 0;
                    padic_reduce(rop, ctx);
                }
                ans = 1;
            }
            else
            {
                ans = 0;
            }
        }

        fmpz_clear(y);
        return ans;
    }
}

