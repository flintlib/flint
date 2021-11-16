/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arith.h"

static void
_fmpz_ui_pow_ui(fmpz_t x, ulong b, ulong e)
{
    if (e <= 1)
    {
        fmpz_set_ui(x, e == 0 ? 1 : b);
    }
    else if (e == 2)
    {
        mp_limb_t t[2];
        umul_ppmm(t[1], t[0], b, b);
        fmpz_set_uiui(x, t[1], t[0]);
    }
    else if (b <= 1)
    {
        fmpz_set_ui(x, b);
    }
    else
    {
        ulong bits = FLINT_BIT_COUNT(b);

        if (e * bits <= FLINT_BITS)
        {
            fmpz_set_ui(x, n_pow(b, e));
        }
        else
        {
            __mpz_struct * z = _fmpz_promote(x);
            flint_mpz_set_ui(z, b);
            flint_mpz_pow_ui(z, z, e);
            _fmpz_demote_val(x);
        }
    }
}

void
arith_bell_number_dobinski(fmpz_t res, ulong n)
{
    fmpz_t P, Q, t;
    fmpz * pow;
    slong N, k, kodd, shift;

    if (n <= 1)
    {
        fmpz_one(res);
        return;
    }

    N = n * (1 + 1.2 / log(n)) + 2;

    fmpz_init(P);
    fmpz_init(Q);
    fmpz_init(t);

    pow = _fmpz_vec_init((N + 2) / 4);

    fmpz_one(P);
    fmpz_mul_2exp(P, P, n);
    fmpz_add_ui(P, P, 2);
    fmpz_set_ui(Q, 5);

    for (k = 3; k <= N; k++)
    {
        fmpz_mul_ui(P, P, k);

        if (k % 2 == 1)
        {
            if (2 * k <= N)
            {
                _fmpz_ui_pow_ui(pow + k / 2, k, n);
                fmpz_add(P, P, pow + k / 2);
            }
            else
            {
                _fmpz_ui_pow_ui(t, k, n);
                fmpz_add(P, P, t);
            }
        }
        else
        {
            kodd = k / 2;
            shift = n;
            while (kodd % 2 == 0)
            {
                kodd = kodd / 2;
                shift += n;
            }

            if (kodd == 1)
                fmpz_one_2exp(t, shift);
            else
                fmpz_mul_2exp(t, pow + kodd / 2, shift);

            fmpz_add(P, P, t);
        }

        fmpz_mul_ui(Q, Q, k);
        fmpz_add_ui(Q, Q, 1);
    }

    fmpz_cdiv_q(res, P, Q);

    _fmpz_vec_clear(pow, (N + 2) / 4);

    fmpz_clear(P);
    fmpz_clear(Q);
    fmpz_clear(t);
}
