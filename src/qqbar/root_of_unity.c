/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpq.h"
#include "qqbar.h"

int
qqbar_is_root_of_unity(slong * p, ulong * q, const qqbar_t x)
{
    ulong n;

    n = fmpz_poly_is_cyclotomic(QQBAR_POLY(x));

    if (n == 0)
        return 0;

    if (q != NULL)
        *q = n;

    if (n == 1)
    {
        if (p != NULL) *p = 0;
    }
    else if (n == 2)
    {
        if (p != NULL) *p = 1;
    }
    else if (n == 3)
    {
        if (p != NULL) *p = (qqbar_sgn_im(x) > 0) ? 1 : 2;
    }
    else if (n == 4)
    {
        if (p != NULL) *p = (qqbar_sgn_im(x) > 0) ? 1 : 3;
    }
    else
    {
        if (p != NULL)
        {
            arb_t t, u;
            acb_t z;
            fmpz_t k;
            slong prec;

            acb_init(z);
            arb_init(t);
            arb_init(u);
            fmpz_init(k);

            prec = 64;  /* more than enough */

            qqbar_get_acb(z, x, prec);
            acb_arg(t, z, prec);
            arb_const_pi(u, prec);
            arb_div(t, t, u, prec);
            arb_mul_2exp_si(t, t, -1);
            arb_mul_ui(t, t, n, prec);

            if (!arb_get_unique_fmpz(k, t))
            {
                flint_throw(FLINT_ERROR, "qqbar_is_root_of_unity: unexpected precision issue\n");
            }

            if (fmpz_sgn(k) < 0)
                fmpz_add_ui(k, k, n);

            *p = fmpz_get_si(k);

            acb_clear(z);
            arb_clear(t);
            arb_clear(u);
            fmpz_clear(k);
        }
    }

    return 1;
}

void
qqbar_root_of_unity(qqbar_t res, slong p, ulong q)
{
    fmpq_t t;
    ulong a, b;
    slong prec;

    fmpq_init(t);

    if (q == 0)
    {
        flint_throw(FLINT_ERROR, "qqbar_root_of_unity: q = 0\n");
    }

    fmpq_set_si(t, p, q);
    fmpz_fdiv_r(fmpq_numref(t), fmpq_numref(t), fmpq_denref(t));

    a = fmpz_get_ui(fmpq_numref(t));
    b = fmpz_get_ui(fmpq_denref(t));

    if (a == 0)
    {
        qqbar_one(res);
    }
    else if (a == 1 && b == 2)
    {
        qqbar_set_si(res, -1);
    }
    else if (a == 1 && b == 4)
    {
        qqbar_i(res);
    }
    else if (a == 3 && b == 4)
    {
        qqbar_i(res);
        qqbar_conj(res, res);
    }
    else
    {
        fmpz_poly_cyclotomic(QQBAR_POLY(res), b);
        fmpq_mul_2exp(t, t, 1);

        for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            arb_sin_cos_pi_fmpq(acb_imagref(QQBAR_ENCLOSURE(res)),
                                acb_realref(QQBAR_ENCLOSURE(res)),
                                t, prec);

            /* todo: this is really unnecessary... */
            if (_qqbar_validate_uniqueness(QQBAR_ENCLOSURE(res),
                    QQBAR_POLY(res), QQBAR_ENCLOSURE(res), prec * 2))
            {
                break;
            }
        }
    }

    fmpq_clear(t);
}

