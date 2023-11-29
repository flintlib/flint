/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "acb_dirichlet.h"

void
acb_dirichlet_root(acb_t z, const acb_dirichlet_roots_t t, ulong k, slong prec)
{
    ulong n = t->order;
    int swap, flip, conjugate;
    slong wp;

    if (k > n)
        k %= n;

    swap = flip = conjugate = 0;

    if (k > n / 2)
    {
        conjugate = 1;
        k = n - k;
    }

    if (n % 2 == 0 && k > n / 4)
    {
        flip = 1;
        k = n / 2 - k;
    }

    if (n % 4 == 0 && k > n / 8)
    {
        swap = 1;
        k = n / 4 - k;
    }

    wp = prec + 6 + 2 * FLINT_BIT_COUNT(t->reduced_order);

    if (k == 0)
    {
        acb_one(z);
    }
    else if (t->depth == 0)
    {
        if (t->use_pow)
        {
            acb_pow_ui(z, t->z, k, wp);
            acb_set_round(z, z, prec);
        }
        else
        {
            ulong r;
            fmpq_t t;
            fmpq_init(t);
            r = n_gcd(n, 2 * k); /* no overflow since since k <= n / 2 */
            fmpz_set_ui(fmpq_numref(t), (2 * k) / r);
            fmpz_set_ui(fmpq_denref(t), n / r);
            arb_sin_cos_pi_fmpq(acb_imagref(z), acb_realref(z), t, prec);
            fmpq_clear(t);
        }
    }
    else if (t->depth == 1)
    {
        acb_set_round(z, t->Z[0] + k, prec);
    }
    else
    {
        slong j;
        ulong r;

        r = k % t->size;
        k = k / t->size;
        acb_set(z, t->Z[0] + r);

        for (j = 1; j < t->depth && k != 0; j++)
        {
            r = k % t->size;
            k = k / t->size;
            acb_mul(z, z, t->Z[j] + r, wp);
        }

        if (k != 0)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        acb_set_round(z, z, prec);
    }

    if (swap)
        arb_swap(acb_realref(z), acb_imagref(z));
    if (flip)
        arb_neg(acb_realref(z), acb_realref(z));
    if (conjugate)
        arb_neg(acb_imagref(z), acb_imagref(z));
}

