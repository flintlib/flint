/*
    Copyright (C) 2013-2014 Fredrik Johansson
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "arb.h"

/* Assumption: p/q <= 2 */
static void
atanh_bs(arb_t s, ulong p, ulong q, slong prec)
{
    fmpz_t pp, qq;

    fmpz_init_set_ui(pp, p);
    fmpz_init_set_ui(qq, q);

    arb_atan_frac_bsplit(s, pp, qq, 1, prec);

    fmpz_clear(pp);
    fmpz_clear(qq);
}

static int n_width(ulong k)
{
    int a, b;
    a = flint_clz(k);
    b = flint_ctz(k);
    return FLINT_BITS - a - b;
}

void
arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, slong prec)
{
    if (prev < 2 || prec < 600 ||
        (prec < ARB_LOG_TAB2_PREC - 64 && n_width(k) <= ARB_LOG_TAB21_BITS + 1)
        || k < prev || (k + prev) < prev ||
        (k - prev) >= 0.25 * (k + prev))
    {
        arf_t t;
        arf_init_set_ui(t, k);
        arb_log_arf(s, t, prec);
        /* no need to clear t */
    }
    else
    {
        arb_t t;
        ulong p, q;

        arb_init(t);

        p = k - prev;
        q = k + prev;

        if ((p % 2 == 0) && (q % 2 == 0))
        {
            p >>= 1;
            q >>= 1;
        }

        atanh_bs(t, p, q, prec);
        arb_mul_2exp_si(t, t, 1);
        arb_add(s, log_prev, t, prec);

        arb_clear(t);
    }
}

void
arb_log_ui(arb_t z, ulong x, slong prec)
{
    if (x == 2)
    {
        arb_const_log2(z, prec);
    }
    else if (x == 10)
    {
        arb_const_log10(z, prec);
    }
    else
    {
        arf_t t;
        arf_init(t);
        arf_set_ui(t, x);
        arb_log_arf(z, t, prec);
        arf_clear(t);
    }
}

void
arb_log_fmpz(arb_t z, const fmpz_t x, slong prec)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    arb_log_arf(z, t, prec);
    arf_clear(t);
}
