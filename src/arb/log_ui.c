/*
    Copyright (C) 2013-2014 Fredrik Johansson
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
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

/* log k from log(prev): the kernels of mp_real (arb_log_arf) at most
   precisions, else log k = log prev + 2 atanh(p/q), p/q = (k - prev) /
   (k + prev), by binary splitting, which wins only at a high precision
   growing with the ratio: measured against the kernels, at prec above
   about 4 10^4 (16 / lg(q/p))^2 bits */
void
arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, slong prec)
{
    ulong p, q;
    double lg;

    if (prev < 2 || k <= prev || (k + prev) < prev)
    {
        arf_t t;
        arf_init_set_ui(t, k);
        arb_log_arf(s, t, prec);
        return;
    }

    p = k - prev;
    q = k + prev;
    lg = log2((double) q / (double) p);

    if (prec < 4e4 * (16.0 / lg) * (16.0 / lg))
    {
        arf_t t;
        arf_init_set_ui(t, k);
        arb_log_arf(s, t, prec);
    }
    else
    {
        arb_t t;

        arb_init(t);

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
