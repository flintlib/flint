/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_powsum_term(acb_ptr res, arb_t log_prev, ulong * prev,
    const acb_t s, ulong k, int integer, int critical_line, slong len, slong prec)
{
    slong i;

    if (integer)
    {
        arb_neg(acb_realref(res), acb_realref(s));
        arb_set_ui(acb_imagref(res), k);
        arb_pow(acb_realref(res), acb_imagref(res), acb_realref(res), prec);
        arb_zero(acb_imagref(res));

        if (len != 1)
        {
            arb_log_ui_from_prev(log_prev, k, log_prev, *prev, prec);
            *prev = k;
        }
    }
    else
    {
        arb_t w;
        arb_init(w);

        arb_log_ui_from_prev(log_prev, k, log_prev, *prev, prec);
        *prev = k;
        arb_mul(w, log_prev, acb_imagref(s), prec);
        arb_sin_cos(acb_imagref(res), acb_realref(res), w, prec);
        arb_neg(acb_imagref(res), acb_imagref(res));

        if (critical_line)
        {
            arb_rsqrt_ui(w, k, prec);
            acb_mul_arb(res, res, w, prec);
        }
        else
        {
            arb_mul(w, acb_realref(s), log_prev, prec);
            arb_neg(w, w);
            arb_exp(w, w, prec);
            acb_mul_arb(res, res, w, prec);
        }

        arb_clear(w);
    }

    if (len > 1)
    {
        arb_neg(log_prev, log_prev);

        for (i = 1; i < len; i++)
        {
            acb_mul_arb(res + i, res + i - 1, log_prev, prec);
            acb_div_ui(res + i, res + i, i, prec);
        }

        arb_neg(log_prev, log_prev);
    }
}

