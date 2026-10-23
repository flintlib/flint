/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_dft.h"

/* Thin wrappers around gr_dft (see gr_dft_acb in gr_dft/acb.c), which
   use fixed-point arithmetic with rigorous error bounds by default
   and fall back to ball arithmetic when fixed point does not
   apply. */

void
acb_dft(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    if (len <= 0)
        return;

    gr_dft_acb(w, v, len, prec);
}

void
acb_dft_inverse(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    if (len <= 0)
        return;

    gr_dft_acb_inverse(w, v, len, prec);
}

void
acb_dft_precomp_init(acb_dft_pre_t pre, slong len, slong prec)
{
    if (len <= 0)
    {
        /* empty plan: transforms of length 0 are no-ops, and clearing
           is safe */
        pre->n = 0;
        pre->which = 0;
        return;
    }

    if (gr_dft_acb_precomp_init(pre, len, prec) != GR_SUCCESS)
        flint_throw(FLINT_ERROR,
                "acb_dft_precomp_init: len = %wd, prec = %wd\n", len, prec);
}

void
acb_dft_precomp_clear(acb_dft_pre_t pre)
{
    gr_dft_acb_precomp_clear(pre);
}

void
acb_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)
{
    if (pre->n <= 0)
        return;

    gr_dft_acb_precomp(w, v, pre, prec);
}

void
acb_dft_inverse_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre,
        slong prec)
{
    if (pre->n <= 0)
        return;

    gr_dft_acb_inverse_precomp(w, v, pre, prec);
}
