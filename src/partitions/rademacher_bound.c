/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "partitions.h"

static void
_arf_sinh(arf_t y, const arf_t x, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_arf(t, x);
    arb_sinh(t, t, prec);
    arb_get_abs_ubound_arf(y, t, prec);
    arb_clear(t);
}

/* Equation (1.8) in the paper */
void
partitions_rademacher_bound(arf_t b, const fmpz_t n, ulong N)
{
    arf_t A, B, C, t, u;
    fmpz_t n1;

    arf_init(A);
    arf_init(B);
    arf_init(C);
    arf_init(t);
    arf_init(u);
    fmpz_init(n1);

    /* bound for 44*pi^2/(225*sqrt(3)) */
    arf_set_si_2exp_si(A, 18695160, -24);

    /* bound for pi*sqrt(2)/75 */
    arf_set_si_2exp_si(B, 993857, -24);

    /* bound for pi*sqrt(2/3) */
    arf_set_si_2exp_si(C, 43035232, -24);

    /* first term: A / sqrt(N) */
    arf_sqrt_ui(t, N, MAG_BITS, ARF_RND_DOWN);
    arf_div(b, A, t, MAG_BITS, ARF_RND_UP);

    /* B * sqrt(N/(n-1)) */
    arf_set_ui(t, N);
    fmpz_sub_ui(n1, n, 1);
    arf_div_fmpz(t, t, n1, MAG_BITS, ARF_RND_UP);
    arf_sqrt(t, t, MAG_BITS, ARF_RND_UP);
    arf_mul(t, B, t, MAG_BITS, ARF_RND_UP);

    /* sinh(C*sqrt(n)/N) */
    arf_sqrt_fmpz(u, n, MAG_BITS, ARF_RND_UP);
    arf_div_ui(u, u, N, MAG_BITS, ARF_RND_UP);
    arf_mul(u, C, u, MAG_BITS, ARF_RND_UP);

    _arf_sinh(u, u, MAG_BITS);

    /* second term: B * ... * sinh... */
    arf_mul(t, t, u, MAG_BITS, ARF_RND_UP);
    arf_add(b, b, t, MAG_BITS, ARF_RND_UP);

    arf_clear(A);
    arf_clear(B);
    arf_clear(C);
    arf_clear(t);
    arf_clear(u);
    fmpz_clear(n1);
}

