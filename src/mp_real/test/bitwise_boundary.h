/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MP_REAL_TEST_BITWISE_BOUNDARY_H
#define MP_REAL_TEST_BITWISE_BOUNDARY_H

#include "arb.h"

#define BOUNDARY_LOG1P 0
#define BOUNDARY_ATAN 1
#define BOUNDARY_TRIG 2

/* (x, n), in [0, 1), at the decision boundaries of the bitwise
   reductions up to index r: with signed digits d_i (in sparse runs, from
   a random start) and, half the time, one signed-digit threshold
   e = +-0.75 2^-j of the vector reductions, the argument with

       1 + x = prod (1 + d_i 2^-i) exp(e)         (log1p),
       atan(x) = sum d_i atan(2^-i) + e           (atan),
       x / 2 = sum d_i atan(2^-i) + e             (sin, cos, tan),

   floored to n limbs and moved by a few ulps; with only positive digits
   (every other call) these are the boundaries of the greedy reductions
   as well */
static void
_test_bitwise_boundary_arg(nn_ptr x, slong n, slong r, int kind,
    flint_rand_t state)
{
    arb_t a, u;
    fmpz_t f;
    slong i, prec = FLINT_BITS * n + 128;
    int sgn = n_randint(state, 2);

    arb_init(a);
    arb_init(u);
    fmpz_init(f);

    if (kind == BOUNDARY_LOG1P)
        arb_one(a);
    else
        arb_zero(a);

    for (i = 1 + (n_randint(state, 3) == 0 ? n_randint(state, r) : 0);
            i <= r + 1; i += 1 + n_randint(state, 4))
    {
        arb_one(u);
        arb_mul_2exp_si(u, u, -i);
        if (sgn && n_randint(state, 2))
            arb_neg(u, u);
        if (kind == BOUNDARY_LOG1P)
        {
            arb_add_ui(u, u, 1, prec);
            arb_mul(a, a, u, prec);
        }
        else
        {
            arb_atan(u, u, prec);
            arb_add(a, a, u, prec);
        }
    }

    if (n_randint(state, 2))
    {
        arb_set_ui(u, 3);
        arb_mul_2exp_si(u, u, -(slong) (1 + n_randint(state, r)) - 2);
        if (n_randint(state, 2))
            arb_neg(u, u);
        if (kind == BOUNDARY_LOG1P)
        {
            arb_exp(u, u, prec);
            arb_mul(a, a, u, prec);
        }
        else
            arb_add(a, a, u, prec);
    }

    if (kind == BOUNDARY_LOG1P)
        arb_sub_ui(a, a, 1, prec);
    else if (kind == BOUNDARY_ATAN)
        arb_tan(a, a, prec);
    else
        arb_mul_2exp_si(a, a, 1);

    flint_mpn_zero(x, n);
    if (arb_is_positive(a) && arf_cmpabs_2exp_si(arb_midref(a), 0) < 0)
    {
        arb_mul_2exp_si(a, a, FLINT_BITS * n);
        arf_get_fmpz(f, arb_midref(a), ARF_RND_FLOOR);
        if (fmpz_sgn(f) > 0 && fmpz_bits(f) <= FLINT_BITS * n)
            fmpz_get_ui_array(x, n, f);
    }

    /* a few ulps either way, staying in [0, 1) */
    if (n_randint(state, 2))
    {
        ulong c = n_randint(state, 3);
        if (mpn_add_1(x, x, n, c))
            mpn_sub_1(x, x, n, c);
    }
    else
    {
        ulong c = n_randint(state, 3);
        if (mpn_sub_1(x, x, n, c))
            mpn_add_1(x, x, n, c);
    }

    arb_clear(a);
    arb_clear(u);
    fmpz_clear(f);
}

/* the measured error of (y, yn) against the reference e, in ulps of
   2^(-FLINT_BITS n) */
static double
_test_bitwise_ulps(nn_srcptr y, slong yn, slong n, const arb_t e)
{
    arb_t va, dd;
    fmpz_t f;
    mag_t mb;
    double u;

    arb_init(va);
    arb_init(dd);
    fmpz_init(f);
    mag_init(mb);
    fmpz_set_ui_array(f, y, yn);
    arb_set_fmpz(va, f);
    arb_mul_2exp_si(va, va, -FLINT_BITS * n);
    arb_sub(dd, va, e, ARF_PREC_EXACT);
    arb_get_mag(mb, dd);
    mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
    u = mag_get_d(mb);
    arb_clear(va);
    arb_clear(dd);
    fmpz_clear(f);
    mag_clear(mb);
    return u;
}

/* x as an arb */
static void
_test_bitwise_arb(arb_t xa, nn_srcptr x, slong n)
{
    fmpz_t f;
    fmpz_init(f);
    fmpz_set_ui_array(f, x, n);
    arb_set_fmpz(xa, f);
    arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
    fmpz_clear(f);
}

/* n and r for a boundary test iteration */
static void
_test_bitwise_boundary_nr(slong * n, int * r, slong iter, flint_rand_t state)
{
    *n = 4 + n_randint(state, (iter % 4 == 0) ? 160 : 40);
    *r = (iter % 3 == 0) ? 0 : 32 + (int) n_randint(state, 737);
    *r = FLINT_MIN(*r, FLINT_BITS * *n - 16);
}

#endif
