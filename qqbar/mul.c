/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_mul(qqbar_t res, const qqbar_t x, const qqbar_t y)
{
    if (qqbar_is_zero(x) || qqbar_is_zero(y))
    {
        qqbar_zero(res);
    }
    else if (qqbar_is_one(x))
    {
        qqbar_set(res, y);
    }
    else if (qqbar_is_one(y))
    {
        qqbar_set(res, x);
    }
    else if (qqbar_is_neg_one(x))
    {
        qqbar_neg(res, y);
    }
    else if (qqbar_is_neg_one(y))
    {
        qqbar_neg(res, x);
    }
    else if (qqbar_is_rational(y))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_neg(a, QQBAR_COEFFS(y));
        fmpz_set(c, QQBAR_COEFFS(y) + 1);

        qqbar_scalar_op(res, x, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else if (qqbar_is_rational(x))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_neg(a, QQBAR_COEFFS(x));
        fmpz_set(c, QQBAR_COEFFS(x) + 1);

        qqbar_scalar_op(res, y, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else if (qqbar_equal(x, y))
    {
        /* This may detect exact square roots and other special cases. */
        qqbar_pow_ui(res, x, 2);
    }
    else if (_fmpz_vec_is_zero(QQBAR_COEFFS(x) + 1, qqbar_degree(x) - 1) &&
             _fmpz_vec_is_zero(QQBAR_COEFFS(y) + 1, qqbar_degree(y) - 1) &&
            fmpz_sgn(QQBAR_COEFFS(x)) < 0 &&
            fmpz_sgn(QQBAR_COEFFS(y)) < 0 &&
            arb_contains_zero(acb_imagref(QQBAR_ENCLOSURE(x))) && arb_is_positive(acb_realref(QQBAR_ENCLOSURE(x))) &&
            arb_contains_zero(acb_imagref(QQBAR_ENCLOSURE(y))) && arb_is_positive(acb_realref(QQBAR_ENCLOSURE(y))))
    {
        /* (p/q)^(1/d) * (r/s)^(1/e) */

        fmpq_t t, u;
        ulong d, e, f, g;

        d = qqbar_degree(x);
        e = qqbar_degree(y);
        g = n_gcd(d, e);
        f = (d / g) * e;

        fmpq_init(t);
        fmpq_init(u);

        fmpz_neg(fmpq_numref(t), QQBAR_COEFFS(x));
        fmpz_set(fmpq_denref(t), QQBAR_COEFFS(x) + d);
        fmpz_neg(fmpq_numref(u), QQBAR_COEFFS(y));
        fmpz_set(fmpq_denref(u), QQBAR_COEFFS(y) + e);

        fmpq_pow_si(t, t, e / g);
        fmpq_pow_si(u, u, d / g);
        fmpq_mul(t, t, u);

        /* todo: recycle the existing enclosures instead of computing
           a numerical f-th root from scratch */
        qqbar_fmpq_root_ui(res, t, f);

        fmpq_clear(t);
        fmpq_clear(u);
    }
    else
    {
        qqbar_binary_op(res, x, y, 2);
    }
}

void
qqbar_mul_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpq(t, y);
    qqbar_mul(res, x, t);
    qqbar_clear(t);
}

void
qqbar_mul_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpz(t, y);
    qqbar_mul(res, x, t);
    qqbar_clear(t);
}

void
qqbar_mul_ui(qqbar_t res, const qqbar_t x, ulong y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_ui(t, y);
    qqbar_mul(res, x, t);
    qqbar_clear(t);
}

void
qqbar_mul_si(qqbar_t res, const qqbar_t x, slong y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_si(t, y);
    qqbar_mul(res, x, t);
    qqbar_clear(t);
}

