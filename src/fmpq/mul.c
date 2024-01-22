/*
    Copyright (C) 2011, 2020 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpq.h"

static ulong _fmpz_gcd_ui(const fmpz_t g, ulong h)
{
    if (!COEFF_IS_MPZ(*g))
        return n_gcd(FLINT_ABS(*g), h);
    else
        return n_gcd(flint_mpz_fdiv_ui(COEFF_TO_PTR(*g), h), h);
}

void
_fmpq_mul_small(fmpz_t rnum, fmpz_t rden, slong op1num, ulong op1den, slong op2num, ulong op2den)
{
    mp_limb_t hi, lo, denhi, denlo;
    int neg;

    if (op1num == 0 || op2num == 0)
    {
        fmpz_zero(rnum);
        fmpz_one(rden);
        return;
    }

    neg = 0;
    if (op1num < 0)
    {
        op1num = -op1num;
        neg = 1;
    }

    if (op2num < 0)
    {
        op2num = -op2num;
        neg = !neg;
    }

    if (op1den == op2den)
    {
        umul_ppmm(hi, lo, op1num, op2num);
        umul_ppmm(denhi, denlo, op1den, op2den);
    }
    else if (op1den == 1)
    {
        ulong t, x;
        t = n_gcd(op1num, op2den);
        x = op1num / t;
        t = op2den / t;
        umul_ppmm(hi, lo, x, op2num);
        umul_ppmm(denhi, denlo, op1den, t);
    }
    else if (op2den == 1)
    {
        ulong t, x;
        t = n_gcd(op2num, op1den);
        x = op2num / t;
        t = op1den / t;
        umul_ppmm(hi, lo, x, op1num);
        umul_ppmm(denhi, denlo, op2den, t);
    }
    else
    {
        ulong t, u, x, y;
        t = n_gcd(op1num, op2den);
        u = n_gcd(op1den, op2num);
        x = op1num / t;
        y = op2num / u;
        umul_ppmm(hi, lo, x, y);
        x = op1den / u;
        y = op2den / t;
        umul_ppmm(denhi, denlo, x, y);
    }

    if (neg)
        fmpz_neg_uiui(rnum, hi, lo);
    else
        fmpz_set_uiui(rnum, hi, lo);
    fmpz_set_uiui(rden, denhi, denlo);
}

void
_fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den,
            const fmpz_t op2num, const fmpz_t op2den)
{
    if (!COEFF_IS_MPZ(*op1num) && !COEFF_IS_MPZ(*op1den) && !COEFF_IS_MPZ(*op2num) && !COEFF_IS_MPZ(*op2den))
    {
        _fmpq_mul_small(rnum, rden, *op1num, *op1den, *op2num, *op2den);
        return;
    }

    /* Common special cases: squaring, same denominator (e.g. both integers) */
    if (((op1num == op2num) && (op1den == op2den)) ||
         fmpz_equal(op1den, op2den))
    {
        fmpz_mul(rnum, op1num, op2num);
        fmpz_mul(rden, op1den, op2den);
    }
    /* Exactly one argument is an integer */
    else if (fmpz_is_one(op1den))
    {
        fmpz_t t, x;
        fmpz_init(t);
        fmpz_init(x);

        fmpz_gcd(t, op1num, op2den);
        fmpz_divexact(x, op1num, t);
        fmpz_mul(rnum, x, op2num);
        fmpz_divexact(t, op2den, t);
        fmpz_mul(rden, op1den, t);

        fmpz_clear(t);
        fmpz_clear(x);
    }
    else if (fmpz_is_one(op2den))
    {
        fmpz_t t, x;
        fmpz_init(t);
        fmpz_init(x);

        fmpz_gcd(t, op2num, op1den);
        fmpz_divexact(x, op2num, t);
        fmpz_mul(rnum, x, op1num);
        fmpz_divexact(t, op1den, t);
        fmpz_mul(rden, op2den, t);

        fmpz_clear(t);
        fmpz_clear(x);
    }
    else
    {
        fmpz_t t, u, x, y;

        fmpz_init(t);
        fmpz_init(u);
        fmpz_init(x);
        fmpz_init(y);

        fmpz_gcd(t, op1num, op2den);
        fmpz_gcd(u, op1den, op2num);

        fmpz_divexact(x, op1num, t);
        fmpz_divexact(y, op2num, u);

        fmpz_mul(rnum, x, y);

        fmpz_divexact(x, op1den, u);
        fmpz_divexact(y, op2den, t);

        fmpz_mul(rden, x, y);

        fmpz_clear(t);
        fmpz_clear(u);
        fmpz_clear(x);
        fmpz_clear(y);
    }
}


void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    _fmpq_mul(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1),
              fmpq_numref(op2), fmpq_denref(op2));
}

void fmpq_mul_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
{
    fmpz_t y;
    *y = 1;
    _fmpq_mul(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op), fmpq_denref(op), x, y);
}

void
_fmpq_mul_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            ulong r)
{
    if (r == 0 || fmpz_is_zero(p))
    {
        fmpz_zero(rnum);
        fmpz_one(rden);
    }
    else if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && r <= COEFF_MAX)
    {
        _fmpq_mul_small(rnum, rden, *p, *q, r, 1);
    }
    else if (r == 1)
    {
        fmpz_set(rnum, p);
        fmpz_set(rden, q);
    }
    else
    {
        ulong g = _fmpz_gcd_ui(q, r);

        if (g == 1)
        {
            fmpz_set(rden, q);
            fmpz_mul_ui(rnum, p, r);
        }
        else
        {
            fmpz_mul_ui(rnum, p, r / g);
            fmpz_divexact_ui(rden, q, g);
        }
    }
}

void fmpq_mul_ui(fmpq_t res, const fmpq_t op1, ulong c)
{
    _fmpq_mul_ui(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}

void
_fmpq_mul_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            slong r)
{
    if (r == 0 || fmpz_is_zero(p))
    {
        fmpz_zero(rnum);
        fmpz_one(rden);
    }
    else if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && r >= COEFF_MIN && r <= COEFF_MAX)
    {
        _fmpq_mul_small(rnum, rden, *p, *q, r, 1);
    }
    else if (r == 1)
    {
        fmpz_set(rnum, p);
        fmpz_set(rden, q);
    }
    else if (r == -WORD(1))
    {
        fmpz_neg(rnum, p);
        fmpz_set(rden, q);
    }
    else
    {
        ulong a, g;

        a = FLINT_ABS(r);
        g = _fmpz_gcd_ui(q, a);

        if (g == 1)
        {
            fmpz_set(rden, q);
            fmpz_mul_si(rnum, p, r);
        }
        else
        {
            /* not using fmpz_mul_si(...)  because of the special case g = -WORD_MIN */
            fmpz_mul_ui(rnum, p, a / g);
            if (r < 0)
                fmpz_inplace_neg(rnum);
            fmpz_divexact_ui(rden, q, g);
        }
    }
}

void fmpq_mul_si(fmpq_t res, const fmpq_t op1, slong c)
{
    _fmpq_mul_si(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}
