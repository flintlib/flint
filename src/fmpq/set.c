/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "fmpq.h"

void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, slong p, ulong q)
{
    if (q == 1 || p == 0)
    {
        fmpz_set_si(rnum, p);
        fmpz_one(rden);
    }
    else
    {
        ulong r = n_gcd(p < 0 ? (-(ulong) p) : (ulong) p, q);

        if (p < 0)
        {
            fmpz_set_ui(rnum, (-(ulong) p) / r);
            fmpz_neg(rnum, rnum);
        }
        else
            fmpz_set_si(rnum, p / r);

        fmpz_set_ui(rden, q / r);
    }
}

void fmpq_set_si(fmpq_t res, slong p, ulong q)
{
    _fmpq_set_si(fmpq_numref(res), fmpq_denref(res), p, q);
}

void _fmpq_set_ui(fmpz_t rnum, fmpz_t rden, ulong p, ulong q)
{
    if (q == 1 || p == 0)
    {
        fmpz_set_ui(rnum, p);
        fmpz_one(rden);
    }
    else
    {
        ulong r = n_gcd(p, q);

        fmpz_set_ui(rnum, p / r);
        fmpz_set_ui(rden, q / r);
    }
}

void fmpq_set_ui(fmpq_t res, ulong p, ulong q)
{
    _fmpq_set_ui(fmpq_numref(res), fmpq_denref(res), p, q);
}

void
fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q)
{
    if (fmpz_is_zero(p))
    {
        fmpq_zero(res);
    }
    else if (fmpz_is_pm1(q) || fmpz_is_pm1(p))
    {
        if (fmpz_sgn(q) < 0)
        {
            fmpz_neg(fmpq_numref(res), p);
            fmpz_neg(fmpq_denref(res), q);
        }
        else
        {
            fmpz_set(fmpq_numref(res), p);
            fmpz_set(fmpq_denref(res), q);
        }
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_gcd(t, p, q);

        if (fmpz_is_one(t))
        {
            fmpz_set(fmpq_numref(res), p);
            fmpz_set(fmpq_denref(res), q);
        }
        else
        {
            fmpz_divexact(fmpq_numref(res), p, t);
            fmpz_divexact(fmpq_denref(res), q, t);
        }

        if (fmpz_sgn(fmpq_denref(res)) < 0)
        {
            fmpz_neg(fmpq_numref(res), fmpq_numref(res));
            fmpz_neg(fmpq_denref(res), fmpq_denref(res));
        }

        fmpz_clear(t);
    }
}

void fmpq_set_mpq(fmpq_t dest, const mpq_t src)
{
    fmpz_set_mpz(fmpq_numref(dest), mpq_numref(src));
    fmpz_set_mpz(fmpq_denref(dest), mpq_denref(src));
}
