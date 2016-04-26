/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)
{
    if (fmpz_cmp_ui(k, UWORD(2)) <= 0 || fmpz_is_zero(h) || fmpz_equal(h, k))
    {
        fmpq_zero(s);
    }
    else if (fmpz_sgn(h) < 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_neg(t, h);
        fmpq_dedekind_sum(s, t, k);
        fmpq_neg(s, s);
        fmpz_clear(t);
    }
    else
    {
        fmpz_t t, u, q;

        fmpz_init(t);
        fmpz_init(u);
        fmpz_init(q);

        fmpz_gcd(q, h, k);
        fmpz_divexact(t, h, q);
        fmpz_divexact(u, k, q);

        if (fmpz_cmp(t, u) > 0)
        {
            fmpq_t r;
            fmpq_init(r);

            /* r = (1 + h(h-3k) + k^2) / (12hk) */
            fmpz_mul_ui(fmpq_numref(r), u, UWORD(3));
            fmpz_sub(fmpq_numref(r), t, fmpq_numref(r));
            fmpz_mul(fmpq_numref(r), fmpq_numref(r), t);
            fmpz_addmul(fmpq_numref(r), u, u);
            fmpz_add_ui(fmpq_numref(r), fmpq_numref(r), UWORD(1));
            fmpz_mul(fmpq_denref(r), t, u);
            fmpz_mul_ui(fmpq_denref(r), fmpq_denref(r), UWORD(12));
            fmpq_canonicalise(r);
            fmpq_dedekind_sum_coprime(s, u, t);
            fmpq_sub(s, r, s);

            fmpq_clear(r);
        }
        else
        {
            fmpq_dedekind_sum_coprime(s, t, u);
        }

        fmpz_clear(t);
        fmpz_clear(u);
        fmpz_clear(q);
    }
}
