/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_mat.h"

void
_fmpq_set_cfrac_basecase(fmpz_t p, fmpz_t t, fmpz_t q, fmpz_t u,
    const fmpz * c, slong n)
{
    slong i;

    fmpz_set(p, c);
    fmpz_one(q);
    fmpz_one(t);
    fmpz_zero(u);

    for (i = 1; i < n; i++)
    {
        fmpz_addmul(t, c + i, p);
        fmpz_addmul(u, c + i, q);
        fmpz_swap(t, p);
        fmpz_swap(u, q);
    }
}

void
_fmpq_set_cfrac_divconquer(fmpz_mat_t P, const fmpz * c, slong n)
{
    if (n < 32)
    {
        _fmpq_set_cfrac_basecase(
            fmpz_mat_entry(P, 0, 0), fmpz_mat_entry(P, 0, 1),
            fmpz_mat_entry(P, 1, 0), fmpz_mat_entry(P, 1, 1), c, n);
    }
    else
    {
        fmpz_mat_t L, R;
        slong m = n / 2;

        fmpz_mat_init(L, 2, 2);
        fmpz_mat_init(R, 2, 2);

        _fmpq_set_cfrac_divconquer(L, c, m);
        _fmpq_set_cfrac_divconquer(R, c + m, n - m);
        fmpz_mat_mul_classical(P, L, R);  /* Should be Strassen */

        fmpz_mat_clear(L);
        fmpz_mat_clear(R);
    }
}

void
fmpq_set_cfrac(fmpq_t x, const fmpz * c, slong n)
{
    if (n <= 64)
    {
        fmpz_t t, u;
        fmpz_init(t);
        fmpz_init(u);
        _fmpq_set_cfrac_basecase(fmpq_numref(x), t, fmpq_denref(x), u, c, n);
        fmpz_clear(t);
        fmpz_clear(u);
    }
    else
    {
        fmpz_mat_t P;
        fmpz_mat_init(P, 2, 2);
        _fmpq_set_cfrac_divconquer(P, c, n);
        fmpz_set(fmpq_numref(x), fmpz_mat_entry(P, 0, 0));
        fmpz_set(fmpq_denref(x), fmpz_mat_entry(P, 1, 0));
        fmpz_mat_clear(P);
    }
}
