/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_vec.h"

/*
    num[i nstride] = num(a_i) (den / den(a_i)) for den a positive common
    multiple of the denominators. The factor den / den(a_i) is 1 when the
    denominator is den (the common case of a common denominator) and den
    when it is 1, and it is computed once for a run of equal denominators.
*/
/* equality of canonical fmpz: small values are equal iff their words are */
#define FMPZ_EQUAL(a, b) (*(a) == *(b) || (COEFF_IS_MPZ(*(a)) && COEFF_IS_MPZ(*(b)) && fmpz_equal(a, b)))

void
_fmpq_vec_get_fmpz_vec_scaled(fmpz * num, slong nstride,
                const fmpq * a, slong astride, slong len, const fmpz_t den)
{
    fmpz_t t;
    const fmpz * lastd = NULL;
    slong i;

    fmpz_init(t);

    for (i = 0; i < len; i++)
    {
        const fmpz * n = fmpq_numref(a + i * astride);
        const fmpz * d = fmpq_denref(a + i * astride);
        fmpz * res = num + i * nstride;

        if (fmpz_is_one(d))
        {
            fmpz_mul(res, n, den);
        }
        else if (FMPZ_EQUAL(d, den))
        {
            fmpz_set(res, n);
        }
        else
        {
            /* the quotient goes to a temporary rather than to res, which
               would otherwise switch between a small value and an mpz */
            if (lastd == NULL || !FMPZ_EQUAL(d, lastd))
            {
                lastd = d;
                fmpz_divexact(t, den, d);
            }
            fmpz_mul(res, n, t);
        }
    }

    fmpz_clear(t);
}
