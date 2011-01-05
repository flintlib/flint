/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"

#if FLINT64
#define SMALL_BERNOULLI_LIMIT 35
#else
#define SMALL_BERNOULLI_LIMIT 27
#endif

static const long bernoulli_numer_small[] = {
    1L, 1L, -1L, 1L, -1L, 5L, -691L, 7L, -3617L, 43867L, -174611L, 854513L,
    -236364091L, 8553103L,
#if FLINT64
    -23749461029L, 8615841276005L, -7709321041217L, 2577687858367L
#endif
};

static const unsigned int bernoulli_denom_small[] = {
    1, 6, 30, 42, 30, 66, 2730, 6, 510, 798, 330, 138,
    2730, 6, 870, 14322, 510, 6
};

void _fmpz_bernoulli_vec_recursive(fmpz_t den, fmpz * b, long n)
{
    fmpz_t t, c, d;
    long j, k, m, mcase;

    fmpz_init(t);
    fmpz_init(c);
    fmpz_init(d);

    fmpz_primorial(den, n + 1);

    /* Initial values */
    for (k = 0; k < FLINT_MIN(SMALL_BERNOULLI_LIMIT, n); k += 2)
    {
        fmpz_mul_si(b + k, den, bernoulli_numer_small[k/2]);
        fmpz_divexact_ui(b + k, b + k, bernoulli_denom_small[k/2]);
    }

    /* Odd values */
    if (n > 1)
        fmpz_divexact_si(b + 1, den, -2L);
    for (k = 3; k < n; k += 2)
        fmpz_zero(b + k);

    /* Ramanujan's recursive formula */
    for (m = SMALL_BERNOULLI_LIMIT + 1; m < n; m += 2)
    {
        mcase = m % 6;

        fmpz_mul_ui(b + m, den, m + 3UL);
        fmpz_divexact_ui(b + m, b + m, 3UL);

        if (mcase == 4)
        {
            fmpz_neg(b + m, b + m);
            fmpz_divexact_ui(b + m, b + m, 2UL);
        }

        /* c = t = binomial(m+3, m) */
        fmpz_set_ui(t, m + 1UL);
        fmpz_mul_ui(t, t, m + 2UL);
        fmpz_mul_ui(t, t, m + 3UL);
        fmpz_divexact_ui(t, t, 6UL);
        fmpz_set(c, t);

        for (j = 1; j <= m / 6; j++)
        {
            /* c = binomial(m+3, m-6*j); */
            fmpz_mul_ui(c, c, m + 6 - 6*j);
            fmpz_mul_ui(c, c, m + 5 - 6*j);
            fmpz_mul_ui(c, c, m + 4 - 6*j);
            fmpz_mul_ui(c, c, m + 3 - 6*j);
            fmpz_mul_ui(c, c, m + 2 - 6*j);
            fmpz_mul_ui(c, c, m + 1 - 6*j);
            fmpz_set_ui(d, 72*j);
            fmpz_mul_ui(d, d, 2*j + 1);
            fmpz_mul_ui(d, d, 3*j - 1);
            fmpz_mul_ui(d, d, 3*j + 1);
            fmpz_mul_ui(d, d, 6*j - 1);
            fmpz_mul_ui(d, d, 6*j + 1);
            fmpz_divexact(c, c, d);
            fmpz_submul(b + m, c, b + (m - 6*j));
        }

        fmpz_divexact(b + m, b + m, t);
    }

    fmpz_clear(t);
    fmpz_clear(c);
    fmpz_clear(d);
}
