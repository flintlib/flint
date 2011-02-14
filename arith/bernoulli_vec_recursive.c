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
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"


void _fmpz_bernoulli_vec_recursive(fmpz_t den, fmpz * b, long n)
{
    fmpz_t t, c, d;
    long j, k, m, mcase;
    int prodsize;

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

        /* All factors are strictly smaller than m + 4; choose prodsize such
           that (m + 4)^prodsize fits in a signed long. */
        {
#if FLINT64
            if      (m < 1444L)       prodsize = 6;
            else if (m < 2097148L)    prodsize = 3;
            else if (m < 3037000495L) prodsize = 2;  /* not very likely... */
            else abort();
#else
            if      (m < 32L)    prodsize = 6;
            else if (m < 1286L)  prodsize = 3;
            else if (m < 46336L) prodsize = 2;
            else abort();
#endif
        }

        /* c = t = binomial(m+3, m) */
        fmpz_set_ui(t, m + 1UL);
        fmpz_mul_ui(t, t, m + 2UL);
        fmpz_mul_ui(t, t, m + 3UL);
        fmpz_divexact_ui(t, t, 6UL);
        fmpz_set(c, t);

        for (j = 6; j <= m; j += 6)
        {
            long r = m - j;

            /* c = binomial(m+3, m-j); */
            switch (prodsize)
            {
                case 2:
                fmpz_mul_ui(c, c, (r+6)*(r+5));
                fmpz_mul_ui(c, c, (r+4)*(r+3));
                fmpz_mul_ui(c, c, (r+2)*(r+1));
                fmpz_set_ui(d,    (j+0)*(j+3));
                fmpz_mul_ui(d, d, (j-2)*(j+2));
                fmpz_mul_ui(d, d, (j-1)*(j+1));
                fmpz_divexact(c, c, d);
                break;

                case 3:
                fmpz_mul_ui(c, c, (r+6)*(r+5)*(r+4));
                fmpz_mul_ui(c, c, (r+3)*(r+2)*(r+1));
                fmpz_set_ui(d,    (j+0)*(j+3)*(j-2));
                fmpz_mul_ui(d, d, (j+2)*(j-1)*(j+1));
                fmpz_divexact(c, c, d);
                break;

                case 6:
                fmpz_mul_ui(c, c,      (r+6)*(r+5)*(r+4)*(r+3)*(r+2)*(r+1));
                fmpz_divexact_ui(c, c, (j+0)*(j+3)*(j-2)*(j+2)*(j-1)*(j+1));
                break;
            }

            fmpz_submul(b + m, c, b + (m - j));
        }
        fmpz_divexact(b + m, b + m, t);
    }

    fmpz_clear(t);
    fmpz_clear(c);
    fmpz_clear(d);
}
