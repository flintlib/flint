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

#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "math.h"
#include "arith.h"
#include "ulong_extras.h"


static __inline__ void
update_count(int * count, long p, ulong k)
{
    if (p < 0)
        p = -p;

    p = p % (12 * k);

    if (p >= 6 * k)
        p = 12 * k - p;

    if (p > 3 * k)
        count[6*k - p]--;
    else
        count[p]++;
}

void
dedekind_cosine_sum_mpfr_fast(mpfr_t sum, ulong k, ulong n)
{
    mpfr_t t, u;
    ulong h, q;
    long i, p;
    int * count;
    double s;

    if (k > 1000000)
    {
        /* TODO: just don't use doubles for large k ... */
        printf("Exception: dedekind_cosine_sum: large k not yet supported\n");
        abort();
    }

    if (k <= 2)
    {
        if (k == 0)
            mpfr_set_ui(sum, 0UL, GMP_RNDN);
        else if (k == 1)
            mpfr_set_ui(sum, 1UL, GMP_RNDN);
        else if (n % 2 == 1)
            mpfr_set_si(sum, -1L, GMP_RNDN);
        else
            mpfr_set_si(sum, 1L, GMP_RNDN);
        return;
    }

    count = calloc((3*k + 1), sizeof(int));
    q = 6 * k;

    for (h = 0; h < (k + 1) / 2; h++)
    {
        if (n_gcd(k, h) == 1UL)
        {
            s = dedekind_sum_coprime_d(h, k);
            p = floor((s * q) + 0.5);

            /* XXX: 32-bit overflow */
            update_count(count, p - 12*h*n, k);
            update_count(count, -p - 12*(k-h)*n, k);
        }
    }

    mpfr_set_si(sum, count[0], MPFR_RNDN);
    mpfr_init2(t, mpfr_get_prec(sum) + 5);
    mpfr_init2(u, mpfr_get_prec(sum) + 5);

    mpfr_const_pi(u, MPFR_RNDN);
    mpfr_div_ui(u, u, 6 * k, MPFR_RNDN);

    for (i = 1; i < 3 * k; i++)
    {
        if (count[i] != 0)
        {
            mpfr_mul_si(t, u, i, MPFR_RNDN);
            mpfr_cos(t, t, MPFR_RNDN);
            if (count[i] != 1)
                mpfr_mul_si(t, t, count[i], MPFR_RNDN);
            mpfr_add(sum, sum, t, MPFR_RNDN);
        }
    }

    mpfr_clear(t);
    mpfr_clear(u);
    free(count);
}
