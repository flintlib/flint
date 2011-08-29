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
#include <string.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "math.h"
#include "arith.h"
#include "ulong_extras.h"


static __inline__ void
update_count(int * count, double p, double k)
{
    if (p < 0)
        p = -p;

    p = fmod(p, 12.0 * k);

    if (p >= 6.0 * k)
        p = 12.0 * k - p;

    if (p > 3.0 * k)
        count[(int)(6.0*k - p)]--;
    else
        count[(int)(p)]++;
}

#define TAB_SIZE 128

void
dedekind_cosine_sum_mpfr_fast(mpfr_t sum, ulong k, ulong n)
{
    mpfr_t t, u, cosa, sina, cc, ss, cs, sc;
    ulong h;
    long i, prev;
    int * count;
    mpfr_t cos_tab[TAB_SIZE];
    mpfr_t sin_tab[TAB_SIZE];
    char have_diff[TAB_SIZE];
    double s, q;
    mpfr_prec_t wp;

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
    q = fmod(12.0 * ((double) n), 12.0 * k);

    for (h = 0; h < (k + 1) / 2; h++)
    {
        if (n_gcd(k, h) == 1UL)
        {
            s = dedekind_sum_coprime_d(h, k);
            s = floor(s * (6.0 * k) + 0.5);
            update_count(count, s - q*h, k);
            update_count(count, -s - q*(k-h), k);
        }
    }

    mpfr_set_si(sum, count[0], MPFR_RNDN);

    wp = mpfr_get_prec(sum) + 5;

    mpfr_init2(t, wp);
    mpfr_init2(u, wp);
    mpfr_init2(cosa, wp);
    mpfr_init2(sina, wp);
    mpfr_init2(cc, wp);
    mpfr_init2(ss, wp);
    mpfr_init2(cs, wp);
    mpfr_init2(sc, wp);

    mpfr_const_pi(u, MPFR_RNDN);
    mpfr_div_ui(u, u, 6 * k, MPFR_RNDN);

    prev = 0;
    memset(have_diff, 0, TAB_SIZE);

    for (i = 1; i < 3 * k; i++)
    {
        if (count[i] != 0)
        {
            if (prev > 2 && i - prev < TAB_SIZE)
            {
                if (!have_diff[i - prev])
                {
                    mpfr_init2(sin_tab[i - prev], wp);
                    mpfr_init2(cos_tab[i - prev], wp);
                    mpfr_mul_si(t, u, i - prev, MPFR_RNDN);
                    mpfr_sin_cos(sin_tab[i - prev],
                                 cos_tab[i - prev], t, GMP_RNDN);
                    have_diff[i - prev] = 1;
                }

                mpfr_mul(cc, cosa, cos_tab[i - prev], MPFR_RNDN);
                mpfr_mul(ss, sina, sin_tab[i - prev], MPFR_RNDN);
                mpfr_mul(cs, cosa, sin_tab[i - prev], MPFR_RNDN);
                mpfr_mul(sc, sina, cos_tab[i - prev], MPFR_RNDN);

                mpfr_sub(cosa, cc, ss, MPFR_RNDN);
                mpfr_add(sina, sc, cs, MPFR_RNDN);
            }
            else
            {
                mpfr_mul_si(t, u, i, MPFR_RNDN);
                mpfr_sin_cos(sina, cosa, t, GMP_RNDN);
            }

            if (count[i] != 1)
                mpfr_mul_si(t, cosa, count[i], MPFR_RNDN);

            mpfr_add(sum, sum, t, MPFR_RNDN);
            prev = i;
        }
    }

    for (i = 0; i < TAB_SIZE; i++)
    {
        if (have_diff[i])
        {
            mpfr_clear(sin_tab[i]);
            mpfr_clear(cos_tab[i]);
        }
    }

    mpfr_clear(t);
    mpfr_clear(u);
    mpfr_clear(cosa);
    mpfr_clear(sina);
    mpfr_clear(cc);
    mpfr_clear(ss);
    mpfr_clear(cs);
    mpfr_clear(sc);
    free(count);
}
