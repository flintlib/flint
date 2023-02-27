/*
    Copyright (C) 2013-2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/thread_support.h"
#include "partitions.h"

/* defined in flint*/
#define NUMBER_OF_SMALL_PARTITIONS 128
FLINT_DLL extern const unsigned int partitions_lookup[NUMBER_OF_SMALL_PARTITIONS];

slong partitions_hrr_needed_terms(double n);

void
partitions_fmpz_fmpz_hrr(fmpz_t p, const fmpz_t n, int use_doubles)
{
    arb_t x;
    arf_t bound;
    slong N;

    arb_init(x);
    arf_init(bound);

    N = partitions_hrr_needed_terms(fmpz_get_d(n));

    partitions_hrr_sum_arb(x, n, 1, N, use_doubles);

    partitions_rademacher_bound(bound, n, N);
    arb_add_error_arf(x, bound);

    if (!arb_get_unique_fmpz(p, x))
    {
        flint_printf("not unique!\n");
        arb_printd(x, 50);
        flint_printf("\n");
        flint_abort();
    }

    arb_clear(x);
    arf_clear(bound);
}

/* To compute p(n) mod 2^64. */
static void
partitions_vec(mp_ptr v, slong len)
{
    slong i, j, n;
    mp_limb_t p;

    for (n = 0; n < FLINT_MIN(len, NUMBER_OF_SMALL_PARTITIONS); n++)
        v[n] = partitions_lookup[n];

    for (n = NUMBER_OF_SMALL_PARTITIONS; n < len; n++)
    {
        p = 0;
        for (i = 1, j = 1; j <= n; j += 3 * i + 1, i++)
            p = v[n - j] - p;
        if ((i & 1) == 0)
            p = -p;
        for (i = 1, j = 2; j <= n; j += 3 * i + 2, i++)
            p = v[n - j] - p;
        if ((i & 1) != 0)
            p = -p;
        v[n] = p;
    }
}

/* The floor+vec method *requires* n <= 1498 for floor(p(n)/2^64)
   to be equal to floor(T/2^64). It is faster up to n ~= 1200.
   With doubles, it is faster up to n ~= 500. */
void
_partitions_fmpz_ui(fmpz_t res, ulong n, int use_doubles)
{
    if (n < NUMBER_OF_SMALL_PARTITIONS)
    {
        fmpz_set_ui(res, partitions_lookup[n]);
    }
    else if (FLINT_BITS == 64 && (n < 500 || (!use_doubles && n < 1200)))
    {
        mp_ptr tmp = flint_malloc((n + 1) * sizeof(mp_limb_t));

        if (n < 417)  /* p(n) < 2^64 */
        {
            partitions_vec(tmp, n + 1);
            fmpz_set_ui(res, tmp[n]);
        }
        else
        {
            arb_t x;
            arb_init(x);
            fmpz_set_ui(res, n);
            partitions_leading_fmpz(x, res, 4 * sqrt(n) - 50);
            arb_mul_2exp_si(x, x, -64);
            arb_floor(x, x, 4 * sqrt(n) - 50);

            if (arb_get_unique_fmpz(res, x))
            {
                fmpz_mul_2exp(res, res, 64);
                partitions_vec(tmp, n + 1);
                fmpz_add_ui(res, res, tmp[n]);
            }
            else
            {
                flint_printf("warning: failed at %wu\n", n);
                fmpz_set_ui(res, n);
                partitions_fmpz_fmpz_hrr(res, res, use_doubles);
            }
            arb_clear(x);
        }
        flint_free(tmp);
    }
    else
    {
        fmpz_set_ui(res, n);
        partitions_fmpz_fmpz_hrr(res, res, use_doubles);
    }
}

void
partitions_fmpz_fmpz(fmpz_t res, const fmpz_t n, int use_doubles)
{
    if (fmpz_cmp_ui(n, 2000) < 0)
    {
        if (fmpz_sgn(n) < 0)
            fmpz_zero(res);
        else
            _partitions_fmpz_ui(res, *n, use_doubles);
    }
    else
    {
        partitions_fmpz_fmpz_hrr(res, n, use_doubles);
    }
}

void
partitions_fmpz_ui(fmpz_t res, ulong n)
{
    _partitions_fmpz_ui(res, n, 0);
}
