/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "thread_support.h"
#include "mpn_extras.h"
#include "impl.h"

/* Shell handlers for sums of three squares (F = I_3), streaming the
   counts from r3 without a table of C:

   r3:     integral shifts (S = 1, g = 1, r0 = 0) with no character or the
           character (1/2,1/2,1/2): C(M) = r3(M), at the residue
           M + sum s_i mod 2 for the character;
   r3_odd: shifts in (1/2,1/2,1/2) + Z^3 with no character (S = 2,
           g = 8, r0 = 3): C(M) = r3(8M + 3), from a table or, for long
           sums (or variant 1), from the t2 divisor sieve with one sparse
           convolution per block. */

typedef struct
{
    ulong roff;
    uint32_t * tab;
    uint16_t * t2;
    ulong len;
}
r3_struct;

static void
r3_clear(void * data)
{
    r3_struct * R = data;
    flint_free(R->tab);
    flint_free(R->t2);
    flint_free(R);
}

/* r3(V) <= 2 (2 sqrt(V) + 1)^2 <= 18 V for V >= 1 (choose two
   coordinates, the third up to sign); V = g M + r0 < g len + r0 */
static ulong
r3_max_count(const void * data, const _acb_mzeta_shells_struct * T)
{
    return 18 * (T->g * T->len + T->r0) + 6;
}

static double
r3_cost_common(ulong len, ulong maxlen)
{
    return (len > maxlen) ? INFINITY : 20.0 * len;
}

/* integral shifts */

static void *
r3_init(const _acb_mzeta_shells_struct * T)
{
    r3_struct * R;
    int ell_ones = 1;
    ulong ssum = 0;
    slong i;

    if (T->d != 3 || !fmpz_mat_is_one(T->F) || T->S != 1 || T->g != 1 || T->r0 != 0)
        return NULL;
    for (i = 0; i < 3; i++)
    {
        ell_ones = ell_ones && fmpz_is_one(T->ell + i);
        ssum += fmpz_fdiv_ui(T->s + i, 2);
    }
    if (!(T->N == 1 || (T->N == 2 && ell_ones)))
        return NULL;

    R = flint_calloc(1, sizeof(r3_struct));
    R->roff = ssum % 2;
    return R;
}

static double
r3_cost(const void * data, const _acb_mzeta_shells_struct * T, ulong len)
{
    return r3_cost_common(len, UWORD(1) << 26);
}

static int
r3_compute(void * data, const _acb_mzeta_shells_struct * T, ulong len)
{
    r3_struct * R = data;
    flint_free(R->tab);
    R->tab = _acb_mzeta_r3_all_table(len);
    R->len = len;
    return 1;
}

static void
r3_block(ulong * out, const void * data, const _acb_mzeta_shells_struct * T, ulong M0, ulong n)
{
    const r3_struct * R = data;
    ulong N = T->N, i;

    flint_mpn_zero(out, n * N);
    for (i = 0; i < n; i++)
        out[i * N + ((N == 2) ? (M0 + i + R->roff) % 2 : 0)] = R->tab[M0 + i];
}

const _acb_mzeta_shells_handler_struct _acb_mzeta_shells_r3 =
{
    "r3", r3_init, r3_clear, r3_cost, r3_compute, r3_block, r3_max_count, NULL
};

/* shifts in (1/2,1/2,1/2) + Z^3 */

static void *
r3_odd_init(const _acb_mzeta_shells_struct * T)
{
    slong i;

    if (T->d != 3 || !fmpz_mat_is_one(T->F) || T->S != 2 || T->N != 1
        || T->g != 8 || T->r0 != 3)
        return NULL;
    for (i = 0; i < 3; i++)
        if (!fmpz_is_odd(T->s + i))
            return NULL;

    return flint_calloc(1, sizeof(r3_struct));
}

static double
r3_odd_cost(const void * data, const _acb_mzeta_shells_struct * T, ulong len)
{
    /* 8 len + 3 and 4 len + 1 (the t2 sieve) fit in a ulong */
    return r3_cost_common(len, (UWORD(1) << (FLINT_BITS - 4)) - 1);
}

static int
r3_odd_compute(void * data, const _acb_mzeta_shells_struct * T, ulong len)
{
    r3_struct * R = data;

    flint_free(R->tab);
    flint_free(R->t2);
    R->tab = NULL;
    R->t2 = NULL;
    if (len <= (UWORD(1) << 22) && T->variant != 1)
        R->tab = _acb_mzeta_r3_3mod8_table(len);
    else
        R->t2 = _acb_mzeta_t2_table(len, flint_get_num_threads());
    R->len = len;
    return 1;
}

static void
r3_odd_block(ulong * out, const void * data, const _acb_mzeta_shells_struct * T, ulong M0, ulong n)
{
    const r3_struct * R = data;
    ulong i;

    if (R->tab != NULL)
    {
        for (i = 0; i < n; i++)
            out[i] = R->tab[M0 + i];
    }
    else
    {
        uint32_t * tmp = flint_malloc(sizeof(uint32_t) * FLINT_MAX(n, 1));
        _acb_mzeta_r3_3mod8_block(tmp, R->t2, M0, n);
        for (i = 0; i < n; i++)
            out[i] = tmp[i];
        flint_free(tmp);
    }
}

const _acb_mzeta_shells_handler_struct _acb_mzeta_shells_r3_odd =
{
    "r3_odd", r3_odd_init, r3_clear, r3_odd_cost, r3_odd_compute, r3_odd_block, r3_max_count, NULL
};
