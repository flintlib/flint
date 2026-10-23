/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "thread_support.h"
#include "impl.h"

/* Shared cache (thread_support.h) of the angles theta_0 = pi/2,
   theta_j = 2 arg(pi_j) of the first num nonreal Gaussian primes, for
   the diophantine (multi-prime) trigonometric reduction
   (sin_cos_diophantine.c), laid out exactly like the logarithm cache
   of log_primes.c: nc limbs per entry, a unit limb
   on top (theta < pi < 4), the lowest fraction limb a guard, each
   entry the exact floor of theta B^(n-1) so that shorter reads are
   exact floors too.  The values come from _mp_real_atan_gauss_vec
   (machin_bsplit.c: the Machin-type sets of machin_tab.c by binary
   splitting in mp_real arithmetic, single arctangents of small quotients
   between neighbouring primes beyond, arb's static table of the first
   13 at small precisions). */

/* real and imaginary parts of the first 64 nonreal Gaussian primes */
const signed char _mp_real_gaussian_primes[2 * MP_REAL_ATAN_GAUSS_MAX] = {
    1, 1, 1, 2, 2, 3, 1, 4, 2, 5, 1, 6, 4, 5, 2, 7, 5, 6, 3, 8, 5, 8, 4, 9,
    1, 10, 3, 10, 7, 8, 4, 11, 7, 10, 6, 11, 2, 13, 9, 10, 7, 12, 1, 14, 2, 15, 8, 13,
    4, 15, 1, 16, 10, 13, 9, 14, 5, 16, 2, 17, 12, 13, 11, 14, 9, 16, 5, 18, 8, 17, 7, 18,
    10, 17, 6, 19, 1, 20, 3, 20, 14, 15, 12, 17, 7, 20, 4, 21, 10, 19, 5, 22, 11, 20, 10, 21,
    14, 19, 13, 20, 1, 24, 8, 23, 5, 24, 17, 18, 16, 19, 4, 25, 13, 22, 6, 25, 12, 23, 1, 26,
    5, 26, 15, 22, 2, 27, 9, 26
};

/* the shared cache of the table (thread_support.h, impl.h) */
static flint_shared_cache_struct _mp_real_atan_gauss_cache =
    FLINT_SHARED_CACHE_INIT(_mp_real_tab_free);

/* res_j = theta_j = 2 arg(pi_j), j < num, as arb balls at prec bits */
void
_mp_real_atan_gauss_vec_arb(arb_ptr res, slong num, slong prec)
{
    mp_real_struct * v;
    slong j;

    v = flint_malloc(num * sizeof(mp_real_struct));
    for (j = 0; j < num; j++)
        mp_real_init(v + j);

    _mp_real_atan_gauss_vec(v, num, (prec + FLINT_BITS - 1) / FLINT_BITS + 2);

    for (j = 0; j < num; j++)
    {
        mp_real_get_arb(res + j, v + j);
        arb_set_round(res + j, res + j, prec);
        mp_real_clear(v + j);
    }
    flint_free(v);
}

/* the shared cache protocol: a version grown from the current one,
   recomputed as a whole (the Machin-type sets evaluate all entries
   together) */
static void *
_mp_real_atan_gauss_build(const void * old, void * args)
{
    _mp_real_tab_struct * T;
    slong nc, num, j, guard, i_from;
    mp_real_struct * v;

    T = _mp_real_tab_new(old, args, &i_from);
    nc = T->nc;
    num = T->num;
    if (i_from >= num)
        return T;

    v = flint_malloc(num * sizeof(mp_real_struct));
    for (j = 0; j < num; j++)
        mp_real_init(v + j);

    /* undetermined floors (a value within the radius of a grid point:
       astronomically rare) are retried at higher precision */
    for (guard = 2; ; guard *= 2)
    {
        _mp_real_atan_gauss_vec(v, num, nc + guard);
        if (_mp_real_store_floors(T->tab, nc, v, num))
            break;
    }

    for (j = 0; j < num; j++)
        mp_real_clear(v + j);
    flint_free(v);

    return T;
}

/* the version covering num entries at nv value fraction limbs (one
   guard limb below, one unit limb above) */
const _mp_real_tab_struct *
_mp_real_atan_gauss_ensure(slong num, slong nv)
{
    _mp_real_tab_req_struct R;

    FLINT_ASSERT(num >= 1 && num <= MP_REAL_ATAN_GAUSS_MAX);

    R.nc = nv + 2;
    R.num = num;
    return flint_shared_cache_ensure(&_mp_real_atan_gauss_cache,
        _mp_real_tab_sufficient, _mp_real_atan_gauss_build, &R);
}


/* free the cache (no other thread may be using it: tuning programs) */
void
_mp_real_atan_gauss_clear(void)
{
    flint_shared_cache_clear(&_mp_real_atan_gauss_cache);
}
