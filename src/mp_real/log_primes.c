/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "thread_support.h"
#include "impl.h"

/* Shared cache (thread_support.h) of fixed-point logarithms of the
   first num primes 2, 3, 5, ..., for the diophantine (multi-prime)
   argument reductions (exp_diophantine.c).

   Each entry occupies nc limbs: a unit limb on top (log 41 < 4, so it
   holds two bits), then nc - 1 fraction limbs, the lowest of which is
   a guard.  The stored value is the EXACT floor of
   log(p) 2^(FLINT_BITS (nc - 1)), so the top nv + 1 limbs of an entry
   are exactly the floor at nv fraction limbs: floors nest, and
   consumers read the prefix they need (_mp_real_machin_entry, impl.h),
   so one version serves every smaller precision.  The reductions rely on the entries
   being one-sided (at most the true value, short by less than one
   ulp of the limbs read).

   The values come from _mp_real_log_primes_vec (machin_bsplit.c:
   the Machin-type sets of machin_tab.c evaluated by binary splitting
   in mp_real arithmetic, each atanh term by Zuniga's series [Zun2025]
   for the logarithm of the same ratio beyond a crossover precision,
   or arb's static 4608-bit table of the first 13 logarithms at small
   precisions), floored through the balls. */

/* the shared cache of the table (thread_support.h, impl.h) */
static flint_shared_cache_struct _mp_real_log_primes_cache =
    FLINT_SHARED_CACHE_INIT(_mp_real_tab_free);

/* e_j = floor(v_j B^(nc - 1)) into nc limbs for j < num, if the balls
   determine these floors; returns 0 otherwise */
int
_mp_real_store_floors(nn_ptr e, slong nc, mp_real_struct * v, slong num)
{
    slong j;

    for (j = 0; j < num; j++)
    {
        /* v_j < B: the floor of (v_j / B) B^nc */
        mp_real_mul_2exp_si(v + j, v + j, -FLINT_BITS);
        if (!_mp_real_get_fixed_floor(e + j * nc, nc, v + j))
            return 0;
    }
    return 1;
}

/* the shared cache protocol: a version grown from the current one,
   recomputed as a whole (the Machin-type sets evaluate all entries
   together) */
static void *
_mp_real_log_primes_build(const void * old, void * args)
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
        _mp_real_log_primes_vec(v, num, nc + guard);
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
_mp_real_log_primes_ensure(slong num, slong nv)
{
    _mp_real_tab_req_struct R;

    FLINT_ASSERT(num >= 1 && num <= MP_REAL_LOG_PRIMES_MAX);

    R.nc = nv + 2;
    R.num = num;
    return flint_shared_cache_ensure(&_mp_real_log_primes_cache,
        _mp_real_tab_sufficient, _mp_real_log_primes_build, &R);
}


/* free the cache (no other thread may be using it: tuning programs) */
void
_mp_real_log_primes_clear(void)
{
    flint_shared_cache_clear(&_mp_real_log_primes_cache);
}

/* the number of value fraction limbs the cached entries currently
   carry (0 when empty) */
slong
_mp_real_log_primes_max_limbs(void)
{
    const _mp_real_tab_struct * T =
        flint_shared_cache_get(&_mp_real_log_primes_cache);

    return (T == NULL) ? 0 : T->nc - 2;
}
