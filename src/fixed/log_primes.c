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
#include "fixed.h"

/* Thread-local cache of fixed-point logarithms of the first
   _fixed_log_primes_num primes 2, 3, 5, ..., for the diophantine (multi-prime)
   argument reductions (exp_diophantine.c).

   Each entry occupies _fixed_log_primes_n limbs: a unit limb on top
   (log 41 < 4, so it holds two bits), then _fixed_log_primes_n - 1
   fraction limbs, the lowest of which is a guard.  The stored value
   is the EXACT floor of log(p) 2^(FLINT_BITS (n - 1)), so the top
   nv + 1 limbs of an entry are exactly the floor at nv fraction
   limbs: floors nest, and consumers read the prefix they need
   (_fixed_log_primes_entry).  The reductions rely on the entries
   being one-sided (at most the true value, short by less than one
   ulp of the limbs read).

   The values come from _fixed_log_primes_vec_fball (machin_bsplit.c:
   the Machin-type sets of machin_tab.c evaluated by binary splitting
   in fball arithmetic, each atanh term by Zuniga's series [Zun2025]
   for the logarithm of the same ratio beyond a crossover precision,
   or arb's static 4608-bit table of the first 13 logarithms at small
   precisions), floored through the balls. */

FLINT_TLS_PREFIX nn_ptr _fixed_log_primes = NULL;
FLINT_TLS_PREFIX slong _fixed_log_primes_n = 0;
FLINT_TLS_PREFIX slong _fixed_log_primes_num = 0;
static FLINT_TLS_PREFIX int _fixed_log_primes_cleanup_registered = 0;

static void
_fixed_log_primes_cleanup(void)
{
    flint_free(_fixed_log_primes);
    _fixed_log_primes = NULL;
    _fixed_log_primes_n = 0;
    _fixed_log_primes_num = 0;
    _fixed_log_primes_cleanup_registered = 0;
}

/* e_j = floor(v_j B^(nc - 1)) into nc limbs for j < num, if the balls
   determine these floors; returns 0 otherwise */
int
_fixed_store_floors(nn_ptr e, slong nc, fball_struct * v, slong num)
{
    slong j;

    for (j = 0; j < num; j++)
    {
        /* v_j < B: the floor of (v_j / B) B^nc */
        fball_mul_2exp_si(v + j, -FLINT_BITS);
        if (!fball_get_fixed_floor(e + j * nc, nc, v + j))
            return 0;
    }
    return 1;
}

void
_fixed_log_primes_ensure(slong num, slong nv)
{
    slong nc, j, guard;
    fball_struct * v;

    FLINT_ASSERT(num >= 1 && num <= FIXED_LOG_PRIMES_MAX);

    /* nv value fraction limbs, one guard limb, one unit limb */
    if (nv + 2 <= _fixed_log_primes_n && num <= _fixed_log_primes_num)
        return;

    nc = FLINT_MAX(nv + 2, _fixed_log_primes_n);
    num = FLINT_MAX(num, _fixed_log_primes_num);

    flint_free(_fixed_log_primes);
    _fixed_log_primes = flint_malloc(num * nc * sizeof(ulong));
    _fixed_log_primes_n = nc;
    _fixed_log_primes_num = num;

    v = flint_malloc(num * sizeof(fball_struct));
    for (j = 0; j < num; j++)
        fball_init(v + j);

    /* an undetermined floor (a value within the radius of a grid
       point: astronomically rare) is retried at higher precision */
    for (guard = 2; ; guard *= 2)
    {
        _fixed_log_primes_vec_fball(v, num, nc + guard);
        if (_fixed_store_floors(_fixed_log_primes, nc, v, num))
            break;
    }

    for (j = 0; j < num; j++)
        fball_clear(v + j);
    flint_free(v);

    if (!_fixed_log_primes_cleanup_registered)
    {
        flint_register_cleanup_function(_fixed_log_primes_cleanup);
        _fixed_log_primes_cleanup_registered = 1;
    }
}

/* the top nv + 1 limbs of entry j: nv fraction limbs with the unit
   limb at index nv; valid until the next ensure call on this thread */
nn_srcptr
_fixed_log_primes_entry(slong j, slong nv)
{
    FLINT_ASSERT(j >= 0 && j < _fixed_log_primes_num);
    FLINT_ASSERT(nv >= 1 && nv + 2 <= _fixed_log_primes_n);
    return _fixed_log_primes + j * _fixed_log_primes_n
        + (_fixed_log_primes_n - (nv + 1));
}

void
_fixed_log_primes_clear(void)
{
    flint_free(_fixed_log_primes);
    _fixed_log_primes = NULL;
    _fixed_log_primes_n = 0;
    _fixed_log_primes_num = 0;
}

/* the number of value fraction limbs the cached entries currently
   carry (0 when empty) */
slong
_fixed_log_primes_max_limbs(void)
{
    return FLINT_MAX(_fixed_log_primes_n - 2, 0);
}
