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
#include "arb.h"
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

   The values are computed by arb_log_primes_vec_bsplit (through
   arb's own log(p) cache for up to 13 primes, which also serves the
   small precisions from its static table; the Machin-type sets of
   machin_tab.c cover up to 48 primes) and floored through arb.
   An mpn binary splitting of the same Machin formulas was
   prototyped (see the fixed-log-primes-bsplit patch): 5-13% faster
   than arb's, mostly through a better product dispatch for the
   short-by-long products of the tree -- a change that belongs in
   flint_mpn_mul itself, after which the wrapper here inherits it. */

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

/* e = floor(x 2^(64 (nc - 1))) into nc limbs if x's radius determines
   that floor; returns 0 otherwise */
static int
_store_floor(nn_ptr e, const arb_t x, slong nc, slong prec)
{
    arb_t y;
    fmpz_t f;
    int ok;

    arb_init(y);
    fmpz_init(f);
    arb_mul_2exp_si(y, x, FLINT_BITS * (nc - 1));
    arb_floor(y, y, FLINT_MAX(prec, FLINT_BITS * nc) + 64);
    ok = arb_get_unique_fmpz(f, y);
    if (ok)
    {
        FLINT_ASSERT(fmpz_sgn(f) > 0 && fmpz_bits(f) <= FLINT_BITS * nc);
        fmpz_get_ui_array(e, nc, f);
    }
    arb_clear(y);
    fmpz_clear(f);
    return ok;
}

void
_fixed_log_primes_ensure(slong num, slong nv)
{
    slong nc, j, wp;
    arb_ptr lp;
    n_primes_t iter;

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

    wp = FLINT_BITS * nc + 64;
    if (num <= ARB_LOG_PRIME_CACHE_NUM)
    {
        _arb_log_p_ensure_cached(wp);
        lp = (arb_ptr) _arb_log_p_cache_vec();
    }
    else
    {
        lp = _arb_vec_init(num);
        arb_log_primes_vec_bsplit(lp, num, wp);
    }

    n_primes_init(iter);
    for (j = 0; j < num; j++)
    {
        nn_ptr e = _fixed_log_primes + j * nc;
        ulong p = n_primes_next(iter);

        if (!_store_floor(e, lp + j, nc, wp))
        {
            /* the floor sits within the radius of an integer
               (astronomically rare): recompute at increasing precision
               until it is determined */
            arb_t x;
            slong p2;

            arb_init(x);
            for (p2 = 2 * wp; ; p2 *= 2)
            {
                arb_log_ui(x, p, p2);
                if (_store_floor(e, x, nc, p2))
                    break;
            }
            arb_clear(x);
        }
    }
    n_primes_clear(iter);

    if (num > ARB_LOG_PRIME_CACHE_NUM)
        _arb_vec_clear(lp, num);

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
