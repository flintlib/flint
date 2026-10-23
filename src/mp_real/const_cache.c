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
#include "mp_real.h"
#include "thread_support.h"
#include "impl.h"

/* The constants of the module, as balls and as limb arrays.

   Every constant c has a from-scratch ball evaluation
   _mp_real_const_X_compute (in its own file).  The limb version writes
   floor(c B^n): n fraction limbs for c in [0, 1), n fraction limbs and
   a units limb for c in [1, B).  Such a floor is computed from a ball
   at a few guard limbs more, retried with more guard limbs until the
   radius determines it uniquely, so it is exact and *err = 1 bounds
   c - floor(c B^n) B^-n in ulps.

   With cache = 1 the floors come from a shared cache
   (thread_support.h), one per constant, read by every thread.  Floors
   nest: the top limbs of floor(c B^N) are floor(c B^n) for n <= N, so
   any request up to the cached precision is a copy.  A longer request
   recomputes the entry at N = max(n + 5, 1.5 N_old) limbs, so that a
   sequence of slowly growing requests costs a constant factor over the
   last one (and the superseded entries, which stay allocated until
   flint_cleanup_master, a bounded multiple of it).  The ball version
   with cache = 1 reads the same entry and encloses the floor with a
   radius of one ulp. */

typedef void (* _const_compute_func)(mp_real_t, slong);

typedef struct
{
    _const_compute_func compute;
    int units;                  /* 1: c in [1, B), a units limb */
}
_const_info_t;

enum
{
    C_PI4 = MP_REAL_CONST_ID_PI4, C_LOG2 = MP_REAL_CONST_ID_LOG2, C_EULER,
    C_E, C_LOG10, C_CATALAN, C_ZETA3, C_ZETA5, C_GAMMA_1_3, C_GAMMA_1_4,
    C_2_DIV_PI = MP_REAL_CONST_ID_2_DIV_PI,
    C_INV_LOG2_FRAC = MP_REAL_CONST_ID_INV_LOG2_FRAC, C_NUM
};

static const _const_info_t _const_info[C_NUM] = {
    { _mp_real_const_pi4_compute, 0 },
    { _mp_real_const_log2_compute, 0 },
    { _mp_real_const_euler_compute, 0 },
    { _mp_real_const_e_compute, 1 },
    { _mp_real_const_log10_compute, 1 },
    { _mp_real_const_catalan_compute, 0 },
    { _mp_real_const_zeta3_compute, 1 },
    { _mp_real_const_zeta5_compute, 1 },
    { _mp_real_const_gamma_1_3_compute, 1 },
    { _mp_real_const_gamma_1_4_compute, 1 },
    { _mp_real_const_2_div_pi_compute, 0 },
    { _mp_real_const_inv_log2_frac_compute, 0 },
};

/* a version of the entry of a constant: floor(c B^n) in n + units limbs */
typedef struct
{
    slong n;
    ulong d[];
}
_const_entry_struct;

#define CONST_CACHE_INIT FLINT_SHARED_CACHE_INIT(flint_free)
static flint_shared_cache_struct _const_cache[C_NUM] = {
    CONST_CACHE_INIT, CONST_CACHE_INIT, CONST_CACHE_INIT, CONST_CACHE_INIT,
    CONST_CACHE_INIT, CONST_CACHE_INIT, CONST_CACHE_INIT, CONST_CACHE_INIT,
    CONST_CACHE_INIT, CONST_CACHE_INIT, CONST_CACHE_INIT, CONST_CACHE_INIT
};

/* free the caches: not while another thread may be using them */
void
_mp_real_const_clear_cache(void)
{
    slong i;

    for (i = 0; i < C_NUM; i++)
        flint_shared_cache_clear(_const_cache + i);
}

/* floor(c B^n) into (e, n + units) */
static void
_const_floor(nn_ptr e, int id, slong n)
{
    const _const_info_t * c = _const_info + id;
    slong guard = 3;
    mp_real_t v;

    mp_real_init(v);
    for (;;)
    {
        c->compute(v, n + guard + c->units);
        /* for c in [1, B): floor(c B^n) = floor((c / B) B^(n + 1)) */
        if (c->units)
            mp_real_mul_2exp_si(v, v, -FLINT_BITS);
        if (_mp_real_get_fixed_floor(e, n + c->units, v))
            break;
        guard += 2 + guard / 2;
    }
    mp_real_clear(v);
}

/* the shared cache protocol: a request is (id, n) */
typedef struct
{
    int id;
    slong n;
}
_const_req_struct;

static int
_const_sufficient(const void * version, void * args)
{
    return ((const _const_entry_struct *) version)->n
        >= ((const _const_req_struct *) args)->n;
}

static void *
_const_build(const void * old, void * args)
{
    const _const_req_struct * R = args;
    int units = _const_info[R->id].units;
    slong N = R->n + 5;
    _const_entry_struct * e;

    if (old != NULL)
    {
        slong N_old = ((const _const_entry_struct *) old)->n;
        N = FLINT_MAX(N, N_old + N_old / 2);
    }

    e = flint_malloc(sizeof(_const_entry_struct) + (N + units) * sizeof(ulong));
    e->n = N;
    _const_floor(e->d, R->id, N);
    return e;
}

/* the cache entry of id covering at least n limbs */
static const _const_entry_struct *
_const_ensure(int id, slong n)
{
    _const_req_struct R;

    R.id = id;
    R.n = n;
    return flint_shared_cache_ensure(_const_cache + id, _const_sufficient,
        _const_build, &R);
}

static void
_const_limbs(nn_ptr res, ulong * err, int id, slong n, int cache)
{
    int units = _const_info[id].units;

    FLINT_ASSERT(n >= 1);

    if (!cache)
    {
        _const_floor(res, id, n);
    }
    else if (n <= MP_REAL_CONST_STATIC_N && (id == MP_REAL_CONST_ID_PI4
        || id == MP_REAL_CONST_ID_LOG2 || id == MP_REAL_CONST_ID_2_DIV_PI
        || id == MP_REAL_CONST_ID_INV_LOG2_FRAC))
    {
        /* the static tables (impl.h), so that the few-limb regimes
           touch no shared cache */
        flint_mpn_copyi(res, _mp_real_const_ptr(id, n), n);
    }
    else
    {
        const _const_entry_struct * e = _const_ensure(id, n);
        flint_mpn_copyi(res, e->d + (e->n - n), n + units);
    }

    if (err != NULL)
        *err = 1;
}

/* the top n limbs of the cached floor, extending the entry if needed
   (the identifiers of impl.h are the enumeration's) */
nn_srcptr
_mp_real_const_cached_ptr(int id, slong n)
{
    const _const_entry_struct * e = _const_ensure(id, n);
    return e->d + (e->n - n);
}

/* the number of limbs the shared cache holds for id (0 if none; the
   static tables cover MP_REAL_CONST_STATIC_N limbs of pi/4, log 2,
   2/pi and 1/log 2 - 1 regardless): lets a caller choose an algorithm
   that avoids a constant not yet computed */
slong
_mp_real_const_cached_limbs(int id)
{
    const _const_entry_struct * e = flint_shared_cache_get(_const_cache + id);
    return (e == NULL) ? 0 : e->n;
}

static void
_const_ball(mp_real_t res, int id, slong n, int cache)
{
    if (!cache)
    {
        _const_info[id].compute(res, n);
    }
    else
    {
        int units = _const_info[id].units;
        nn_ptr t;
        TMP_INIT;

        TMP_START;
        t = TMP_ALLOC((n + units) * sizeof(ulong));
        _const_limbs(t, NULL, id, n, 1);
        _mp_real_set_mpn_2exp(res, t, n + units, -FLINT_BITS * n);
        _mp_real_add_error_ulps_at(res, 1.0, -n);
        TMP_END;
    }
}

#define DEF_CONST(name, id)                                             \
void                                                                    \
mp_real_const_ ## name(mp_real_t res, slong n, int cache)               \
{                                                                       \
    _const_ball(res, id, n, cache);                                     \
}                                                                       \
                                                                        \
void                                                                    \
_mp_real_const_ ## name(nn_ptr res, ulong * err, slong n, int cache)    \
{                                                                       \
    _const_limbs(res, err, id, n, cache);                               \
}

DEF_CONST(pi4, C_PI4)
DEF_CONST(log2, C_LOG2)
DEF_CONST(euler, C_EULER)
DEF_CONST(e, C_E)
DEF_CONST(log10, C_LOG10)
DEF_CONST(catalan, C_CATALAN)
DEF_CONST(zeta3, C_ZETA3)
DEF_CONST(zeta5, C_ZETA5)
DEF_CONST(gamma_1_3, C_GAMMA_1_3)
DEF_CONST(gamma_1_4, C_GAMMA_1_4)
DEF_CONST(2_div_pi, C_2_DIV_PI)
DEF_CONST(inv_log2_frac, C_INV_LOG2_FRAC)
