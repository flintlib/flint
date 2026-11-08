/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"
#include "dense_rt_tables.h"

#if FLINT_BITS == 64

#if FLINT_USES_PTHREAD
#include <stdatomic.h>
#endif
#include "thread_support.h"

/*
    Dense-table kernels beyond the static per-size kernels, at
    MP_REAL_DENSE_W + 1 to MP_REAL_DENSE_RT_MAX fraction limbs, with
    tables partly computed at run time.  The reductions are those of the
    static per-size kernels (see dev/gen_mp_real_dense_tables.py and
    dev/gen_mp_real_dense_bodies.py): nine levels of seven bits leave a
    remainder below (1 + 2^-8) 2^-63, and a tenth level takes the rest of
    the top limb (index j <= 2, multiplier j / B, four table entries),
    leaving it below 2^-64 except in rare cases (a lowered table value or
    index estimate at a boundary; then below 2^-63).  So the hardcoded
    series of the 2^-64 tier (exp_rs.c, trig_rs.c, up to 21, 22 and 35
    limbs for exp, sin/cos and atan/atanh) apply directly, the rare cases
    taking the 2^-32 tier, and beyond them the generic series of
    series_rs.c.

    Every multiplier c is one limb (j 2^-7k for the levels of seven bits,
    as in the static kernels).  For log1p, D_k is the one-limb value of
    the static tables up to level nine and c - c^2 (exact, two limbs) at
    the tenth.  The index bounds are those of the static tables.

    The tables hold floors (or lower bounds within one ulp); a kernel at
    n limbs reads the top n limbs of an entry.  The first two levels of
    every family (and the first level of sin and cos) are static, at
    RT_S = MP_REAL_DENSE_RT_W = 32 limbs (dense_rt_tables.c,
    dev/gen_mp_real_dense_rt_tables.py, 215 KB): their arguments are the
    largest, so they are the most expensive to compute.  The levels from
    the third on are filled as they are needed, in two tiers: at RT_WN =
    21 limbs for the kernels up to 21 limbs and at RT_WW = 32 limbs
    beyond (an entry costs about the square of its width, so a process
    working at the smaller sizes fills only the cheaper ones).  The
    storage of a family in a tier is a shared cache (thread_support.h)
    with a single version, allocated on its first use (afterwards an
    acquire load, a plain load on x86-64, finds it); each entry has a
    state byte (empty, being written, ready).  A kernel reading an entry
    that is not ready computes it into a buffer of its own (at one more
    limb, by the Horner scheme of _dense_horner_entry, lowered by its
    error bound: a lower bound within one ulp), uses that, and stores it
    in the table if it wins the compare-and-swap of the state from empty
    to being written, publishing it with a release store of the state
    ready.  So a process pays only for the entries its arguments touch
    (about 0.5-3 us per level and new index at 21 limbs, 1-5 us at 32), and
    no reader ever waits.  The storage lives until flint_cleanup_master:
    163, 152 and 158 KB for exp, log1p and atan (with sin and cos) in
    the narrow tier, 249, 231 and 241 KB in the wide one, committed by
    the operating system as entries are written.

    At RT_SL_MIN .. RT_SL_MAX limbs the kernels run as straight-line
    instances (see below): unrolled carry chains and, with the x86-64
    ADX assembly, rows in one pass of mulx, adcx and adox
    (dense_rt_asm.inc); log1p ends with the hand series of
    dense_rt_series.inc at 11 to 14 limbs.  Every variant gives the same
    result bits.
*/

#define LV MP_REAL_DENSE_RT_LEVELS
#define BITS 7

#if LV != 10
#error "the index bounds assume ten levels"
#endif

/* index bounds per level (the largest j), as in dense_tables.h */
static const short _exp_jmax[LV + 1] = { 127, 127, 128, 191, 129, 128, 128, 128, 128, 128, 3 };
static const short _log_jmax[LV + 1] = { 127, 127, 127, 127, 127, 127, 128, 128, 128, 128, 3 };
static const short _atan_jmax[LV + 1] = { 127, 127, 133, 133, 133, 133, 133, 133, 133, 133, 3 };

/* the weight 2^-e of the index unit at level k: seven bits per level,
   the last level the rest of the top limb */
#define LEVEL_E(k) (((k) < LV) ? BITS * (k) : FLINT_BITS)

#define RT_EXP MP_REAL_DENSE_RT_EXP
#define RT_LOG MP_REAL_DENSE_RT_LOG1P
#define RT_ATAN MP_REAL_DENSE_RT_ATAN

/* The static levels (dense_rt_tables.c) hold floors at RT_S limbs; the
   run-time levels 3 .. LV are kept in two tiers, the narrow one at
   RT_WN limbs for n <= RT_WN and the wide one at RT_WW limbs beyond (so
   that a process working at the smaller sizes fills only the cheaper
   entries) */
#define RT_S MP_REAL_DENSE_RT_W
#define RT_WN 21
#define RT_WW MP_REAL_DENSE_RT_MAX
#define RT_WMAX FLINT_MAX(RT_WN, RT_WW)

#if RT_S < RT_WMAX
#error "the static tables must be generated at MP_REAL_DENSE_RT_MAX limbs at least"
#endif

/* the static levels 0 .. 2: level 0 is the first level of sin and cos
   (T_1, kept with the atan family), level 1 the first level of every
   family (exp: L_1) */
static const ulong * const _rt_static[3][3] = {
    { NULL, _mp_real_dense_rt_exp_l1, _mp_real_dense_rt_exp_g2 },
    { NULL, _mp_real_dense_rt_log_m1, _mp_real_dense_rt_log_m2 },
    { _mp_real_dense_rt_sin_cos_t1, _mp_real_dense_rt_atan_a1,
      _mp_real_dense_rt_atan_a2 },
};

static const short * const _rt_jmax[3] = { _exp_jmax, _log_jmax, _atan_jmax };

/* entry states */
#define RT_EMPTY 0
#define RT_BUSY 1
#define RT_READY 2

#if FLINT_USES_PTHREAD
typedef _Atomic(unsigned char) _rt_state_t;
#else
typedef unsigned char _rt_state_t;
#endif

/* the run-time levels 3 .. LV of a family in a tier: entry j of level
   k at the limbs tab + (off[k] + j) w, with the state state[off[k] + j] */
typedef struct
{
    slong off[LV + 2];
    slong w;
    ulong * tab;
    _rt_state_t * state;
}
_dense_rt_struct;

typedef _dense_rt_struct * _dense_rt_ptr;

static void
_dense_rt_free(void * version)
{
    _dense_rt_ptr T = version;

    flint_free(T->tab);
    flint_free((void *) T->state);
    flint_free(T);
}

#define DENSE_RT_CACHE_INIT FLINT_SHARED_CACHE_INIT(_dense_rt_free)
static flint_shared_cache_struct _dense_rt_cache[2][3] = {
    { DENSE_RT_CACHE_INIT, DENSE_RT_CACHE_INIT, DENSE_RT_CACHE_INIT },
    { DENSE_RT_CACHE_INIT, DENSE_RT_CACHE_INIT, DENSE_RT_CACHE_INIT }
};

/* entry computation ********************************************************/

/* D_k(j) as two limbs (d1, d0), the value (d1 + d0 / B) / B */
FLINT_FORCE_INLINE void
_dense_log_d(ulong * d1, ulong * d0, ulong j, slong k)
{
    if (k == 1)
    {
        *d1 = _mp_real_dense_log_d1[j];
        *d0 = 0;
    }
    else if (k <= 9)
    {
        slong s = BITS * k;
        ulong c2 = (2 * s <= FLINT_BITS) ? ((j * j) << (FLINT_BITS - 2 * s))
            : ((j * j + (UWORD(1) << (2 * s - FLINT_BITS)) - 1) >> (2 * s - FLINT_BITS));
        *d1 = (j << (FLINT_BITS - s)) - c2;
        *d0 = 0;
    }
    else
    {
        /* the last level: c = j / B, c^2 = j^2 / B^2 exactly */
        sub_ddmmss(*d1, *d0, j, 0, 0, j * j);
    }
}

/* Y, an approximation at m = w + 1 limbs within e units, lowered by
   e + 1 units (to zero if below) and truncated to w limbs at r: a lower
   bound within one ulp (and 2 e + 1 units of the extra limb) */
static void
_dense_round_entry(nn_ptr r, nn_ptr Y, slong w, ulong e)
{
    if (mpn_sub_1(Y, Y, w + 1, e + 1))
        flint_mpn_zero(Y, w + 1);
    flint_mpn_copyi(r, Y + 1, w);
}

/*
    An entry f(x) for a short argument x = (m1 + m0 / B) / B < 2^-12
    (one limb, or two at the last level of log1p; below 2^-13 from the
    third level on), by the Horner scheme in
    the multiplier x, every step a pass of mpn_mul_1 or two:

        log(1 + x)  = x V,  V = sum_(i >= 0) (-x)^i / (i + 1)
        -log(1 - x) = x V,  V = sum_(i >= 0) x^i / (i + 1)
        atan(x)     = x V,  V = sum_(i >= 0) (-x^2)^i / (2 i + 1)

    (z = x, or x^2 applied as x twice).  With q_i the denominators,
    V_i = 1/q_i + s z V_(i+1) >= 0 (for s = -1 as z V_(i+1) < 2^-12
    < 1/q_i).  The denominators run in groups of consecutive i, formed
    greedily from the top, with product Q < 2^62: within a group the sum
    is kept as U = Q V, so that a step adds the integer Q / q_i at the
    units limb (or subtracts z U from it); between groups U is rescaled
    by Q' / Q (one mpn_mul_1, one mpn_divrem_1).  At L = w + 1 fraction
    limbs, step i only needs the limbs above floor(i zb / 64), z < 2^-zb
    damping its errors by z^i; so each step errs by at most 2 units of
    the result (a floor per pass), each rescale and the final division
    by 1, the final product by 2, the tail by less than 1: below 3 N + 4
    units in all, N the number of terms.  Costs about N L / 2 limb
    products: 1.3 to 2.5 times less than the series of series_rs.c on
    these arguments (whose powers are full products).
*/

#define DH_LOG1P 0
#define DH_NLOG1M 1
#define DH_ATAN 2

/* (D, len) = (S, len) x / B for x = (m1 + m0 / B) / B (D != S),
   below by at most one unit (two with m0); D[-1] is scratch */
static void
_dh_mulx(nn_ptr D, nn_srcptr S, slong len, ulong m1, ulong m0)
{
    D[len - 1] = mpn_mul_1(D - 1, S, len, m1);
    if (m0 != 0 && len >= 2)
    {
        /* + S m0 / B^2, the low limb of S only reaching the dropped
           limbs */
        ulong c = mpn_addmul_1(D - 1, S + 1, len - 1, m0);
        add_ssaaaa(D[len - 1], D[len - 2], D[len - 1], D[len - 2], 0, c);
    }
}

static void
_dense_horner_entry(nn_ptr r, slong w, int kind, ulong m1, ulong m0)
{
    ulong A[RT_WMAX + 4], C[RT_WMAX + 4];
    nn_ptr U = A + 1, T = C + 1;
    slong L = w + 1, N, zb, i, lo, lo2, len, hi, glo;
    ulong Q, Qn, q, err;
    int alt = (kind != DH_NLOG1M), sq = (kind == DH_ATAN);

    if (m1 == 0 && m0 == 0)
    {
        flint_mpn_zero(r, w);
        return;
    }

    FLINT_ASSERT(m0 == 0 || !sq);

    /* x < 2^-zb, z < 2^-(zb (1 + sq)); terms while z^N > 2^-(64 L + 4).
       m1 = 0 for one entry, D = 2^-64 - 2^-128 at the last level of
       log1p (clz of zero is undefined: bsr leaves its destination
       alone) */
    zb = (m1 != 0) ? flint_clz(m1) : FLINT_BITS + flint_clz(m0);
    FLINT_ASSERT(zb >= 12);
    zb *= (1 + sq);
    N = (FLINT_BITS * L + 4 + zb - 1) / zb;

#define DH_Q(ii) (sq ? 2 * (ulong) (ii) + 1 : (ulong) (ii) + 1)
#define DH_LO(ii) FLINT_MIN((slong) (((ii) * zb) / FLINT_BITS), L)

    /* the top group, hi down to glo */
    hi = N - 1;
    Q = 1;
    for (glo = hi; glo >= 0 && Q <= (UWORD(1) << 62) / DH_Q(glo); glo--)
        Q *= DH_Q(glo);
    glo++;

    /* U = Q / q_(N-1) over the window of step N - 1 */
    lo = DH_LO(hi);
    flint_mpn_zero(U + lo, L - lo);
    U[L] = Q / DH_Q(hi);
    err = 0;

    for (i = hi - 1; i >= 0; i--)
    {
        q = DH_Q(i);

        /* a new group: rescale by Q' / Q */
        if (i < glo)
        {
            Qn = 1;
            for (glo = i; glo >= 0 && Qn <= (UWORD(1) << 62) / DH_Q(glo); glo--)
                Qn *= DH_Q(glo);
            glo++;
            len = L + 1 - lo;
            T[L + 1] = mpn_mul_1(T + lo, U + lo, len, Qn);
            mpn_divrem_1(U + lo, 0, T + lo, len + 1, Q);
            Q = Qn;
            err++;
        }

        /* widen the window to that of step i */
        lo2 = DH_LO(i);
        if (lo2 < lo)
        {
            flint_mpn_zero(U + lo2, lo - lo2);
            lo = lo2;
        }
        len = L + 1 - lo;

        /* T = z U */
        if (sq)
        {
            _dh_mulx(T + lo, U + lo, len, m1, 0);
            _dh_mulx(U + lo, T + lo, len, m1, 0);
            FLINT_SWAP(nn_ptr, U, T);
        }
        else
            _dh_mulx(T + lo, U + lo, len, m1, m0);

        /* U = Q / q -+ z U */
        if (alt)
        {
            ulong b = (len > 1) ? (mpn_neg(T + lo, T + lo, len - 1) != 0) : 0;
            T[L] = Q / q - T[L] - b;
        }
        else
            T[L] += Q / q;
        FLINT_SWAP(nn_ptr, U, T);
        err += 2;
    }

    /* V = U / Q, then x V */
    T[L + 1] = 0;
    mpn_divrem_1(T, 0, U, L + 1, Q);
    _dh_mulx(U, T, L + 1, m1, m0);
    err += 4;

#undef DH_Q
#undef DH_LO

    FLINT_ASSERT(U[L] == 0);
    _dense_round_entry(r, U, w, err);
}

/* entry j of level k >= 3 of a family at w limbs, by the Horner
   scheme, with c = j 2^-e:
   exp: log(1 + c);
   log1p: -log(1 - D), D = (d1 B + d0) / B^2;
   atan: atan(c) */
static void
_dense_entry(nn_ptr r, slong w, int which, slong k, ulong j)
{
    slong e = LEVEL_E(k);

    FLINT_ASSERT(k >= 3 && k <= LV);

    if (which == RT_EXP)
        _dense_horner_entry(r, w, DH_LOG1P, j << (FLINT_BITS - e), 0);
    else if (which == RT_LOG)
    {
        ulong d1, d0;
        _dense_log_d(&d1, &d0, j, k);
        _dense_horner_entry(r, w, DH_NLOG1M, d1, d0);
    }
    else
        _dense_horner_entry(r, w, DH_ATAN, j << (FLINT_BITS - e), 0);
}

/* the cache ****************************************************************/

/* the single version of a family in a tier (never replaced) */
static int
_dense_rt_sufficient(const void * FLINT_UNUSED(version), void * FLINT_UNUSED(args))
{
    return 1;
}

/* the storage of a family in a tier, allocated on its first use: the
   request is (tier, which) */
static void *
_dense_rt_build(const void * FLINT_UNUSED(old), void * args)
{
    const int * R = args;
    int tier = R[0], which = R[1];
    _dense_rt_ptr T;
    slong k, num = 0;

    T = flint_malloc(sizeof(_dense_rt_struct));
    T->w = tier ? RT_WW : RT_WN;
    for (k = 0; k <= LV + 1; k++)
    {
        T->off[k] = num;
        if (k >= 3 && k <= LV)
            num += _rt_jmax[which][k] + 1;
    }
    T->tab = flint_malloc(num * T->w * sizeof(ulong));
    T->state = flint_calloc(num, sizeof(_rt_state_t));
    return T;
}

FLINT_FORCE_INLINE const _dense_rt_struct *
_dense_rt(int tier, int which)
{
    const _dense_rt_struct * T = flint_shared_cache_get(&_dense_rt_cache[tier][which]);
    int R[2];

    if (FLINT_LIKELY(T != NULL))
        return T;
    R[0] = tier;
    R[1] = which;
    return flint_shared_cache_ensure(&_dense_rt_cache[tier][which],
        _dense_rt_sufficient, _dense_rt_build, R);
}

#if defined(__GNUC__)
#define RT_NOINLINE __attribute__((noinline))
#else
#define RT_NOINLINE
#endif

/* an entry that is not ready: computed into buf, and stored by the
   thread that claims it */
static RT_NOINLINE nn_srcptr
_dense_rt_fill(const _dense_rt_struct * T, int which, slong k, ulong j,
    nn_ptr buf)
{
    slong i = T->off[k] + j, w = T->w;

    _dense_entry(buf, w, which, k, j);
#if FLINT_USES_PTHREAD
    {
        unsigned char expected = RT_EMPTY;

        if (atomic_compare_exchange_strong_explicit(&T->state[i], &expected,
                RT_BUSY, memory_order_relaxed, memory_order_relaxed))
        {
            flint_mpn_copyi(T->tab + i * w, buf, w);
            atomic_store_explicit(&T->state[i], RT_READY, memory_order_release);
        }
    }
#else
    flint_mpn_copyi(T->tab + i * w, buf, w);
    T->state[i] = RT_READY;
#endif
    return buf;
}

/* the top n limbs of entry j of level k of a family: static for
   k <= 2, else from the cache of the tier (wide: n > RT_WN) or computed
   into buf (RT_WMAX limbs) */
FLINT_FORCE_INLINE nn_srcptr
_rt_entry(const _dense_rt_struct * T, int wide, int which, slong k,
    ulong j, slong n, nn_ptr buf)
{
    slong i, w = wide ? RT_WW : RT_WN;

    FLINT_ASSERT(k == 0 ? (which == RT_ATAN && j < 64)
                        : j <= (ulong) _rt_jmax[which][k]);

    if (k <= 2)
        return _rt_static[which][k] + j * RT_S + (RT_S - n);

    i = T->off[k] + j;
#if FLINT_USES_PTHREAD
    if (FLINT_LIKELY(atomic_load_explicit(&T->state[i], memory_order_acquire) == RT_READY))
#else
    if (FLINT_LIKELY(T->state[i] == RT_READY))
#endif
        return T->tab + i * w + (w - n);
    return _dense_rt_fill(T, which, k, j, buf) + (w - n);
}

/* the number of (tier, family) storages allocated so far (tests of
   which regimes touch the run-time tables) */
slong
_mp_real_dense_rt_allocated(void)
{
    slong count = 0;
    int tier, which;

    for (tier = 0; tier < 2; tier++)
        for (which = 0; which < 3; which++)
            count += (flint_shared_cache_get(&_dense_rt_cache[tier][which]) != NULL);
    return count;
}

/* entry j of level k (3 <= k <= the last level) of a family in the tier
   of the kernels at n limbs, filled if needed, into (r, w), with its
   argument x = (a[1] + a[0] / B) / B: the entry is a lower bound for
   log(1 + x) (exp), -log(1 - x) (log1p) or atan(x) within one ulp of
   B^-w (tests); returns w */
slong
_mp_real_dense_rt_entry(nn_ptr r, ulong * a, int which, slong n, slong k,
    ulong j)
{
    int wide = (n > RT_WN);
    const _dense_rt_struct * T = _dense_rt(wide, which);
    ulong buf[RT_WMAX];

    FLINT_ASSERT(k >= 3 && k <= LV && j <= (ulong) _rt_jmax[which][k]);

    flint_mpn_copyi(r, _rt_entry(T, wide, which, k, j, T->w, buf), T->w);
    if (which == RT_LOG)
        _dense_log_d(&a[1], &a[0], j, k);
    else
    {
        a[1] = j << (FLINT_BITS - LEVEL_E(k));
        a[0] = 0;
    }
    return T->w;
}

/* the largest index at level k of a family */
slong
_mp_real_dense_rt_jmax(int which, slong k)
{
    FLINT_ASSERT(k >= 0 && k <= LV);
    return _rt_jmax[which][k];
}

slong
_mp_real_dense_rt_prepare(int which, slong n)
{
    int wide = (n > RT_WN);
    const _dense_rt_struct * T = _dense_rt(wide, which);
    ulong buf[RT_WMAX];
    slong k, j;

    for (k = 3; k <= LV; k++)
        for (j = 0; j <= _rt_jmax[which][k]; j++)
            _rt_entry(T, wide, which, k, j, T->w, buf);
    return T->w;
}

/* kernel helpers ***********************************************************/

/* the top n limbs of an entry */

/* the level multiplier c = j 2^-e (e = LEVEL_E(k) <= 64) as one limb
   m, c = m / B */
FLINT_FORCE_INLINE ulong
_dense_mult(ulong j, slong k)
{
    return j << (FLINT_BITS - LEVEL_E(k));
}

/* floor(t 2^e) for the n-limb fraction t (t below 2^-(e - 8)) */
FLINT_FORCE_INLINE ulong
_dense_index(nn_srcptr t, slong n, slong k)
{
    return t[n - 1] >> (FLINT_BITS - LEVEL_E(k));
}

/*
    The kernels are written once, as inline bodies taking the size n and
    a flag sl, and instantiated (where add-with-carry intrinsics exist)
    for each n from RT_SL_MIN to RT_SL_MAX with sl = 1: with n a constant
    the limb loops of the helpers below unroll into straight-line
    products and carry chains, with no calls and no loop control.  The
    other sizes (and other machines) take sl = 0, the mpn functions.
*/

#define RT_SL_MIN 11
#define RT_SL_MAX 16

#if defined(__x86_64__) || defined(_M_X64)
# define RT_SL 1
# if defined(_MSC_VER) && !defined(__clang__)
#  include <intrin.h>
# else
#  include <x86intrin.h>
# endif
FLINT_FORCE_INLINE unsigned char
_rt_adc(unsigned char c, ulong a, ulong b, ulong * r)
{
    unsigned long long t;
    c = _addcarry_u64(c, a, b, &t);
    *r = t;
    return c;
}
FLINT_FORCE_INLINE unsigned char
_rt_sbb(unsigned char c, ulong a, ulong b, ulong * r)
{
    unsigned long long t;
    c = _subborrow_u64(c, a, b, &t);
    *r = t;
    return c;
}
#elif defined(__has_builtin)
# if __has_builtin(__builtin_addcll) && __has_builtin(__builtin_subcll)
#  define RT_SL 1
FLINT_FORCE_INLINE unsigned char
_rt_adc(unsigned char c, ulong a, ulong b, ulong * r)
{
    unsigned long long co;
    *r = __builtin_addcll(a, b, c, &co);
    return (unsigned char) co;
}
FLINT_FORCE_INLINE unsigned char
_rt_sbb(unsigned char c, ulong a, ulong b, ulong * r)
{
    unsigned long long co;
    *r = __builtin_subcll(a, b, c, &co);
    return (unsigned char) co;
}
# endif
#endif

#ifndef RT_SL
# define RT_SL 0
#endif

#if defined(__GNUC__)
# define RT_UNROLL _Pragma("GCC unroll 24")
#else
# define RT_UNROLL
#endif

/* the products as plain C where the compiler has a double-limb type (so
   that it schedules them freely), else the longlong.h macro; the guard
   and the __uint128_t spelling (a predefined typedef, accepted by
   -pedantic) as in crt_helpers.h */
#if defined(__GNUC__) && defined(__SIZEOF_INT128__) && FLINT_BITS == 64
# define RT_UMUL(hi, lo, a, b) \
    do { \
        __uint128_t __p = (__uint128_t) (a) * (b); \
        (lo) = (ulong) __p; \
        (hi) = (ulong) (__p >> 64); \
    } while (0)
#else
# define RT_UMUL umul_ppmm
#endif

/* (r, n) = (a, n) + (b, n), returning the carry */
FLINT_FORCE_INLINE ulong
_rt_add_n(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n, int sl)
{
#if RT_SL
    if (sl)
    {
        unsigned char c = 0;
        slong i;
        RT_UNROLL
        for (i = 0; i < n; i++)
            c = _rt_adc(c, a[i], b[i], r + i);
        return c;
    }
#endif
    return mpn_add_n(r, a, b, n);
}

/* (r, n) = (a, n) - (b, n), returning the borrow */
FLINT_FORCE_INLINE ulong
_rt_sub_n(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n, int sl)
{
#if RT_SL
    if (sl)
    {
        unsigned char c = 0;
        slong i;
        RT_UNROLL
        for (i = 0; i < n; i++)
            c = _rt_sbb(c, a[i], b[i], r + i);
        return c;
    }
#endif
    return mpn_sub_n(r, a, b, n);
}

/* (T, n + 1) = (F, n) m: the products, then one carry chain */
FLINT_FORCE_INLINE void
_rt_mul_1(nn_ptr T, nn_srcptr F, slong n, ulong m, int sl)
{
#if RT_SL
    if (sl)
    {
        ulong h[RT_SL_MAX], l[RT_SL_MAX];
        unsigned char c = 0;
        slong i;
        RT_UNROLL
        for (i = 0; i < n; i++)
            RT_UMUL(h[i], l[i], F[i], m);
        T[0] = l[0];
        RT_UNROLL
        for (i = 1; i < n; i++)
            c = _rt_adc(c, l[i], h[i - 1], T + i);
        T[n] = h[n - 1] + c;
        return;
    }
#endif
    T[n] = mpn_mul_1(T, F, n, m);
}

/* the rows by a one-limb c: with P = floor(s c / B) (the low limb of the
   product dropped), d += P (returning the carry out), d -= P (d >= P),
   P itself (returned at T or T + 1, T with room for n + 1 limbs); one
   pass of mulx, adcx and adox at RT_ASM_MIN .. RT_ASM_MAX limbs on
   builds with the x86-64 ADX assembly (dense_rt_asm.inc), else a
   product row and a carry chain; d may equal s */

#if RT_SL && defined(FLINT_HAVE_ASSEMBLY_x86_64_adx) && FLINT_HAVE_ASSEMBLY_x86_64_adx && defined(__GNUC__)
# define RT_ASM 1
# include "dense_rt_asm.inc"
#else
# define RT_ASM 0
#endif

FLINT_FORCE_INLINE ulong
_rt_addmul_shift(nn_ptr d, nn_srcptr s, slong n, ulong c, int sl)
{
    ulong T[RT_WMAX + 1];
#if RT_ASM
    if (sl && n >= RT_ASM_MIN && n <= RT_ASM_MAX)
        return _rt_asm_addmul_shift(d, s, n, c);
#endif
    _rt_mul_1(T, s, n, c, sl);
    return _rt_add_n(d, d, T + 1, n, sl);
}

FLINT_FORCE_INLINE void
_rt_submul_shift(nn_ptr d, nn_srcptr s, slong n, ulong c, int sl)
{
    ulong T[RT_WMAX + 1];
#if RT_ASM
    if (sl && n >= RT_ASM_MIN && n <= RT_ASM_MAX)
    {
        _rt_asm_submul_shift(d, s, n, c);
        return;
    }
#endif
    _rt_mul_1(T, s, n, c, sl);
    _rt_sub_n(d, d, T + 1, n, sl);
}

FLINT_FORCE_INLINE nn_srcptr
_rt_mul_shift(nn_ptr T, nn_srcptr s, slong n, ulong c, int sl)
{
#if RT_ASM
    if (sl && n >= RT_ASM_MIN && n <= RT_ASM_MAX)
    {
        _rt_asm_mul_shift(T, s, n, c);
        return T;
    }
#endif
    _rt_mul_1(T, s, n, c, sl);
    return T + 1;
}

/* F c at the scale of F, c = (m1 + m0 / B) / B, truncated (at most two
   ulps low), as the n limbs T + 1 of (T, n + 1) */
FLINT_FORCE_INLINE void
_dense_mulc(nn_ptr T, nn_srcptr F, slong n, ulong m1, ulong m0, int sl)
{
    _rt_mul_1(T, F, n, m1, sl);
    if (m0 != 0)
    {
        /* + (F m0) / B, the low limb of F only reaching the dropped
           limbs */
        ulong U[RT_WMAX + 1];
        _rt_mul_1(U, F + 1, n - 1, m0, sl);
        T[n] += _rt_add_n(T, T, U, n, sl);
    }
}

/* the kernels **************************************************************/

#include "hand_mulhi.inc"
#include "dense_rt_series.inc"

/* exp(v), v in [0, 1) at n fraction limbs: (y, n + 1) */
FLINT_FORCE_INLINE void
_rt_exp(nn_ptr y, nn_srcptr v, slong n, int sl)
{
    int wide = (n > RT_WN);
    const _dense_rt_struct * T = _dense_rt(wide, RT_EXP);
    slong k;
    ulong t[RT_WMAX], ys[RT_WMAX + 1], eb[RT_WMAX];
    ulong j[LV + 1], m, c, e, cy;

    FLINT_ASSERT(n >= 2 && n <= RT_WW);

    /* t = v - L_1[j_1], then t <- t - log(1 + c_k) */
    j[1] = v[n - 1] >> (FLINT_BITS - BITS);
    _rt_sub_n(t, v, _rt_entry(T, wide, RT_EXP, 1, j[1], n, eb), n, sl);
    /* (the second level, static, peeled off the loop) */
#define RT_STEP(kk) \
    do { \
        j[kk] = _dense_index(t, n, kk); \
        FLINT_ASSERT(j[kk] <= (ulong) _exp_jmax[kk]); \
        if (j[kk] != 0) \
            _rt_sub_n(t, t, _rt_entry(T, wide, RT_EXP, kk, j[kk], n, eb), n, sl); \
    } while (0)
    RT_STEP(2);
    for (k = 3; k <= LV; k++)
        RT_STEP(k);
#undef RT_STEP

    /* 1 + F = exp(t), t < 2^-64 (but for rare cases below 2^-63) */
    _mp_real_exp_rs(ys, NULL, t, n);

    /* 1 + F <- (1 + F)(1 + c_k): F + c + F c, F c truncated */
    for (k = LV; k >= 2; k--)
    {
        if (j[k] == 0)
            continue;
        m = _dense_mult(j[k], k);
        _rt_addmul_shift(ys, ys, n, m, sl);
        ys[n - 1] += m;
    }

    /* (e, F) = (1 + F)(1 + C_1[j_1]) < 2.02, then times 2^s */
    c = _mp_real_dense_exp_c1[j[1]];
    cy = _rt_addmul_shift(ys, ys, n, c, sl);
    add_ssaaaa(e, ys[n - 1], 0, ys[n - 1], 0, c);
    e += 1 + cy;

    if (j[1] >= MP_REAL_DENSE_EXP_S1)
    {
        y[n] = (e << 1) | (ys[n - 1] >> (FLINT_BITS - 1));
        mpn_lshift(y, ys, n, 1);
    }
    else
    {
        flint_mpn_copyi(y, ys, n);
        y[n] = e;
    }
}

/* log(1 + v), v in [0, 1) at n fraction limbs: (y, n) */
FLINT_FORCE_INLINE void
_rt_log1p(nn_ptr y, nn_srcptr v, slong n, int sl)
{
    int wide = (n > RT_WN);
    const _dense_rt_struct * T = _dense_rt(wide, RT_LOG);
    slong k, un;
    ulong eb[RT_WMAX];
    ulong u[RT_WMAX], acc[RT_WMAX];
    ulong P[RT_WMAX + 1], S[RT_WMAX + 1];
    ulong t[RT_WMAX], ys[RT_WMAX];
    ulong j, d1, d0;

    FLINT_ASSERT(n >= 2 && n <= RT_WW);

    flint_mpn_copyi(u, v, n);

    /* u <- u - D - u D (u D truncated), acc += M_k = -log(1 - D) */
    j = u[n - 1] >> (FLINT_BITS - BITS);
    flint_mpn_copyi(acc, _rt_entry(T, wide, RT_LOG, 1, j, n, eb), n);
    d1 = _mp_real_dense_log_d1[j];
    _rt_submul_shift(u, u, n, d1, sl);
    u[n - 1] -= d1;
    /* (the second level, static, peeled off the loop) */
#define RT_STEP(kk) \
    do { \
        j = _dense_index(u, n, kk); \
        FLINT_ASSERT(j <= (ulong) _log_jmax[kk]); \
        if (j != 0) \
        { \
            _rt_add_n(acc, acc, _rt_entry(T, wide, RT_LOG, kk, j, n, eb), n, sl); \
            _dense_log_d(&d1, &d0, j, kk); \
            _rt_submul_shift(u, u, n, d1, sl); \
            u[n - 1] -= d1; \
        } \
    } while (0)
    RT_STEP(2);
    for (k = 3; k < LV; k++)
        RT_STEP(k);
#undef RT_STEP
    /* the top level, D = j / B - j^2 / B^2 */
    j = u[n - 1];
    FLINT_ASSERT(j <= (ulong) _log_jmax[LV]);
    if (j != 0)
    {
        _rt_add_n(acc, acc, _rt_entry(T, wide, RT_LOG, LV, j, n, eb), n, sl);
        _dense_log_d(&d1, &d0, j, LV);
        _dense_mulc(P, u, n, d1, d0, sl);
        _rt_sub_n(u, u, P + 1, n, sl);
        sub_ddmmss(u[n - 1], u[n - 2], u[n - 1], u[n - 2], d1, d0);
    }

    /* + log(1 + u), u < 2^-64 (but for rare cases below 2^-63) */
    if ((sl || !RT_SL) && n >= RT_LOG1P_SERIES_MIN && n <= RT_LOG1P_SERIES_MAX)
    {
        /* the hand series (no division) */
        switch (n)
        {
            case 11: _rt_log1p_series_11(ys, u); break;
            case 12: _rt_log1p_series_12(ys, u); break;
            case 13: _rt_log1p_series_13(ys, u); break;
            default: _rt_log1p_series_14(ys, u); break;
        }
        _rt_add_n(y, acc, ys, n, sl);
        return;
    }

    /* 2 atanh(u / (2 + u)) */
    un = n;
    while (un > 0 && u[un - 1] == 0)
        un--;
    if (un == 0)
    {
        flint_mpn_copyi(y, acc, n);
        return;
    }
    flint_mpn_copyi(S, u, n);
    S[n] = 2;
    flint_mpn_zero(t, n);
    flint_mpn_divapprox_fraction(t, u, un, S, n + 1, n);
    _mp_real_atanh_rs(ys, NULL, t, n);
    mpn_lshift(ys, ys, n, 1);
    _rt_add_n(y, acc, ys, n, sl);
}

/* atan(v), v in [0, 1) at n fraction limbs: (y, n) */
FLINT_FORCE_INLINE void
_rt_atan(nn_ptr y, nn_srcptr v, slong n, int sl)
{
    int wide = (n > RT_WN);
    const _dense_rt_struct * T = _dense_rt(wide, RT_ATAN);
    slong k, yn;
    ulong eb[RT_WMAX];
    ulong X[RT_WMAX + 1], Y[RT_WMAX];
    ulong PX[RT_WMAX + 1];
    ulong acc[RT_WMAX], t[RT_WMAX], ys[RT_WMAX];
    ulong j, c, m, R, rr, q, Ys;
    nn_srcptr Q;

    FLINT_ASSERT(n >= 2 && n <= RT_WW);

    /* level 1: c = j_1 / 2^7 <= v, X = 1 + c v, Y = v - c */
    j = v[n - 1] >> (FLINT_BITS - BITS);
    c = j << (FLINT_BITS - BITS);
    flint_mpn_copyi(X, _rt_mul_shift(PX, v, n, c, sl), n);
    X[n] = 1;
    flint_mpn_copyi(Y, v, n);
    Y[n - 1] &= (UWORD(1) << (FLINT_BITS - BITS)) - 1;
    flint_mpn_copyi(acc, _rt_entry(T, wide, RT_ATAN, 1, j, n, eb), n);

    /* R ~ 2^64 / X */
    udiv_qrnnd(R, rr, (UWORD(1) << (FLINT_BITS - 1)) - 1, ~UWORD(0),
        (UWORD(1) << (FLINT_BITS - 1)) | (X[n - 1] >> 1));
    (void) rr;

    /* (the second level, static, peeled off the loop); j = floor(2^e Y
       / X) from the estimate, lowered by 2^-12, then (X, Y) <- (X + c Y,
       Y - c X) on the old values */
#define RT_STEP(kk) \
    do { \
        slong s = LEVEL_E(kk) - 9; \
        Ys = (Y[n - 1] << s) | (Y[n - 2] >> (FLINT_BITS - s)); \
        q = n_mulhi(Ys, R); \
        q -= q >> 12; \
        j = q >> (FLINT_BITS - 9); \
        FLINT_ASSERT(j <= (ulong) _atan_jmax[kk]); \
        if (j != 0) \
        { \
            m = _dense_mult(j, kk); \
            Q = _rt_mul_shift(PX, X, n, m, sl); \
            _rt_addmul_shift(X, Y, n, m, sl); \
            _rt_sub_n(Y, Y, Q, n, sl); \
            Y[n - 1] -= m; \
            _rt_add_n(acc, acc, _rt_entry(T, wide, RT_ATAN, kk, j, n, eb), n, sl); \
        } \
    } while (0)
    RT_STEP(2);
    for (k = 3; k <= LV; k++)
        RT_STEP(k);
#undef RT_STEP

    /* + atan(Y / X), Y / X < 2^-64 (but for rare cases below 2^-63) */
    yn = n;
    while (yn > 0 && Y[yn - 1] == 0)
        yn--;
    if (yn == 0)
    {
        flint_mpn_copyi(y, acc, n);
        return;
    }
    flint_mpn_zero(t, n);
    flint_mpn_divapprox_fraction(t, Y, yn, X, n + 1, n);
    _mp_real_atan_rs(ys, NULL, t, n);
    _rt_add_n(y, acc, ys, n, sl);
}

/* (sin v, cos v), v in [0, 1) at n fraction limbs: (ysin, n + 1),
   (ycos, n + 1) */
FLINT_FORCE_INLINE void
_rt_sin_cos(nn_ptr ysin, nn_ptr ycos, nn_srcptr v, slong n, int sl)
{
    int wide = (n > RT_WN);
    const _dense_rt_struct * T = _dense_rt(wide, RT_ATAN);
    slong k;
    ulong eb[RT_WMAX];
    ulong t[RT_WMAX], X[RT_WMAX + 1], Y[RT_WMAX + 1];
    ulong PY[RT_WMAX + 1];
    ulong X2[RT_WMAX + 1], Y2[RT_WMAX], XY[RT_WMAX];
    ulong NF[RT_WMAX + 1], I[RT_WMAX + 2];
    ulong j[LV + 1], m, di, si, ni;
    nn_srcptr Q;

    FLINT_ASSERT(n >= 2 && n <= RT_WW);

    /* t = w = v / 2 (the lowest bit of v dropped) */
    mpn_rshift(t, v, n, 1);

    /* t <- t - atan(c_k), c_1 = C_1[j_1] / B */
    j[1] = t[n - 1] >> (FLINT_BITS - BITS);
    _rt_sub_n(t, t, _rt_entry(T, wide, RT_ATAN, 0, j[1], n, eb), n, sl);
    /* (the second level, static, peeled off the loop) */
#define RT_STEP(kk) \
    do { \
        j[kk] = _dense_index(t, n, kk); \
        FLINT_ASSERT(j[kk] <= (ulong) _atan_jmax[kk]); \
        if (j[kk] != 0) \
            _rt_sub_n(t, t, _rt_entry(T, wide, RT_ATAN, kk, j[kk], n, eb), n, sl); \
    } while (0)
    RT_STEP(2);
    for (k = 3; k <= LV; k++)
        RT_STEP(k);
#undef RT_STEP

    /* (X, Y) = (cos t, sin t), X with its integral limb */
    _mp_real_sin_cos_rs(Y, X, NULL, t, n);

    /* (X, Y) <- (X - c Y, Y + c X) */
    for (k = LV; k >= 1; k--)
    {
        if (k == 1)
            m = _mp_real_dense_sin_cos_c1[j[1]];
        else if (j[k] == 0)
            continue;
        else
            m = _dense_mult(j[k], k);
        Q = _rt_mul_shift(PY, Y, n, m, sl);
        _rt_addmul_shift(Y, X, n, m, sl);
        if (X[n] != 0)
            Y[n - 1] += m;
        X[n] -= _rt_sub_n(X, X, Q, n, sl);
    }

    /* X^2, Y^2, XY, N = X^2 + Y^2 (in [1, 1.3] up to rounding) */
    flint_mpn_sqrhigh_n(X2, X, n);
    flint_mpn_sqrhigh_n(Y2, Y, n);
    flint_mpn_mulhigh_n(XY, X, Y, n);
    di = X[n];
    if (di)
    {
        /* X^2 = 1 + 2 F + F^2, XY = F Y + Y */
        di += _rt_add_n(X2, X2, X, n, sl);
        di += _rt_add_n(X2, X2, X, n, sl);
        _rt_add_n(XY, XY, Y, n, sl);
    }
    si = XY[n - 1] >> (FLINT_BITS - 1);
    mpn_lshift(XY, XY, n, 1);
    ni = di + _rt_add_n(NF, X2, Y2, n, sl);
    di -= _rt_sub_n(X2, X2, Y2, n, sl);

    /* 1/N; N >= 1 exactly, so a rounded N below 1 is taken as 1 */
    if (ni == 0)
        flint_mpn_zero(NF, n);
    if ((sl || !RT_SL) && n <= FLINT_MPN_RECIP_SMALL_MAX)
    {
        /* inline (with run-time n, a switch over its sizes): only in the
           straight-line instance at 11 limbs, the generic body (from 17
           limbs on) never reaching it */
        _flint_mpn_recip_small(I, NF, n, FLINT_BITS * n - 3);
    }
    else
    {
        /* floor(B^n / N) or one more, lowered by two ulps (below B^n) */
        NF[n] = 1;
        flint_mpn_invapprox(I, NF, n + 1, 2 * n);
        mpn_sub_1(I, I, n + 1, 2);
        if (I[n] != 0)
            flint_mpn_store(I, n, ~UWORD(0));
    }

    /* cos v = D / N, sin v = 2 XY / N (D, 2 XY with integral limbs di, si) */
    flint_mpn_mulhigh_n(ycos, X2, I, n);
    flint_mpn_mulhigh_n(ysin, XY, I, n);
    ycos[n] = di ? _rt_add_n(ycos, ycos, I, n, sl) : 0;
    ysin[n] = si ? _rt_add_n(ysin, ysin, I, n, sl) : 0;
}

/* the instances: straight-line at RT_SL_MIN .. RT_SL_MAX limbs */

#if RT_SL
#define RT_INSTANCES(N) \
static void _rt_exp_##N(nn_ptr y, nn_srcptr v) { _rt_exp(y, v, N, 1); } \
static void _rt_log1p_##N(nn_ptr y, nn_srcptr v) { _rt_log1p(y, v, N, 1); } \
static void _rt_atan_##N(nn_ptr y, nn_srcptr v) { _rt_atan(y, v, N, 1); } \
static void _rt_sin_cos_##N(nn_ptr ys, nn_ptr yc, nn_srcptr v) { _rt_sin_cos(ys, yc, v, N, 1); }
RT_INSTANCES(11)
RT_INSTANCES(12)
RT_INSTANCES(13)
RT_INSTANCES(14)
RT_INSTANCES(15)
RT_INSTANCES(16)

typedef void (*_rt_fn)(nn_ptr, nn_srcptr);
typedef void (*_rt_fn2)(nn_ptr, nn_ptr, nn_srcptr);
static const _rt_fn _rt_exp_tab[] = { _rt_exp_11, _rt_exp_12, _rt_exp_13,
    _rt_exp_14, _rt_exp_15, _rt_exp_16 };
static const _rt_fn _rt_log1p_tab[] = { _rt_log1p_11, _rt_log1p_12,
    _rt_log1p_13, _rt_log1p_14, _rt_log1p_15, _rt_log1p_16 };
static const _rt_fn _rt_atan_tab[] = { _rt_atan_11, _rt_atan_12,
    _rt_atan_13, _rt_atan_14, _rt_atan_15, _rt_atan_16 };
static const _rt_fn2 _rt_sin_cos_tab[] = { _rt_sin_cos_11, _rt_sin_cos_12,
    _rt_sin_cos_13, _rt_sin_cos_14, _rt_sin_cos_15, _rt_sin_cos_16 };
#define RT_IS_SL(n) ((n) >= RT_SL_MIN && (n) <= RT_SL_MAX)
#else
#define RT_IS_SL(n) 0
#endif

void
_mp_real_exp_dense_rt(nn_ptr y, nn_srcptr v, slong n)
{
#if RT_SL
    if (RT_IS_SL(n))
    {
        _rt_exp_tab[n - RT_SL_MIN](y, v);
        return;
    }
#endif
    _rt_exp(y, v, n, 0);
}

void
_mp_real_log1p_dense_rt(nn_ptr y, nn_srcptr v, slong n)
{
#if RT_SL
    if (RT_IS_SL(n))
    {
        _rt_log1p_tab[n - RT_SL_MIN](y, v);
        return;
    }
#endif
    _rt_log1p(y, v, n, 0);
}

void
_mp_real_atan_dense_rt(nn_ptr y, nn_srcptr v, slong n)
{
#if RT_SL
    if (RT_IS_SL(n))
    {
        _rt_atan_tab[n - RT_SL_MIN](y, v);
        return;
    }
#endif
    _rt_atan(y, v, n, 0);
}

void
_mp_real_sin_cos_dense_rt(nn_ptr ysin, nn_ptr ycos, nn_srcptr v, slong n)
{
#if RT_SL
    if (RT_IS_SL(n))
    {
        _rt_sin_cos_tab[n - RT_SL_MIN](ysin, ycos, v);
        return;
    }
#endif
    _rt_sin_cos(ysin, ycos, v, n, 0);
}

#endif
