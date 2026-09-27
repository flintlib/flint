/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "thread_support.h"
#include "fmpz_poly/impl.h"

#if FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#endif

/*
    Multimodular GCD.

    Let A, B be primitive with len(A) >= len(B), G = gcd(A, B) with
    lc(G) > 0 and gamma = gcd(lc(A), lc(B)). For a prime p not dividing
    lc(A) lc(B), deg gcd(A mod p, B mod p) >= deg G with equality except
    for finitely many (unlucky) p, in which case the monic gcd modulo p
    scaled by gamma is the image of H = (gamma / lc(G)) G. We collect such
    images for a set of primes of the minimal degree seen and reconstruct H
    by the CRT; the candidate G is the primitive part of H. The result is
    certified by checking that G divides A and B.

    Compared to the classical algorithm which handles one prime at a time:

    * Primes are processed in batches: A and B are reduced modulo all primes
      of a batch at once (using fmpz_multi_mod_ui where this is faster than
      reducing modulo each prime separately), and the modular GCDs of a
      batch are independent and are computed in parallel.

    * The residues are kept and combined with a single fast multi-CRT when a
      reconstruction is attempted, instead of an incremental CRT per prime
      (which costs quadratic time in the number of primes).

    * Output-sensitive termination: only a random linear combination
      s = sum +/- H_i is maintained by incremental CRT (cheap: a single
      integer). A full reconstruction is attempted only when s has
      stabilized, i.e. when its symmetric residue is much smaller than the
      modulus. The batch size grows geometrically, so few primes are wasted.

    * Cheaper certification: exact division over Z is expensive, so we
      instead compute the quotients A / G and B / G multimodularly (their
      images modulo the primes used for the GCD are cheap to compute as
      divisions of nmod_polys, and more primes are added if needed, again
      with early termination) and certify by checking G * (A / G) = A
      and G * (B / G) = B. When the modulus is large enough, this follows
      from the congruences and coefficient bounds; otherwise we multiply.

    * The primes are chosen for the nmod_poly arithmetic. Unless the
      coefficients are small, the first prime is a 20-bit prime (the modular
      gcd is cheapest with small primes, and the first prime usually
      suffices to prove that the gcd is trivial).
      Then follow 62-bit primes (which allow Montgomery arithmetic in the
      Euclidean GCD), or for long polynomials 50-bit FFT primes for which
      fft_small multiplication needs a single transform (making half-GCD
      up to twice as fast per bit).

    * If gamma = gcd(lc(A), lc(B)) is larger than gcd(A(0), B(0)), we
      normalise the images by their constant coefficients instead, which
      avoids needlessly large images.

    * If deg gcd(A mod p, B mod p) = deg B, we directly test whether B
      divides A.
*/

/* Stop condition for the probes: the symmetric residue must be
   this many bits smaller than the modulus. If the value is not yet
   determined, this happens with probability about 2^(1 - PROBE_MARGIN);
   the cost of a spurious reconstruction attempt is small since most wrong
   candidates are rejected quickly. */
#define PROBE_MARGIN 8

/* Verification is cheap for short candidate gcds; for these we attempt
   reconstruction with a smaller margin after the first primes (typically,
   when the first small prime shows that the gcd is nontrivial). */
#define SHORT_GCD_LEN 8
#define SHORT_GCD_PROBE_MARGIN 2

/* Minimum number of (limbs x primes) to reduce per thread */
#define MULTI_MOD_MIN_WORK_PER_THREAD 20000

/* Use FFT primes when len2 is at least this. */
#define FFT_PRIMES_CUTOFF 2500

/* Primes of the default sequence (largest primes < 2^62); tabulated
   to avoid primality tests for small problems. */
#if FLINT_BITS == 64
#define NUM_TAB_PRIMES 63
static const ulong _gcd_primes_tab[NUM_TAB_PRIMES] = {
    UWORD(0x3fffffffffffffc7), UWORD(0x3fffffffffffffa9), UWORD(0x3fffffffffffff8b),
    UWORD(0x3fffffffffffff71), UWORD(0x3fffffffffffff67), UWORD(0x3fffffffffffff59),
    UWORD(0x3fffffffffffff55), UWORD(0x3fffffffffffff3d), UWORD(0x3fffffffffffff35),
    UWORD(0x3ffffffffffffeef), UWORD(0x3ffffffffffffee1), UWORD(0x3ffffffffffffec3),
    UWORD(0x3ffffffffffffe45), UWORD(0x3ffffffffffffe1d), UWORD(0x3ffffffffffffe11),
    UWORD(0x3ffffffffffffdc1), UWORD(0x3ffffffffffffdbb), UWORD(0x3ffffffffffffda5),
    UWORD(0x3ffffffffffffd87), UWORD(0x3ffffffffffffd69), UWORD(0x3ffffffffffffd03),
    UWORD(0x3ffffffffffffcfb), UWORD(0x3ffffffffffffcf7), UWORD(0x3ffffffffffffce9),
    UWORD(0x3ffffffffffffcd3), UWORD(0x3ffffffffffffcc1), UWORD(0x3ffffffffffffc65),
    UWORD(0x3ffffffffffffc2b), UWORD(0x3ffffffffffffc1f), UWORD(0x3ffffffffffffc17),
    UWORD(0x3ffffffffffffc11), UWORD(0x3ffffffffffffc07), UWORD(0x3ffffffffffffb53),
    UWORD(0x3ffffffffffffb27), UWORD(0x3ffffffffffffaf3), UWORD(0x3ffffffffffffab7),
    UWORD(0x3ffffffffffffa67), UWORD(0x3ffffffffffffa15), UWORD(0x3ffffffffffff9ef),
    UWORD(0x3ffffffffffff9d9), UWORD(0x3ffffffffffff9d3), UWORD(0x3ffffffffffff9c5),
    UWORD(0x3ffffffffffff9af), UWORD(0x3ffffffffffff977), UWORD(0x3ffffffffffff95f),
    UWORD(0x3ffffffffffff95b), UWORD(0x3ffffffffffff959), UWORD(0x3ffffffffffff8e1),
    UWORD(0x3ffffffffffff8a7), UWORD(0x3ffffffffffff889), UWORD(0x3ffffffffffff87d),
    UWORD(0x3ffffffffffff805), UWORD(0x3ffffffffffff7e7), UWORD(0x3ffffffffffff7c9),
    UWORD(0x3ffffffffffff7a3), UWORD(0x3ffffffffffff775), UWORD(0x3ffffffffffff757),
    UWORD(0x3ffffffffffff739), UWORD(0x3ffffffffffff713), UWORD(0x3ffffffffffff6d1),
    UWORD(0x3ffffffffffff6c1), UWORD(0x3ffffffffffff6b9), UWORD(0x3ffffffffffff6a3),
};
#else
#define NUM_TAB_PRIMES 0
#endif

/* Prime sequence *************************************************************/

typedef struct
{
    ulong p;            /* last prime returned (in the counting phase) */
    slong index;        /* number of primes returned from the tables */
    slong fft_count;    /* number of fft_small context primes to use first */
    ulong fft_c;        /* if nonzero, next FFT prime candidate c 2^32 + 1 */
    int up;             /* count upwards (testing mode) */
    int first;          /* first prime not yet returned */
}
_prime_iter_t;

static void
_prime_iter_init(_prime_iter_t * it, slong len2, ulong start, int small_first)
{
    (void) len2;    /* only used with fft_small */

    it->index = 0;
    it->fft_count = 0;
    it->fft_c = 0;
    it->first = small_first;

    if (start != 0)
    {
        /* testing mode: primes > start, counting upwards (start must leave
           room for enough primes below 2^FLINT_BITS) */
        FLINT_ASSERT(start < (UWORD(1) << (FLINT_BITS - 1)));
        it->p = start;
        it->up = 1;
        return;
    }

    it->up = 0;

#if FLINT_HAVE_FFT_SMALL && FLINT_BITS == 64
    if (len2 >= FFT_PRIMES_CUTOFF)
    {
        it->fft_count = MPN_CTX_NCRTS;
        it->fft_c = (UWORD(1) << 18) - 1;
    }
#endif

#if NUM_TAB_PRIMES
    it->p = _gcd_primes_tab[NUM_TAB_PRIMES - 1];
#else
    it->p = (UWORD(1) << (FLINT_BITS - 2)) + 1;
#endif
}

static ulong
_prime_iter_next(_prime_iter_t * it)
{
    if (it->up)
    {
        it->p = n_nextprime(it->p, 1);
        return it->p;
    }

#if FLINT_BITS == 64
    if (it->first)
    {
        it->first = 0;
        return UWORD(1048573);
    }
#endif

#if FLINT_HAVE_FFT_SMALL && FLINT_BITS == 64
    if (it->fft_c != 0)
    {
        mpn_ctx_struct * R = get_default_mpn_ctx();
        ulong q;
        slong i;

        /* First the primes of the default fft_small context (for which
           multiplication needs a single transform), then other 50-bit
           primes q = c 2^32 + 1 (for which fft_small also computes
           products with a single transform). */
        if (it->index < it->fft_count)
            return R->ffts[it->index++].mod.n;

        while (it->fft_c >= (UWORD(1) << 17))
        {
            q = it->fft_c * (UWORD(1) << 32) + 1;
            it->fft_c--;

            for (i = 0; i < it->fft_count; i++)
                if (q == R->ffts[i].mod.n)
                    break;

            if (i == it->fft_count && n_is_prime(q))
                return q;
        }

        /* Exhausted (not reached in practice); continue with the
           default primes. */
        it->fft_c = 0;
        it->index = 0;
        it->fft_count = 0;
    }
#endif

#if NUM_TAB_PRIMES
    if (it->index < it->fft_count + NUM_TAB_PRIMES)
    {
        FLINT_ASSERT(n_is_prime(_gcd_primes_tab[it->index - it->fft_count]));
        return _gcd_primes_tab[(it->index++) - it->fft_count];
    }
#endif

    do {
        it->p -= 2;
    } while (!n_is_prime(it->p));

    return it->p;
}

/* Get the next batch of n primes and set up the reduction */
static void
_get_primes(nn_ptr primes, nmod_t * mods, slong n, _prime_iter_t * it)
{
    slong k;

    for (k = 0; k < n; k++)
    {
        primes[k] = _prime_iter_next(it);
        nmod_init(mods + k, primes[k]);
    }
}

/* Multimodular reduction ****************************************************/

/* Reduce the vector (A, len) modulo all primes; the residues modulo prime k
   go to res + k * stride.

   Entries that are not much larger than the product P of the primes are
   reduced using the comb (subproduct tree). Much larger entries are first
   reduced modulo P when there are many primes (this is what makes the
   whole algorithm quasilinear in the bit size: reducing separately modulo
   each prime costs O(bits) per prime), and otherwise reduced separately
   modulo each prime (which is faster for few primes). */
#define MULTI_MOD_REDUCE_PROD_MIN_PRIMES 16

static int
_use_comb(const fmpz_t c, slong num_primes)
{
    return num_primes > 1 && COEFF_IS_MPZ(*c) &&
        ((slong) fmpz_size(c) <= 2 * num_primes || num_primes >= MULTI_MOD_REDUCE_PROD_MIN_PRIMES);
}

static void
_fmpz_vec_multi_mod_ui_strided(nn_ptr res, slong stride, const fmpz * A, slong len,
    const fmpz_comb_struct * comb, const fmpz_t prod,
    const nmod_t * mods, slong num_primes)
{
    slong j, k;
    fmpz_comb_temp_t temp;
    nn_ptr tmp = NULL;
    fmpz_t r;

    fmpz_init(r);

    if (comb != NULL)
    {
        fmpz_comb_temp_init(temp, comb);
        tmp = flint_malloc(sizeof(ulong) * num_primes);
    }

    for (j = 0; j < len; j++)
    {
        fmpz c = A[j];

        if (!COEFF_IS_MPZ(c))
        {
            if (c >= 0)
            {
                for (k = 0; k < num_primes; k++)
                    res[k * stride + j] = ((ulong) c < mods[k].n) ? (ulong) c : nmod_set_ui(c, mods[k]);
            }
            else
            {
                ulong v = -(ulong) c;
                for (k = 0; k < num_primes; k++)
                    res[k * stride + j] = nmod_neg((v < mods[k].n) ? v : nmod_set_ui(v, mods[k]), mods[k]);
            }
        }
        else if (comb == NULL || !_use_comb(A + j, num_primes))
        {
            for (k = 0; k < num_primes; k++)
                res[k * stride + j] = fmpz_get_nmod(A + j, mods[k]);
        }
        else
        {
            if ((slong) fmpz_size(A + j) > 2 * num_primes)
            {
                fmpz_fdiv_r(r, A + j, prod);
                fmpz_multi_mod_ui(tmp, r, comb, temp);
            }
            else
            {
                fmpz_multi_mod_ui(tmp, A + j, comb, temp);
            }

            for (k = 0; k < num_primes; k++)
                res[k * stride + j] = tmp[k];
        }
    }

    if (comb != NULL)
    {
        fmpz_comb_temp_clear(temp);
        flint_free(tmp);
    }

    fmpz_clear(r);
}

typedef struct
{
    nn_ptr res[3];
    slong stride[3];
    const fmpz * vec[3];
    slong len[3];
    slong chunk;
    const fmpz_comb_struct * comb;
    const fmpz * prod;
    nn_srcptr primes;
    const nmod_t * mods;
    slong num_primes;
}
_multi_mod_arg_t;

static void
_multi_mod_worker(slong i, void * varg)
{
    _multi_mod_arg_t * arg = (_multi_mod_arg_t *) varg;
    slong v, start, stop;

    /* the chunks cover the concatenation of the three vectors */
    start = i * arg->chunk;
    stop = start + arg->chunk;

    for (v = 0; v < 3; v++)
    {
        slong a = FLINT_MAX(start, 0);
        slong b = FLINT_MIN(stop, arg->len[v]);

        if (a < b)
            _fmpz_vec_multi_mod_ui_strided(arg->res[v] + a, arg->stride[v],
                arg->vec[v] + a, b - a, arg->comb, arg->prod, arg->mods, arg->num_primes);

        start -= arg->len[v];
        stop -= arg->len[v];
    }
}

/* Reduce up to three vectors modulo a batch of primes */
static void
_multi_mod_batch(nn_ptr Ares, slong Astride, const fmpz * A, slong lenA,
                 nn_ptr Bres, slong Bstride, const fmpz * B, slong lenB,
                 nn_ptr Cres, slong Cstride, const fmpz * C, slong lenC,
                 nn_srcptr primes, const nmod_t * mods, slong num_primes)
{
    fmpz_comb_t comb;
    fmpz_t prod;
    int use_comb = 0;
    slong j, total_len, limbs, num_chunks, num_threads;
    _multi_mod_arg_t arg;

    total_len = lenA + lenB + lenC;

    /* Estimate the work and decide whether the comb is needed */
    limbs = 0;
    for (j = 0; j < lenA; j++)
    {
        limbs += fmpz_size(A + j);
        use_comb = use_comb || _use_comb(A + j, num_primes);
    }
    for (j = 0; j < lenB; j++)
    {
        limbs += fmpz_size(B + j);
        use_comb = use_comb || _use_comb(B + j, num_primes);
    }
    for (j = 0; j < lenC; j++)
    {
        limbs += fmpz_size(C + j);
        use_comb = use_comb || _use_comb(C + j, num_primes);
    }

    fmpz_init(prod);

    if (use_comb)
    {
        fmpz_comb_init2(comb, primes, num_primes, FMPZ_COMB_MOD);
        _fmpz_ui_vec_prod(prod, primes, num_primes);
    }

    arg.res[0] = Ares; arg.stride[0] = Astride; arg.vec[0] = A; arg.len[0] = lenA;
    arg.res[1] = Bres; arg.stride[1] = Bstride; arg.vec[1] = B; arg.len[1] = lenB;
    arg.res[2] = Cres; arg.stride[2] = Cstride; arg.vec[2] = C; arg.len[2] = lenC;
    arg.comb = use_comb ? comb : NULL;
    arg.prod = prod;
    arg.primes = primes;
    arg.mods = mods;
    arg.num_primes = num_primes;

    num_threads = ((limbs + total_len) * num_primes) / MULTI_MOD_MIN_WORK_PER_THREAD;
    num_threads = FLINT_MAX(num_threads, 1);
    num_threads = FLINT_MIN(num_threads, flint_get_num_available_threads() + 1);
    num_chunks = (num_threads == 1) ? 1 : 4 * num_threads;
    arg.chunk = (total_len + num_chunks - 1) / num_chunks;
    num_chunks = (total_len + arg.chunk - 1) / arg.chunk;

    if (num_threads == 1)
        _multi_mod_worker(0, &arg);
    else
        flint_parallel_do(_multi_mod_worker, &arg, num_chunks, num_threads,
            FLINT_PARALLEL_DYNAMIC);

    if (use_comb)
        fmpz_comb_clear(comb);

    fmpz_clear(prod);
}

/* Residue cache ****************************************************************/

/*
    Supplies the residues of up to three vectors modulo successive primes.
    The residues are computed in chunks whose size doubles (reducing modulo
    a large batch of primes at once is much cheaper per prime than reducing
    modulo a few primes at a time), and handed out in whatever smaller
    batches the caller asks for. For each prime, the residues of the three
    vectors are stored consecutively (stride len[0] + len[1] + len[2]).
*/
typedef struct
{
    const fmpz * vec[3];
    slong len[3];
    slong stride;
    nn_ptr res;
    nn_ptr primes;
    nmod_t * mods;
    slong start;        /* first unused cached prime */
    slong num;          /* number of cached primes */
    slong alloc;
    slong fetched;      /* total number of primes fetched */
    _prime_iter_t * it;
}
_residue_cache_t;

static void
_residue_cache_init(_residue_cache_t * rc, const fmpz * A, slong lenA,
    const fmpz * B, slong lenB, const fmpz * C, slong lenC, _prime_iter_t * it)
{
    rc->vec[0] = A; rc->len[0] = lenA;
    rc->vec[1] = B; rc->len[1] = lenB;
    rc->vec[2] = C; rc->len[2] = lenC;
    rc->stride = lenA + lenB + lenC;
    rc->res = NULL;
    rc->primes = NULL;
    rc->mods = NULL;
    rc->start = rc->num = rc->alloc = rc->fetched = 0;
    rc->it = it;
}

static void
_residue_cache_clear(_residue_cache_t * rc)
{
    flint_free(rc->res);
    flint_free(rc->primes);
    flint_free(rc->mods);
}

/* Makes n primes available at rc->start (the caller advances rc->start
   by the number it uses). At most max_fetch primes (but at least n) are
   fetched at once. */
static void
_residue_cache_require(_residue_cache_t * rc, slong n, slong max_fetch)
{
    slong avail = rc->num - rc->start;
    slong m;

    if (avail >= n)
        return;

    /* move the unused residues to the front */
    if (rc->start != 0)
    {
        memmove(rc->res, rc->res + rc->start * rc->stride, sizeof(ulong) * avail * rc->stride);
        memmove(rc->primes, rc->primes + rc->start, sizeof(ulong) * avail);
        memmove(rc->mods, rc->mods + rc->start, sizeof(nmod_t) * avail);
        rc->num = avail;
        rc->start = 0;
    }

    m = FLINT_MAX(n - avail, FLINT_MIN(rc->fetched, max_fetch));
    m = FLINT_MAX(m, 1);

    if (rc->num + m > rc->alloc)
    {
        rc->alloc = rc->num + m;
        rc->res = flint_realloc(rc->res, sizeof(ulong) * rc->alloc * rc->stride);
        rc->primes = flint_realloc(rc->primes, sizeof(ulong) * rc->alloc);
        rc->mods = flint_realloc(rc->mods, sizeof(nmod_t) * rc->alloc);
    }

    _get_primes(rc->primes + rc->num, rc->mods + rc->num, m, rc->it);

    {
        nn_ptr r = rc->res + rc->num * rc->stride;
        _multi_mod_batch(r, rc->stride, rc->vec[0], rc->len[0],
                         r + rc->len[0], rc->stride, rc->vec[1], rc->len[1],
                         r + rc->len[0] + rc->len[1], rc->stride, rc->vec[2], rc->len[2],
                         rc->primes + rc->num, rc->mods + rc->num, m);
    }

    rc->num += m;
    rc->fetched += m;
}

/* Probes *********************************************************************/

/* The probe for a vector (h_0, ..., h_{len-1}) is s = sum_j w_j h_j with
   deterministic pseudorandom signs w_j = +/- 1. (Small weights keep s
   small, so that stabilization is detected as early as possible.) */
#if FLINT_BITS == 64
#define PROBE_HASH_MULTIPLIER UWORD(0x9e3779b97f4a7c15)
#else
#define PROBE_HASH_MULTIPLIER UWORD(0x9e3779b9)
#endif

static int
_weight_is_negative(slong j)
{
    return ((((ulong) j + 1) * PROBE_HASH_MULTIPLIER) >> (FLINT_BITS - 1)) != 0;
}

static ulong
_probe_nmod(nn_srcptr h, slong len, nmod_t mod)
{
    ulong s = 0;
    slong j;

    for (j = 0; j < len; j++)
    {
        if (_weight_is_negative(j))
            s = nmod_sub(s, h[j], mod);
        else
            s = nmod_add(s, h[j], mod);
    }

    return s;
}

/* Probe value s modulo sprod. New residues are collected and combined
   with s a batch at a time (using a fast multi-CRT for the batch), so that
   the total cost is quasilinear in the number of primes. */
typedef struct
{
    fmpz_t s;
    fmpz_t sprod;
    nn_ptr pend_res;        /* pending residues */
    nn_ptr pend_primes;     /* and their primes */
    slong pend_len;
    slong pend_alloc;
}
_probe_t;

static void
_probe_init(_probe_t * P)
{
    fmpz_init(P->s);
    fmpz_init_set_ui(P->sprod, 1);
    P->pend_res = P->pend_primes = NULL;
    P->pend_len = P->pend_alloc = 0;
}

static void
_probe_clear(_probe_t * P)
{
    fmpz_clear(P->s);
    fmpz_clear(P->sprod);
    flint_free(P->pend_res);
    flint_free(P->pend_primes);
}

static void
_probe_reset(_probe_t * P)
{
    fmpz_zero(P->s);
    fmpz_one(P->sprod);
    P->pend_len = 0;
}

static void
_probe_add(_probe_t * P, nn_srcptr h, slong len, nmod_t mod)
{
    if (P->pend_len == P->pend_alloc)
    {
        P->pend_alloc = FLINT_MAX(2 * P->pend_alloc, 16);
        P->pend_res = flint_realloc(P->pend_res, sizeof(ulong) * P->pend_alloc);
        P->pend_primes = flint_realloc(P->pend_primes, sizeof(ulong) * P->pend_alloc);
    }

    P->pend_res[P->pend_len] = _probe_nmod(h, len, mod);
    P->pend_primes[P->pend_len] = mod.n;
    P->pend_len++;
}

static void
_probe_flush(_probe_t * P)
{
    fmpz_t t, tprod;

    if (P->pend_len == 0)
        return;

    if (P->pend_len == 1)
    {
        fmpz_CRT_ui(P->s, P->s, P->sprod, P->pend_res[0], P->pend_primes[0], 1);
        fmpz_mul_ui(P->sprod, P->sprod, P->pend_primes[0]);
    }
    else
    {
        fmpz_init(t);
        fmpz_init(tprod);

        fmpz_multi_CRT_ui_once(t, tprod, P->pend_res, P->pend_primes, P->pend_len, 1);

        if (fmpz_is_one(P->sprod))
            fmpz_swap(P->s, t);
        else
            fmpz_CRT(P->s, P->s, P->sprod, t, tprod, 1);

        fmpz_mul(P->sprod, P->sprod, tprod);

        fmpz_clear(t);
        fmpz_clear(tprod);
    }

    P->pend_len = 0;
}

static int
_probe_stable_margin(_probe_t * P, slong margin)
{
    _probe_flush(P);
    return (slong) fmpz_bits(P->s) + margin <= (slong) fmpz_bits(P->sprod);
}

static int
_probe_stable(_probe_t * P)
{
    return _probe_stable_margin(P, PROBE_MARGIN);
}

/* CRT a set of residue vectors stored with the given stride */
static void
_crt_strided(fmpz * res, nn_srcptr R, slong stride, slong len,
                nn_srcptr primes, slong num_primes)
{
    nn_srcptr * residues;
    slong k;

    residues = flint_malloc(sizeof(nn_srcptr) * num_primes);
    for (k = 0; k < num_primes; k++)
        residues[k] = R + k * stride;

    _fmpz_vec_multi_CRT_ui(res, residues, len, primes, num_primes, 1);
    flint_free(residues);
}

/* Modular GCDs ***************************************************************/

typedef struct
{
    nn_srcptr ABmod;    /* images of A (len1), then B (len2); stride ABstride */
    slong ABstride;
    nn_ptr Hmod;        /* stride len2 */
    slong * hlen;
    const nmod_t * mods;
    slong len1;
    slong len2;
}
_gcd_batch_arg_t;

static void
_gcd_batch_worker(slong k, void * varg)
{
    _gcd_batch_arg_t * arg = (_gcd_batch_arg_t *) varg;
    slong len1 = arg->len1, len2 = arg->len2;
    nn_srcptr a = arg->ABmod + k * arg->ABstride;
    nn_srcptr b = a + len1;

    /* The prime divides a leading coefficient */
    if (a[len1 - 1] == 0 || b[len2 - 1] == 0)
        arg->hlen[k] = -1;
    else
        arg->hlen[k] = _nmod_poly_gcd(arg->Hmod + k * len2, a, len1, b, len2, arg->mods[k]);
}

/* Modular quotients ***********************************************************/

/*
    Computes the quotient P / G modulo a batch of primes. If Hres is
    nonNULL, this is for stored primes where the images of P and of
    H = (gamma / lc(G)) G are available (and the division is known to be
    exact); otherwise, P and G have been reduced from scratch and exactness
    is checked.
*/
typedef struct
{
    /* stored primes */
    nn_srcptr Hres;     /* stride Hstride: image of H at offset 0, of P at Poff */
    slong Hstride;
    slong Poff;
    /* new primes */
    nn_srcptr PGmod;    /* images of P (lenP), then G (lenG); stride lenP + lenG */
    nn_ptr R;           /* scratch, lenG - 1 per prime */

    nn_ptr Qmod;        /* stride lenP - lenG + 1 */
    int * status;       /* 1: ok, 0: not divisible, -1: bad prime */
    const nmod_t * mods;
    slong lenP;
    slong lenG;
}
_div_batch_arg_t;

static void
_div_batch_worker(slong k, void * varg)
{
    _div_batch_arg_t * arg = (_div_batch_arg_t *) varg;
    slong lenP = arg->lenP, lenG = arg->lenG;
    nn_ptr q = arg->Qmod + k * (lenP - lenG + 1);
    nmod_t mod = arg->mods[k];

    if (arg->Hres != NULL)
    {
        nn_srcptr h = arg->Hres + k * arg->Hstride;

        _nmod_poly_div(q, h + arg->Poff, lenP, h, lenG, mod);
        arg->status[k] = 1;
    }
    else
    {
        nn_srcptr p = arg->PGmod + k * (lenP + lenG);
        nn_srcptr g = p + lenP;
        nn_ptr r = arg->R + k * (lenG - 1);

        if (g[lenG - 1] == 0)
        {
            arg->status[k] = -1;
            return;
        }

        _nmod_poly_divrem(q, r, p, lenP, g, lenG, mod);
        arg->status[k] = _nmod_vec_is_zero(r, lenG - 1);
    }
}

/* Given that G * Q = A mod M where M >= 2^mbits, check whether
   G * Q == A. If M is large enough, this follows from bounds for the
   coefficients; otherwise, we multiply. */
static int
_check_product(const fmpz * A, slong lenA, const fmpz * G, slong lenG,
    const fmpz * Q, slong lenQ, fmpz * T);

static int
_check_product_mod(const fmpz * A, slong lenA, const fmpz * G, slong lenG,
    const fmpz * Q, slong lenQ, fmpz * T, slong mbits)
{
    slong bA, bG, bQ;

    bA = _fmpz_vec_max_bits(A, lenA);
    bG = _fmpz_vec_max_bits(G, lenG);
    bQ = _fmpz_vec_max_bits(Q, lenQ);
    bA = FLINT_ABS(bA);
    bG = FLINT_ABS(bG);
    bQ = FLINT_ABS(bQ);

    /* |coefficients of G Q| <= min(lenG, lenQ) |G| |Q|; need both sides < M / 2 */
    if (bA + 1 < mbits && bG + bQ + (slong) FLINT_BIT_COUNT(FLINT_MIN(lenG, lenQ)) + 1 < mbits)
        return 1;

    return _check_product(A, lenA, G, lenG, Q, lenQ, T);
}

/* Check whether G * Q == A */
static int
_check_product(const fmpz * A, slong lenA, const fmpz * G, slong lenG,
    const fmpz * Q, slong lenQ, fmpz * T)
{
    if (lenG >= lenQ)
        _fmpz_poly_mul(T, G, lenG, Q, lenQ);
    else
        _fmpz_poly_mul(T, Q, lenQ, G, lenG);

    return _fmpz_vec_equal(T, A, lenA);
}

/* Number of threads to use for a batch of modular gcds or divisions of
   polynomials of length about len (not worth it for short polynomials) */
static int
_num_threads_for_batch(slong batch_size, slong len)
{
    if (batch_size < 2 || len < 200)
        return 1;
    return FLINT_MIN(batch_size, flint_get_num_available_threads() + 1);
}

/*
    Certify that the candidate G (primitive, positive leading coefficient,
    length lenG > 1) divides P.

    Hres contains, for num_stored primes (with product >= 2^stored_bits),
    the images of H = (gamma / lc(G)) G (or (gamma / G(0)) G), and at
    offset Poff the images of P. The quotient P / G is computed modulo
    these primes (as needed), then modulo more primes if needed. When the
    reconstructed quotient has stabilized, it is verified by bounds or by
    multiplication.

    If the quotient is expected to need many more primes than those
    stored (P has much larger coefficients than G), or if G or the
    quotient is short, exact division over Z is cheaper.
*/
static int
_verify_divides(const fmpz * P, slong lenP, const fmpz * G, slong lenG,
    nn_srcptr Hres, slong Hstride, slong Poff, nn_srcptr Hprimes,
    slong num_stored, slong stored_bits,
    const fmpz_t gamma, int trailing, _prime_iter_t * it)
{
    slong lenQ = lenP - lenG + 1;
    slong k, num, alloc, next_crt, num_used_stored, mbits, qbits;
    nn_ptr Qres, primes, scale = NULL;
    _residue_cache_t rc;
    fmpz * Q, * T;
    _probe_t PQ;
    int result = -1;

    Q = _fmpz_vec_init(lenQ);
    T = _fmpz_vec_init(lenP);

    /* rough estimate of the size of the quotient */
    {
        slong bP = _fmpz_vec_max_bits(P, lenP);
        slong bG = _fmpz_vec_max_bits(G, lenG);
        qbits = FLINT_ABS(bP) - FLINT_ABS(bG) + FLINT_BIT_COUNT(lenP) + 1;
    }

    /* If the quotient needs many more primes than those stored (P has much
       larger coefficients than G), reducing P modulo all the extra primes
       costs more than dividing exactly over Z. */
    if (lenG <= 8 || lenQ <= 8 || qbits > 3 * stored_bits + 32 * (FLINT_BITS - 4))
    {
        /* For short Q, quotient + multiplication is fastest (note that
           _fmpz_poly_div with exact = 1 checks the coefficient divisions,
           unlike _fmpz_poly_divexact, which requires knowing that the
           division is exact); otherwise, division with remainder. */
        if (lenQ <= 8)
            result = _fmpz_poly_div(Q, P, lenP, G, lenG, 1)
                && _check_product(P, lenP, G, lenG, Q, lenQ, T);
        else
            result = _fmpz_poly_divides(Q, P, lenP, G, lenG);

        _fmpz_vec_clear(Q, lenQ);
        _fmpz_vec_clear(T, lenP);

        return result;
    }

    _probe_init(&PQ);
    _residue_cache_init(&rc, P, lenP, G, lenG, NULL, 0, it);

    /* residues of gamma and of the normalising coefficient of G modulo
       the stored primes */
    if (num_stored > 0)
    {
        fmpz * v = _fmpz_vec_init(2);
        nmod_t * smods = flint_malloc(sizeof(nmod_t) * num_stored);

        for (k = 0; k < num_stored; k++)
            nmod_init(smods + k, Hprimes[k]);

        scale = flint_malloc(sizeof(ulong) * 2 * num_stored);
        fmpz_set(v, gamma);
        fmpz_set(v + 1, G + (trailing ? 0 : lenG - 1));
        _multi_mod_batch(scale, 2, v, 2, NULL, 0, NULL, 0, NULL, 0, NULL, 0,
                         Hprimes, smods, num_stored);
        _fmpz_vec_clear(v, 2);
        flint_free(smods);
    }

    alloc = FLINT_MAX(num_stored, 4);
    Qres = flint_malloc(sizeof(ulong) * alloc * lenQ);
    primes = flint_malloc(sizeof(ulong) * alloc);

    num = 0;
    num_used_stored = 0;
    next_crt = 0;
    mbits = 0;

    while (result == -1)
    {
        if (num >= next_crt && num > 0 && _probe_stable(&PQ))
        {
            _crt_strided(Q, Qres, lenQ, lenQ, primes, num);

            if (_check_product_mod(P, lenP, G, lenG, Q, lenQ, T, mbits))
            {
                result = 1;
                break;
            }

            /* The probe was fooled (unlikely); try again later */
            next_crt = num + FLINT_MAX(1, num / 4);
        }

        /* More primes: first the stored ones, then new ones */
        {
            slong bs, first_stored = num_used_stored;
            nn_ptr R = NULL, Qmod, bprimes;
            nmod_t * mods;
            int * status;
            _div_batch_arg_t arg;
            int stored = (num_used_stored < num_stored);

            bs = FLINT_MAX(1, num / 4);
            bs = FLINT_MAX(bs, FLINT_MIN(num, flint_get_num_available_threads() + 1));

            if (stored)
                bs = FLINT_MIN(bs, num_stored - num_used_stored);

            Qmod = flint_malloc(sizeof(ulong) * bs * lenQ);
            bprimes = flint_malloc(sizeof(ulong) * bs);
            mods = flint_malloc(sizeof(nmod_t) * bs);
            status = flint_malloc(sizeof(int) * bs);

            if (stored)
            {
                for (k = 0; k < bs; k++)
                {
                    bprimes[k] = Hprimes[num_used_stored + k];
                    nmod_init(mods + k, bprimes[k]);
                }

                arg.Hres = Hres + num_used_stored * Hstride;
                arg.Hstride = Hstride;
                arg.Poff = Poff;
                arg.PGmod = NULL;
                num_used_stored += bs;
            }
            else
            {
                /* the reduction dominates for new primes; fetch in
                   doubling chunks */
                _residue_cache_require(&rc, bs, WORD_MAX);
                R = flint_malloc(sizeof(ulong) * bs * (lenG - 1));

                for (k = 0; k < bs; k++)
                {
                    bprimes[k] = rc.primes[rc.start + k];
                    mods[k] = rc.mods[rc.start + k];
                }

                arg.Hres = NULL;
                arg.PGmod = rc.res + rc.start * rc.stride;
                rc.start += bs;
            }

            arg.R = R;
            arg.Qmod = Qmod;
            arg.status = status;
            arg.mods = mods;
            arg.lenP = lenP;
            arg.lenG = lenG;

            flint_parallel_do(_div_batch_worker, &arg, bs,
                _num_threads_for_batch(bs, lenP), FLINT_PARALLEL_STRIDED);

            for (k = 0; k < bs; k++)
            {
                nn_ptr q = Qmod + k * lenQ;
                nmod_t mod = mods[k];

                if (status[k] == 0)
                {
                    result = 0;
                    break;
                }

                if (status[k] == 1)
                {
                    /* The quotient by H is (lc(G) / gamma) times the
                       quotient by G (or (G(0) / gamma) times) */
                    if (stored)
                    {
                        slong i = first_stored + k;
                        ulong c = nmod_mul(scale[2 * i], nmod_inv(scale[2 * i + 1], mod), mod);
                        _nmod_vec_scalar_mul_nmod(q, q, lenQ, c, mod);
                    }

                    if (num == alloc)
                    {
                        alloc = 2 * alloc;
                        Qres = flint_realloc(Qres, sizeof(ulong) * alloc * lenQ);
                        primes = flint_realloc(primes, sizeof(ulong) * alloc);
                    }

                    _nmod_vec_set(Qres + num * lenQ, q, lenQ);
                    primes[num] = mod.n;
                    mbits += FLINT_BIT_COUNT(mod.n) - 1;
                    _probe_add(&PQ, q, lenQ, mod);
                    num++;
                }
            }

            flint_free(R);
            flint_free(Qmod);
            flint_free(bprimes);
            flint_free(mods);
            flint_free(status);
        }
    }

    _probe_clear(&PQ);
    _residue_cache_clear(&rc);
    flint_free(scale);
    flint_free(Qres);
    flint_free(primes);
    _fmpz_vec_clear(Q, lenQ);
    _fmpz_vec_clear(T, lenP);

    return result;
}

/*
    Certify that the candidate G (primitive, positive leading coefficient,
    length lenG > 1) divides A and B. Hres contains, for num_stored primes,
    the images of H followed by the images of A and B.
*/
static int
_verify_candidate(const fmpz * G, slong lenG,
    const fmpz * A, slong len1, const fmpz * B, slong len2,
    nn_srcptr Hres, slong Hstride, nn_srcptr Hprimes, slong num_stored,
    const fmpz_t gamma, int trailing, _prime_iter_t * it)
{
    slong k, stored_bits = 0;

    /* Quick rejection */
    if (!fmpz_divisible(B, G) || !fmpz_divisible(A, G))
        return 0;

    for (k = 0; k < num_stored; k++)
        stored_bits += FLINT_BIT_COUNT(Hprimes[k]) - 1;

    /* B first: it is usually cheaper */
    return _verify_divides(B, len2, G, lenG, Hres, Hstride, lenG + len1,
                Hprimes, num_stored, stored_bits, gamma, trailing, it)
        && _verify_divides(A, len1, G, lenG, Hres, Hstride, lenG,
                Hprimes, num_stored, stored_bits, gamma, trailing, it);
}

/* Upper bound for log2 of the 2-norm of (A, len) */
static slong
_fmpz_vec_log2_norm_upper(const fmpz * A, slong len)
{
    slong bits = _fmpz_vec_max_bits(A, len);
    bits = FLINT_ABS(bits);
    return bits + (FLINT_BIT_COUNT(len) + 1) / 2;
}

/* If first_prime is nonzero, use the primes greater than first_prime
   (in increasing order) instead of the default sequence (for testing).
   Requires first_prime < 2^(FLINT_BITS - 1). */
void _fmpz_poly_gcd_modular_primes(fmpz * res, const fmpz * poly1, slong len1,
                    const fmpz * poly2, slong len2, ulong first_prime)
{
    fmpz_t ac, bc, d, gamma;
    const fmpz * A, * B;
    fmpz * Acopy = NULL, * Bcopy = NULL;
    fmpz * Gc;              /* candidate gcd */
    slong k;
    slong hlen_cur;         /* length of the candidate gcd images */
    slong num_stored, alloc_stored, Hstride;
    nn_ptr Hres;            /* stored images of H, A, B */
    nn_ptr primes;          /* primes of the stored images */
    slong batch_size, next_attempt, bound_bits, total_bits;
    int done = 0, trailing;
    nn_ptr Hmod;
    _residue_cache_t rc;
    slong * hlens;
    slong res_len = 0;
    _prime_iter_t it;
    _probe_t P;

    fmpz_init(ac);
    fmpz_init(bc);
    fmpz_init(d);

    /* compute gcd of content of poly1 and poly2 */
    _fmpz_vec_content(ac, poly1, len1);
    _fmpz_vec_content(bc, poly2, len2);
    fmpz_gcd(d, ac, bc);

    /* special case, one of the polys is a constant */
    if (len2 == 1) /* if len1 == 1 then so does len2 */
    {
        fmpz_set(res, d);

        fmpz_clear(ac);
        fmpz_clear(bc);
        fmpz_clear(d);
        return;
    }

    /* We work with A and B primitive and ensure that B has positive
       leading coefficient, so that B itself is a candidate gcd. The input
       is only copied if needed. (The output is written only at the end,
       so aliasing is allowed.) */
    if (fmpz_is_one(ac))
    {
        A = poly1;
    }
    else
    {
        Acopy = _fmpz_vec_init(len1);
        _fmpz_vec_scalar_divexact_fmpz(Acopy, poly1, len1, ac);
        A = Acopy;
    }

    if (fmpz_sgn(poly2 + len2 - 1) < 0)
        fmpz_neg(bc, bc);

    if (fmpz_is_one(bc))
    {
        B = poly2;
    }
    else
    {
        Bcopy = _fmpz_vec_init(len2);
        _fmpz_vec_scalar_divexact_fmpz(Bcopy, poly2, len2, bc);
        B = Bcopy;
    }

    fmpz_clear(ac);
    fmpz_clear(bc);

    Gc = _fmpz_vec_init(len2);

    /* The images of the gcd modulo p are normalised to have leading
       coefficient gamma = gcd(lc(A), lc(B)), making them images of
       H = (gamma / lc(G)) G. If gcd(A(0), B(0)) is smaller, we normalise
       to have that constant coefficient instead (then H = (gamma / G(0)) G).
       This avoids coefficient growth if A and B have a large common
       factor in their leading coefficients. */
    fmpz_init(gamma);
    fmpz_gcd(gamma, A + len1 - 1, B + len2 - 1);
    trailing = 0;

    if (!fmpz_is_one(gamma) && !fmpz_is_zero(A) && !fmpz_is_zero(B))
    {
        fmpz_t gamma0;
        fmpz_init(gamma0);
        fmpz_gcd(gamma0, A, B);
        if (fmpz_bits(gamma0) < fmpz_bits(gamma))
        {
            fmpz_swap(gamma, gamma0);
            trailing = 1;
        }
        fmpz_clear(gamma0);
    }

    /* Mignotte bound for the height of H = (gamma / lc(G)) G:
       ||H||_1 <= 2^deg(G) (|gamma| / |lc(A)|) ||A||_2, and similarly with B
       (and for the reversed polynomials in the case of trailing
       normalisation). A modulus larger than twice this allows correct
       reconstruction if the degree is correct. We only use this to limit the
       batch sizes and force a reconstruction attempt. */
    {
        slong ba, bb;
        ba = _fmpz_vec_log2_norm_upper(A, len1) - (fmpz_bits(A + (trailing ? 0 : len1 - 1)) - 1);
        bb = _fmpz_vec_log2_norm_upper(B, len2) - (fmpz_bits(B + (trailing ? 0 : len2 - 1)) - 1);
        bound_bits = (len2 - 1) + FLINT_MIN(ba, bb) + fmpz_bits(gamma) + 2;
    }

    _probe_init(&P);

    Hstride = len1 + 2 * len2;
    alloc_stored = 0;
    num_stored = 0;
    Hres = NULL;
    primes = NULL;
    hlen_cur = len2;
    next_attempt = 1;
    total_bits = 0;

    Hmod = NULL;
    hlens = NULL;

    /* The first prime serves mainly to detect a trivial gcd (the common
       case), for which a small prime is the cheapest choice: the modular
       gcd is up to twice as fast with a 20-bit prime as with a 62-bit
       prime. The probability that a random 20-bit prime is unlucky
       is still negligible for typical input. If the coefficients are
       small, however, a nontrivial gcd may be recovered from a single
       large prime, and a small first prime would likely be wasted. */
    {
        slong ba = _fmpz_vec_max_bits(A, len1);
        slong bb = _fmpz_vec_max_bits(B, len2);
        int small_first = (FLINT_ABS(ba) + FLINT_ABS(bb) >= 2 * FLINT_BITS);
        _prime_iter_init(&it, len2, first_prime, small_first);
    }

    _residue_cache_init(&rc, A, len1, B, len2, gamma, 1, &it);

    /* The first prime alone: frequently decides that the gcd is trivial
       or that B divides A. */
    batch_size = 1;

    while (!done)
    {
        _gcd_batch_arg_t arg;
        nn_srcptr ABmod;
        const nmod_t * mods;

        Hmod = flint_realloc(Hmod, sizeof(ulong) * batch_size * len2);
        hlens = flint_realloc(hlens, sizeof(slong) * batch_size);

        /* residues of A, B, gamma; fetch at most up to the Mignotte bound */
        {
            slong max_fetch = WORD_MAX;
            if (total_bits <= bound_bits)
                max_fetch = (bound_bits - total_bits) / (FLINT_BITS - 14) + 1;
            _residue_cache_require(&rc, batch_size, max_fetch);
        }

        ABmod = rc.res + rc.start * rc.stride;
        mods = rc.mods + rc.start;
        rc.start += batch_size;

        arg.ABmod = ABmod;
        arg.ABstride = rc.stride;
        arg.Hmod = Hmod;
        arg.hlen = hlens;
        arg.mods = mods;
        arg.len1 = len1;
        arg.len2 = len2;

        flint_parallel_do(_gcd_batch_worker, &arg, batch_size,
            _num_threads_for_batch(batch_size, len2), FLINT_PARALLEL_STRIDED);

        for (k = 0; k < batch_size; k++)
        {
            slong hlen = hlens[k];
            nn_ptr h = Hmod + k * len2;
            nmod_t mod = mods[k];
            nn_srcptr ab = ABmod + k * rc.stride;
            ulong c;

            if (hlen < 0)   /* bad prime */
                continue;

            if (hlen == 1)  /* coprime */
            {
                fmpz_one(Gc);
                res_len = 1;
                done = 1;
                break;
            }

            if (hlen > hlen_cur)    /* unlucky prime */
                continue;

            if (hlen < hlen_cur)    /* all previous primes were unlucky */
            {
                hlen_cur = hlen;
                num_stored = 0;
                total_bits = 0;
                _probe_reset(&P);
                next_attempt = 1;
            }

            /* Scale to leading (or constant) coefficient gamma */
            {
                ulong g_mod = ab[len1 + len2];

                /* With trailing normalisation, primes dividing G(0) or
                   gamma are useless. (With leading normalisation, both are
                   nonzero since the prime does not divide lc(A) lc(B).) */
                c = h[trailing ? 0 : hlen - 1];
                if (c == 0 || g_mod == 0)
                    continue;

                c = nmod_mul(g_mod, nmod_inv(c, mod), mod);
            }
            _nmod_vec_scalar_mul_nmod(h, h, hlen, c, mod);

            if (num_stored == alloc_stored)
            {
                alloc_stored = FLINT_MAX(2 * alloc_stored, 4);
                Hres = flint_realloc(Hres, sizeof(ulong) * alloc_stored * Hstride);
                primes = flint_realloc(primes, sizeof(ulong) * alloc_stored);
            }

            /* Store the images of H, A, B */
            _nmod_vec_set(Hres + num_stored * Hstride, h, hlen);
            _nmod_vec_set(Hres + num_stored * Hstride + hlen, ab, len1 + len2);
            primes[num_stored] = mod.n;
            num_stored++;
            total_bits += FLINT_BIT_COUNT(mod.n) - 1;

            if (hlen == len2)
            {
                /* Degree of the image equals deg(B): either B | A, or the
                   prime is unlucky. */
                if (_verify_candidate(B, len2, A, len1, B, len2,
                        Hres, Hstride, primes, num_stored, gamma, trailing, &it))
                {
                    _fmpz_vec_set(Gc, B, len2);
                    res_len = len2;
                    done = 1;
                    break;
                }

                hlen_cur = len2 - 1;
                num_stored = 0;
                total_bits = 0;
                _probe_reset(&P);
                next_attempt = 1;
                continue;
            }

            _probe_add(&P, h, hlen, mod);
        }

        if (done)
            break;

        if (num_stored >= next_attempt && (_probe_stable(&P) || total_bits > bound_bits
                || (hlen_cur <= SHORT_GCD_LEN && num_stored <= 2
                    && _probe_stable_margin(&P, SHORT_GCD_PROBE_MARGIN))))
        {
            fmpz_t cont;
            int ok;

            _crt_strided(Gc, Hres, Hstride, hlen_cur, primes, num_stored);

            fmpz_init(cont);
            _fmpz_vec_content(cont, Gc, hlen_cur);
            if (fmpz_sgn(Gc + hlen_cur - 1) < 0)
                fmpz_neg(cont, cont);
            if (!fmpz_is_one(cont))
                _fmpz_vec_scalar_divexact_fmpz(Gc, Gc, hlen_cur, cont);
            fmpz_clear(cont);

            ok = _verify_candidate(Gc, hlen_cur, A, len1, B, len2,
                    Hres, Hstride, primes, num_stored, gamma, trailing, &it);

            if (ok)
            {
                res_len = hlen_cur;
                break;
            }

            /* Either the probe was fooled, or all primes so far were
               unlucky; require more primes before the next attempt. */
            next_attempt = num_stored + FLINT_MAX(1, num_stored / 4);
        }

        /* Choose the next batch size: grow geometrically, so that at most
           a constant fraction of the work is wasted, but don't overshoot
           the Mignotte bound by much. With multiple threads, use at least
           one prime per thread. */
        batch_size = FLINT_MAX(1, num_stored / 4);
        batch_size = FLINT_MAX(batch_size, FLINT_MIN(num_stored, flint_get_num_available_threads() + 1));

        if (total_bits <= bound_bits)
        {
            slong needed = (bound_bits - total_bits) / (FLINT_BITS - 14) + 1;
            batch_size = FLINT_MIN(batch_size, needed);
        }
    }

    /* finally multiply by content */
    if (fmpz_is_one(d))
        _fmpz_vec_swap(res, Gc, res_len);
    else
        _fmpz_vec_scalar_mul_fmpz(res, Gc, res_len, d);
    _fmpz_vec_zero(res + res_len, len2 - res_len);

    flint_free(Hmod);
    flint_free(hlens);
    _residue_cache_clear(&rc);
    flint_free(Hres);
    flint_free(primes);

    _probe_clear(&P);
    fmpz_clear(gamma);
    fmpz_clear(d);
    if (Acopy != NULL)
        _fmpz_vec_clear(Acopy, len1);
    if (Bcopy != NULL)
        _fmpz_vec_clear(Bcopy, len2);
    _fmpz_vec_clear(Gc, len2);
}

void _fmpz_poly_gcd_modular(fmpz * res, const fmpz * poly1, slong len1,
                                        const fmpz * poly2, slong len2)
{
    _fmpz_poly_gcd_modular_primes(res, poly1, len1, poly2, len2, 0);
}

void
fmpz_poly_gcd_modular(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    if (poly1->length < poly2->length)
    {
        fmpz_poly_gcd_modular(res, poly2, poly1);
    }
    else /* len1 >= len2 >= 0 */
    {
        const slong len1 = poly1->length;
        const slong len2 = poly2->length;

        if (len1 == 0) /* len1 = len2 = 0 */
        {
            fmpz_poly_zero(res);
        }
        else if (len2 == 0) /* len1 > len2 = 0 */
        {
            if (fmpz_sgn(poly1->coeffs + (len1 - 1)) > 0)
                fmpz_poly_set(res, poly1);
            else
                fmpz_poly_neg(res, poly1);
        }
        else /* len1 >= len2 >= 1 */
        {
            /* underscore function automatically aliases */
            fmpz_poly_fit_length(res, len2);

            _fmpz_poly_gcd_modular(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);

            _fmpz_poly_set_length(res, len2);
            _fmpz_poly_normalise(res);
        }
    }
}
