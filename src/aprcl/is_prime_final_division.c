/*
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2026 FLINT contributors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#if FLINT_USES_PTHREAD
# include <pthread.h>
#endif
#include "mpn_extras.h"
#include "fmpz.h"
#include "aprcl.h"

/* is {x, xn} divisible by {d, dn}? (both normalised, dn >= 1) */
static int
_aprcl_mpn_divisible(nn_srcptr x, slong xn, nn_srcptr d, slong dn)
{
    __mpz_struct xz, dz;

    xz._mp_d = (nn_ptr) x;
    xz._mp_size = xz._mp_alloc = xn;
    dz._mp_d = (nn_ptr) d;
    dz._mp_size = dz._mp_alloc = dn;

    return mpz_divisible_p(&xz, &dz);
}

/*
    Final trial division step of APRCL (step 7 of algorithm 9.1.28 in
    Cohen's "A Course in Computational Algebraic Number Theory").

    Given n coprime to s with s^2 > n (the walk itself is correct for any
    s, but the conclusion that n is prime needs s^2 > n), the earlier steps
    have shown that
    every divisor of n is congruent to n^i mod s for some 0 <= i < r.
    If n is composite it has a prime divisor a with 1 < a <= sqrt(n) < s,
    so a is exactly the least residue of n^i mod s for some 1 <= i < r.
    We therefore walk the residues n^i mod s, i = 1, 2, ..., and test
    divisibility of n by those residues a with 1 < a <= sqrt(n) (residues
    that are even, for odd n, are skipped too). Since n^r = 1 mod s the
    walk stops when the residue 1 is reached, or after r - 1 steps.

    The walk uses mpn arithmetic with a preconditioned multiplication by
    the fixed multiplier n mod s. The range of exponents is split into
    chunks which are processed in parallel; each chunk computes its first
    power with a modular exponentiation. If a chunk reaches the residue 1
    at n^i then the residues that follow are n^1, n^2, ..., which are
    covered by the first chunk, so the chunk may stop early.
*/

/* Reference implementation with fmpz arithmetic (used for tiny s) */
static int
_aprcl_is_prime_final_division_fmpz(const fmpz_t n, const fmpz_t s, ulong r,
                                                            const fmpz_t bound)
{
    int result = 1, n_odd = fmpz_is_odd(n);
    ulong i;
    fmpz_t npow, nmul, rem;

    fmpz_init(rem);
    fmpz_init(nmul);
    fmpz_init(npow);
    fmpz_mod(nmul, n, s);       /* nmul = n mod s */
    fmpz_set(npow, nmul);       /* npow = n^i mod s */

    for (i = 1; i < r; i++)
    {
        if (fmpz_is_one(npow))
            break;

        if (!(n_odd && fmpz_is_even(npow)) && fmpz_cmp(npow, bound) <= 0)
        {
            fmpz_mod(rem, n, npow);
            if (fmpz_is_zero(rem) && !fmpz_is_one(npow))
            {
                result = 0;
                break;
            }
        }

        fmpz_mul(npow, npow, nmul);
        fmpz_mod(npow, npow, s);
    }

    fmpz_clear(npow);
    fmpz_clear(nmul);
    fmpz_clear(rem);

    return result;
}

typedef struct
{
    nn_srcptr n;        /* n, nn limbs */
    slong nn;
    nn_srcptr s;        /* s, sn limbs */
    nn_srcptr snormed;  /* s shifted left by norm bits */
    nn_srcptr sinv;     /* precomputed inverse of snormed */
    slong sn;
    ulong norm;
    nn_srcptr nmul;     /* n mod s, unshifted */
    nn_srcptr nmul_pre; /* precomputed data for multiplication by nmul, or NULL */
    int precond;        /* MPN_MULMOD_PRECOND_{NONE,SHOUP,MATRIX} */
    nn_srcptr bound;    /* floor(sqrt(n)) padded to sn limbs */
    int n_odd;
    ulong r;
    ulong chunk;        /* exponents per chunk */
    slong num_chunks;
    int found_divisor;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
} _aprcl_fd_args;

#define FD_CHECK_INTERVAL 1024

static void
_aprcl_fd_worker(slong c, void * varg)
{
    _aprcl_fd_args * args = (_aprcl_fd_args *) varg;
    slong sn = args->sn, nn = args->nn;
    ulong norm = args->norm;
    ulong i, i0, i1;
    /* the preconditioned methods work on unshifted operands, mulmod_preinvn on shifted ones */
    int unshifted = (args->precond != MPN_MULMOD_PRECOND_NONE);
    int stop = 0;
    nn_ptr npow, t, nmul_sh;
    TMP_INIT;

    /* this chunk covers exponents i0 <= i < i1 */
    i0 = 1 + c * args->chunk;
    i1 = FLINT_MIN(i0 + args->chunk, args->r);
    if (i0 >= i1)
        return;

    TMP_START;
    npow = TMP_ALLOC(3 * sn * sizeof(ulong));
    t = npow + sn;
    nmul_sh = t + sn;

    /* npow = n^{i0} mod s */
    {
        fmpz_t f, e, m;
        fmpz_init(f);
        fmpz_init(m);
        fmpz_init_set_ui(e, i0);
        fmpz_set_ui_array(f, args->n, nn);
        fmpz_set_ui_array(m, args->s, sn);
        fmpz_powm(f, f, e, m);
        fmpz_get_ui_array(npow, sn, f);
        fmpz_clear(f);
        fmpz_clear(e);
        fmpz_clear(m);
    }

    if (!unshifted)
    {
        /* mulmod_preinvn works with operands shifted left by norm bits */
        if (norm)
        {
            mpn_lshift(nmul_sh, args->nmul, sn, norm);
            mpn_lshift(npow, npow, sn, norm);
        }
        else
            flint_mpn_copyi(nmul_sh, args->nmul, sn);
    }

    for (i = i0; i < i1 && !stop; i++)
    {
        nn_srcptr a;
        slong an;

        /* a = unshifted residue n^i mod s */
        if (unshifted || norm == 0)
            a = npow;
        else
        {
            mpn_rshift(t, npow, sn, norm);
            a = t;
        }

        an = sn;
        MPN_NORM(a, an);

        /* residue 1: the powers cycle from here on */
        if (an == 1 && a[0] == 1)
            break;

        /* candidate divisor: 1 < a <= sqrt(n), and a odd if n is odd */
        if (an > 0 && !(args->n_odd && (a[0] & 1) == 0)
                   && mpn_cmp(a, args->bound, sn) <= 0
                   && _aprcl_mpn_divisible(args->n, nn, a, an))
        {
#if FLINT_USES_PTHREAD
            pthread_mutex_lock(&args->mutex);
#endif
            args->found_divisor = 1;
#if FLINT_USES_PTHREAD
            pthread_mutex_unlock(&args->mutex);
#endif
            break;
        }

        /* npow = npow * nmul mod s */
        if (args->precond == MPN_MULMOD_PRECOND_MATRIX)
            flint_mpn_mulmod_precond_matrix(npow, args->nmul_pre, npow, sn,
                                            args->snormed, args->sinv, norm);
        else if (args->precond == MPN_MULMOD_PRECOND_SHOUP)
            flint_mpn_mulmod_precond_shoup(npow, args->nmul, args->nmul_pre,
                                                        npow, sn, args->s, norm);
        else
            flint_mpn_mulmod_preinvn(npow, npow, nmul_sh, sn,
                                            args->snormed, args->sinv, norm);

        /* periodically check whether another chunk found a divisor */
        if (args->num_chunks > 1 && (i & (FD_CHECK_INTERVAL - 1)) == 0)
        {
#if FLINT_USES_PTHREAD
            pthread_mutex_lock(&args->mutex);
#endif
            stop = args->found_divisor;
#if FLINT_USES_PTHREAD
            pthread_mutex_unlock(&args->mutex);
#endif
        }
    }

    TMP_END;
}

int
aprcl_is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r)
{
    int result;
    slong sn, nn, num_threads;
    ulong norm;
    fmpz_t bound, nmul;
    nn_ptr nlimbs, slimbs, snormed, sinv, nmul_limbs, nmul_pre, bound_limbs;
    _aprcl_fd_args args;

    if (fmpz_cmp_ui(n, 3) <= 0 || r <= 1)
        return 1;

    /*
        bound = min(floor(sqrt(n)), s - 1): only residues up to sqrt(n)
        need testing, and every residue is below s (so the bound fits in
        the limbs of s even if s^2 < n).
    */
    fmpz_init(bound);
    fmpz_sqrt(bound, n);
    if (fmpz_cmp(bound, s) >= 0)
    {
        fmpz_set(bound, s);
        fmpz_sub_ui(bound, bound, 1);
    }

    sn = fmpz_size(s);

    if (sn < 2 || fmpz_sgn(s) < 0)
    {
        result = _aprcl_is_prime_final_division_fmpz(n, s, r, bound);
        fmpz_clear(bound);
        return result;
    }

    nn = fmpz_size(n);

    nlimbs = flint_malloc((nn + 6 * sn) * sizeof(ulong));
    slimbs = nlimbs + nn;
    snormed = slimbs + sn;
    sinv = snormed + sn;
    nmul_limbs = sinv + sn;
    nmul_pre = nmul_limbs + sn;
    bound_limbs = nmul_pre + sn;

    fmpz_get_ui_array(nlimbs, nn, n);
    fmpz_get_ui_array(slimbs, sn, s);
    fmpz_get_ui_array(bound_limbs, sn, bound); /* bound < s fits in sn limbs */

    norm = flint_clz(slimbs[sn - 1]);
    if (norm)
        mpn_lshift(snormed, slimbs, sn, norm);
    else
        flint_mpn_copyi(snormed, slimbs, sn);
    flint_mpn_preinvn(sinv, snormed, sn);

    fmpz_init(nmul);
    fmpz_mod(nmul, n, s);
    fmpz_get_ui_array(nmul_limbs, sn, nmul);
    fmpz_clear(nmul);

    args.n = nlimbs;
    args.nn = nn;
    args.s = slimbs;
    args.snormed = snormed;
    args.sinv = sinv;
    args.sn = sn;
    args.norm = norm;
    args.nmul = nmul_limbs;
    args.bound = bound_limbs;
    args.n_odd = fmpz_is_odd(n);
    args.r = r;
    args.found_divisor = 0;

    /* preconditioned multiplication by the fixed multiplier nmul */
    args.precond = flint_mpn_mulmod_want_precond(sn, r, norm);
    args.nmul_pre = NULL;
    if (args.precond == MPN_MULMOD_PRECOND_MATRIX)
    {
        nmul_pre = flint_malloc(flint_mpn_mulmod_precond_matrix_alloc(sn) * sizeof(ulong));
        flint_mpn_mulmod_precond_matrix_precompute(nmul_pre, nmul_limbs, sn,
                                                            snormed, sinv, norm);
        args.nmul_pre = nmul_pre;
    }
    else if (args.precond == MPN_MULMOD_PRECOND_SHOUP)
    {
        flint_mpn_mulmod_precond_shoup_precompute(nmul_pre, nmul_limbs, sn,
                                                            snormed, sinv, norm);
        args.nmul_pre = nmul_pre;
    }

    /* split the exponents 1 <= i < r into chunks, one per thread */
    num_threads = flint_get_num_threads();
    if (r < 4096)
        num_threads = 1;
    args.num_chunks = FLINT_MIN(num_threads, (slong) (r - 1));
    args.chunk = (r - 1 + args.num_chunks - 1) / args.num_chunks;

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&args.mutex, NULL);
#endif

    flint_parallel_do(_aprcl_fd_worker, &args, args.num_chunks, 0, 0);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&args.mutex);
#endif

    result = !args.found_divisor;

    if (args.precond == MPN_MULMOD_PRECOND_MATRIX)
        flint_free(nmul_pre);
    flint_free(nlimbs);
    fmpz_clear(bound);

    return result;
}
