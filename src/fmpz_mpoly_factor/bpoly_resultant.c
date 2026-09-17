/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "ulong_extras.h"
#include "n_poly.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "thread_support.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

/* Bits of the primes used for the modular images. Staying below
   FLINT_BITS - 1 keeps NMOD_CAN_USE_SHOUP true for the images. */
#define MODULAR_PRIME_BITS (FLINT_BITS - 2)

/* The images may instead be taken at primes p = m 2^k + 1, at which the
   multipoint algorithm evaluates with a DFT rather than with the Bluestein
   products of the geometric method. The transforms of fft_small work with
   doubles and accept no modulus of more than this many bits, so such a prime
   carries fewer bits and more of them are needed; the trade is worth making
   only under the conditions in _use_fft_primes below. Their product is
   bounded below using one bit less, since the primes are only known to be
   below 2^FFT_PRIME_BITS and above half of it. */
#define FFT_PRIME_BITS 50

/* Beyond this exponent, too few m remain below 2^FFT_PRIME_BITS for the
   primes to be found. */
#define FFT_PRIME_MAX_DEPTH 32

/* Number of points from which the DFT evaluation gains more than the extra
   primes cost. */
#define MODULAR_FFT_PRIME_POINTS 2048

/* When the reconstruction is not required to be proved, it stops once it has
   been left unchanged by primes whose product carries this many bits. */
#define MODULAR_STABLE_BITS 100

/* Bound for the coefficients of res_y(A, B), where A and B are seen as
   polynomials in y with coefficients in Z[x]:

      ||res||_oo <= (sum_i ||A_i||_1^2)^(n/2) (sum_j ||B_j||_1^2)^(m/2)

   with m = deg_y(A) and n = deg_y(B). The coefficients of the resultant are
   bounded by its maximum modulus on the unit circle |x| = 1, where the
   univariate bound |res(F, G)| <= ||F||_2^deg(G) ||G||_2^deg(F) applies with
   ||A(x, .)||_2^2 = sum_i |A_i(x)|^2 <= sum_i ||A_i||_1^2. */
static void
_bivariate_resultant_bound(fmpz_t bound,
                           const fmpz_poly_struct * A, slong lenA,
                           const fmpz_poly_struct * B, slong lenB)
{
    fmpz_t sa, sb, t, u;
    slong i, k;

    fmpz_init(sa);
    fmpz_init(sb);
    fmpz_init(t);
    fmpz_init(u);

    for (i = 0; i < lenA; i++)
    {
        fmpz_zero(t);
        for (k = 0; k < A[i].length; k++)
        {
            fmpz_abs(u, A[i].coeffs + k);
            fmpz_add(t, t, u);
        }
        fmpz_addmul(sa, t, t);
    }

    for (i = 0; i < lenB; i++)
    {
        fmpz_zero(t);
        for (k = 0; k < B[i].length; k++)
        {
            fmpz_abs(u, B[i].coeffs + k);
            fmpz_add(t, t, u);
        }
        fmpz_addmul(sb, t, t);
    }

    fmpz_pow_ui(sa, sa, lenB - 1);
    fmpz_pow_ui(sb, sb, lenA - 1);
    fmpz_mul(bound, sa, sb);
    fmpz_sqrt(bound, bound);
    fmpz_add_ui(bound, bound, 1);

    fmpz_clear(sa);
    fmpz_clear(sb);
    fmpz_clear(t);
    fmpz_clear(u);
}

/* res_y(A, B) mod p, written to (out, outlen) with zero padding. A and B must
   have nonzero leading coefficients mod p, so that the degrees in y are
   preserved and the reduction commutes with the resultant. The n_poly
   structures are supplied by the caller and reused across primes. */
static int
_bivariate_resultant_nmod(nn_ptr out, slong outlen,
                          const fmpz_poly_struct * A, slong lenA,
                          const fmpz_poly_struct * B, slong lenB, ulong p,
                          n_poly_struct * Ap, n_poly_struct * Bp, n_poly_t rp)
{
    nmod_t mod;
    slong i, k, len;

    nmod_init(&mod, p);

    for (i = 0; i < lenA; i++)
    {
        for (k = 0; k < A[i].length; k++)
            Ap[i].coeffs[k] = fmpz_fdiv_ui(A[i].coeffs + k, p);
        len = A[i].length;
        while (len > 0 && Ap[i].coeffs[len - 1] == 0)
            len--;
        Ap[i].length = len;
    }

    for (i = 0; i < lenB; i++)
    {
        for (k = 0; k < B[i].length; k++)
            Bp[i].coeffs[k] = fmpz_fdiv_ui(B[i].coeffs + k, p);
        len = B[i].length;
        while (len > 0 && Bp[i].coeffs[len - 1] == 0)
            len--;
        Bp[i].length = len;
    }

    if (!_n_bpoly_mod_resultant(rp, Ap, lenA, Bp, lenB, mod))
        return 0;

    /* outlen comes from the degree bound on the resultant, so the image
       cannot be longer; truncating silently would hide a wrong bound */
    FLINT_ASSERT(rp->length <= outlen);

    for (k = 0; k < outlen; k++)
        out[k] = (k < rp->length) ? rp->coeffs[k] : 0;

    return 1;
}

/* The images at the several primes are independent, and are split over the
   thread pool. Each worker reduces the inputs and runs the multipoint
   algorithm in its own buffers; the multipoint algorithm is itself threaded,
   but the pool hands out each thread once, so an inner call simply finds
   none left and runs serially rather than oversubscribing. */
typedef struct
{
    const fmpz_poly_struct * A;
    const fmpz_poly_struct * B;
    slong lenA;
    slong lenB;
    slong outlen;
    slong lo;
    slong hi;
    slong nworkers;
    nn_srcptr primes;
    nn_ptr residues;
    n_poly_struct ** Ap;
    n_poly_struct ** Bp;
    n_poly_struct * rp;
    int * status;
}
_modular_args_t;

static void
_modular_worker(slong widx, void * argsv)
{
    _modular_args_t * W = argsv;
    slong i;

    for (i = W->lo + widx; i < W->hi; i += W->nworkers)
    {
        if (!_bivariate_resultant_nmod(W->residues + i * W->outlen,
                W->outlen, W->A, W->lenA, W->B, W->lenB, W->primes[i],
                W->Ap[widx], W->Bp[widx], W->rp + widx))
            W->status[widx] = 0;
    }
}

/* Lists num primes not dividing l, which the images are then taken at. With
   fft set they are the largest primes p = m 2^depth + 1 below 2^FFT_PRIME_BITS,
   at which the multipoint algorithm can evaluate with a transform of depth up
   to depth; otherwise they are the primes just above 2^MODULAR_PRIME_BITS.
   Returns the number listed, which is less than num only in the first case,
   when the exponent leaves too few candidates. */
static slong
_modular_primes(nn_ptr primes, slong num, int fft, flint_bitcnt_t depth,
                const fmpz_t l)
{
    slong n = 0;

    if (fft)
    {
#if FLINT_HAVE_FFT_SMALL
        ulong m;

        /* stopping halfway keeps every prime above 2^(FFT_PRIME_BITS - 1) */
        for (m = (UWORD(1) << (FFT_PRIME_BITS - depth)) - 1;
             m > (UWORD(1) << (FFT_PRIME_BITS - 1 - depth)) && n < num; m--)
        {
            ulong p = (m << depth) + 1;

            if (!n_is_prime(p))
                continue;
            if (!fft_small_prime_supports_depth(p, depth))
                continue;
            if (fmpz_fdiv_ui(l, p) == 0)
                continue;

            primes[n++] = p;
        }
#endif
    }
    else
    {
        ulong p = UWORD(1) << MODULAR_PRIME_BITS;

        while (n < num)
        {
            p = n_nextprime(p, 0);

            if (fmpz_fdiv_ui(l, p) == 0)
                continue;

            primes[n++] = p;
        }
    }

    return n;
}

/* Whether to spend the extra primes on a modulus admitting a DFT. The images
   have to be taken by the multipoint algorithm for a DFT to be reachable at
   all, and the evaluation has to be a large enough share of one: it costs
   about (lenA + lenB) npoints log(npoints) word operations against the
   lenA lenB npoints of the resultants at the points, and only the first of
   the two is made faster. */
static int
_use_fft_primes(slong lenA, slong lenB, slong npoints, flint_bitcnt_t depth)
{
#if FLINT_HAVE_FFT_SMALL
    if (depth > FFT_PRIME_MAX_DEPTH)
        return 0;

    if (!n_bpoly_mod_resultant_cutoff(lenA, lenB, npoints))
        return 0;

    return npoints >= MODULAR_FFT_PRIME_POINTS;
#else
    return 0;
#endif
}

/* CRT the residues of `len` coefficients, laid out prime-major with the given
   stride, over the first `num` primes. */
static void
_crt_reconstruct(fmpz * out, slong len, nn_srcptr residues, slong stride,
                 nn_srcptr primes, slong num)
{
    fmpz_comb_t comb;
    fmpz_comb_temp_t temp;
    nn_ptr r;
    slong i, k;

    r = flint_malloc(num * sizeof(ulong));

    fmpz_comb_init(comb, primes, num);
    fmpz_comb_temp_init(temp, comb);

    for (k = 0; k < len; k++)
    {
        for (i = 0; i < num; i++)
            r[i] = residues[i * stride + k];

        fmpz_multi_CRT_ui(out + k, r, comb, temp, 1);
    }

    fmpz_comb_temp_clear(temp);
    fmpz_comb_clear(comb);
    flint_free(r);
}

/* res_y(A, B) over Z, where A and B are seen as polynomials in y of lengths
   lenA and lenB with coefficients in Z[x]. A and B are working copies, which
   are modified. See fmpz_mpoly_factor.h for what proved means. */
int
_fmpz_bpoly_resultant(fmpz_poly_t res,
                      fmpz_poly_struct * A, slong lenA,
                      fmpz_poly_struct * B, slong lenB, int proved)
{
    fmpz_poly_t pa, pb, t;
    fmpz_t ca, cb, l, bound, u;
    nn_ptr residues, primes, Abuf, Bbuf;
    n_poly_struct ** Ap, ** Bp;
    n_poly_struct * Apbuf, * Bpbuf, * rp;
    fmpz * cand, * prev;
    slong totA, totB, totAall, totBall;
    slong i, j, k, dxA, dxB, outlen, num, nprimes, nworkers, checkpoint;
    flint_bitcnt_t bound_bits, depth, pbits, stable_bits;
    int fft;
    int status = 1;
    int * wstatus;
    _modular_args_t wargs;

    fmpz_init(ca);
    fmpz_init(cb);
    fmpz_init(l);
    fmpz_init(u);
    fmpz_poly_init(pa);
    fmpz_poly_init(pb);
    fmpz_poly_init(t);

    /* integer contents */
    for (i = 0; i < lenA; i++)
        _fmpz_vec_content_chained(ca, A[i].coeffs, A[i].length, ca);
    if (!fmpz_is_one(ca))
        for (i = 0; i < lenA; i++)
            fmpz_poly_scalar_divexact_fmpz(A + i, A + i, ca);

    for (i = 0; i < lenB; i++)
        _fmpz_vec_content_chained(cb, B[i].coeffs, B[i].length, cb);
    if (!fmpz_is_one(cb))
        for (i = 0; i < lenB; i++)
            fmpz_poly_scalar_divexact_fmpz(B + i, B + i, cb);

    /* contents in x; res(c(x) A, B) = c(x)^deg_y(B) res(A, B) */
    for (i = 0; i < lenA; i++)
    {
        fmpz_poly_gcd(pa, pa, A + i);
        if (pa->length == 1)
            break;
    }
    if (pa->length > 1)
        for (i = 0; i < lenA; i++)
            fmpz_poly_divexact(A + i, A + i, pa);

    for (i = 0; i < lenB; i++)
    {
        fmpz_poly_gcd(pb, pb, B + i);
        if (pb->length == 1)
            break;
    }
    if (pb->length > 1)
        for (i = 0; i < lenB; i++)
            fmpz_poly_divexact(B + i, B + i, pb);

    /* a prime must not divide the leading coefficient in y of either input,
       or the degree in y would drop in the modular image */
    _fmpz_vec_content(l, A[lenA - 1].coeffs, A[lenA - 1].length);
    _fmpz_vec_content(u, B[lenB - 1].coeffs, B[lenB - 1].length);
    fmpz_mul(l, l, u);

    dxA = 0;
    for (i = 0; i < lenA; i++)
        dxA = FLINT_MAX(dxA, A[i].length);
    dxB = 0;
    for (i = 0; i < lenB; i++)
        dxB = FLINT_MAX(dxB, B[i].length);

    /* deg_x(res) <= deg_y(B) deg_x(A) + deg_y(A) deg_x(B) */
    outlen = (lenB - 1) * (dxA - 1) + (lenA - 1) * (dxB - 1) + 1;

    fmpz_init(bound);
    _bivariate_resultant_bound(bound, A, lenA, B, lenB);
    bound_bits = fmpz_bits(bound) + 2;

    /* modular images, reused across primes */
    totA = 0;
    for (i = 0; i < lenA; i++)
        totA += A[i].length;
    totB = 0;
    for (i = 0; i < lenB; i++)
        totB += B[i].length;

    totAall = totA;
    totBall = totB;

    /* The number of primes that reaching the bound takes is known in advance,
       and so is the list itself; the images at them are independent and are
       computed in parallel. */
    depth = FLINT_BIT_COUNT(outlen - 1) + 1;
    fft = _use_fft_primes(lenA, lenB, outlen, depth);
    pbits = fft ? FFT_PRIME_BITS - 1 : MODULAR_PRIME_BITS;

    nprimes = (bound_bits + pbits - 1) / pbits;
    nworkers = FLINT_MIN(nprimes, flint_get_num_threads());

    Abuf = _nmod_vec_init(FLINT_MAX(nworkers * totA, 1));
    Bbuf = _nmod_vec_init(FLINT_MAX(nworkers * totB, 1));
    Ap = flint_malloc(nworkers * sizeof(n_poly_struct *));
    Bp = flint_malloc(nworkers * sizeof(n_poly_struct *));
    Apbuf = flint_malloc(nworkers * lenA * sizeof(n_poly_struct));
    Bpbuf = flint_malloc(nworkers * lenB * sizeof(n_poly_struct));

    rp = flint_malloc(nworkers * sizeof(n_poly_struct));
    wstatus = flint_malloc(nworkers * sizeof(int));

    for (j = 0; j < nworkers; j++)
    {
        Ap[j] = Apbuf + j * lenA;
        Bp[j] = Bpbuf + j * lenB;

        totA = 0;
        for (i = 0; i < lenA; i++)
        {
            Ap[j][i].coeffs = Abuf + j * totAall + totA;
            Ap[j][i].alloc = A[i].length;
            Ap[j][i].length = 0;
            totA += A[i].length;
        }
        totB = 0;
        for (i = 0; i < lenB; i++)
        {
            Bp[j][i].coeffs = Bbuf + j * totBall + totB;
            Bp[j][i].alloc = B[i].length;
            Bp[j][i].length = 0;
            totB += B[i].length;
        }

        n_poly_init(rp + j);
        wstatus[j] = 1;
    }

    residues = flint_malloc(nprimes * outlen * sizeof(ulong));
    primes = flint_malloc(nprimes * sizeof(ulong));
    cand = _fmpz_vec_init(outlen);
    prev = _fmpz_vec_init(outlen);

    /* At depth d there are 2^(FFT_PRIME_BITS - 1 - d) candidates m, about a
       seventeenth of which give a prime, so running out takes a resultant of
       2^(FFT_PRIME_BITS - 1 - FFT_PRIME_MAX_DEPTH) / 17 primes at the largest
       exponent, far more than any input reaches. Falling back rather than
       failing costs a comparison and keeps the constants above free to
       change. */
    if (fft && _modular_primes(primes, nprimes, 1, depth, l) < nprimes)
    {
        fft = 0;
        nprimes = (bound_bits + MODULAR_PRIME_BITS - 1) / MODULAR_PRIME_BITS;
    }

    if (!fft)
        _modular_primes(primes, nprimes, 0, 0, l);

    /* Proved, primes are used until their product exceeds twice the bound
       above, which makes the reconstruction exact. Unproved, the images are
       computed in growing rounds and the reconstruction is stopped once it
       has been left unchanged by primes whose product carries
       MODULAR_STABLE_BITS bits, which is a good deal faster whenever the
       resultant is much smaller than the bound. That test is a heuristic and
       not a proof: it only fails to detect a wrong candidate when the primes
       of the last round divide the difference, and since they are picked
       deterministically an input can be built for which they do. */
    num = 0;
    stable_bits = 0;
    checkpoint = proved ? nprimes : nworkers;

    while (num < nprimes)
    {
        slong batch = FLINT_MIN(checkpoint, nprimes) - num;

        wargs.A = A; wargs.B = B;
        wargs.lenA = lenA; wargs.lenB = lenB;
        wargs.outlen = outlen;
        wargs.lo = num; wargs.hi = num + batch;
        wargs.nworkers = FLINT_MIN(batch, nworkers);
        wargs.primes = primes; wargs.residues = residues;
        wargs.Ap = Ap; wargs.Bp = Bp; wargs.rp = rp;
        wargs.status = wstatus;

        flint_parallel_do(_modular_worker, &wargs, wargs.nworkers, 0, 0);

        for (j = 0; j < nworkers; j++)
            if (!wstatus[j])
                status = 0;

        if (!status)
            break;

        num += batch;

        if (proved)
            break;

        /* prev keeps the reconstruction of the round before, so that cand is
           the newest one whichever way the loop is left */
        _fmpz_vec_swap(cand, prev, outlen);
        _crt_reconstruct(cand, outlen, residues, outlen, primes, num);

        if (num > batch && _fmpz_vec_equal(cand, prev, outlen))
        {
            for (j = num - batch; j < num; j++)
                stable_bits += FLINT_BIT_COUNT(primes[j]);

            if (stable_bits > MODULAR_STABLE_BITS)
                break;
        }
        else
        {
            stable_bits = 0;
        }

        /* rounds grow geometrically, so that the reconstructions the test
           costs amount to a bounded multiple of the last of them */
        checkpoint = FLINT_MAX(2 * checkpoint, num + nworkers);
    }

    if (status)
    {
        if (proved)
            _crt_reconstruct(cand, outlen, residues, outlen, primes, num);

        fmpz_poly_fit_length(res, outlen);
        for (k = 0; k < outlen; k++)
            fmpz_swap(res->coeffs + k, cand + k);
        _fmpz_poly_set_length(res, outlen);
        _fmpz_poly_normalise(res);

        /* put back the contents removed above */
        fmpz_pow_ui(ca, ca, lenB - 1);
        fmpz_pow_ui(cb, cb, lenA - 1);
        fmpz_mul(ca, ca, cb);
        fmpz_poly_scalar_mul_fmpz(res, res, ca);

        if (pa->length > 1)
        {
            fmpz_poly_pow(t, pa, lenB - 1);
            fmpz_poly_mul(res, res, t);
        }

        if (pb->length > 1)
        {
            fmpz_poly_pow(t, pb, lenA - 1);
            fmpz_poly_mul(res, res, t);
        }
    }

    for (j = 0; j < nworkers; j++)
        n_poly_clear(rp + j);

    _nmod_vec_clear(Abuf);
    _nmod_vec_clear(Bbuf);
    flint_free(Ap);
    flint_free(Bp);
    flint_free(Apbuf);
    flint_free(Bpbuf);
    flint_free(rp);
    flint_free(wstatus);

    _fmpz_vec_clear(cand, outlen);
    _fmpz_vec_clear(prev, outlen);
    flint_free(residues);
    flint_free(primes);

    fmpz_clear(ca);
    fmpz_clear(cb);
    fmpz_clear(l);
    fmpz_clear(u);
    fmpz_clear(bound);
    fmpz_poly_clear(pa);
    fmpz_poly_clear(pb);
    fmpz_poly_clear(t);

    return status;
}

int
fmpz_bpoly_resultant(fmpz_poly_t res, const fmpz_bpoly_t A,
                     const fmpz_bpoly_t B, int proved)
{
    slong lenA = A->length;
    slong lenB = B->length;
    fmpz_poly_struct * Ac, * Bc;
    slong i;
    int success;

    if (lenA == 0 || lenB == 0)
    {
        fmpz_poly_zero(res);
        return 1;
    }

    if (lenB > lenA)
        return fmpz_bpoly_resultant(res, B, A, proved) &&
            (((lenA | lenB) & 1) != 0 ||
                (fmpz_poly_neg(res, res), 1));

    if (lenB <= 1)
        return 0;

    /* the driver removes the contents of its inputs, so it works on copies */
    Ac = flint_malloc(lenA * sizeof(fmpz_poly_struct));
    Bc = flint_malloc(lenB * sizeof(fmpz_poly_struct));

    for (i = 0; i < lenA; i++)
    {
        fmpz_poly_init(Ac + i);
        fmpz_poly_set(Ac + i, A->coeffs + i);
    }
    for (i = 0; i < lenB; i++)
    {
        fmpz_poly_init(Bc + i);
        fmpz_poly_set(Bc + i, B->coeffs + i);
    }

    success = _fmpz_bpoly_resultant(res, Ac, lenA, Bc, lenB, proved);

    for (i = 0; i < lenA; i++)
        fmpz_poly_clear(Ac + i);
    for (i = 0; i < lenB; i++)
        fmpz_poly_clear(Bc + i);

    flint_free(Ac);
    flint_free(Bc);

    return success;
}
