/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz.h"
#include "fmpz_poly.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly_factor.h"

/* The variables other than the one eliminated are packed into one variable t
   by Kronecker substitution (mpoly/resultant_kronecker.c), and the resultant
   is taken over (t, y) by the bivariate resultant _fmpz_bpoly_resultant_len,
   which reduces modulo primes and leaves each image to the bivariate
   resultant over Z/pZ and whichever algorithm it picks. The dense form holds
   lenA times tlenA entries for A, and the resultant needs npoints evaluation
   points, one per coefficient of its dense box; these bound what the
   conversion can cost. They are also what keeps the size arithmetic from
   overflowing, which is why the comparisons are made in double. Lower than
   over a word-size prime field: the multimodular algorithm holds the
   residues of every prime at once, so its working space is a multiple of
   the number of points rather than of the input alone. */
#define BIVARIATE_MAX_AREA 1073741824.0     /* 2^30 */
#define BIVARIATE_MAX_POINTS 268435456.0    /* 2^28 */

/* A sparse input expands into a box far larger than itself. Below
   BIVARIATE_FREE_AREA entries that does not matter, the box being small
   whatever the input; above it the box may exceed the input by at most
   BIVARIATE_MAX_EXPANSION, so that a few terms of high degree are left to the
   sparse algorithm rather than materialised. */
#define BIVARIATE_FREE_AREA 8388608.0       /* 2^23 */
#define BIVARIATE_MAX_EXPANSION 64.0

/* Two polynomials of only a few terms each have a resultant the subresultant
   PRS finds almost at once, and at low degree the dense algorithm cannot make
   back the cost of the whole box. Measured over Z at degree 8, four terms
   each: three to a hundred times slower, at every coefficient size from 32 to
   10000 bits. From BIVARIATE_MIN_AREA in the degrees up it wins again, the
   PRS by then paying for the degrees whether or not the coefficients are
   there. */
#define BIVARIATE_MIN_TERMS 5
#define BIVARIATE_MIN_AREA 256

/* With two or more variables packed into t, the subresultant PRS runs over a
   polynomial ring in several variables, where it pays for the terms the
   pseudo-remainders fill in and for the growth of their coefficients.
   Measured on 3 to 5 variables at 32 and 1000 bits, the dense algorithm is
   ahead when both inputs have at least KRONECKER_MIN_TERMS terms and

      (T^2/KRONECKER_TERMS_SCALE)^(L - 2) >= exp(KRONECKER_LOG_BASE),

   T being the smaller number of terms and L the smaller length in y, the
   number of points hardly mattering against the coefficient growth of the
   PRS. Linear inputs in y (L = 2) are a single step and stay with the PRS.

   That was measured on inputs dense in a box, whose resultant fills its own
   box. The resultant of inputs dense in total degree only fills the part of
   its box within the bound on its total degree, about half of it with two
   variables packed and a sixth with three, which the PRS does not pay for
   and the dense algorithm does; as over nmod, the left-hand side is lowered
   by KRONECKER_FILL_POWER log(nterms/npoints), where nterms counts the
   monomials of that part. Over box inputs, and over inputs dense in total
   degree or structured as in implicitisation, integration and composed
   sums, this rule loses at most a factor 4.6 on an input taking a tenth of
   a millisecond and 2.7 above that, where the rule without that term lost
   up to a factor 50. */
#define KRONECKER_MIN_TERMS 4
#define KRONECKER_TERMS_SCALE 6.0
#define KRONECKER_LOG_BASE 8.0
#define KRONECKER_FILL_POWER 16.0

/* Whether the dense algorithm is expected to be faster than the subresultant
   PRS on the univariate forms, for inputs of Alen and Blen terms packed as
   described by K. */
int
_fmpz_mpoly_resultant_kronecker_cutoff(const mpoly_resultant_kronecker_t K,
                                       slong Alen, slong Blen)
{
    double areaA = (double) K->lenA * K->tlenA;
    double areaB = (double) K->lenB * K->tlenB;

    if (areaA > BIVARIATE_FREE_AREA && areaA > BIVARIATE_MAX_EXPANSION * Alen)
        return 0;

    if (areaB > BIVARIATE_FREE_AREA && areaB > BIVARIATE_MAX_EXPANSION * Blen)
        return 0;

    if (K->nused >= 2)
    {
        double T = FLINT_MIN(Alen, Blen);
        slong L = FLINT_MIN(K->lenA, K->lenB);

        if (T < KRONECKER_MIN_TERMS)
            return 0;

        return (L - 2)*log(T*T/KRONECKER_TERMS_SCALE)
                + KRONECKER_FILL_POWER*log(K->nterms/K->npoints) >=
                                                     KRONECKER_LOG_BASE;
    }

    if (FLINT_MIN(Alen, Blen) < BIVARIATE_MIN_TERMS &&
            K->lenA * K->lenB < BIVARIATE_MIN_AREA)
        return 0;

    return n_bpoly_mod_resultant_cutoff(FLINT_MAX(K->lenA, K->lenB),
                            FLINT_MIN(K->lenA, K->lenB), (slong) K->npoints);
}

/* the coefficients of A in the dense form of K, which are polynomials in t
   indexed by the degree in y */
static void
_fmpz_mpoly_kronecker_pack(fmpz_bpoly_t Ab, const fmpz_mpoly_t A, slong len,
          int isB, const mpoly_resultant_kronecker_t K,
          const fmpz_mpoly_ctx_t ctx)
{
    slong i, * ey, * et, * tlen;

    ey = (slong *) flint_malloc((2*A->length + len)*sizeof(slong));
    et = ey + A->length;
    tlen = et + A->length;

    mpoly_resultant_kronecker_pack(ey, et, A->exps, A->length, A->bits, isB,
                                                               K, ctx->minfo);
    for (i = 0; i < len; i++)
        tlen[i] = 0;
    for (i = 0; i < A->length; i++)
        tlen[ey[i]] = FLINT_MAX(tlen[ey[i]], et[i] + 1);

    fmpz_bpoly_fit_length(Ab, len);
    for (i = 0; i < len; i++)
    {
        fmpz_poly_zero(Ab->coeffs + i);
        fmpz_poly_fit_length(Ab->coeffs + i, tlen[i]);
    }
    Ab->length = len;

    for (i = 0; i < A->length; i++)
        fmpz_set(Ab->coeffs[ey[i]].coeffs + et[i], A->coeffs + i);

    for (i = 0; i < len; i++)
        _fmpz_poly_set_length(Ab->coeffs + i, tlen[i]);

    flint_free(ey);
}

/* moves the coefficients of r into R */
static void
_fmpz_mpoly_kronecker_unpack(fmpz_mpoly_t R, fmpz_poly_t r,
          const mpoly_resultant_kronecker_t K, const fmpz_mpoly_ctx_t ctx)
{
    slong i, k, N, len;

    len = 0;
    for (i = 0; i < r->length; i++)
        len += !fmpz_is_zero(r->coeffs + i);

    fmpz_mpoly_fit_length_reset_bits(R, len, K->Rbits, ctx);
    N = mpoly_words_per_exp(K->Rbits, ctx->minfo);

    k = 0;
    for (i = 0; i < r->length; i++)
    {
        if (fmpz_is_zero(r->coeffs + i))
            continue;

        fmpz_swap(R->coeffs + k, r->coeffs + i);
        mpoly_resultant_kronecker_unpack(R->exps + N*k, i, K, ctx->minfo);
        k++;
    }

    _fmpz_mpoly_set_length(R, len, ctx);
    fmpz_mpoly_sort_terms(R, ctx);
}

/* res_var(A, B) through Kronecker substitution and the dense bivariate
   multimodular algorithm. Returns 0, leaving R untouched, when the
   substitution does not apply, when the sizes are beyond the limits above,
   or when check_cutoff is set and the dense algorithm is not expected to be
   faster. */
static int
_fmpz_mpoly_resultant_kronecker(fmpz_mpoly_t R, const fmpz_mpoly_t A,
           const fmpz_mpoly_t B, slong var, int check_cutoff,
           const fmpz_mpoly_ctx_t ctx)
{
    mpoly_resultant_kronecker_t K;
    fmpz_bpoly_t Ab, Bb;
    fmpz_poly_t res;
    int success = 0;

    if (A->length == 0 || B->length == 0)
        return 0;

    mpoly_resultant_kronecker_init(K, ctx->minfo->nvars);

    if (!mpoly_resultant_kronecker_set(K, A->exps, A->length, A->bits,
                       B->exps, B->length, B->bits, var,
                       BIVARIATE_MAX_POINTS, ctx->minfo))
        goto cleanup_K;

    if ((double) K->lenA * K->tlenA > BIVARIATE_MAX_AREA ||
        (double) K->lenB * K->tlenB > BIVARIATE_MAX_AREA)
        goto cleanup_K;

    if (check_cutoff &&
        !_fmpz_mpoly_resultant_kronecker_cutoff(K, A->length, B->length))
        goto cleanup_K;

    fmpz_bpoly_init(Ab);
    fmpz_bpoly_init(Bb);
    fmpz_poly_init(res);

    _fmpz_mpoly_kronecker_pack(Ab, A, K->lenA, 0, K, ctx);
    _fmpz_mpoly_kronecker_pack(Bb, B, K->lenB, 1, K, ctx);

    /* the packed forms are working copies, which the algorithm may modify;
       the box of the resultant is smaller than the bound the algorithm reads
       off the degrees in t when the bounds of K are of weight above 0 */
    if (K->lenA >= K->lenB)
    {
        success = _fmpz_bpoly_resultant_len(res, Ab->coeffs, K->lenA,
                          Bb->coeffs, K->lenB, (slong) K->npoints, 1);
    }
    else
    {
        success = _fmpz_bpoly_resultant_len(res, Bb->coeffs, K->lenB,
                          Ab->coeffs, K->lenA, (slong) K->npoints, 1);

        /* res(A, B) = (-1)^(deg A deg B) res(B, A) */
        if (success && ((K->lenA | K->lenB) & 1) == 0)
            fmpz_poly_neg(res, res);
    }

    if (success)
        _fmpz_mpoly_kronecker_unpack(R, res, K, ctx);

    fmpz_poly_clear(res);
    fmpz_bpoly_clear(Ab);
    fmpz_bpoly_clear(Bb);

cleanup_K:

    mpoly_resultant_kronecker_clear(K);

    return success;
}

int fmpz_mpoly_resultant_kronecker(fmpz_mpoly_t R, const fmpz_mpoly_t A,
           const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)
{
    return _fmpz_mpoly_resultant_kronecker(R, A, B, var, 0, ctx);
}

int fmpz_mpoly_resultant(fmpz_mpoly_t R, const fmpz_mpoly_t A,
                   const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpoly_univar_t Ax, Bx;

    if (_fmpz_mpoly_resultant_kronecker(R, A, B, var, 1, ctx))
        return 1;

    fmpz_mpoly_univar_init(Ax, ctx);
    fmpz_mpoly_univar_init(Bx, ctx);

    fmpz_mpoly_to_univar(Ax, A, var, ctx);
    fmpz_mpoly_to_univar(Bx, B, var, ctx);

    success = fmpz_mpoly_univar_resultant(R, Ax, Bx, ctx);

    fmpz_mpoly_univar_clear(Ax, ctx);
    fmpz_mpoly_univar_clear(Bx, ctx);

    return success;
}
