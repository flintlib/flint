/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "mpoly.h"

/* The exponent side of res_y(A, B) through Kronecker substitution, shared by
   the coefficient rings.

   Write x_j for the variables other than y, n_A = deg_y A and n_B = deg_y B.
   Three ring homomorphisms of the coefficient ring commute with the
   resultant, the leading coefficients in y staying nonzero:

   - a monomial factor: res(x^a A, B) = x^(a n_B) res(A, B), so A and B are
     shifted down to their minimum exponent in each x_j;
   - deflation x_j -> x_j^g, injective, where g divides every exponent of
     x_j left in A and B;
   - Kronecker substitution x_j -> t^(S_j), which is injective on the
     polynomials of degree less than D_j in each x_j when S_j is the product
     of the D_i before it.

   With D_j one more than a bound on the degree in x_j of the resultant of
   the shifted and deflated polynomials, the resultant is recovered from its
   image in t digit by digit, and the bivariate algorithm needs D_1 ... D_k
   points, the number of coefficients of the dense box of the resultant.

   The bound comes from the Sylvester matrix. Give y the weight w >= 0 and let
   a(w) be the largest of e + w m over the terms x_j^e y^m of A, b(w) that of
   B. Row i of A holds the coefficient of y^(c - i) in column c, of degree at
   most a(w) + w i - w c in x_j, and the rows of B are alike, so that every
   term of the determinant has degree at most

       f(w) = n_B a(w) + n_A b(w) - w n_A n_B

   in x_j. The weight 0 gives n_B deg_j A + n_A deg_j B, which is reached by
   inputs dense in a box, and the weight 1 the Bezout bound d_A d_B for
   inputs of total degrees d_A and d_B in (x_j, y), half as much when those
   are their degrees in y. Since f(w) >= n_B deg_j A and n_A deg_j B, the
   substitution is also injective on A and B, which keeps their leading
   coefficients in y. */

/* Every weight gives a bound, so that _degree_bound tries the slopes of
   about WEIGHT_CANDIDATES vertices of the hulls at most. */
#define WEIGHT_CANDIDATES 32

/* The tables of E[m] for each variable and for the sum of the exponents
   cost (n_A + n_B + 2)(nvars + 1) words; above WEIGHT_TABLE_FACTOR times the
   number of terms plus WEIGHT_TABLE_FREE, the bounds of weight 0 are kept. */
#define WEIGHT_TABLE_FACTOR 4
#define WEIGHT_TABLE_FREE 1024

/* The monomials of the dense box within the bound on the total degree are
   counted by a sum of 2^k binomials of alternating signs over the k packed
   variables. Above BOX_SIMPLEX_MAX_VARS variables they are not counted and
   the whole box is taken instead, which leaves out the fill term of the
   cutoffs. */
#define BOX_SIMPLEX_MAX_VARS 12

void
mpoly_resultant_kronecker_init(mpoly_resultant_kronecker_t K, slong nvars)
{
    ulong * d;

    FLINT_ASSERT(nvars > 0);

    d = (ulong *) flint_malloc(7*nvars*sizeof(ulong));

    K->Amin    = d; d += nvars;
    K->Bmin    = d; d += nvars;
    K->stride  = d; d += nvars;
    K->bound   = d; d += nvars;
    K->kstride = d; d += nvars;
    K->rshift  = d; d += nvars;
    K->exps    = d;
}

void
mpoly_resultant_kronecker_clear(mpoly_resultant_kronecker_t K)
{
    flint_free(K->Amin);
}

/* minimum and maximum exponents of each variable, and the gcd of the
   differences of the exponents to those of the first term */
static void
_limits(ulong * min, ulong * max, ulong * g, ulong * first, ulong * e,
        const ulong * exps, slong len, flint_bitcnt_t bits,
        const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong N = mpoly_words_per_exp_sp(bits, mctx);
    slong i, j;

    mpoly_get_monomial_ui_sp(first, exps, bits, mctx);

    for (j = 0; j < nvars; j++)
    {
        min[j] = max[j] = first[j];
        g[j] = 0;
    }

    for (i = 1; i < len; i++)
    {
        mpoly_get_monomial_ui_sp(e, exps + N*i, bits, mctx);

        for (j = 0; j < nvars; j++)
        {
            min[j] = FLINT_MIN(min[j], e[j]);
            max[j] = FLINT_MAX(max[j], e[j]);
            g[j] = n_gcd(g[j], e[j] >= first[j] ? e[j] - first[j]
                                                : first[j] - e[j]);
        }
    }
}

/* Only the largest exponent E[m] among the terms of each degree m in y
   matters, and a(w) is the largest E[m] + w m over the vertices of the upper
   hull of the points (m, E[m]). Here E[m] is stored plus one, 0 marking a
   degree without terms, and the vertices are written to v from the last
   highest one on, the only ones a weight w >= 0 can pick. */
static slong
_upper_hull(slong * v, const ulong * E, slong n)
{
    slong m, h, top;

    h = 0;
    for (m = 0; m <= n; m++)
    {
        if (E[m] == 0)
            continue;

        /* drop the last vertex while it is not strictly above the segment
           from the one before it to m */
        while (h >= 2 && ((slong) E[v[h - 1]] - (slong) E[v[h - 2]])*(m - v[h - 2])
                      <= ((slong) E[m] - (slong) E[v[h - 2]])*(v[h - 1] - v[h - 2]))
            h--;

        v[h++] = m;
    }

    top = 0;
    for (m = 1; m < h; m++)
        if (E[v[m]] >= E[v[top]])
            top = m;

    for (m = top; m < h; m++)
        v[m - top] = v[m];

    return h - top;
}

/* q f(p/q) */
static ulong
_weighted_bound(const slong * vA, slong hA, const ulong * EA, ulong nA,
                const slong * vB, slong hB, const ulong * EB, ulong nB,
                ulong p, ulong q)
{
    ulong a = 0, b = 0;
    slong i;

    for (i = 0; i < hA; i++)
        a = FLINT_MAX(a, q*(EA[vA[i]] - 1) + p*vA[i]);

    for (i = 0; i < hB; i++)
        b = FLINT_MAX(b, q*(EB[vB[i]] - 1) + p*vB[i]);

    return nB*a + nA*b - p*nA*nB;
}

/* The minimum of f, rounded down. It is convex and piecewise linear, so that
   the minimum is at w = 0 or where a or b changes slope, at
   w = (E[m] - E[m'])/(m' - m) for consecutive vertices m < m' of either hull.
   Every w gives a bound, so that past WEIGHT_CANDIDATES vertices only some of
   them are tried. The caller makes sure that q f(p/q) fits in a word. */
static ulong
_degree_bound(slong * vA, const ulong * EA, ulong nA,
              slong * vB, const ulong * EB, ulong nB)
{
    slong hA = _upper_hull(vA, EA, nA);
    slong hB = _upper_hull(vB, EB, nB);
    slong i, step;
    ulong best, p, q;

    best = _weighted_bound(vA, hA, EA, nA, vB, hB, EB, nB, 0, 1);

    step = 1 + (hA + hB)/WEIGHT_CANDIDATES;

    for (i = 0; i + 1 < hA; i += step)
    {
        p = EA[vA[i]] - EA[vA[i + 1]];
        q = vA[i + 1] - vA[i];
        best = FLINT_MIN(best,
                _weighted_bound(vA, hA, EA, nA, vB, hB, EB, nB, p, q)/q);
    }

    for (i = 0; i + 1 < hB; i += step)
    {
        p = EB[vB[i]] - EB[vB[i + 1]];
        q = vB[i + 1] - vB[i];
        best = FLINT_MIN(best,
                _weighted_bound(vA, hA, EA, nA, vB, hB, EB, nB, p, q)/q);
    }

    return best;
}

/* the number of points of [0, D_0) x ... x [0, D_(k-1)) of coordinate sum at
   most T, by inclusion and exclusion */
static double
_box_simplex(const double * D, slong k, double T)
{
    double total = 0.0, s, c;
    slong S, j;
    int neg;

    for (S = 0; S < (WORD(1) << k); S++)
    {
        s = T;
        neg = 0;
        for (j = 0; j < k; j++)
        {
            if (S & (WORD(1) << j))
            {
                s -= D[j];
                neg = !neg;
            }
        }

        if (s < 0.0)
            continue;

        c = 1.0;
        for (j = 1; j <= k; j++)
            c = c*(s + j)/j;

        total += neg ? -c : c;
    }

    return total;
}

/* the tables of E[m] of the shifted and deflated exponents: column j for
   each variable j other than var, and column nvars for the sum of their
   exponents */
static void
_weight_tables(ulong * E, slong n, const ulong * exps, slong len,
               flint_bitcnt_t bits, const ulong * min,
               const mpoly_resultant_kronecker_t K, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong N = mpoly_words_per_exp_sp(bits, mctx);
    ulong * e = K->exps;
    ulong d, sum;
    slong i, j, m;

    for (i = 0; i < len; i++)
    {
        mpoly_get_monomial_ui_sp(e, exps + N*i, bits, mctx);
        m = e[K->var];

        sum = 0;
        for (j = 0; j < nvars; j++)
        {
            if (j == K->var)
                continue;

            d = (e[j] - min[j])/K->stride[j];
            sum += d;
            E[j*(n + 1) + m] = FLINT_MAX(E[j*(n + 1) + m], d + 1);
        }

        E[nvars*(n + 1) + m] = FLINT_MAX(E[nvars*(n + 1) + m], sum + 1);
    }
}

/* Sets up the substitution for res_var(A, B). Returns 0 when it does not
   apply: an input with exponents of more than one word, a degree below one in
   var, no other variable left after shifting and deflating, more than
   max_points coefficients in the dense box of the resultant, or exponents of
   the resultant not fitting in one word. */
int
mpoly_resultant_kronecker_set(mpoly_resultant_kronecker_t K,
                   const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
                   const ulong * Bexps, slong Blen, flint_bitcnt_t Bbits,
                   slong var, double max_points, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong j;
    ulong nA, nB, S, tA, tB, dA, dB, sumA, sumB, tot;
    ulong * Amax, * Bmax, * Ag, * Bg, * first, * EA, * EB;
    slong * vA, * vB;
    double npoints, rmax, rsum, limit, nmax;
    double * D;
    int weighted, success = 0;
    TMP_INIT;

    if (Alen < 1 || Blen < 1 || Abits > FLINT_BITS || Bbits > FLINT_BITS)
        return 0;

    TMP_START;

    Amax  = (ulong *) TMP_ALLOC(5*nvars*sizeof(ulong));
    Bmax  = Amax + nvars;
    Ag    = Bmax + nvars;
    Bg    = Ag + nvars;
    first = Bg + nvars;

    _limits(K->Amin, Amax, Ag, first, K->exps, Aexps, Alen, Abits, mctx);
    _limits(K->Bmin, Bmax, Bg, first, K->exps, Bexps, Blen, Bbits, mctx);

    nA = Amax[var];
    nB = Bmax[var];

    /* every length and degree below is at most the number of points, so
       bounding that keeps the arithmetic in slong */
    max_points = FLINT_MIN(max_points, (double) (WORD_MAX/4));

    if (nA < 1 || nB < 1 || (double) nA >= max_points ||
                            (double) nB >= max_points)
        goto cleanup;

    /* the exponents of the resultant with the bounds of weight 0, in double
       first, which also keeps those bounds and their sum in a word */
    limit = ldexp(1.0, FLINT_BITS - 3);
    rsum = 0.0;
    for (j = 0; j < nvars; j++)
    {
        if (j == var)
            continue;

        K->stride[j] = n_gcd(Ag[j], Bg[j]);
        if (K->stride[j] == 0)
            K->stride[j] = 1;

        dA = (Amax[j] - K->Amin[j])/K->stride[j];
        dB = (Bmax[j] - K->Bmin[j])/K->stride[j];

        rmax = ((double) nB*dA + (double) nA*dB)*K->stride[j]
             + (double) K->Amin[j]*nB + (double) K->Bmin[j]*nA;

        if (rmax >= limit)
            goto cleanup;

        rsum += rmax;
    }

    if (rsum >= limit)
        goto cleanup;

    /* now exactly */
    K->var = var;
    K->lenA = nA + 1;
    K->lenB = nB + 1;

    sumA = sumB = 0;
    for (j = 0; j < nvars; j++)
    {
        if (j == var)
        {
            K->stride[j] = 1;
            K->bound[j] = 1;
            continue;
        }

        dA = (Amax[j] - K->Amin[j])/K->stride[j];
        dB = (Bmax[j] - K->Bmin[j])/K->stride[j];
        sumA += dA;
        sumB += dB;

        K->bound[j] = nB*dA + nA*dB;
    }

    tot = nB*sumA + nA*sumB;

    /* the bounds of other weights, when the tables are affordable and q f(p/q)
       fits in a word: its terms are at most n f(0) with n = max(n_A, n_B) */
    nmax = (double) FLINT_MAX(nA, nB);
    weighted = (double) (nA + nB + 2)*(nvars + 1) <=
                     (double) WEIGHT_TABLE_FACTOR*(Alen + Blen)*(nvars + 1)
                                                         + WEIGHT_TABLE_FREE
            && 8.0*nmax*((double) tot + 1.0) < ldexp(1.0, FLINT_BITS - 1);

    if (weighted)
    {
        EA = (ulong *) flint_calloc((nvars + 1)*(nA + nB + 2), sizeof(ulong));
        EB = EA + (nvars + 1)*(nA + 1);
        vA = (slong *) flint_malloc((nA + nB + 2)*sizeof(slong));
        vB = vA + nA + 1;

        _weight_tables(EA, nA, Aexps, Alen, Abits, K->Amin, K, mctx);
        _weight_tables(EB, nB, Bexps, Blen, Bbits, K->Bmin, K, mctx);

        for (j = 0; j < nvars; j++)
            if (j != var && K->bound[j] > 0)
                K->bound[j] = _degree_bound(vA, EA + j*(nA + 1), nA,
                                            vB, EB + j*(nB + 1), nB);

        tot = _degree_bound(vA, EA + nvars*(nA + 1), nA,
                            vB, EB + nvars*(nB + 1), nB);

        flint_free(EA);
        flint_free(vA);
    }

    /* the dense box of the resultant */
    npoints = 1.0;
    for (j = 0; j < nvars; j++)
    {
        if (j == var)
            continue;

        K->bound[j] += 1;
        npoints *= (double) K->bound[j];
    }

    if (npoints > max_points)
        goto cleanup;

    K->npoints = npoints;

    S = 1;
    tA = tB = 0;
    K->nused = 0;
    D = (double *) TMP_ALLOC(nvars*sizeof(double));
    for (j = 0; j < nvars; j++)
    {
        if (j == var)
        {
            K->kstride[j] = 0;
            K->rshift[j] = 0;
            K->exps[j] = 0;
            continue;
        }

        dA = (Amax[j] - K->Amin[j])/K->stride[j];
        dB = (Bmax[j] - K->Bmin[j])/K->stride[j];

        K->rshift[j] = K->Amin[j]*nB + K->Bmin[j]*nA;
        K->exps[j] = (K->bound[j] - 1)*K->stride[j] + K->rshift[j];

        if (K->bound[j] > 1)
        {
            K->kstride[j] = S;
            tA += dA*S;
            tB += dB*S;
            S *= K->bound[j];
            D[K->nused] = (double) K->bound[j];
            K->nused++;
        }
        else
        {
            K->kstride[j] = 0;
        }
    }

    if (K->nused == 0)
        goto cleanup;

    K->tlenA = tA + 1;
    K->tlenB = tB + 1;

    /* the monomials of the box within the bound on the total degree */
    K->nterms = npoints;
    if (weighted && K->nused <= BOX_SIMPLEX_MAX_VARS)
    {
        double c = _box_simplex(D, K->nused, (double) tot);
        K->nterms = FLINT_MAX(1.0, FLINT_MIN(npoints, c));
    }

    K->Rbits = mpoly_exp_bits_required_ui(K->exps, mctx);
    K->Rbits = mpoly_fix_bits(K->Rbits, mctx);
    FLINT_ASSERT(K->Rbits <= FLINT_BITS);

    success = 1;

cleanup:

    TMP_END;

    return success;
}

/* the degree in y and the exponent of t of each term of A, or of B if B is
   set */
void
mpoly_resultant_kronecker_pack(slong * ey, slong * et,
                   const ulong * exps, slong len, flint_bitcnt_t bits, int B,
                   const mpoly_resultant_kronecker_t K, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong N = mpoly_words_per_exp_sp(bits, mctx);
    const ulong * min = B ? K->Bmin : K->Amin;
    ulong * e = K->exps;
    slong i, j;
    ulong t;

    for (i = 0; i < len; i++)
    {
        mpoly_get_monomial_ui_sp(e, exps + N*i, bits, mctx);

        t = 0;
        for (j = 0; j < nvars; j++)
            if (K->kstride[j] != 0)
                t += (e[j] - min[j])/K->stride[j]*K->kstride[j];

        ey[i] = e[K->var];
        et[i] = t;
    }
}

/* the exponent, packed in K->Rbits, of the monomial of the resultant whose
   image is t^t */
void
mpoly_resultant_kronecker_unpack(ulong * exp, ulong t,
                   const mpoly_resultant_kronecker_t K, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    ulong * e = K->exps;
    slong j;

    for (j = 0; j < nvars; j++)
    {
        if (K->kstride[j] == 0)
            e[j] = K->rshift[j];
        else
            e[j] = (t/K->kstride[j]) % K->bound[j]*K->stride[j] + K->rshift[j];
    }

    mpoly_set_monomial_ui(exp, e, K->Rbits, mctx);
}
