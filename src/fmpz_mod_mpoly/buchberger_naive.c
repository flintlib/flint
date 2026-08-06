/*
    Copyright (C) 2020, 2026 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/* Implements a version of GROEBNERNEW2 (Becker and Weispfenning, p. 232).
   Parts of the implementation were inspired by SymPy's groebner().
   The current version of the code has been substantially rewritten using
   Claude to incorporate optimizations from BW that were omitted or
   done incorrectly in the original hand-written version.

   Note: the fmpz_mpoly and fmpz_mod_mpoly implementations are essentially
   identical. Any future improvements should be applied to both. */

/* Helper functions for index pairs. */

typedef struct
{
    slong a;
    slong b;
} pair_t;

typedef struct
{
    pair_t * pairs;
    slong length;
    slong alloc;
}
pairs_struct;

typedef pairs_struct pairs_t[1];

static void
pairs_init(pairs_t vec)
{
    vec->pairs = NULL;
    vec->length = 0;
    vec->alloc = 0;
}

static void
pairs_fit_length(pairs_t vec, slong len)
{
    if (len > vec->alloc)
    {
        if (len < 2 * vec->alloc)
            len = 2 * vec->alloc;

        vec->pairs = flint_realloc(vec->pairs, len * sizeof(pair_t));
        vec->alloc = len;
    }
}

static void
pairs_clear(pairs_t vec)
{
    flint_free(vec->pairs);
}

static void
pairs_append(pairs_t vec, slong i, slong j)
{
    pairs_fit_length(vec, vec->length + 1);
    vec->pairs[vec->length].a = i;
    vec->pairs[vec->length].b = j;
    vec->length++;
}

/* Two selection strategies are available.

   BUCHBERGER_DEGREE_SELECTION 1 (default): Degree strategy.
   Select the pair whose lcm(LM(f), LM(g)) has the smallest total degree,
   for ALL monomial orderings. For a lex ordering, equal-degree pairs are
   broken by lex comparison. For degree orderings (degrevlex, deglex),
   equal-degree pairs are broken by position in the pairs vector.

   Motivation: pure lex comparison (BW's strategy, below) causes
   pathological blowup on certain lex systems that we encounter in
   practice in the Calcium test suite.

   BUCHBERGER_DEGREE_SELECTION 0: BW Normal strategy.
   Select the pair whose lcm(LM(f), LM(g)) is smallest in the monomial order
   being used, for all orderings. For lex this is a pure lexicographic
   comparison. For degree orderings it coincides with minimum total degree
   (same as the degree strategy), so the two variants differ only for lex.
*/
#define BUCHBERGER_DEGREE_SELECTION 1

static pair_t
fmpz_mod_mpoly_select_pop_pair(pairs_t pairs, const fmpz_mod_mpoly_vec_t G,
                           const fmpz_mod_mpoly_ctx_t ctx)
{
    slong len, choice, nvars;
    pair_t result;
 
    nvars = ctx->minfo->nvars;
    len = pairs->length;
    choice = 0;
 
    if (len > 1)
    {
        slong i, j, a, b;
        ulong * exp;
        ulong * lcm;
        ulong * best_lcm;
        ulong l, total;
        int best;
 
        exp = flint_malloc(sizeof(ulong) * G->length * nvars);
        lcm = flint_malloc(sizeof(ulong) * (nvars + 1));
        best_lcm = flint_malloc(sizeof(ulong) * (nvars + 1));
 
        for (i = 0; i <= nvars; i++)
            best_lcm[i] = UWORD_MAX;
 
        for (i = 0; i < G->length; i++)
            fmpz_mod_mpoly_get_term_exp_ui(exp + i * nvars, G->p + i, 0, ctx);
 
        for (i = 0; i < len; i++)
        {
            a = pairs->pairs[i].a;
            b = pairs->pairs[i].b;
            total = 0;
            best = 1;
 
            /* Compute the full lcm vector and total degree unconditionally.
             * Both strategies need lcm[] fully populated before comparing
             * so that best_lcm[] is always written with complete data. */
            for (j = 0; j < nvars; j++)
            {
                l = FLINT_MAX(exp[a * nvars + j], exp[b * nvars + j]);
                lcm[j] = l;
                total += l;
            }
 
#if BUCHBERGER_DEGREE_SELECTION
            /* Degree strategy  */
            if (total > best_lcm[nvars])
            {
                best = 0;
            }
            else if (total == best_lcm[nvars] && ctx->minfo->ord == ORD_LEX)
            {
                for (j = 0; j < nvars; j++)
                {
                    if (lcm[j] > best_lcm[j]) { best = 0; break; }
                    if (lcm[j] < best_lcm[j]) break;
                }
            }
            else if (total == best_lcm[nvars])
            {
                best = 0;
            }

#else
            /* BW Normal strategy */
            if (ctx->minfo->ord == ORD_LEX)
            {
                for (j = 0; j < nvars; j++)
                {
                    if (lcm[j] > best_lcm[j]) { best = 0; break; }
                    if (lcm[j] < best_lcm[j]) break;
                }
            }
            else
            {
                if (total > best_lcm[nvars])
                    best = 0;
                else if (total == best_lcm[nvars])
                    best = 0;
            }
#endif
            if (best)
            {
                for (j = 0; j < nvars; j++)
                    best_lcm[j] = lcm[j];
 
                best_lcm[nvars] = total;
                choice = i;
            }
        }
 
        flint_free(exp);
        flint_free(lcm);
        flint_free(best_lcm);
    }
 
    result = pairs->pairs[choice];
    pairs->pairs[choice] = pairs->pairs[pairs->length - 1];
    pairs->length--;
 
    return result;
}

/* Evaluation limits to allow aborting a basis that spirals out of control. */
static int
within_limits(const fmpz_mod_mpoly_t poly, slong poly_len_limit,
              const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_length(poly, ctx) > poly_len_limit)
        return 0;

    return 1;
}

/* Helper functions for exponent vectors, packed to 1 ulong per exponent.
   TODO: better packing for small exponents? If we ensure that all
   polynomials pack to the same number of bits (resizing if necessary
   during the algorithm), we could maybe work directly with the
   exponent vectors stored in the polynomials? */

static int
monomial_divides(const ulong * a, const ulong * b, slong nvars)
{
    slong i;

    for (i = 0; i < nvars; i++)
        if (a[i] > b[i])
            return 0;

    return 1;
}

static int
monomial_disjoint(const ulong * a, const ulong * b, slong nvars)
{
    slong i;

    for (i = 0; i < nvars; i++)
        if (a[i] && b[i])
            return 0;

    return 1;
}

static void
monomial_lcm(ulong * out, const ulong * a, const ulong * b, slong nvars)
{
    slong i;

    for (i = 0; i < nvars; i++)
        out[i] = FLINT_MAX(a[i], b[i]);
}

static int
monomial_equal(const ulong * a, const ulong * b, slong nvars)
{
    slong i;

    for (i = 0; i < nvars; i++)
        if (a[i] != b[i])
            return 0;

    return 1;
}

/* Basis-element filtering. When a new polynomial h at index ih is added to
   the basis, any existing basis element g whose leading monomial is
   divisible by LM(h) is redundant: every future S-polynomial involving g as
   a basis candidate is already subsumed by pairs involving h. Removing such
   g from G_active shrinks the candidate set for future calls to update_pairs.

   Note: pairs already in the pair set that mention a removed index are left
   in place. Their S-polynomials are not guaranteed to reduce to zero
   immediately (they may still yield new generators), and update_pairs's own
   pair-filtering step already handles the subset of redundant old pairs.
   The polynomial store G is also left untouched; removed indices remain valid
   for spoly and reduction computations. */
static void
filter_redundant(slong * G_active, slong * G_active_len,
                 slong ih, const ulong * exp, slong nvars)
{
    const ulong * mh;
    slong i, new_active_len;

    mh = exp + ih * nvars;

    new_active_len = 0;
    for (i = 0; i < *G_active_len; i++)
    {
        slong ig = G_active[i];

        if (!monomial_divides(mh, exp + ig * nvars, nvars))
            G_active[new_active_len++] = ig;
    }
    *G_active_len = new_active_len;
}

/* Gebauer-Moller pair-set update (Becker & Weispfenning, p. 230).

   Given the current pair set B and a new polynomial at index ih (with leading
   monomial exp[ih * nvars ...]), update B as follows:

   1. Compute the set E of new pairs (ih, ig) to add, using the chain
      criterion to discard candidates dominated by a pair with smaller lcm,
      and the product criterion to discard pairs with disjoint leading
      monomials (the C -> D -> E reduction in BW's description).

   2. Remove from B any existing pair (ig1, ig2) made redundant by ih:
      those where LM(ih) | lcm(LM(ig1), LM(ig2)), subject to two
      correctness guards that prevent over-aggressive removal.

   3. Append E to B.
*/
static void
update_pairs(const slong * G_active, slong G_active_len,
             pairs_t B, slong ih, const ulong * exp, slong nvars)
{
    const ulong * mh;
    ulong * lcm_hg;
    ulong * lcm_tmp;
    slong * D_ig;
    slong * E_ig;
    slong D_alloc, D_len, E_len;
    slong i, k;
    int dominated;

    mh = exp + ih * nvars;
    lcm_hg  = flint_malloc(nvars * sizeof(ulong));
    lcm_tmp = flint_malloc(nvars * sizeof(ulong));

    /* Step 1 (C -> D): decide which candidate pairs (ih, ig) to accept.

       Process G_active left-to-right.  At step i, the "C" set of BW's
       algorithm is G_active[i + 1 ..] (not yet visited) and "D" is the
       pairs already accepted.  Accept (ih, ig) unless some p already in D
       or still in C satisfies lcm(LM(h), LM(p)) | lcm(LM(h), LM(g)) --
       the pair (h, p) then dominates (h, g) by the chain criterion.
       Coprime pairs (product criterion) are also accepted here; Step 2
       removes them again. */
    D_alloc = G_active_len + 1;
    D_len   = 0;
    D_ig    = flint_malloc(D_alloc * sizeof(slong));

    for (i = 0; i < G_active_len; i++)
    {
        slong ig = G_active[i];
        const ulong * mg = exp + ig * nvars;

        if (!monomial_disjoint(mh, mg, nvars))
        {
            monomial_lcm(lcm_hg, mh, mg, nvars);
            dominated = 0;

            for (k = 0; k < D_len && !dominated; k++)
            {
                monomial_lcm(lcm_tmp, mh, exp + D_ig[k] * nvars, nvars);
                if (monomial_divides(lcm_tmp, lcm_hg, nvars))
                    dominated = 1;
            }

            for (k = i + 1; k < G_active_len && !dominated; k++)
            {
                monomial_lcm(lcm_tmp, mh, exp + G_active[k] * nvars, nvars);
                if (monomial_divides(lcm_tmp, lcm_hg, nvars))
                    dominated = 1;
            }

            if (dominated)
                continue;
        }

        if (D_len == D_alloc)
        {
            D_alloc *= 2;
            D_ig = flint_realloc(D_ig, D_alloc * sizeof(slong));
        }

        D_ig[D_len++] = ig;
    }

    /* Step 2 (D -> E): apply the product criterion to D. */
    E_ig  = flint_malloc((D_len + 1) * sizeof(slong));
    E_len = 0;

    for (i = 0; i < D_len; i++)
    {
        slong ig = D_ig[i];

        if (!monomial_disjoint(mh, exp + ig * nvars, nvars))
            E_ig[E_len++] = ig;
    }

    flint_free(D_ig);

    /* Step 3: remove stale pairs from B (Gebauer-Moller criterion).

       Discard (ig1, ig2) if all three conditions hold:
         (a) LM(h) | lcm(LM(ig1), LM(ig2))
         (b) lcm(LM(ig1), LM(h)) != lcm(LM(ig1), LM(ig2))
         (c) lcm(LM(ig2), LM(h)) != lcm(LM(ig1), LM(ig2))
      
       The equality guards (b) and (c) prevent incorrectly removing a pair
       for which one element's lcm with h already equals the pair's lcm. */
    {
        slong new_len = 0;

        for (k = 0; k < B->length; k++)
        {
            slong ig1 = B->pairs[k].a;
            slong ig2 = B->pairs[k].b;
            const ulong * mg1 = exp + ig1 * nvars;
            const ulong * mg2 = exp + ig2 * nvars;

            monomial_lcm(lcm_hg, mg1, mg2, nvars);

            if (!monomial_divides(mh, lcm_hg, nvars))
                goto keep;

            monomial_lcm(lcm_tmp, mg1, mh, nvars);
            if (monomial_equal(lcm_tmp, lcm_hg, nvars))
                goto keep;

            monomial_lcm(lcm_tmp, mg2, mh, nvars);
            if (monomial_equal(lcm_tmp, lcm_hg, nvars))
                goto keep;

            continue;   /* discard */

        keep:

            B->pairs[new_len++] = B->pairs[k];
        }

        B->length = new_len;
    }

    /* Step 4: append the surviving new pairs E to B. */
    for (i = 0; i < E_len; i++)
        pairs_append(B, ih, E_ig[i]);

    flint_free(E_ig);
    flint_free(lcm_hg);
    flint_free(lcm_tmp);
}

int
fmpz_mod_mpoly_buchberger_naive_with_limits(fmpz_mod_mpoly_vec_t G,
    const fmpz_mod_mpoly_vec_t F,
    slong ideal_len_limit, slong poly_len_limit,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    pairs_t B;
    fmpz_mod_mpoly_t h;
    slong * G_active;
    ulong * exp;
    slong i, ih, G_active_len, G_active_alloc, exp_alloc, nvars;
    pair_t pair;
    int success;

    fmpz_mod_mpoly_vec_set_monic_unique(G, F, ctx);

    if (G->length <= 1)
        return 1;

    if (G->length >= ideal_len_limit)
        return 0;

    for (i = 0; i < G->length; i++)
        if (!within_limits(fmpz_mod_mpoly_vec_entry(G, i), poly_len_limit, ctx))
            return 0;

    nvars = ctx->minfo->nvars;

    /* Cache leading monomial exponent vectors for fast comparisons. */
    exp_alloc = G->length + 16;
    exp = flint_malloc(exp_alloc * nvars * sizeof(ulong));

    for (i = 0; i < G->length; i++)
        fmpz_mod_mpoly_get_term_exp_ui(exp + i * nvars,
                                       fmpz_mod_mpoly_vec_entry(G, i), 0, ctx);

    /* G_active tracks which indices in G belong to the current basis,
       allowing update_pairs to enumerate candidates without disturbing the
       storage layout that spoly and reduction use. */
    G_active_alloc = exp_alloc;
    G_active = flint_malloc(G_active_alloc * sizeof(slong));
    G_active_len = 0;

    pairs_init(B);
    fmpz_mod_mpoly_init(h, ctx);

    /* Build the initial pair set by simulating GROEBNERNEW2's initialisation
       loop: start with G_active = {0} and call update_pairs for
       each subsequent initial polynomial in turn. The original code only
       applied the product criterion here; update_pairs applies the full
       Gebauer-Moller criteria, eliminating more redundant pairs up front. */
    G_active[G_active_len++] = 0;

    for (i = 1; i < G->length; i++)
    {
        update_pairs(G_active, G_active_len, B, i, exp, nvars);
        filter_redundant(G_active, &G_active_len, i, exp, nvars);
        G_active[G_active_len++] = i;
    }

    /* Main loop -- GROEBNERNEW2 */
    success = 1;

    while (B->length != 0)
    {
        pair = fmpz_mod_mpoly_select_pop_pair(B, G, ctx);

        fmpz_mod_mpoly_spoly(h, fmpz_mod_mpoly_vec_entry(G, pair.a),
                                fmpz_mod_mpoly_vec_entry(G, pair.b), ctx);
        fmpz_mod_mpoly_reduction_monic_part(h, h, G, ctx);

        if (!fmpz_mod_mpoly_is_zero(h, ctx))
        {
            if (G->length >= ideal_len_limit ||
                !within_limits(h, poly_len_limit, ctx))
            {
                success = 0;
                break;
            }

            ih = G->length;
            fmpz_mod_mpoly_vec_append(G, h, ctx);

            if (ih >= exp_alloc)
            {
                exp_alloc = FLINT_MAX(exp_alloc * 2, ih + 1);
                exp = flint_realloc(exp, exp_alloc * nvars * sizeof(ulong));
            }

            fmpz_mod_mpoly_get_term_exp_ui(exp + ih * nvars,
                                           fmpz_mod_mpoly_vec_entry(G, ih), 0,
                                           ctx);

            if (G_active_len + 1 > G_active_alloc)
            {
                G_active_alloc = FLINT_MAX(G_active_alloc * 2,
                                           G_active_len + 2);
                G_active = flint_realloc(G_active,
                                         G_active_alloc * sizeof(slong));
            }

            update_pairs(G_active, G_active_len, B, ih, exp, nvars);
            filter_redundant(G_active, &G_active_len, ih, exp, nvars);
            G_active[G_active_len++] = ih;
        }
    }

    /* Compact G to the active basis in-place. G has accumulated every
       polynomial ever added, including those later made redundant by
       basis-element filtering.

       G_active is strictly increasing throughout the algorithm: indices are
       appended in the order they are added to the store (always ih = old
       G->length, strictly larger than everything already there), and
       filter_redundant only removes elements while preserving the order of
       survivors. Therefore G_active[i] >= i for all i. We exploit this to
       compact G with a single forward pass. */
    for (i = 0; i < G_active_len; i++)
    {
        FLINT_ASSERT(G_active[i] >= i);
        fmpz_mod_mpoly_swap(fmpz_mod_mpoly_vec_entry(G, i),
                        fmpz_mod_mpoly_vec_entry(G, G_active[i]), ctx);
    }

    fmpz_mod_mpoly_clear(h, ctx);
    pairs_clear(B);
    flint_free(G_active);
    flint_free(exp);

    return success;
}

void
fmpz_mod_mpoly_buchberger_naive(fmpz_mod_mpoly_vec_t G,
    const fmpz_mod_mpoly_vec_t F, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_buchberger_naive_with_limits(G, F, WORD_MAX, WORD_MAX, ctx);
}
