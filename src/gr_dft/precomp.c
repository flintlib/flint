/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "gr_dft.h"

/* Tuning: with AUTO, switch from plain Cooley-Tukey to the four-step
   algorithm at this depth. */

/* Tuning: primes up to this bound use the direct O(p^2) kernel inside
   the mixed-radix algorithm; larger primes attempt Bluestein. */
#define GR_DFT_PRIME_CUTOFF 19

/* Tuning: in complex Karatsuba mode over a ring with a working
   precision (e.g. acb over arb), the 3-real-multiplication products
   are only enabled at this precision and above; below it, a native
   complex multiplication is faster than three multiplications plus
   extra additions in the real ring. The special rotations (the classes
   in wclass) are used at any precision. Over exact rings, where the
   number of multiplications is what matters, Karatsuba products are
   always enabled. */
#define GR_DFT_KARATSUBA_CUTOFF_PREC 1024

int
gr_dft_default_root(gr_ptr w, ulong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (n <= 1)
    {
        status |= gr_one(w, ctx);
    }
    else if (n == 2)
    {
        status |= gr_neg_one(w, ctx);
    }
    else
    {
        /* w = exp(-2 pi i / n), for rings which support exp_pi_i */
        gr_ptr t, u;
        fmpq_t q;
        GR_TMP_INIT2(t, u, ctx);
        fmpq_init(q);
        fmpq_set_si(q, -2, n);
        status |= gr_set_fmpq(t, q, ctx);
        status |= gr_exp_pi_i(w, t, ctx);

        if (status != GR_SUCCESS)
        {
            /* fallback for rings with pi, i and exp but no exp_pi_i */
            status = GR_SUCCESS;
            status |= gr_pi(t, ctx);
            status |= gr_i(u, ctx);
            status |= gr_mul(t, t, u, ctx);
            status |= gr_mul_si(t, t, -2, ctx);
            status |= gr_div_ui(t, t, n, ctx);
            status |= gr_exp(w, t, ctx);
        }

        fmpq_clear(q);
        GR_TMP_CLEAR2(t, u, ctx);
    }

    return status;
}

/* For the default root w = exp(-2 pi i / n), attempt to compute the table
   entries directly as w^j = exp_pi_i(-2j/n). Besides giving tighter
   enclosures than a power chain over inexact rings, this typically yields
   the entries at the quarter points (1, -i, -1, i) exactly, allowing them
   to be classified as free rotations in complex Karatsuba mode. */
static int
_gr_dft_roots_canonical(gr_ptr roots, ulong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong j;
    gr_ptr t;
    fmpq_t q;

    /* Fast path for acb: one primitive root computed at elevated
       precision, the remaining entries by a running product rounded
       back to the working precision, with exact values (and error
       resets) at the quarter points. The running product accumulates
       a radius of order j ulp at the elevated precision, so the
       rounded entries have radii within a small factor of one ulp at
       the working precision -- comparable to evaluating each entry
       with exp, at a fraction of the cost. */
    if (ctx->which_ring == GR_CTX_CC_ACB && sz == sizeof(acb_struct))
    {
        acb_ptr tab = (acb_ptr) roots;
        acb_t r, w;
        slong prec, eprec;
        ulong L;

        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        eprec = prec + FLINT_BIT_COUNT(n) + 16;

        acb_init(r);
        acb_init(w);

        fmpq_init(q);
        fmpq_set_si(q, -2, n);
        arb_sin_cos_pi_fmpq(acb_imagref(w), acb_realref(w), q, eprec);
        fmpq_clear(q);

        /* entries up to an eighth (4 | n), quarter (2 | n) or half
           (odd n) of the circle by a running product; the rest by
           the eight-fold symmetry */
        if (n % 4 == 0)
            L = (n + 7) / 8;
        else if (n % 2 == 0)
            L = (n + 3) / 4;
        else
            L = n / 2;

        acb_one(r);
        for (j = 0; j <= L; j++)
        {
            acb_set_round(tab + j, r, prec);
            acb_mul(r, r, w, eprec);
        }

        /* w^(n/4 - j) = -y - x i, w^(n/2 - j) = -x + y i,
           w^(n - j) = x - y i, for w^j = x + y i (exact copies) */
        if (n % 4 == 0)
        {
            for (j = L + 1; j <= n / 4; j++)
            {
                arb_neg(acb_realref(tab + j), acb_imagref(tab + n / 4 - j));
                arb_neg(acb_imagref(tab + j), acb_realref(tab + n / 4 - j));
            }
        }
        if (n % 2 == 0)
        {
            for (j = FLINT_MAX(L, n / 4) + 1; j <= n / 2; j++)
            {
                arb_neg(acb_realref(tab + j), acb_realref(tab + n / 2 - j));
                arb_set(acb_imagref(tab + j), acb_imagref(tab + n / 2 - j));
            }
        }
        for (j = n / 2 + 1; j < n; j++)
            acb_conj(tab + j, tab + n - j);

        acb_clear(r);
        acb_clear(w);

        return GR_SUCCESS;
    }

    GR_TMP_INIT(t, ctx);
    fmpq_init(q);

    for (j = 0; j < n && status == GR_SUCCESS; j++)
    {
        fmpz_set_ui(fmpq_numref(q), 2 * j);
        fmpz_neg(fmpq_numref(q), fmpq_numref(q));
        fmpz_set_ui(fmpq_denref(q), n);
        fmpq_canonicalise(q);

        status |= gr_set_fmpq(t, q, ctx);
        status |= gr_exp_pi_i(GR_ENTRY(roots, j, sz), t, ctx);
    }

    fmpq_clear(q);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

static void
_gr_dft_precomp_zero(gr_dft_pre_t P)
{
    P->n = 0;
    P->depth = 0;
    P->alg = GR_DFT_ALG_NAIVE;
    P->flags = 0;
    P->ctx = NULL;
    P->real_ctx = NULL;
    P->roots = NULL;
    P->wtab = NULL;
    P->wclass = NULL;
    P->n1 = P->n2 = 0;
    P->pfa_a = P->pfa_b = 0;
    P->radices = NULL;
    P->num_radices = 0;
    P->conv_len = 0;
    P->bl_kern = NULL;
    P->bl_wtab = NULL;
    P->P1 = P->P2 = NULL;
    P->threads = NULL;
    P->num_threads = 0;
    P->serial_block = 0;
    P->nfixed_root_err = 0.0;
    P->bl_shifted = 0;
    P->stage_tab = NULL;
    P->stage_len = 0;
}



/* Phase 1 of plan construction: the layout. Resolves the algorithm
   and computes the complete decomposition (sub-plan layouts included)
   without touching any ring elements, so that cost and error bounds
   (gr_dft_precomp_nfixed_bound) can be evaluated before an internal
   working precision is chosen and before the root tables are
   computed. complex_mode indicates whether the plan will be realized
   in complex mode (it stands in for real_ctx != NULL in the automatic
   algorithm selection). The nfixed_root_err field is filled with a
   worst-case estimate for the doubling root-table construction, to be
   replaced by the actual bound when the plan is realized over a
   fixed-point context. */
static int
_gr_dft_layout(gr_dft_pre_t P, ulong n, int alg, int flags, int complex_mode)
{
    int status = GR_SUCCESS;
    int pow2, i;
    n_factor_t fac;

    _gr_dft_precomp_zero(P);

    if (n == 0)
        return GR_DOMAIN;

    pow2 = ((n & (n - 1)) == 0);

    n_factor_init(&fac);
    if (n >= 2 && !pow2)
        n_factor(&fac, n, 0);

    P->n = n;
    P->depth = pow2 ? (int) FLINT_BIT_COUNT(n) - 1 : -1;

    /* Resolve and validate the algorithm */
    if (alg == GR_DFT_ALG_AUTO)
    {
        if (pow2)
        {
            /* split-radix saves multiplications in complex mode and,
               like plain Cooley-Tukey, threads with full parallel
               width; it does not support scrambled ordering. (The
               Bailey algorithm remains available explicitly; in
               measurements up to n = 2^20 its cache blocking did not
               beat the plain transforms, serially or threaded.) */
            if (complex_mode && !(flags & GR_DFT_SCRAMBLED))
                alg = GR_DFT_ALG_SPLIT;
            else
                alg = GR_DFT_ALG_CT;
        }
        else
        {
            alg = (fac.num >= 2) ? GR_DFT_ALG_PFA : GR_DFT_ALG_MIXED;
        }
    }
    else
    {
        if (!pow2 && (alg == GR_DFT_ALG_CT || alg == GR_DFT_ALG_BAILEY ||
                      alg == GR_DFT_ALG_SPLIT))
            return GR_DOMAIN;
        if (alg == GR_DFT_ALG_PFA && (pow2 || fac.num < 2))
            return GR_DOMAIN;
        if (alg == GR_DFT_ALG_BLUESTEIN && (n < 3 || n % 2 == 0))
            return GR_DOMAIN;
        if (alg == GR_DFT_ALG_MIXED && pow2 && n >= 2)
        {
            /* supported, but requires the radix list */
            fac.num = 1;
            fac.p[0] = 2;
            fac.exp[0] = P->depth;
        }
        if (alg == GR_DFT_ALG_MIXED && n < 2)
            alg = GR_DFT_ALG_NAIVE;
    }
    if (alg == GR_DFT_ALG_BAILEY && P->depth < 2)
        alg = GR_DFT_ALG_CT;
    if (alg == GR_DFT_ALG_SPLIT && P->depth < 2)
        alg = GR_DFT_ALG_CT;
    P->alg = alg;

    /* scrambled ordering is only supported by NAIVE, CT and BAILEY
       (all power-of-two only, except NAIVE which supports it for
       power-of-two lengths) */
    if (alg != GR_DFT_ALG_NAIVE && alg != GR_DFT_ALG_CT &&
        alg != GR_DFT_ALG_BAILEY)
        flags &= ~GR_DFT_SCRAMBLED;
    if (!pow2)
        flags &= ~GR_DFT_SCRAMBLED;
    P->flags = flags;

    /* worst-case error estimate for a doubling-constructed root table
       (see _gr_dft_nfixed_root_err_bound; the constant is the modulus
       rounding of the most expensive complex multiplication) */
    P->nfixed_root_err = _gr_dft_nfixed_root_err_bound(n, 6.0 * 1.4142135623730951);

    /* Four-step decomposition n = n1 * n2 with sub-plans */
    if (alg == GR_DFT_ALG_BAILEY)
    {
        ulong n1, n2;

        n1 = UWORD(1) << (P->depth / 2);
        n2 = n / n1;
        P->n1 = n1;
        P->n2 = n2;

        P->P1 = flint_malloc(sizeof(gr_dft_pre_struct));
        P->P2 = flint_malloc(sizeof(gr_dft_pre_struct));

        /* The sub-plans always use plain Cooley-Tukey and natural
           ordering. */
        status |= _gr_dft_layout(P->P1, n1, GR_DFT_ALG_CT, 0, complex_mode);
        status |= _gr_dft_layout(P->P2, n2, GR_DFT_ALG_CT, 0, complex_mode);
    }

    /* Good-Thomas decomposition n = n1 * n2 into coprime parts, with
       n1 the full power of the smallest prime factor. There are no
       twiddle factors; the input is gathered by j = (n2 j1 + n1 j2)
       mod n and the output entry (k1, k2) is scattered to
       k = (pfa_a k1 + pfa_b k2) mod n where pfa_a = 1 mod n1,
       0 mod n2 and pfa_b = 0 mod n1, 1 mod n2. */
    if (status == GR_SUCCESS && alg == GR_DFT_ALG_PFA)
    {
        ulong n1, n2;

        n1 = fac.p[0];
        for (i = 1; i < fac.exp[0]; i++)
            n1 *= fac.p[0];
        n2 = n / n1;
        P->n1 = n1;
        P->n2 = n2;
        P->pfa_a = n2 * n_invmod(n2 % n1, n1);
        P->pfa_b = n1 * n_invmod(n1 % n2, n2);

        P->P1 = flint_malloc(sizeof(gr_dft_pre_struct));
        P->P2 = flint_malloc(sizeof(gr_dft_pre_struct));

        status |= _gr_dft_layout(P->P1, n1, GR_DFT_ALG_AUTO, 0, complex_mode);
        status |= _gr_dft_layout(P->P2, n2, GR_DFT_ALG_AUTO, 0, complex_mode);
    }

    /* Mixed-radix: flatten the factorization into a list of prime
       radices, one per recursion level. If n = p^e for a single prime
       p exceeding the direct-kernel cutoff, lay out a Bluestein plan
       for the length-p transforms; if it cannot be realized later
       (e.g. when the ring cannot provide a root of unity for the
       required power-of-two length), the realization falls back to
       the direct kernel. */
    if (status == GR_SUCCESS && alg == GR_DFT_ALG_MIXED)
    {
        slong num = 0, e;

        for (i = 0; i < fac.num; i++)
            num += fac.exp[i];

        P->radices = flint_malloc(num * sizeof(ulong));
        P->num_radices = num;

        num = 0;
        for (i = 0; i < fac.num; i++)
            for (e = 0; e < fac.exp[i]; e++)
                P->radices[num++] = fac.p[i];

        if (fac.num == 1 && fac.p[0] > GR_DFT_PRIME_CUTOFF)
        {
            P->P1 = flint_malloc(sizeof(gr_dft_pre_struct));
            status |= _gr_dft_layout(P->P1, fac.p[0],
                    GR_DFT_ALG_BLUESTEIN, 0, complex_mode);
        }
    }

    /* Bluestein: reduce to a cyclic convolution of power-of-two length
       conv_len >= 2n - 1, computed with a power-of-two sub-plan (with
       scrambled ordering; since the kernel is transformed by the same
       sub-plan, the orderings are consistent and the reordering passes
       are saved). */
    if (status == GR_SUCCESS && alg == GR_DFT_ALG_BLUESTEIN)
    {
        ulong N = UWORD(1) << FLINT_BIT_COUNT(2 * n - 2);
        P->conv_len = N;

        P->P1 = flint_malloc(sizeof(gr_dft_pre_struct));
        status |= _gr_dft_layout(P->P1, N, GR_DFT_ALG_AUTO,
                GR_DFT_SCRAMBLED, complex_mode);
    }

    return status;
}

static int _gr_dft_bluestein_realize_kernel(gr_dft_pre_t P, gr_ctx_t ctx);

/* Phase 2 of plan construction: bind the ring and compute the root
   tables, the complex mode tables, and the convolution kernels, over
   the layout produced by _gr_dft_layout. w is a principal n-th root
   of unity, or NULL when canonical is set and the table can be built
   directly (standard root over rings with exp, or the fixed-point
   builder). Recurses into the sub-plan layouts, feeding them the
   appropriate powers of w. */
static int
_gr_dft_realize(gr_dft_pre_t P, gr_srcptr w, int canonical,
        gr_ctx_struct * real_ctx, gr_ctx_t ctx,
        gr_srcptr rsrc, ulong rstride, double rsrc_err)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong n = P->n;
    int i;

    P->ctx = ctx;
    P->real_ctx = real_ctx;

    /* Canonical PFA plans need no root table: the algorithm performs
       no multiplications by roots (no twiddle factors), and the
       sub-plans, being canonical themselves, build their own (small)
       tables directly. (Plans with user-supplied roots keep the
       table, from which the sub-plan roots are derived.) */
    if (canonical && P->alg == GR_DFT_ALG_PFA)
    {
        status |= _gr_dft_realize(P->P1, NULL, 1, real_ctx, ctx,
                NULL, 0, 0.0);
        status |= _gr_dft_realize(P->P2, NULL, 1, real_ctx, ctx,
                NULL, 0, 0.0);
        return status;
    }

    /* Root table w^0, ..., w^(n-1) */
    P->roots = flint_malloc(n * sz);
    _gr_vec_init(P->roots, n, ctx);

    if (rsrc != NULL)
    {
        /* the roots of a sub-transform of a length-m parent, using
           the root w^(m/n), are strided entries of the parent table:
           (w^(m/n))^j = w^(j m/n); copy them instead of recomputing */
        ulong j;

        for (j = 0; j < n; j++)
            status |= gr_set(GR_ENTRY(P->roots, j, sz),
                    GR_ENTRY(rsrc, j * rstride, sz), ctx);

        P->nfixed_root_err = rsrc_err;
    }
    else if (_gr_dft_ctx_is_nfixed_complex(ctx))
    {
        /* fixed-point contexts: table built by binary-splitting
           products from a high-precision primitive root, with a
           computed ulp error bound; w is ignored (and may be NULL) */
        status = _gr_dft_nfixed_roots(P->roots, n, &P->nfixed_root_err, ctx);
    }
    else if (canonical && n > 2)
    {
        status = _gr_dft_roots_canonical(P->roots, n, ctx);
        if (status != GR_SUCCESS && w != NULL)
            status = _gr_vec_set_powers(P->roots, w, n, ctx);
    }
    else if (w != NULL)
    {
        status |= _gr_vec_set_powers(P->roots, w, n, ctx);
    }
    else if (n > 2)
    {
        status = GR_UNABLE;
    }
    else
    {
        /* n <= 2: the table is 1 (and -1), independent of w */
        if (n == 2)
            status |= gr_neg_one(GR_ENTRY(P->roots, 1, sz), ctx);
        status |= gr_one(P->roots, ctx);
    }

    /* Sanity checks for user-supplied roots: w should be a principal
       n-th root of unity. In particular w^(n/2) = -1 for even n, and
       w^(n/p) != 1 for every prime p dividing n. These may be
       inconclusive (T_UNKNOWN) over inexact rings, in which case we
       proceed. Canonical tables are correct by construction and are
       not checked (over approximate rings such as the fixed-point
       contexts, the checks could not certify anything anyway). */
    if (status == GR_SUCCESS && !canonical && n % 2 == 0)
    {
        if (gr_is_neg_one(GR_ENTRY(P->roots, n / 2, sz), ctx) == T_FALSE)
            status = GR_DOMAIN;
    }
    if (status == GR_SUCCESS && !canonical && n >= 2)
    {
        n_factor_t fac;
        n_factor_init(&fac);
        if ((n & (n - 1)) != 0)
            n_factor(&fac, n, 0);
        for (i = 0; i < fac.num; i++)
        {
            if (fac.p[i] != 2 &&
                gr_is_one(GR_ENTRY(P->roots, n / fac.p[i], sz), ctx) == T_TRUE)
                status = GR_DOMAIN;
        }
    }

    /* Complex mode: classify the special rotations (wclass), and, when
       profitable, precompute for each root c + d*i the table entries
       (c, d - c, d + c) as elements of the real ring, allowing
       multiplication by a generic root to be done with 3 instead of 4
       real multiplications (complex Karatsuba). */
    if (status == GR_SUCCESS && real_ctx != NULL)
    {
        slong rsz = real_ctx->sizeof_elem;
        slong rprec;
        int use_karatsuba;
        ulong j;

        /* Karatsuba trades one multiplication for three additions,
           which loses at low precision where a native complex
           multiplication is cheap. */
        use_karatsuba = !(gr_ctx_get_real_prec(&rprec, ctx) == GR_SUCCESS
                && rprec < GR_DFT_KARATSUBA_CUTOFF_PREC);

        /* fixed-point contexts always use the schoolbook complex
           product: the Karatsuba table entries d - c, d + c and the
           intermediate sum x_r + x_i can reach sqrt(2) in magnitude,
           which does not fit the representation */
        if (_gr_dft_ctx_is_nfixed(real_ctx))
            use_karatsuba = 0;

        /* Elements of ctx must consist of a real and an imaginary part
           stored contiguously as two elements of real_ctx. */
        if (sz != 2 * rsz)
        {
            status = GR_UNABLE;
        }
        else
        {
            P->wclass = flint_malloc(n);
            for (j = 0; j < n; j++)
                P->wclass[j] = GR_DFT_ROOT_GENERIC;

            if (use_karatsuba)
            {
                P->wtab = flint_malloc(3 * n * rsz);
                _gr_vec_init(P->wtab, 3 * n, real_ctx);

                for (j = 0; j < n; j++)
                {
                    gr_srcptr c = GR_ENTRY(P->roots, j, sz);
                    gr_srcptr d = GR_ENTRY(c, 1, rsz);
                    gr_ptr t = GR_ENTRY(P->wtab, 3 * j, rsz);

                    status |= gr_set(t, c, real_ctx);
                    status |= gr_sub(GR_ENTRY(t, 1, rsz), d, c, real_ctx);
                    status |= gr_add(GR_ENTRY(t, 2, rsz), d, c, real_ctx);
                }
            }

            /* Since complex mode always uses the standard root
               w = exp(-2 pi i / n), the classes of the powers at the
               quarter and eighth points (when these are integers) are
               known structurally: w^(n/2) = -1, w^(n/4) = -i,
               w^(3n/4) = i, and the odd multiples of n/8 have d = -c
               or d = c exactly. No value inspection is needed (nor
               could it certify these forms over inexact rings); the
               cheap multiplication formulas use at most the enclosure
               of c, so enclosures of the exact products are still
               obtained over ball rings. */
            if (n % 2 == 0)
                P->wclass[n / 2] = GR_DFT_ROOT_NEG_ONE;
            if (n % 4 == 0)
            {
                P->wclass[n / 4] = GR_DFT_ROOT_NEG_I;
                P->wclass[3 * n / 4] = GR_DFT_ROOT_I;
            }
            if (n % 8 == 0)
            {
                P->wclass[n / 8] = GR_DFT_ROOT_DPC_ZERO;
                P->wclass[3 * n / 8] = GR_DFT_ROOT_DMC_ZERO;
                P->wclass[5 * n / 8] = GR_DFT_ROOT_DPC_ZERO;
                P->wclass[7 * n / 8] = GR_DFT_ROOT_DMC_ZERO;
            }
        }
    }

    /* Pack the per-stage twiddle tables and drop the serial table:
       plain Cooley-Tukey and (complex-mode) split-radix read their
       twiddles in fixed strided sequences per stage, and packing
       these contiguously makes every table access during a transform
       a sequential walk, where the strided reads into a length-n
       table miss the cache once per entry as soon as the stride
       exceeds a line. The total size is the same as the serial table
       (sum 2^(s-1) = n - 1). Inverse transforms use the identities
       w^(-j rstep) = -stage_s[hm - j] (a reversed walk with the
       negation folded into the butterfly) and, for split,
       w^(-k r) = i w^((m/4-k) r), w^(-3k r) = -i w^(3(m/4-k) r)
       (reversed pair walks with free quarter rotations), so no
       inverse tables are needed. Plans using the complex Karatsuba
       tables keep the serial layout: they only arise at high
       precision, where transforms are compute-bound. */
    if (status == GR_SUCCESS && P->wtab == NULL && P->depth >= 2 &&
        (P->alg == GR_DFT_ALG_CT ||
         (P->alg == GR_DFT_ALG_SPLIT && real_ctx != NULL)))
    {
        ulong j, m;
        int s;

        if (P->alg == GR_DFT_ALG_CT)
        {
            P->stage_len = n - 1;
            P->stage_tab = flint_malloc(P->stage_len * sz);
            _gr_vec_init(P->stage_tab, P->stage_len, ctx);

            for (s = P->depth; s >= 1; s--)
            {
                gr_ptr base = GR_ENTRY(P->stage_tab, n - (UWORD(1) << s), sz);
                ulong hm = UWORD(1) << (s - 1), rstep = n >> s;

                for (j = 0; j < hm; j++)
                    status |= gr_set(GR_ENTRY(base, j, sz),
                            GR_ENTRY(P->roots, j * rstep, sz), ctx);
            }
        }
        else
        {
            P->stage_len = n - 2;
            P->stage_tab = flint_malloc(P->stage_len * sz);
            _gr_vec_init(P->stage_tab, P->stage_len, ctx);

            for (m = 4; m <= n; m *= 2)
            {
                gr_ptr base = GR_ENTRY(P->stage_tab, n - m, sz);
                ulong q = m / 4, rstep = n / m, k, e;

                for (k = 0; k < q; k++)
                {
                    status |= gr_set(GR_ENTRY(base, 2 * k, sz),
                            GR_ENTRY(P->roots, k * rstep, sz), ctx);
                    e = 3 * k * rstep;
                    if (e >= n)
                        e -= n;
                    status |= gr_set(GR_ENTRY(base, 2 * k + 1, sz),
                            GR_ENTRY(P->roots, e, sz), ctx);
                }
            }
        }

        _gr_vec_clear(P->roots, n, ctx);
        flint_free(P->roots);
        P->roots = NULL;

        if (P->wclass != NULL)
        {
            flint_free(P->wclass);
            P->wclass = NULL;
        }
    }

    /* Sub-plans (their root tables are strided copies of ours) */
    if (status == GR_SUCCESS &&
        (P->alg == GR_DFT_ALG_BAILEY || P->alg == GR_DFT_ALG_PFA))
    {
        status |= _gr_dft_realize(P->P1, NULL, canonical, real_ctx, ctx,
                P->roots, P->n2, P->nfixed_root_err);
        status |= _gr_dft_realize(P->P2, NULL, canonical, real_ctx, ctx,
                P->roots, P->n1, P->nfixed_root_err);
    }

    if (status == GR_SUCCESS && P->alg == GR_DFT_ALG_MIXED && P->P1 != NULL)
    {
        ulong p = P->P1->n;
        int bstatus = _gr_dft_realize(P->P1, NULL, canonical, real_ctx, ctx,
                P->roots, n / p, P->nfixed_root_err);

        if (bstatus != GR_SUCCESS)
        {
            /* fall back to the direct kernel */
            gr_dft_precomp_clear(P->P1);
            flint_free(P->P1);
            P->P1 = NULL;
        }
    }

    if (status == GR_SUCCESS && P->alg == GR_DFT_ALG_BLUESTEIN)
    {
        /* the power-of-two sub-plan needs its own canonical root */
        if (_gr_dft_ctx_is_nfixed_complex(ctx))
        {
            status |= _gr_dft_realize(P->P1, NULL, 1, real_ctx, ctx,
                    NULL, 0, 0.0);
        }
        else
        {
            gr_ptr w2;
            GR_TMP_INIT(w2, ctx);

            status |= gr_dft_default_root(w2, P->conv_len, ctx);
            if (status == GR_SUCCESS)
                status |= _gr_dft_realize(P->P1, w2, 1, real_ctx, ctx,
                        NULL, 0, 0.0);

            GR_TMP_CLEAR(w2, ctx);
        }

        if (status == GR_SUCCESS)
            status |= _gr_dft_bluestein_realize_kernel(P, ctx);
    }

    return status;
}

static int
_gr_dft_precomp_init_root(gr_dft_pre_t P, gr_srcptr w, ulong n,
        int alg, int flags, int canonical, gr_ctx_struct * real_ctx, gr_ctx_t ctx)
{
    int status = _gr_dft_layout(P, n, alg, flags, real_ctx != NULL);

    if (status == GR_SUCCESS)
        status |= _gr_dft_realize(P, w, canonical, real_ctx, ctx, NULL, 0, 0.0);

    return status;
}

/* Public layout / realize interface (canonical roots), allowing cost
   and error bounds to be inspected before an internal precision is
   chosen and the tables are computed. */
int
_gr_dft_precomp_init_layout(gr_dft_pre_t P, ulong n, int alg, int flags,
        int complex_mode)
{
    int status = _gr_dft_layout(P, n, alg, flags, complex_mode);

    if (status != GR_SUCCESS)
    {
        gr_dft_precomp_clear(P);
        _gr_dft_precomp_zero(P);
    }

    return status;
}

int
_gr_dft_precomp_realize(gr_dft_pre_t P, gr_ctx_struct * real_ctx, gr_ctx_t ctx)
{
    int status = _gr_dft_realize(P, NULL, 1, real_ctx, ctx, NULL, 0, 0.0);

    if (status != GR_SUCCESS)
    {
        gr_dft_precomp_clear(P);
        _gr_dft_precomp_zero(P);
    }

    return status;
}

/* Bluestein kernel: writing t = (n+1)/2 = 1/2 mod n, we have
   jk = t (j^2 + k^2 - (j-k)^2) mod n, so with the chirp
   c_j = w^(t j^2) (entries of the root table),

       X_k = c_k sum_j (c_j x_j) w^(-t (k-j)^2),

   a cyclic convolution of length n of the chirped input with the
   n-periodic even kernel h_j = w^(-t j^2). The convolution is embedded
   into one of power-of-two length conv_len >= 2n - 1 and evaluated by
   a forward transform, a pointwise multiplication by the transformed
   kernel (precomputed here, with the 1/conv_len scaling of the inverse
   transform folded in), and an inverse transform. */
static int
_gr_dft_bluestein_realize_kernel(gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong n = P->n, N = P->conv_len, j, t2, e, ninv;
    gr_ptr h;

    /* kernel h, embedded with wraparound: h'_j = h_j for 0 <= j < n,
       h'_(N-j) = h_j for 1 <= j < n, zero elsewhere */
    h = flint_malloc(N * sz);
    _gr_vec_init(h, N, ctx);
    status |= _gr_vec_zero(h, N, ctx);

    t2 = (n + 1) / 2;
    ninv = n_preinvert_limb(n);

    for (j = 0; j < n; j++)
    {
        e = n_mulmod2_preinv(j, j, n, ninv);
        e = n_mulmod2_preinv(e, t2, n, ninv);
        e = (e == 0) ? 0 : n - e;   /* h_j = w^(-t j^2) */

        status |= gr_set(GR_ENTRY(h, j, sz), GR_ENTRY(P->roots, e, sz), ctx);
        if (j > 0)
            status |= gr_set(GR_ENTRY(h, N - j, sz), GR_ENTRY(h, j, sz), ctx);
    }

    if (_gr_dft_ctx_is_nfixed_complex(ctx))
    {
        /* Over the fixed-point contexts the folded scaling must be
           applied before the transform of the kernel, so that its
           intermediate values stay below 1; the extra factor 1/2
           keeps the sqrt(2) rotation intermediates in range as well
           (|h| <= 1/(2N) gives a transform peak of at most
           sqrt(2)/2), and is undone by doubling the outputs of the
           main transform (P->bl_shifted). The division truncates with
           error below 1 ulp per part. */
        status |= _gr_vec_div_scalar_ui(h, h, N, 2 * N, ctx);
        P->bl_shifted = 1;
    }

    P->bl_kern = flint_malloc(N * sz);
    _gr_vec_init(P->bl_kern, N, ctx);

    status |= _gr_dft_precomp_raw(P->bl_kern, h, 0, P->P1, ctx);

    _gr_vec_clear(h, N, ctx);
    flint_free(h);

    /* fold in the 1/conv_len scaling of the inverse sub-transform
       (already included, together with the extra 1/2, over the
       fixed-point contexts) */
    if (!P->bl_shifted)
        status |= _gr_vec_div_scalar_ui(P->bl_kern, P->bl_kern, N, N, ctx);

    /* complex Karatsuba multiplication tables for the kernel entries,
       matching the enablement decision of the main table */
    if (status == GR_SUCCESS && P->wtab != NULL)
    {
        gr_ctx_struct * rctx = P->real_ctx;
        slong rsz = rctx->sizeof_elem;

        P->bl_wtab = flint_malloc(3 * N * rsz);
        _gr_vec_init(P->bl_wtab, 3 * N, rctx);

        for (j = 0; j < N; j++)
        {
            gr_srcptr c = GR_ENTRY(P->bl_kern, j, sz);
            gr_srcptr d = GR_ENTRY(c, 1, rsz);
            gr_ptr u = GR_ENTRY(P->bl_wtab, 3 * j, rsz);

            status |= gr_set(u, c, rctx);
            status |= gr_sub(GR_ENTRY(u, 1, rsz), d, c, rctx);
            status |= gr_add(GR_ENTRY(u, 2, rsz), d, c, rctx);
        }
    }

    return status;
}

static int
_gr_dft_precomp_init_root_checked(gr_dft_pre_t P, gr_srcptr w, ulong n,
        int alg, int flags, int canonical, gr_ctx_struct * real_ctx, gr_ctx_t ctx)
{
    int status = _gr_dft_precomp_init_root(P, w, n, alg, flags, canonical, real_ctx, ctx);

    if (status != GR_SUCCESS)
    {
        gr_dft_precomp_clear(P);
        _gr_dft_precomp_zero(P);
    }

    return status;
}

int
gr_dft_precomp_init_root(gr_dft_pre_t P, gr_srcptr w, ulong n,
        int alg, int flags, gr_ctx_t ctx)
{
    return _gr_dft_precomp_init_root_checked(P, w, n, alg, flags, 0, NULL, ctx);
}

int
gr_dft_precomp_init(gr_dft_pre_t P, ulong n, int alg, int flags, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr w;

    if (_gr_dft_ctx_is_nfixed_complex(ctx))
        return _gr_dft_precomp_init_root_checked(P, NULL, n, alg, flags, 1, NULL, ctx);

    GR_TMP_INIT(w, ctx);

    status |= gr_dft_default_root(w, n, ctx);

    if (status == GR_SUCCESS)
        status |= _gr_dft_precomp_init_root_checked(P, w, n, alg, flags, 1, NULL, ctx);
    else
        _gr_dft_precomp_zero(P);

    GR_TMP_CLEAR(w, ctx);
    return status;
}

int
gr_dft_precomp_init_karatsuba(gr_dft_pre_t P, ulong n, int alg, int flags,
        gr_ctx_t real_ctx, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr w;

    if (_gr_dft_ctx_is_nfixed_complex(ctx))
        return _gr_dft_precomp_init_root_checked(P, NULL, n, alg, flags, 1, real_ctx, ctx);

    GR_TMP_INIT(w, ctx);

    status |= gr_dft_default_root(w, n, ctx);

    if (status == GR_SUCCESS)
        status |= _gr_dft_precomp_init_root_checked(P, w, n, alg, flags, 1, real_ctx, ctx);
    else
        _gr_dft_precomp_zero(P);

    GR_TMP_CLEAR(w, ctx);
    return status;
}

void
gr_dft_precomp_clear(gr_dft_pre_t P)
{
    if (P->roots != NULL)
    {
        _gr_vec_clear(P->roots, P->n, P->ctx);
        flint_free(P->roots);
        P->roots = NULL;
    }

    if (P->stage_tab != NULL)
    {
        _gr_vec_clear(P->stage_tab, P->stage_len, P->ctx);
        flint_free(P->stage_tab);
        P->stage_tab = NULL;
    }

    if (P->wtab != NULL)
    {
        _gr_vec_clear(P->wtab, 3 * P->n, P->real_ctx);
        flint_free(P->wtab);
        P->wtab = NULL;
    }

    if (P->wclass != NULL)
    {
        flint_free(P->wclass);
        P->wclass = NULL;
    }

    if (P->radices != NULL)
    {
        flint_free(P->radices);
        P->radices = NULL;
    }

    if (P->bl_kern != NULL)
    {
        _gr_vec_clear(P->bl_kern, P->conv_len, P->ctx);
        flint_free(P->bl_kern);
        P->bl_kern = NULL;
    }

    if (P->bl_wtab != NULL)
    {
        _gr_vec_clear(P->bl_wtab, 3 * P->conv_len, P->real_ctx);
        flint_free(P->bl_wtab);
        P->bl_wtab = NULL;
    }

    if (P->P1 != NULL)
    {
        gr_dft_precomp_clear(P->P1);
        flint_free(P->P1);
        P->P1 = NULL;
    }

    if (P->P2 != NULL)
    {
        gr_dft_precomp_clear(P->P2);
        flint_free(P->P2);
        P->P2 = NULL;
    }
}

/* Attach worker threads to the plan (and, recursively, its sub-plans).
   The handles are borrowed for the lifetime of the plan: they are used
   by every transform, are not given back to the thread pool by
   gr_dft_precomp_clear, and remain the responsibility of the caller.
   Passing NULL, 0 detaches, restoring the default behavior of
   requesting workers from the global thread pool during transforms. */
void
gr_dft_precomp_set_threads(gr_dft_pre_t P, thread_pool_handle * threads,
        slong num_threads)
{
    P->threads = threads;
    P->num_threads = threads ? num_threads : 0;

    if (P->P1 != NULL)
        gr_dft_precomp_set_threads(P->P1, threads, num_threads);
    if (P->P2 != NULL)
        gr_dft_precomp_set_threads(P->P2, threads, num_threads);
}

/* Set the threading granularity (and that of the sub-plans): work is
   divided into at most n / serial_block chunks, so in particular
   transforms with n < 2 * serial_block run serially. Zero restores the
   default GR_DFT_SERIAL_BLOCK_DEFAULT. */
void
gr_dft_precomp_set_serial_block(gr_dft_pre_t P, slong serial_block)
{
    P->serial_block = serial_block;

    if (P->P1 != NULL)
        gr_dft_precomp_set_serial_block(P->P1, serial_block);
    if (P->P2 != NULL)
        gr_dft_precomp_set_serial_block(P->P2, serial_block);
}

void
gr_dft_precomp_output_perm(ulong * perm, const gr_dft_pre_t P)
{
    ulong i, n = P->n;

    if (!(P->flags & GR_DFT_SCRAMBLED) || P->alg == GR_DFT_ALG_SPLIT)
    {
        for (i = 0; i < n; i++)
            perm[i] = i;
    }
    else if (P->alg == GR_DFT_ALG_BAILEY)
    {
        /* entry i = k2 * n1 + k1 holds coefficient k2 + n2 * k1 */
        ulong n1 = P->n1, n2 = P->n2;

        for (i = 0; i < n; i++)
            perm[i] = (i / n1) + n2 * (i % n1);
    }
    else
    {
        for (i = 0; i < n; i++)
            perm[i] = n_revbin(i, P->depth);
    }
}
