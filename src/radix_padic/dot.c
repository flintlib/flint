/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"
#include "radix_padic.h"
#include "gr.h"
#include "longlong.h"
#include "mpn_extras.h"

int
radix_padic_dot_strided_naive(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1, slong stride1,
    const radix_padic_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
{
    if (len <= 0)
    {
        if (initial == NULL)
            return radix_padic_zero(res, ctx);
        return radix_padic_set(res, initial, ctx);
    }

    int status = GR_SUCCESS;

    radix_padic_t t;
    radix_padic_init(t, ctx);

    if (initial == NULL)
    {
        status |= radix_padic_mul(res, vec1, vec2, ctx);
    }
    else
    {
        if (subtract)
            status |= radix_padic_neg(res, initial, ctx);
        else
            status |= radix_padic_set(res, initial, ctx);

        status |= radix_padic_mul(t, vec1, vec2, ctx);
        status |= radix_padic_add(res, res, t, ctx);
    }

    for (slong i = 1; i < len; i++)
    {
        status |= radix_padic_mul(t, vec1 + i * stride1, vec2 + i * stride2, ctx);
        status |= radix_padic_add(res, res, t, ctx);
    }

    if (subtract)
        status |= radix_padic_neg(res, res, ctx);

    radix_padic_clear(t, ctx);
    return status;
}

/*
    Optimized strided dot product for radix_padic:

        res = (initial or 0)  +/-  sum_{k=0}^{len-1} a_k * b_k,

    with a_k = vec1[k*stride1], b_k = vec2[k*stride2], and the sign chosen by
    `subtract` (the initial term is always added; only the sum is negated).

    Two passes over the data:

    Pass 1 finds the minimum valuation vmin over the nonzero contributions, the
    smallest absolute precision N implied by the terms / initial and the context's
    relative and absolute precision (so the window of digit positions to compute is
    [vmin, N)), and an upper bound on the result size. If every contribution is
    exact and the context imposes no precision cap, N is infinite and the result is
    returned exactly.

    Pass 2 evaluates the sum into deferred-carry accumulators. A term a_k*b_k sits
    at digit offset d = (v_{a_k}+v_{b_k}) - vmin; write d = o*e + r with e the
    digits per limb. Terms sharing the sub-limb shift r go into one accumulator, at
    limb offset o, so the only per-limb shift by p^r happens once per residue at the
    end rather than once per term. Within an accumulator the limb cross-products are
    summed into 3-word slots with umul_ppmm + add_sssaaaaaa, deferring every divrem
    by B = p^e until finalisation (exactly the slot scheme of radix_mulmid). Positive
    and negative contributions go to separate accumulators so each stays unsigned.

    Pass 3 carry-propagates each accumulator into a normalized integer, shifts it by
    its residue r, and sums the residues; the negative side is subtracted from the
    positive. The resulting unit, valuation vmin and precision N are handed to
    _radix_padic_finalize, which canonicalises (handling low-digit cancellation) and
    reduces / keeps-exact as appropriate.
*/

#define DOT_INF WORD_MAX
#define DOT_MIN2(a, b) ((a) < (b) ? (a) : (b))

#define DOT_MULMID_CUTOFF 140

/* carry-propagate a 3-words-per-slot accumulator into a normalized radix_integer */
static void
_radix_dot_normalize_acc(radix_integer_t T, nn_srcptr acc, slong nslots,
    const radix_t radix)
{
    nmod_t B = radix->B;
    nn_ptr d = radix_integer_fit_limbs(T, nslots + 3, radix);
    ulong cy[3] = { 0, 0, 0 };
    slong s, outn = 0;

    if (B.norm == 0)
    {
        for (s = 0; s < nslots; s++)
        {
            add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0],
                acc[3 * s + 2], acc[3 * s + 1], acc[3 * s]);
            d[outn++] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, B.n, B.ninv);
        }
        while (cy[0] | cy[1] | cy[2])
            d[outn++] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, B.n, B.ninv);
    }
    else
    {
        for (s = 0; s < nslots; s++)
        {
            add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0],
                acc[3 * s + 2], acc[3 * s + 1], acc[3 * s]);
            d[outn++] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, B.n, B.ninv, B.norm);
        }
        while (cy[0] | cy[1] | cy[2])
            d[outn++] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, B.n, B.ninv, B.norm);
    }

    while (outn > 0 && d[outn - 1] == 0)
        outn--;
    T->size = outn;
}

int
radix_padic_dot_strided_delayed(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1, slong stride1,
    const radix_padic_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong prec_abs = RADIX_PADIC_CTX_PREC_ABS(ctx);
    slong prec_rel = RADIX_PADIC_CTX_PREC_REL(ctx);

    slong k, vmin, maxhi, Nterms, N, Wdig, nslots;
    int have_value;
    const radix_padic_struct * a;
    const radix_padic_struct * b;
    int init_nonzero;

    nn_ptr * accp;
    nn_ptr * accn;
    radix_integer_t Mp, Mn, T;
    radix_integer_t big, Ptmp;   /* running total + scratch for the radix_mulmid path */

    if (len <= 0)
    {
        if (initial == NULL)
            return radix_padic_zero(res, ctx);
        return radix_padic_set(res, initial, ctx);
    }

    /* ---------------- Pass 1: window and precision ---------------- */
    vmin = DOT_INF;
    maxhi = -DOT_INF;            /* highest absolute digit position touched + slack */
    Nterms = DOT_INF;            /* min absolute precision from terms/initial */
    have_value = 0;
    init_nonzero = (initial != NULL && initial->u.size != 0);

    if (init_nonzero)
    {
        slong ni = FLINT_ABS(initial->u.size);
        vmin = DOT_MIN2(vmin, initial->v);
        maxhi = FLINT_MAX(maxhi, initial->v + ni * e);
        have_value = 1;
    }
    if (initial != NULL && initial->N != RADIX_PADIC_EXACT)
        Nterms = DOT_MIN2(Nterms, initial->N);

    for (k = 0; k < len; k++)
    {
        slong va, vb, Na, Nb, vk, Nk;

        a = vec1 + k * stride1;
        b = vec2 + k * stride2;

        /* exact zero annihilates exactly: contributes nothing at all */
        if ((a->u.size == 0 && a->N == RADIX_PADIC_EXACT)
         || (b->u.size == 0 && b->N == RADIX_PADIC_EXACT))
            continue;

        va = a->v; vb = b->v; Na = a->N; Nb = b->N;
        vk = va + vb;

        Nk = RADIX_PADIC_EXACT;
        if (Na != RADIX_PADIC_EXACT) Nk = DOT_MIN2(Nk, vb + Na);
        if (Nb != RADIX_PADIC_EXACT) Nk = DOT_MIN2(Nk, va + Nb);
        if (Nk != RADIX_PADIC_EXACT) Nterms = DOT_MIN2(Nterms, Nk);

        if (a->u.size != 0 && b->u.size != 0)
        {
            slong na = FLINT_ABS(a->u.size), nb = FLINT_ABS(b->u.size);
            vmin = DOT_MIN2(vmin, vk);
            maxhi = FLINT_MAX(maxhi, vk + (na + nb) * e);
            have_value = 1;
        }
    }

    /* no actual value anywhere: result is 0 + O(p^Nterms) (or exact 0) */
    if (!have_value)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = Nterms;            /* EXACT if everything cancelled exactly */
        return _radix_padic_finalize(res, ctx);
    }

    /* absolute precision to compute: terms, context absolute, context relative.
       Using vmin (not the post-cancellation valuation) for the relative cap is the
       accepted slightly-suboptimal choice. */
    N = Nterms;
    N = DOT_MIN2(N, prec_abs);
    if (prec_rel != RADIX_PADIC_PREC_INF)
        N = DOT_MIN2(N, vmin + prec_rel);

    if (N != DOT_INF && N <= vmin)
    {
        /* whole result below the horizon */
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = N;
        return _radix_padic_finalize(res, ctx);
    }

    /* number of result limbs (slots) to cover */
    {
        slong nslots_ext = (maxhi - vmin + e - 1) / e + 1;
        if (N == DOT_INF)
            nslots = nslots_ext;
        else
        {
            slong nslots_prec = (N - vmin + e - 1) / e + 2;
            nslots = DOT_MIN2(nslots_prec, nslots_ext);
        }
        if (nslots < 1)
            nslots = 1;
    }
    Wdig = (N == DOT_INF) ? DOT_INF : (N - vmin);

    /* ---------------- Pass 2: accumulate ---------------- */
    accp = flint_calloc(e, sizeof(nn_ptr));
    accn = flint_calloc(e, sizeof(nn_ptr));
    radix_integer_init(big, radix);     /* large products go here (signed) */
    radix_integer_init(Ptmp, radix);

    /* initial term (always added; its own unit sign chooses the side) */
    if (init_nonzero)
    {
        slong ni = FLINT_ABS(initial->u.size);
        nn_srcptr ui = initial->u.d;
        slong d = initial->v - vmin;        /* >= 0 */
        slong r = d % e, o = d / e, i;
        nn_ptr acc;

        if (initial->u.size < 0)
        {
            if (accn[r] == NULL) accn[r] = flint_calloc(3 * nslots, sizeof(ulong));
            acc = accn[r];
        }
        else
        {
            if (accp[r] == NULL) accp[r] = flint_calloc(3 * nslots, sizeof(ulong));
            acc = accp[r];
        }

        for (i = 0; i < ni; i++)
        {
            slong slot = o + i;
            if (slot >= nslots)
                break;
            acc[3 * slot] = ui[i];
        }
    }

    for (k = 0; k < len; k++)
    {
        slong va, vb, vk, d, r, o, na, nb, m, mmax;
        slong win, hi, an_t, bn_t;
        nn_srcptr ua, ub;
        nn_ptr acc;
        int neg;

        a = vec1 + k * stride1;
        b = vec2 + k * stride2;

        if (a->u.size == 0 || b->u.size == 0)
            continue;

        va = a->v; vb = b->v; vk = va + vb;
        d = vk - vmin;                              /* >= 0 */
        if (Wdig != DOT_INF && d >= Wdig)
            continue;                               /* term entirely below horizon */

        r = d % e; o = d / e;
        if (o >= nslots)
            continue;

        na = FLINT_ABS(a->u.size); nb = FLINT_ABS(b->u.size);
        ua = a->u.d; ub = b->u.d;
        neg = (a->u.size < 0) ^ (b->u.size < 0) ^ (subtract != 0);

        /*
            Window for this term: the result keeps limbs in slots [o, nslots), so we
            need product limbs [0, hi) with hi = min(na+nb, nslots - o). Only operand
            limbs below hi affect those, giving the truncated operand sizes an_t, bn_t.
            (Inexact terms may truncate hi < na+nb; the dropped high limbs land at
            slots >= nslots, beyond the precision window. Exact terms are never
            truncated, since nslots covers their full extent.)
        */
        win = nslots - o;
        hi = FLINT_MIN(na + nb, win);
        an_t = FLINT_MIN(na, hi);
        bn_t = FLINT_MIN(nb, hi);

        if (FLINT_MIN(an_t, bn_t) >= DOT_MULMID_CUTOFF && hi >= DOT_MULMID_CUTOFF)
        {
            /*
                Large product: out of the schoolbook regime. Multiply it out with
                radix_mulmid (which picks FFT/KS here) and add the normalized product,
                shifted into place by the full digit offset d, to the running total.
            */
            nn_ptr pd = radix_integer_fit_limbs(Ptmp, hi, radix);
            radix_mulmid(pd, ua, an_t, ub, bn_t, 0, hi, radix);
            Ptmp->size = hi;
            while (Ptmp->size > 0 && pd[Ptmp->size - 1] == 0)
                Ptmp->size--;

            if (neg)
                radix_integer_sublsh(big, big, Ptmp, d, radix);   /* big -= Ptmp << d */
            else
                radix_integer_addlsh(big, big, Ptmp, d, radix);   /* big += Ptmp << d */

            continue;
        }

        /* Small product: deferred per-residue accumulator. */
        if (neg)
        {
            if (accn[r] == NULL) accn[r] = flint_calloc(3 * nslots, sizeof(ulong));
            acc = accn[r];
        }
        else
        {
            if (accp[r] == NULL) accp[r] = flint_calloc(3 * nslots, sizeof(ulong));
            acc = accp[r];
        }

        /*
            Slot-outer / index-inner ordering, as in radix_mulmid_classical: the
            outer loop walks the product output positions m (slot = o + m) and the
            inner loop sums the pairs ua[ii]*ub[m-ii] landing in that slot. The
            three accumulator words for the slot are held in registers (cy0,cy1,cy2)
            across the whole inner loop, with a single load and store per slot.
            Slots at index >= nslots are beyond the precision window and skipped.
        */
        mmax = na + nb - 2;                     /* top product position (na,nb >= 1) */
        if (mmax > nslots - 1 - o)
            mmax = nslots - 1 - o;              /* clamp to the window */

        for (m = 0; m <= mmax; m++)
        {
            slong slot = o + m;
            slong iilo = (m >= nb) ? (m - nb + 1) : 0;     /* max(0, m-(nb-1)) */
            slong iihi = (m < na) ? m : (na - 1);          /* min(na-1, m)     */
            slong ii;
            ulong cy0 = acc[3 * slot];
            ulong cy1 = acc[3 * slot + 1];
            ulong cy2 = acc[3 * slot + 2];

            for (ii = iilo; ii <= iihi; ii++)
            {
                ulong hilo, lo;
                umul_ppmm(hilo, lo, ua[ii], ub[m - ii]);
                add_sssaaaaaa(cy2, cy1, cy0, cy2, cy1, cy0, 0, hilo, lo);
            }

            acc[3 * slot]     = cy0;
            acc[3 * slot + 1] = cy1;
            acc[3 * slot + 2] = cy2;
        }
    }

    /* ---------------- Pass 3: finalise accumulators ---------------- */
    radix_integer_init(Mp, radix);
    radix_integer_init(Mn, radix);
    radix_integer_init(T, radix);

    {
        slong r;
        for (r = 0; r < e; r++)
        {
            if (accp[r] != NULL)
            {
                _radix_dot_normalize_acc(T, accp[r], nslots, radix);
                radix_integer_addlsh(Mp, Mp, T, r, radix);   /* Mp += T << r */
                flint_free(accp[r]);
            }
            if (accn[r] != NULL)
            {
                _radix_dot_normalize_acc(T, accn[r], nslots, radix);
                radix_integer_addlsh(Mn, Mn, T, r, radix);
                flint_free(accn[r]);
            }
        }
    }
    flint_free(accp);
    flint_free(accn);

    radix_integer_sub(&res->u, Mp, Mn, radix);       /* M = Mp - Mn (signed) */
    if (big->size != 0)
        radix_integer_add(&res->u, &res->u, big, radix);   /* + large-product total */
    res->v = vmin;
    res->N = N;

    radix_integer_clear(Mp, radix);
    radix_integer_clear(Mn, radix);
    radix_integer_clear(T, radix);
    radix_integer_clear(big, radix);
    radix_integer_clear(Ptmp, radix);

    return _radix_padic_finalize(res, ctx);
}

int
radix_padic_dot_strided(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1, slong stride1,
    const radix_padic_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
{
    if (len <= 7)
        return radix_padic_dot_strided_naive(res, initial, subtract, vec1, stride1, vec2, stride2, len, ctx);
    else
        return radix_padic_dot_strided_delayed(res, initial, subtract, vec1, stride1, vec2, stride2, len, ctx);
}

int
radix_padic_dot(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1,
    const radix_padic_struct * vec2, slong len, gr_ctx_t ctx)
{
    return radix_padic_dot_strided(res, initial, subtract, vec1, 1, vec2, 1, len, ctx);
}

int
radix_padic_dot_rev(radix_padic_t res, const radix_padic_t initial,
    int subtract, const radix_padic_struct * vec1,
    const radix_padic_struct * vec2, slong len, gr_ctx_t ctx)
{
    return radix_padic_dot_strided(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ctx);
}

