#include <string.h>
/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"
#include "padic_radix.h"
#include "impl.h"
int _padic_radix_addmul_sub(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, int sub, gr_ctx_t ctx);
#include "gr.h"
#include "longlong.h"
#include "mpn_extras.h"

int
padic_radix_dot_strided_naive(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1, slong stride1,
    const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
{
    if (len <= 0)
    {
        if (initial == NULL)
            return padic_radix_zero(res, ctx);
        return padic_radix_set(res, initial, ctx);
    }

    int status = GR_SUCCESS;

    /* mul-then-add steps as addmul (same result) */
    if (initial == NULL)
    {
        status |= padic_radix_mul(res, vec1, vec2, ctx);
    }
    else
    {
        if (subtract)
            status |= padic_radix_neg(res, initial, ctx);
        else
            status |= padic_radix_set(res, initial, ctx);

        status |= _padic_radix_addmul_sub(res, vec1, vec2, 0, ctx);
    }

    for (slong i = 1; i < len; i++)
        status |= _padic_radix_addmul_sub(res, vec1 + i * stride1, vec2 + i * stride2, 0, ctx);

    if (subtract)
        status |= padic_radix_neg(res, res, ctx);

    return status;
}

/*
    Optimized strided dot product for padic_radix:

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
    _padic_radix_finalize, which canonicalises (handling low-digit cancellation) and
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

/* res = P - N for normalised digit vectors (no heap temporaries) */
static void
_dot_set_diff(radix_integer_t res, nn_srcptr P, slong pn, nn_srcptr Nd, slong nn, const radix_t radix)
{
    int sgn;
    nn_srcptr a, b;
    slong an, bn, rn;
    nn_ptr d;

    if (pn != nn)
        sgn = (pn > nn) ? 1 : -1;
    else
        sgn = (pn == 0) ? 0 : mpn_cmp(P, Nd, pn);
    if (sgn == 0)
    {
        res->size = 0;
        return;
    }
    if (sgn > 0) { a = P; an = pn; b = Nd; bn = nn; }
    else         { a = Nd; an = nn; b = P; bn = pn; }
    d = radix_integer_fit_limbs(res, an, radix);
    if (bn == 0)
        flint_mpn_copyi(d, a, an);
    else
        radix_sub(d, a, an, b, bn, radix);
    rn = an;
    while (rn > 0 && d[rn - 1] == 0)
        rn--;
    res->size = (sgn > 0) ? rn : -rn;
}

/* q = (hi B^? ...): one step of dividing by B, hi < B */
static inline ulong
_dot_divrem2(ulong * r, ulong hi, ulong lo, nmod_t B)
{
    ulong q, rr;
    if (B.norm == 0)
    {
        udiv_qrnnd_preinv(q, rr, hi, lo, B.n, B.ninv);
        *r = rr;
    }
    else
    {
        udiv_qrnnd_preinv(q, rr, (hi << B.norm) | (lo >> (FLINT_BITS - B.norm)),
                          lo << B.norm, B.n << B.norm, B.ninv);
        *r = rr >> B.norm;
    }
    return q;
}

/*
    Fused residue combination: comb[5 s .. 5 s + 4] += acc[3 s .. 3 s + 2] * p^r
    for every slot s, then one carry-propagating normalisation in radix B.
    Replaces one normalisation plus one shifted addition per residue.
*/
static void
_dot_comb_add(nn_ptr comb, nn_srcptr acc, ulong pr, slong nslots)
{
    slong s;
    for (s = 0; s < nslots; s++)
    {
        ulong h0, l0, h1, l1, h2, l2, w0, w1, w2, w3, c;
        nn_ptr t = comb + 5 * s;
        umul_ppmm(h0, l0, acc[3 * s], pr);
        umul_ppmm(h1, l1, acc[3 * s + 1], pr);
        umul_ppmm(h2, l2, acc[3 * s + 2], pr);
        /* (h2 h1 h0) shifted one word, plus (l2 l1 l0) */
        w0 = l0;
        add_ssaaaa(c, w1, 0, h0, 0, l1);
        add_ssaaaa(w3, w2, h2, h1, 0, l2);
        add_ssaaaa(w3, w2, w3, w2, 0, c);
        /* t += (w3 w2 w1 w0) */
        add_sssaaaaaa(c, t[1], t[0], 0, t[1], t[0], 0, w1, w0);
        add_sssaaaaaa(t[4], t[3], t[2], t[4], t[3], t[2], 0, w3, w2);
        add_sssaaaaaa(t[4], t[3], t[2], t[4], t[3], t[2], 0, 0, c);
    }
}

static slong
_dot_comb_normalize(nn_ptr d, nn_srcptr comb, slong nslots, const radix_t radix)
{
    nmod_t B = radix->B;
    ulong cy[5] = { 0, 0, 0, 0, 0 }, c;
    slong s, i, outn = 0;
    for (s = 0; ; s++)
    {
        if (s < nslots)
        {
            nn_srcptr t = comb + 5 * s;
            add_sssaaaaaa(c, cy[1], cy[0], 0, cy[1], cy[0], 0, t[1], t[0]);
            add_sssaaaaaa(cy[4], cy[3], cy[2], cy[4], cy[3], cy[2], t[4], t[3], t[2]);
            add_sssaaaaaa(cy[4], cy[3], cy[2], cy[4], cy[3], cy[2], 0, 0, c);
        }
        else if ((cy[0] | cy[1] | cy[2] | cy[3] | cy[4]) == 0)
            break;
        {
            ulong r = 0;
            i = 4;
            while (i > 0 && cy[i] == 0)
                i--;
            if (i == 0 || (i == 1 && cy[1] < B.n))
            {
                /* < B^2: one step */
                if (i == 0)
                    cy[0] = _dot_divrem2(&r, 0, cy[0], B);
                else
                {
                    cy[0] = _dot_divrem2(&r, cy[1], cy[0], B);
                    cy[1] = 0;
                }
            }
            else
                for ( ; i >= 0; i--)
                    cy[i] = _dot_divrem2(&r, r, cy[i], B);
            d[outn++] = r;
        }
    }
    while (outn > 0 && d[outn - 1] == 0)
        outn--;
    return outn;
}

int
padic_radix_dot_strided_delayed(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1, slong stride1,
    const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong prec_abs = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong prec_rel = PADIC_RADIX_CTX_PREC_REL(ctx);

    slong k, vmin, maxhi, Nterms, N, Wdig, nslots;
    int have_value, keep_exact;
    const padic_radix_struct * a;
    const padic_radix_struct * b;
    int init_nonzero;

    nn_ptr accp[FLINT_BITS], accn[FLINT_BITS];   /* e < FLINT_BITS */
    slong nused = 0;
    ulong maskp = 0, maskn = 0;    /* residues with an accumulator */
    TMP_INIT;
    radix_integer_t big, Ptmp;   /* running total + scratch for the radix_mulmid path */

    if (len <= 0)
    {
        if (initial == NULL)
            return padic_radix_zero(res, ctx);
        return padic_radix_set(res, initial, ctx);
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
    if (initial != NULL && initial->N != PADIC_RADIX_EXACT)
        Nterms = DOT_MIN2(Nterms, initial->N);

    for (k = 0; k < len; k++)
    {
        slong va, vb, Na, Nb, vk, Nk;

        a = vec1 + k * stride1;
        b = vec2 + k * stride2;

        /* exact zero annihilates exactly: contributes nothing at all */
        if ((a->u.size == 0 && a->N == PADIC_RADIX_EXACT)
         || (b->u.size == 0 && b->N == PADIC_RADIX_EXACT))
            continue;

        va = a->v; vb = b->v; Na = a->N; Nb = b->N;
        vk = va + vb;

        /* error term of a_k b_k; an inexact zero 0 + O(p^N) has valuation >= N
           (its stored v = 0 is meaningless), cf. padic_radix_mul */
        Nk = PADIC_RADIX_EXACT;
        if (Na != PADIC_RADIX_EXACT) Nk = DOT_MIN2(Nk, ((b->u.size == 0) ? Nb : vb) + Na);
        if (Nb != PADIC_RADIX_EXACT) Nk = DOT_MIN2(Nk, ((a->u.size == 0) ? Na : va) + Nb);
        if (Nk != PADIC_RADIX_EXACT) Nterms = DOT_MIN2(Nterms, Nk);

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
        return _padic_radix_finalize(res, ctx);
    }

    /* absolute precision to compute: terms, context absolute, context relative.
       Using vmin (not the post-cancellation valuation) for the relative cap is the
       accepted slightly-suboptimal choice. */
    N = Nterms;
    N = DOT_MIN2(N, prec_abs);
    if (prec_rel != PADIC_RADIX_PREC_INF)
        N = DOT_MIN2(N, vmin + prec_rel);

    if (N != DOT_INF && N <= vmin)
    {
        /* whole result below the horizon */
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = N;
        return _padic_radix_finalize(res, ctx);
    }

    /* number of result limbs (slots) to cover */
    keep_exact = 0;
    {
        slong nslots_ext = (maxhi - vmin + e - 1) / e + 1;
        if (N == DOT_INF)
            nslots = nslots_ext;
        else
        {
            slong nslots_prec = (N - vmin + e - 1) / e;
            nslots = DOT_MIN2(nslots_prec, nslots_ext);

            /* all terms exact and the exact sum possibly within the context
               precision: compute it fully, and let finalisation keep it
               exact if it fits (as for sums of exact products done with
               mul and add) */
            if (Nterms == DOT_INF && nslots_ext <= nslots_prec + 1)
            {
                nslots = nslots_ext;
                keep_exact = 1;
            }
        }
        if (nslots < 1)
            nslots = 1;
    }

    /* in a context without finite precision, a window (the extent of the
       terms, or the precision they determine) beyond the size limit of the
       context: the addmul chain only builds what the partial sums need
       (getting exact cancellations of terms with huge valuations right),
       and gives up if that is too large */
    if (prec_rel == PADIC_RADIX_PREC_INF && prec_abs == PADIC_RADIX_PREC_INF
        && !_padic_radix_exact_digits_ok(FLINT_MIN(nslots, WORD_MAX / e - 1) * e, ctx))
        return padic_radix_dot_strided_naive(res, initial, subtract, vec1, stride1, vec2, stride2, len, ctx);

    Wdig = (N == DOT_INF || keep_exact) ? DOT_INF : (N - vmin);

    /* ---------------- Pass 2: accumulate ---------------- */
    /* per-residue accumulators from a TMP arena, zeroed when first used */
    TMP_START;
#define DOT_ACC(A, M, r) do { ulong * _m = &(M); \
        if (!((*_m >> (r)) & 1)) { A[r] = TMP_ALLOC(3 * nslots * sizeof(ulong)); \
        flint_mpn_zero(A[r], 3 * nslots); *_m |= UWORD(1) << (r); nused++; } } while (0)
    radix_integer_init(big, radix);     /* large products go here (signed) */
    radix_integer_init(Ptmp, radix);

    /* initial term (always added; its own unit sign chooses the side) */
    if (init_nonzero)
    {
        slong ni = FLINT_ABS(initial->u.size);
        nn_srcptr ui = initial->u.d;
        slong d = initial->v - vmin;        /* >= 0 */
        slong o = (slong) radix_div_exp((ulong) d, radix), r = d - o * e, i;
        nn_ptr acc;

        if (initial->u.size < 0)
        {
            DOT_ACC(accn, maskn, r);
            acc = accn[r];
        }
        else
        {
            DOT_ACC(accp, maskp, r);
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

        o = (slong) radix_div_exp((ulong) d, radix);
        r = d - o * e;
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
            DOT_ACC(accn, maskn, r);
            acc = accn[r];
        }
        else
        {
            DOT_ACC(accp, maskp, r);
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
    {
        slong r, nacc = 0, rlast = 0;
        int lastneg = 0;
        nacc = nused;
        if (nacc == 1)
        {
            lastneg = (maskn != 0);
            rlast = flint_ctz(lastneg ? maskn : maskp);
        }
        if (nacc <= 1 && big->size == 0)
        {
            /* a single accumulator: its residue r just adds to the valuation,
               no digit shift (finalize canonicalises) */
            if (nacc == 0)
                res->u.size = 0;
            else
            {
                _radix_dot_normalize_acc(&res->u, lastneg ? accn[rlast] : accp[rlast], nslots, radix);
                if (lastneg)
                    res->u.size = -res->u.size;
            }
            res->v = vmin + rlast;
        }
        else
        {
            nn_ptr comb = TMP_ALLOC(5 * nslots * sizeof(ulong));
            nn_ptr pd = TMP_ALLOC(2 * (nslots + 7) * sizeof(ulong)), nd = pd + nslots + 7;
            slong pn = 0, nn = 0;
            int side;
            for (side = 0; side < 2; side++)
            {
                nn_ptr * A = side ? accn : accp;
                ulong m = side ? maskn : maskp;
                if (m == 0)
                    continue;
                flint_mpn_zero(comb, 5 * nslots);
                while (m != 0)
                {
                    r = flint_ctz(m);
                    m &= m - 1;
                    _dot_comb_add(comb, A[r], radix->bpow[r], nslots);
                }
                if (side)
                    nn = _dot_comb_normalize(nd, comb, nslots, radix);
                else
                    pn = _dot_comb_normalize(pd, comb, nslots, radix);
            }
            _dot_set_diff(&res->u, pd, pn, nd, nn, radix);       /* M = Mp - Mn (signed) */
            if (big->size != 0)
                radix_integer_add(&res->u, &res->u, big, radix);   /* + large-product total */
            res->v = vmin;
        }
    }
    res->N = keep_exact ? PADIC_RADIX_EXACT : N;

    radix_integer_clear(big, radix);
    radix_integer_clear(Ptmp, radix);
    TMP_END;
#undef DOT_ACC

    return _padic_radix_finalize(res, ctx);
}

/*
    Longest dot product for which the addmul chain beats the delayed dot
    (measured): 4 terms at 1 limb of precision, 3 at 2 limbs, 2 beyond.
*/
static slong
_padic_radix_dot_naive_max(gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong r = PADIC_RADIX_CTX_PREC_REL(ctx), nl;

    if (r == PADIC_RADIX_PREC_INF)
        r = PADIC_RADIX_CTX_PREC_ABS(ctx);
    if (r <= 0 || r == PADIC_RADIX_PREC_INF)
        return 3;
    nl = (slong) radix_div_exp((ulong) r + radix->exp - 1, radix);
    return (nl == 1) ? 4 : (nl == 2) ? 3 : 2;
}

int
padic_radix_dot_strided(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1, slong stride1,
    const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
{
    slong nmax = _padic_radix_dot_naive_max(ctx);
    if (len <= nmax)
        return padic_radix_dot_strided_naive(res, initial, subtract, vec1, stride1, vec2, stride2, len, ctx);
    /* Short dots with differing valuation sums pay a per-residue cost in the
       delayed dot; up to 2 limbs of precision the addmul chain wins below
       about 8 terms in that case. */
    if (len <= 7 && nmax >= 3)
    {
        slong k, v0 = vec1->v + vec2->v;
        for (k = 1; k < len; k++)
            if (vec1[k * stride1].v + vec2[k * stride2].v != v0)
                return padic_radix_dot_strided_naive(res, initial, subtract, vec1, stride1, vec2, stride2, len, ctx);
    }
    return padic_radix_dot_strided_delayed(res, initial, subtract, vec1, stride1, vec2, stride2, len, ctx);
}

int
padic_radix_dot(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1,
    const padic_radix_struct * vec2, slong len, gr_ctx_t ctx)
{
    return padic_radix_dot_strided(res, initial, subtract, vec1, 1, vec2, 1, len, ctx);
}

int
padic_radix_dot_rev(padic_radix_t res, const padic_radix_t initial,
    int subtract, const padic_radix_struct * vec1,
    const padic_radix_struct * vec2, slong len, gr_ctx_t ctx)
{
    return padic_radix_dot_strided(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ctx);
}

