/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz.h"
#include "nmod_poly.h"
#include "mpn_mod.h"
#include "gr.h"
#include "gr_poly.h"
#include "padic_radix.h"
#include "impl.h"

/*
    Polynomial (middle) product through the packed limb-aligned form, as
    padic_radix_mat_mul_packed: coefficient k of the product is the dot
    product sum_i a_i b_(k-i), computed with the semantics of
    padic_radix_dot_strided_delayed; outputs with a single term use
    padic_radix_mul. The coefficients are converted once, and the exponent
    bookkeeping is a set of min-plus convolutions.

    res[0, nhi - nlo) = coefficients nlo, ..., nhi - 1 of poly1 * poly2.
    Returns GR_UNABLE (without touching res) if the kernel does not apply.
*/

typedef struct
{
    slong * vv;
    slong * ve;
    slong * NN;
    slong * ww;
    slong * sz;
    nn_ptr U;
}
_pack_vec_struct;

#define ADD3(c2, c1, c0, ph, pl) add_sssaaaaaa(c2, c1, c0, c2, c1, c0, 0, ph, pl)

static int
_pack_vec_meta(_pack_vec_struct * P, const padic_radix_struct * x, slong len,
    const radix_struct * radix)
{
    padic_radix_pack_meta_struct meta;
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!_padic_radix_pack_meta(&meta, x + i, radix))
            return 0;
        P->vv[i] = meta.vv;
        P->ve[i] = meta.ve;
        P->NN[i] = meta.NN;
        P->ww[i] = meta.ww;
    }

    return 1;
}

static void
_pack_vec_units(_pack_vec_struct * P, const padic_radix_struct * x, slong len,
    slong L, slong stride, const radix_struct * radix)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        nn_ptr U = P->U + stride * i;

        if (P->vv[i] >= PADIC_RADIX_PACK_NONE)
        {
            P->sz[i] = 0;
            flint_mpn_zero(U, stride);
        }
        else
        {
            P->sz[i] = _padic_radix_pack_unit(U, x + i, P->ww[i], L, radix);
            if (stride > L)
                U[1] = 0;
        }
    }
}

/* coefficient k of the product, all of whose terms are exact, as a
   delayed dot product (keeping an exact sum exact when it fits) */
static int
_exact_output(padic_radix_struct * c, const padic_radix_struct * poly1, slong len1,
    const padic_radix_struct * poly2, slong len2, slong k, gr_ctx_t ctx)
{
    slong ilo = FLINT_MAX(0, k - len2 + 1), ihi = FLINT_MIN(k, len1 - 1);
    return padic_radix_dot_strided_delayed(c, NULL, 0, poly1 + ilo, 1,
        poly2 + (k - ilo), -1, ihi - ilo + 1, ctx);
}

int
_padic_radix_poly_mulmid_packed(padic_radix_struct * res,
    const padic_radix_struct * poly1, slong len1,
    const padic_radix_struct * poly2, slong len2,
    slong nlo, slong nhi, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong k, i, L, stride, nout = nhi - nlo;
    _pack_vec_struct P1, P2;
    slong * ex, * PN, * PW, * Pslots;
    nn_ptr acc;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (PADIC_RADIX_CTX_PREC_REL(ctx) == PADIC_RADIX_PREC_INF
        && PADIC_RADIX_CTX_PREC_ABS(ctx) == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;

    if (nout <= 0)
        return GR_SUCCESS;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    TMP_START;

    ex = TMP_ALLOC(sizeof(slong) * (5 * (len1 + len2) + 3 * nout));
    P1.vv = ex;           P1.ve = P1.vv + len1;  P1.NN = P1.ve + len1;
    P1.ww = P1.NN + len1; P1.sz = P1.ww + len1;
    P2.vv = P1.sz + len1; P2.ve = P2.vv + len2;  P2.NN = P2.ve + len2;
    P2.ww = P2.NN + len2; P2.sz = P2.ww + len2;
    PN = P2.sz + len2;
    PW = PN + nout;
    Pslots = PW + nout;

    if (!_pack_vec_meta(&P1, poly1, len1, radix) || !_pack_vec_meta(&P2, poly2, len2, radix))
    {
        TMP_END;
        return GR_UNABLE;
    }

    /* exponent bookkeeping: min-plus convolutions */
    L = 1;
    for (k = nlo; k < nhi; k++)
    {
        slong ilo = FLINT_MAX(0, k - len2 + 1), ihi = FLINT_MIN(k, len1 - 1);
        slong vmin = 2 * PADIC_RADIX_PACK_BIG, Nt = 2 * PADIC_RADIX_PACK_BIG;
        slong W = 2 * PADIC_RADIX_PACK_BIG, N, s;

        for (i = ilo; i <= ihi; i++)
        {
            slong q = k - i;
            slong x1 = P1.vv[i] + P2.vv[q];
            slong e1 = P1.NN[i] + P2.ve[q];
            slong e2 = P1.ve[i] + P2.NN[q];
            slong f = P1.ww[i] + P2.ww[q];
            vmin = FLINT_MIN(vmin, x1);
            e1 = FLINT_MIN(e1, e2);
            Nt = FLINT_MIN(Nt, e1);
            W = FLINT_MIN(W, f);
        }

        if (vmin >= PADIC_RADIX_PACK_NONE)
        {
            N = (Nt >= PADIC_RADIX_PACK_NONE) ? PADIC_RADIX_EXACT : Nt;
            s = 0;
        }
        else if (Nt >= PADIC_RADIX_PACK_NONE)
        {
            N = PADIC_RADIX_PACK_EXACT_SUM;     /* only exact terms */
            s = 0;
        }
        else
        {
            N = _padic_radix_pack_prec(vmin, Nt, ctx);
            s = (N <= vmin) ? 0 : (slong) radix_div_exp((ulong) (N - e * W) + e - 1, radix);
        }

        PN[k - nlo] = N;
        PW[k - nlo] = W;
        Pslots[k - nlo] = s;
        if (ihi > ilo)
            L = FLINT_MAX(L, s);
    }

    stride = FLINT_MAX(L, 2);

    P1.U = TMP_ALLOC(sizeof(ulong) * stride * (len1 + len2));
    P2.U = P1.U + stride * len1;
    _pack_vec_units(&P1, poly1, len1, L, stride, radix);
    _pack_vec_units(&P2, poly2, len2, L, stride, radix);

    acc = TMP_ALLOC(sizeof(ulong) * 3 * FLINT_MAX(L, 2));

    for (k = nlo; k < nhi; k++)
    {
        slong ilo = FLINT_MAX(0, k - len2 + 1), ihi = FLINT_MIN(k, len1 - 1);
        slong N = PN[k - nlo], W = PW[k - nlo], nslots = Pslots[k - nlo];
        padic_radix_struct * c = res + (k - nlo);
        nn_ptr d;

        if (ihi == ilo)
        {
            /* a single term: the product, with its exact semantics */
            status |= padic_radix_mul(c, poly1 + ilo, poly2 + (k - ilo), ctx);
            continue;
        }

        if (N == PADIC_RADIX_PACK_EXACT_SUM)
        {
            status |= _exact_output(c, poly1, len1, poly2, len2, k, ctx);
            continue;
        }

        if (nslots == 0)
        {
            radix_integer_zero(&c->u, radix);
            c->v = 0;
            c->N = N;
            status |= _padic_radix_finalize(c, ctx);
            continue;
        }

        if (L <= 2)
        {
            ulong c0 = 0, c1 = 0, c2 = 0, d0 = 0, d1 = 0, d2 = 0, ph, pl;

            if (nslots == 1)
            {
                for (i = ilo; i <= ihi; i++)
                {
                    slong q = k - i;
                    if (P1.ww[i] + P2.ww[q] == W)
                    {
                        umul_ppmm(ph, pl, P1.U[stride * i], P2.U[stride * q]);
                        ADD3(c2, c1, c0, ph, pl);
                    }
                }
            }
            else
            {
                for (i = ilo; i <= ihi; i++)
                {
                    slong q = k - i;
                    slong o = P1.ww[i] + P2.ww[q] - W;
                    ulong a0 = P1.U[stride * i], a1 = P1.U[stride * i + 1];
                    ulong b0 = P2.U[stride * q], b1 = P2.U[stride * q + 1];

                    if (o == 0)
                    {
                        umul_ppmm(ph, pl, a0, b0);
                        ADD3(c2, c1, c0, ph, pl);
                        umul_ppmm(ph, pl, a0, b1);
                        ADD3(d2, d1, d0, ph, pl);
                        umul_ppmm(ph, pl, a1, b0);
                        ADD3(d2, d1, d0, ph, pl);
                    }
                    else if (o == 1)
                    {
                        umul_ppmm(ph, pl, a0, b0);
                        ADD3(d2, d1, d0, ph, pl);
                    }
                }
            }

            acc[0] = c0; acc[1] = c1; acc[2] = c2;
            acc[3] = d0; acc[4] = d1; acc[5] = d2;
        }
        else
        {
            flint_mpn_zero(acc, 3 * nslots);

            for (i = ilo; i <= ihi; i++)
            {
                slong q = k - i;
                slong an = P1.sz[i], bn = P2.sz[q], o;

                o = P1.ww[i] + P2.ww[q] - W;
                if (o >= nslots || an == 0 || bn == 0)
                    continue;

                _padic_radix_pack_addmul(acc, P1.U + stride * i, an, P2.U + stride * q, bn, o, nslots);
            }
        }

        d = radix_integer_fit_limbs(&c->u, nslots + 3, radix);
        c->u.size = _padic_radix_pack_normalise(d, acc, nslots, radix);
        status |= _padic_radix_pack_finish(c, W, N, ctx);
    }

    TMP_END;
    return status;
}

/* floor(x / e) */
static slong
_fdiv(slong x, slong e)
{
    return (x >= 0) ? x / e : -((-x + e - 1) / e);
}

/*
    Polynomial (middle) product by Kronecker substitution in radix B.
    Relative to the base frames wa = min w(a_i) and wb = min w(b_j) (over the
    nonzero coefficients), coefficient i is the residue U_i B^(w_i - wa)
    truncated to the S_max limbs that can reach any output window, placed in
    a slot of S = 2 S_max + O(1) limbs of one radix integer; one radix_mulmid
    (Kronecker substitution or the small-prime FFT in radix B, no conversion
    to binary) then gives every output in its own slot, without overlap.

    The precisions of the outputs come from min-plus convolutions of the
    exponent sequences (_padic_radix_minplus_conv): exact for short or
    small-range sequences (or convex ones), and otherwise a lower bound by
    convex hulls, so that the result may then claim less precision than
    padic_radix_dot_strided_delayed would (but never more). When they are
    exact the outputs are those of _padic_radix_poly_mulmid_packed.
*/
static int
_poly_mulmid_frames(padic_radix_struct * res,
    const padic_radix_struct * poly1, slong len1,
    const padic_radix_struct * poly2, slong len2,
    slong nlo, slong nhi, int alg, int * spread, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong k, i, nout = nhi - nlo;
    slong wa, wb, wamax, wbmax, Wb, Smax, Sown, S, cl, an, bn;
    _pack_vec_struct P1, P2;
    slong * ex, * V, * E, * E2, * PN, * XW, * SW, nex, nouts;
    nn_ptr A, B, C;
    int status = GR_SUCCESS, exact;
    TMP_INIT;

    if (PADIC_RADIX_CTX_PREC_REL(ctx) == PADIC_RADIX_PREC_INF
        && PADIC_RADIX_CTX_PREC_ABS(ctx) == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;

    if (nout <= 0)
        return GR_SUCCESS;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    TMP_START;

    ex = TMP_ALLOC(sizeof(slong) * (4 * (len1 + len2) + 6 * nout));
    P1.vv = ex;           P1.ve = P1.vv + len1;  P1.NN = P1.ve + len1;
    P1.ww = P1.NN + len1;
    P2.vv = P1.ww + len1; P2.ve = P2.vv + len2;  P2.NN = P2.ve + len2;
    P2.ww = P2.NN + len2;
    V = P2.ww + len2;
    E = V + nout;
    E2 = E + nout;
    PN = E2 + nout;
    XW = PN + nout;
    SW = XW + nout;

    if (!_pack_vec_meta(&P1, poly1, len1, radix) || !_pack_vec_meta(&P2, poly2, len2, radix))
        goto unable;

    /* precisions */
    exact = !(PADIC_RADIX_CTX_FLAGS(ctx) & PADIC_RADIX_FAST_BOUNDS);
    _padic_radix_minplus_conv(V, P1.vv, len1, P2.vv, len2, nlo, nhi, exact);
    _padic_radix_minplus_conv(E, P1.NN, len1, P2.ve, len2, nlo, nhi, exact);
    _padic_radix_minplus_conv(E2, P1.ve, len1, P2.NN, len2, nlo, nhi, exact);

    /* base frames */
    wa = wb = PADIC_RADIX_PACK_BIG;
    wamax = wbmax = -PADIC_RADIX_PACK_BIG;
    for (i = 0; i < len1; i++)
        if (P1.ww[i] < PADIC_RADIX_PACK_NONE)
        {
            wa = FLINT_MIN(wa, P1.ww[i]);
            wamax = FLINT_MAX(wamax, P1.ww[i]);
        }
    for (i = 0; i < len2; i++)
        if (P2.ww[i] < PADIC_RADIX_PACK_NONE)
        {
            wb = FLINT_MIN(wb, P2.ww[i]);
            wbmax = FLINT_MAX(wbmax, P2.ww[i]);
        }

    Wb = (wa < PADIC_RADIX_PACK_NONE && wb < PADIC_RADIX_PACK_NONE) ? wa + wb : 0;

    /*
        Outputs with only exact terms (and more than one): as in
        padic_radix_dot_strided_delayed, the exact sum is kept when its a
        priori extent max (v_i + e |u_i| + v_j + e |u_j|) over the terms is at
        most one limb beyond the context precision; it is then computed
        modulo B^w for a window w covering that extent, one carry limb and
        one sign limb. The extents are a max-plus convolution (over the exact
        entries, which are the only ones such outputs involve).
    */
    nex = 0;
    for (k = 0; k < nout; k++)
    {
        slong kk = k + nlo, ilo = FLINT_MAX(0, kk - len2 + 1), ihi = FLINT_MIN(kk, len1 - 1);
        XW[k] = 0;
        if (V[k] < PADIC_RADIX_PACK_NONE && FLINT_MIN(E[k], E2[k]) >= PADIC_RADIX_PACK_NONE && ihi > ilo)
        {
            XW[k] = 1;
            nex++;
        }
    }

    if (nex != 0)
    {
        slong * ha, * hb, * H;

        ha = TMP_ALLOC(sizeof(slong) * (len1 + len2 + nout));
        hb = ha + len1;
        H = hb + len2;

        for (i = 0; i < len1; i++)
            ha[i] = (poly1[i].u.size != 0 && poly1[i].N == PADIC_RADIX_EXACT)
                ? -(poly1[i].v + e * FLINT_ABS(poly1[i].u.size)) : 2 * PADIC_RADIX_PACK_BIG;
        for (i = 0; i < len2; i++)
            hb[i] = (poly2[i].u.size != 0 && poly2[i].N == PADIC_RADIX_EXACT)
                ? -(poly2[i].v + e * FLINT_ABS(poly2[i].u.size)) : 2 * PADIC_RADIX_PACK_BIG;

        _padic_radix_minplus_conv(H, ha, len1, hb, len2, nlo, nhi, 1);

        for (k = 0; k < nout; k++)
        {
            if (XW[k])
            {
                slong vmin = V[k], hi = -H[k];
                slong Ncap = _padic_radix_pack_prec(vmin, PADIC_RADIX_PACK_NONE, ctx);
                slong ext = (hi - vmin + e - 1) / e + 1, prec = (Ncap - vmin + e - 1) / e;

                if (Ncap != PADIC_RADIX_EXACT && Ncap > vmin && ext <= prec + 1)
                    XW[k] = (hi - e * Wb + e - 1) / e + 2;
                else
                    XW[k] = 0;
            }
        }
    }

    /* windows in the base frame (SW) and in the own frames of the outputs
       (the maximum Sown); the outputs that are single products are done
       separately with padic_radix_mul and need no window */
    Sown = 0;
    for (k = 0; k < nout; k++)
    {
        slong N, Nt = FLINT_MIN(E[k], E2[k]);
        slong kk = k + nlo, ilo = FLINT_MAX(0, kk - len2 + 1), ihi = FLINT_MIN(kk, len1 - 1);

        if (V[k] >= PADIC_RADIX_PACK_NONE)
            N = (Nt >= PADIC_RADIX_PACK_NONE) ? PADIC_RADIX_EXACT : Nt;
        else if (XW[k] != 0)
            N = PADIC_RADIX_PACK_EXACT_SUM;     /* exact sum, window XW[k] */
        else
            N = _padic_radix_pack_prec(V[k], Nt, ctx);

        PN[k] = N;
        SW[k] = 0;

        if (ihi <= ilo)
            continue;

        if (N == PADIC_RADIX_PACK_EXACT_SUM)
            SW[k] = XW[k];
        else if (V[k] < PADIC_RADIX_PACK_NONE && N != PADIC_RADIX_EXACT && N > V[k])
            SW[k] = (slong) radix_div_exp((ulong) (N - e * Wb) + e - 1, radix);

        if (SW[k] != 0)
            Sown = FLINT_MAX(Sown, SW[k] - (_fdiv(V[k], e) - Wb));
    }

    /*
        Outputs whose windows in the base frame are much wider than in their
        own frames (all their terms have valuations far above the base
        frame): when there are few of them, they are computed separately by
        dot products (marked SW = -1); otherwise (valuations or precisions
        spread far beyond the precision along the outputs) the packed
        kernel, which uses the frame of every output, is cheaper, and the
        slots could be arbitrarily large.
    */
    nouts = 0;
    for (k = 0; k < nout; k++)
        nouts += (SW[k] > 4 * Sown + 8);

    if (nouts != 0)
    {
        /* at most about the cost of the transform, and a small part of
           that of the classical algorithm */
        double budget = 4.0 * (len1 + len2) * (4 * Sown + 9) * (FLINT_BIT_COUNT(len1 + len2) + 1);
        budget = FLINT_MIN(budget, (double) len1 * len2 / 16);

        if ((double) nouts * FLINT_MIN(len1, len2) > budget)
        {
            if (spread != NULL)
                *spread = 1;
            goto unable;
        }

        for (k = 0; k < nout; k++)
            if (SW[k] > 4 * Sown + 8)
                SW[k] = -1;
    }

    Smax = 0;
    for (k = 0; k < nout; k++)
        Smax = FLINT_MAX(Smax, SW[k]);

    if (Smax == 0)
    {
        /* every output is zero to its precision (or a single product) */
        for (k = nlo; k < nhi; k++)
        {
            slong ilo = FLINT_MAX(0, k - len2 + 1), ihi = FLINT_MIN(k, len1 - 1);
            padic_radix_struct * c = res + (k - nlo);
            if (ihi == ilo)
                status |= padic_radix_mul(c, poly1 + ilo, poly2 + (k - ilo), ctx);
            else if (SW[k - nlo] == -1)
                status |= padic_radix_dot_strided_delayed(c, NULL, 0, poly1 + ilo, 1,
                    poly2 + (k - ilo), -1, ihi - ilo + 1, ctx);
            else if (PN[k - nlo] == PADIC_RADIX_PACK_EXACT_SUM)
                status |= _exact_output(c, poly1, len1, poly2, len2, k, ctx);
            else
            {
                radix_integer_zero(&c->u, radix);
                c->v = 0;
                c->N = PN[k - nlo];
                status |= _padic_radix_finalize(c, ctx);
            }
        }
        TMP_END;
        return status;
    }

    if (alg == 1)
    {
        /* products modulo B^Smax: nmod_poly (one limb) or mpn_mod_poly
           (up to MPN_MOD_MAX_LIMBS binary limbs); the outputs are then
           converted to radix B, Smax limbs each, in slots of S limbs */
        slong nl;
        fmpz_t M;

        fmpz_init(M);
        fmpz_set_ui(M, radix->B.n);
        fmpz_pow_ui(M, M, Smax);
        nl = fmpz_size(M);

        if (Smax == 1)
        {
            nn_ptr a1, b1, c1;

            S = 1;
            a1 = TMP_ALLOC(sizeof(ulong) * (len1 + len2 + 2 * nout));
            b1 = a1 + len1;
            c1 = b1 + len2;
            C = c1 + nout;

            for (i = 0; i < len1; i++)
                a1[i] = (P1.ww[i] == wa) ? (_padic_radix_pack_unit(a1 + i, poly1 + i, P1.ww[i], 1, radix), a1[i]) : 0;
            for (i = 0; i < len2; i++)
                b1[i] = (P2.ww[i] == wb) ? (_padic_radix_pack_unit(b1 + i, poly2 + i, P2.ww[i], 1, radix), b1[i]) : 0;

            _nmod_poly_mulmid(c1, a1, len1, b1, len2, nlo, nhi, radix->B);
            flint_mpn_copyi(C, c1, nout);
        }
        else if (nl >= MPN_MOD_MIN_LIMBS && nl <= MPN_MOD_MAX_LIMBS)
        {
            gr_ctx_t mctx;
            nn_ptr am, bm, cm, tmp, bin;
            slong bl;

            GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(mctx, M));

            S = Smax + 1;
            am = TMP_ALLOC(sizeof(ulong) * ((len1 + len2 + nout) * nl + nout * S + 2 * Smax + 4
                + radix_set_mpn_need_alloc(nl, radix)));
            bm = am + len1 * nl;
            cm = bm + len2 * nl;
            C = cm + nout * nl;
            tmp = C + nout * S;
            bin = tmp + Smax + 1;
            flint_mpn_zero(am, (len1 + len2) * nl);

            for (i = 0; i < len1; i++)
                if (P1.ww[i] < PADIC_RADIX_PACK_NONE && P1.ww[i] - wa < Smax)
                {
                    slong o = P1.ww[i] - wa;
                    flint_mpn_zero(tmp, o);
                    _padic_radix_pack_unit(tmp + o, poly1 + i, P1.ww[i], Smax - o, radix);
                    bl = radix_get_mpn(bin, tmp, Smax, radix);
                    flint_mpn_copyi(am + i * nl, bin, FLINT_MIN(bl, nl));
                }

            for (i = 0; i < len2; i++)
                if (P2.ww[i] < PADIC_RADIX_PACK_NONE && P2.ww[i] - wb < Smax)
                {
                    slong o = P2.ww[i] - wb;
                    flint_mpn_zero(tmp, o);
                    _padic_radix_pack_unit(tmp + o, poly2 + i, P2.ww[i], Smax - o, radix);
                    bl = radix_get_mpn(bin, tmp, Smax, radix);
                    flint_mpn_copyi(bm + i * nl, bin, FLINT_MIN(bl, nl));
                }

            if (len1 >= len2)
                status |= _mpn_mod_poly_mulmid(cm, am, len1, bm, len2, nlo, nhi, mctx);
            else
                status |= _mpn_mod_poly_mulmid(cm, bm, len2, am, len1, nlo, nhi, mctx);

            for (k = 0; k < nout; k++)
            {
                nn_srcptr z = cm + k * nl;
                slong zn = nl, dn;

                while (zn > 0 && z[zn - 1] == 0)
                    zn--;
                dn = (zn == 0) ? 0 : radix_set_mpn(bin, z, zn, radix);
                dn = FLINT_MIN(dn, Smax);
                flint_mpn_copyi(C + k * S, bin, dn);
                flint_mpn_zero(C + k * S + dn, S - dn);
            }

            gr_ctx_clear(mctx);
        }
        else
        {
            fmpz_clear(M);
            goto unable;
        }

        fmpz_clear(M);
        goto extract;
    }

    /* slot width: a sum of min(len1, len2) products of two S_max-limb
       numbers, below min(len1, len2) B^(2 S_max) */
    cl = 1;
    {
        ulong t = (ulong) FLINT_MIN(len1, len2), Bn = radix->B.n;
        while (t >= Bn)
        {
            t /= Bn;
            cl++;
        }
    }
    S = 2 * Smax + cl;

    an = len1 * S;
    bn = len2 * S;
    A = TMP_ALLOC(sizeof(ulong) * (an + bn + nout * S));
    B = A + an;
    C = B + bn;
    flint_mpn_zero(A, an + bn);

    for (i = 0; i < len1; i++)
        if (P1.ww[i] < PADIC_RADIX_PACK_NONE)
        {
            slong o = P1.ww[i] - wa;
            if (o < Smax)
                _padic_radix_pack_unit(A + i * S + o, poly1 + i, P1.ww[i], Smax - o, radix);
        }

    for (i = 0; i < len2; i++)
        if (P2.ww[i] < PADIC_RADIX_PACK_NONE)
        {
            slong o = P2.ww[i] - wb;
            if (o < Smax)
                _padic_radix_pack_unit(B + i * S + o, poly2 + i, P2.ww[i], Smax - o, radix);
        }

    radix_mulmid(C, A, an, B, bn, nlo * S, nhi * S, radix);

extract:
    for (k = nlo; k < nhi; k++)
    {
        slong ilo = FLINT_MAX(0, k - len2 + 1), ihi = FLINT_MIN(k, len1 - 1);
        padic_radix_struct * c = res + (k - nlo);
        slong N = PN[k - nlo], nl, j0;
        nn_srcptr slot = C + (k - nlo) * S;

        if (ihi == ilo)
        {
            status |= padic_radix_mul(c, poly1 + ilo, poly2 + (k - ilo), ctx);
            continue;
        }

        if (SW[k - nlo] == -1)
        {
            status |= padic_radix_dot_strided_delayed(c, NULL, 0, poly1 + ilo, 1,
                poly2 + (k - ilo), -1, ihi - ilo + 1, ctx);
            continue;
        }

        if (N == PADIC_RADIX_PACK_EXACT_SUM)
        {
            /* the exact sum modulo B^w, |sum| < B^(w - 1): signed digits */
            slong w = XW[k - nlo], sz;
            nn_ptr d = radix_integer_fit_limbs(&c->u, w, radix);
            int neg = (slot[w - 1] != 0);

            if (neg)
                radix_neg(d, slot, w, radix);
            else
                flint_mpn_copyi(d, slot, w);
            sz = w;
            while (sz > 0 && d[sz - 1] == 0)
                sz--;
            c->u.size = neg ? -sz : sz;
            c->v = e * Wb;
            c->N = PADIC_RADIX_EXACT;
            status |= _padic_radix_finalize(c, ctx);
            continue;
        }

        if (V[k - nlo] >= PADIC_RADIX_PACK_NONE || N == PADIC_RADIX_EXACT || N <= V[k - nlo])
            nl = 0;
        else
            nl = (slong) radix_div_exp((ulong) (N - e * Wb) + e - 1, radix);

        j0 = 0;
        while (j0 < nl && slot[j0] == 0)
            j0++;

        if (j0 == nl)
        {
            radix_integer_zero(&c->u, radix);
            c->v = 0;
            c->N = N;
            status |= _padic_radix_finalize(c, ctx);
        }
        else
        {
            nn_ptr d = radix_integer_fit_limbs(&c->u, nl - j0, radix);
            slong sz = nl - j0;
            flint_mpn_copyi(d, slot + j0, sz);
            while (sz > 0 && d[sz - 1] == 0)
                sz--;
            c->u.size = sz;
            status |= _padic_radix_pack_finish(c, Wb + j0, N, ctx);
        }
    }

    TMP_END;
    return status;

unable:
    TMP_END;
    return GR_UNABLE;
}

int
_padic_radix_poly_mulmid_KS(padic_radix_struct * res,
    const padic_radix_struct * poly1, slong len1,
    const padic_radix_struct * poly2, slong len2,
    slong nlo, slong nhi, gr_ctx_t ctx)
{
    return _poly_mulmid_frames(res, poly1, len1, poly2, len2, nlo, nhi, 0, NULL, ctx);
}

int
_padic_radix_poly_mulmid_mod(padic_radix_struct * res,
    const padic_radix_struct * poly1, slong len1,
    const padic_radix_struct * poly2, slong len2,
    slong nlo, slong nhi, gr_ctx_t ctx)
{
    return _poly_mulmid_frames(res, poly1, len1, poly2, len2, nlo, nhi, 1, NULL, ctx);
}

/* below these lengths the classical algorithm (dot products, or addmul
   chains for the shortest ones) is as fast */
#define PACKED_MIN_LEN 4

/*
    Algorithm choice (0: packed classical, 1: radix Kronecker substitution,
    2: products modulo B^S with nmod_poly / mpn_mod_poly) from the number
    of limbs L of the units and an estimate S of the largest output window
    (L, plus one if the valuations are not multiples of e, plus the spread
    of the limb frames, which the packed kernel avoids by using the frame of
    every output). Measured crossovers with the packed kernel, for uniform
    valuations: modulo B (one limb windows) from length 4; two limbs from
    about 64; up to 10 limbs from about 24; for wider windows, the radix KS
    (which avoids conversions and beats mpn_mod there) from about 32.
*/
/*
    Whether several outputs at a sample of positions have windows, in the
    base frame of the Kronecker substitution, much wider than in their own
    frames (valuations spread beyond the precision, e.g. alternating or
    trending valuations), in which case the frames routes would give up
    (after the min-plus precision computation) or be slow. Uses the context
    precision as an upper bound for the precision of every output.
*/
static int
_sampled_spread(const padic_radix_struct * poly1, slong len1, const padic_radix_struct * poly2,
    slong len2, slong nlo, slong nhi, slong Wb, gr_ctx_t ctx)
{
    slong e = PADIC_RADIX_CTX_RADIX(ctx)->exp;
    slong rel = PADIC_RADIX_CTX_PREC_REL(ctx), pabs = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong s, nsamp = FLINT_MIN(nhi - nlo, 16), ninfl = 0;
    double wb[16], maxo = 0;

    for (s = 0; s < nsamp; s++)
    {
        slong k = nlo + (nsamp == 1 ? 0 : (s * (nhi - 1 - nlo)) / (nsamp - 1));
        slong i, ilo = FLINT_MAX(0, k - len2 + 1), ihi = FLINT_MIN(k, len1 - 1);
        double V = HUGE_VAL, N;

        wb[s] = 0;
        if (ihi <= ilo)     /* single products are done separately */
            continue;

        for (i = ilo; i <= ihi; i++)
            if (poly1[i].u.size != 0 && poly2[k - i].u.size != 0)
                V = FLINT_MIN(V, (double) poly1[i].v + (double) poly2[k - i].v);

        if (V == HUGE_VAL)
            continue;

        N = HUGE_VAL;
        if (rel != PADIC_RADIX_PREC_INF)
            N = V + rel;
        if (pabs != PADIC_RADIX_PREC_INF)
            N = FLINT_MIN(N, (double) pabs);
        if (N == HUGE_VAL || N <= V)
            continue;

        wb[s] = ceil((N - e * (double) Wb) / e);
        maxo = FLINT_MAX(maxo, ceil((N - e * floor(V / e)) / e));
    }

    /* a few inflated outputs are computed separately by the frames routes */
    for (s = 0; s < nsamp; s++)
        ninfl += (wb[s] > 4 * maxo + 8);

    return ninfl >= 3;
}

static int
_choose_alg(const padic_radix_struct * poly1, slong len1, const padic_radix_struct * poly2,
    slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong i, L = 1, wlo1 = WORD_MAX, whi1 = WORD_MIN, wlo2 = WORD_MAX, whi2 = WORD_MIN, shift = 0, n, S;
    double thresh;

    n = FLINT_MIN(len1, len2);

    for (i = 0; i < len1; i++)
        if (poly1[i].u.size != 0)
        {
            slong w = poly1[i].v >= 0 ? poly1[i].v / (slong) radix->exp : -1 - (-poly1[i].v) / (slong) radix->exp;
            L = FLINT_MAX(L, FLINT_ABS(poly1[i].u.size));
            shift |= (poly1[i].v - w * (slong) radix->exp);
            wlo1 = FLINT_MIN(wlo1, w);
            whi1 = FLINT_MAX(whi1, w);
        }
    for (i = 0; i < len2; i++)
        if (poly2[i].u.size != 0)
        {
            slong w = poly2[i].v >= 0 ? poly2[i].v / (slong) radix->exp : -1 - (-poly2[i].v) / (slong) radix->exp;
            L = FLINT_MAX(L, FLINT_ABS(poly2[i].u.size));
            shift |= (poly2[i].v - w * (slong) radix->exp);
            wlo2 = FLINT_MIN(wlo2, w);
            whi2 = FLINT_MAX(whi2, w);
        }

    if (wlo1 == WORD_MAX || wlo2 == WORD_MAX)
        return 0;

    /* frames beyond the precision window do not widen the windows */
    S = FLINT_MIN(L + (shift != 0) + (whi1 - wlo1) + (whi2 - wlo2), 4 * L);

    if (S <= 1)
        return 2;

    thresh = (S == 2) ? 64 : (S <= 10) ? 24 : 32;
    if (S > L)
        thresh *= pow((double) S / L, 1.5);

    if (n < thresh)
        return 0;

    /* (only needed when the frames are spread out) */
    if ((whi1 - wlo1) + (whi2 - wlo2) > 2 * L + 2
        && _sampled_spread(poly1, len1, poly2, len2, nlo, nhi, wlo1 + wlo2, ctx))
        return 0;

    return (S <= 10) ? 2 : 1;
}

int
_padic_radix_poly_mulmid(padic_radix_struct * res,
    const padic_radix_struct * poly1, slong len1,
    const padic_radix_struct * poly2, slong len2,
    slong nlo, slong nhi, gr_ctx_t ctx)
{
    int status;

    if (FLINT_MIN(len1, len2) >= PACKED_MIN_LEN)
    {
        int alg = _choose_alg(poly1, FLINT_MIN(len1, nhi), poly2, FLINT_MIN(len2, nhi), nlo, nhi, ctx);

        int spread = 0;

        if (alg == 2)
        {
            status = _poly_mulmid_frames(res, poly1, len1, poly2, len2, nlo, nhi, 1, &spread, ctx);
            if (status != GR_UNABLE)
                return status;
            alg = (FLINT_MIN(len1, len2) >= 32) && !spread;
        }

        if (alg == 1)
        {
            status = _poly_mulmid_frames(res, poly1, len1, poly2, len2, nlo, nhi, 0, &spread, ctx);
            if (status != GR_UNABLE)
                return status;
        }

        status = _padic_radix_poly_mulmid_packed(res, poly1, len1, poly2, len2, nlo, nhi, ctx);
        if (status != GR_UNABLE)
            return status;
    }

    return _gr_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, ctx);
}

int
_padic_radix_poly_mullow(padic_radix_struct * res,
    const padic_radix_struct * poly1, slong len1,
    const padic_radix_struct * poly2, slong len2,
    slong n, gr_ctx_t ctx)
{
    return _padic_radix_poly_mulmid(res, poly1, len1, poly2, len2, 0, n, ctx);
}
