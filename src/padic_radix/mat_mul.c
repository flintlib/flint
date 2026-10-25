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
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "mpn_mod.h"
#include "gr.h"
#include "gr_mat.h"
#include "padic_radix.h"
#include "impl.h"

/*
    Matrix multiplication through the packed limb-aligned form (impl.h).
    Every output C_ij is the dot product of row i of A and column j of B,
    computed with exactly the semantics of padic_radix_dot_strided_delayed
    (same precision, same value modulo it), but:

      - the exponent bookkeeping of all m n outputs is done up front as
        branch-free min-plus products over arrays of exponents;
      - each entry of A and B is converted once (O(mk + kn)) to its
        limb-aligned residue, so that all terms of an output land at
        whole-limb offsets of one accumulator: no per-residue accumulators,
        digit shifts or per-dot allocations;
      - at precisions of up to two limbs (the windows of all outputs fit in
        two limbs) the accumulation runs in registers.

    Returns GR_UNABLE (without touching C) when the kernel does not apply:
    contexts without finite precision (exact results need full products)
    and elements with extreme exponents.
*/

/* exponent arrays (struct of arrays) and residues of a packed matrix */
typedef struct
{
    slong * vv;
    slong * ve;
    slong * NN;
    slong * ww;
    slong * sz;
    nn_ptr U;
}
_pack_mat_struct;

/* c = the exact sum whose residue modulo B^s has the digits d[0, dn)
   (dn <= s), relative to the limb frame W, given that |sum| < B^(s - 1):
   a nonzero top limb means a negative sum */
/* floor(x / e) */
static slong
_fdiv(slong x, slong e)
{
    return (x >= 0) ? x / e : -((-x + e - 1) / e);
}

static int
_exact_sum_finish(padic_radix_struct * c, nn_srcptr d, slong dn, slong s, slong W, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    nn_ptr u = radix_integer_fit_limbs(&c->u, s, radix);
    slong sz;
    int neg;

    flint_mpn_copyi(u, d, dn);
    flint_mpn_zero(u + dn, s - dn);
    neg = (u[s - 1] != 0);
    if (neg)
        radix_neg(u, u, s, radix);
    sz = s;
    while (sz > 0 && u[sz - 1] == 0)
        sz--;
    c->u.size = neg ? -sz : sz;
    c->v = (slong) radix->exp * W;
    c->N = PADIC_RADIX_EXACT;
    return _padic_radix_finalize(c, ctx);
}

#define ADD3(c2, c1, c0, ph, pl) add_sssaaaaaa(c2, c1, c0, c2, c1, c0, 0, ph, pl)

int
padic_radix_mat_mul_packed(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong m = A->r, kk = A->c, n = B->c;
    slong i, j, t, L, stride, nA, nB, * hA, * hB;
    int anyexA = 0, anyexB = 0;
    _pack_mat_struct PA, PB;
    slong * PN, * PW, * Pslots, * ex;
    nn_ptr acc;
    padic_radix_pack_meta_struct meta;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (B->r != kk || C->r != m || C->c != n)
        return GR_DOMAIN;

    if (PADIC_RADIX_CTX_PREC_REL(ctx) == PADIC_RADIX_PREC_INF
        && PADIC_RADIX_CTX_PREC_ABS(ctx) == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;

    if (m == 0 || n == 0)
        return GR_SUCCESS;

    if (kk == 0)
        return gr_mat_zero(C, ctx);

    nA = m * kk;
    nB = kk * n;

    TMP_START;

    ex = TMP_ALLOC(sizeof(slong) * (5 * (nA + nB) + 3 * m * n));
    PA.vv = ex;           PA.ve = PA.vv + nA;   PA.NN = PA.ve + nA;
    PA.ww = PA.NN + nA;   PA.sz = PA.ww + nA;
    PB.vv = PA.sz + nA;   PB.ve = PB.vv + nB;   PB.NN = PB.ve + nB;
    PB.ww = PB.NN + nB;   PB.sz = PB.ww + nB;
    PN = PB.sz + nB;
    PW = PN + m * n;
    Pslots = PW + m * n;

    /* A row-major, B transposed (column j of B contiguous) */
    for (i = 0; i < m; i++)
        for (t = 0; t < kk; t++)
        {
            slong q = i * kk + t;
            if (!_padic_radix_pack_meta(&meta, gr_mat_entry_srcptr(A, i, t, ctx), radix))
                goto unable;
            PA.vv[q] = meta.vv; PA.ve[q] = meta.ve; PA.NN[q] = meta.NN; PA.ww[q] = meta.ww;
        }

    for (t = 0; t < kk; t++)
        for (j = 0; j < n; j++)
        {
            slong q = j * kk + t;
            if (!_padic_radix_pack_meta(&meta, gr_mat_entry_srcptr(B, t, j, ctx), radix))
                goto unable;
            PB.vv[q] = meta.vv; PB.ve[q] = meta.ve; PB.NN[q] = meta.NN; PB.ww[q] = meta.ww;
        }

    /* a priori extents v + e |u| of the exact nonzero entries (for the
       outputs whose terms are all exact, see padic_radix_dot_strided_delayed) */
    hA = TMP_ALLOC(sizeof(slong) * (nA + nB));
    hB = hA + nA;
    for (i = 0; i < m; i++)
        for (t = 0; t < kk; t++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(A, i, t, ctx);
            int isex = (x->u.size != 0 && x->N == PADIC_RADIX_EXACT);
            hA[i * kk + t] = isex ? x->v + e * FLINT_ABS(x->u.size) : -PADIC_RADIX_PACK_BIG;
            anyexA |= isex;
        }
    for (t = 0; t < kk; t++)
        for (j = 0; j < n; j++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(B, t, j, ctx);
            int isex = (x->u.size != 0 && x->N == PADIC_RADIX_EXACT);
            hB[j * kk + t] = isex ? x->v + e * FLINT_ABS(x->u.size) : -PADIC_RADIX_PACK_BIG;
            anyexB |= isex;
        }

    /* exponent bookkeeping for all outputs */
    L = 1;
    for (i = 0; i < m; i++)
    {
        const slong * avv = PA.vv + i * kk, * ave = PA.ve + i * kk;
        const slong * aNN = PA.NN + i * kk, * aww = PA.ww + i * kk;

        for (j = 0; j < n; j++)
        {
            const slong * bvv = PB.vv + j * kk, * bve = PB.ve + j * kk;
            const slong * bNN = PB.NN + j * kk, * bww = PB.ww + j * kk;
            slong vmin = 2 * PADIC_RADIX_PACK_BIG, Nt = 2 * PADIC_RADIX_PACK_BIG;
            slong W = 2 * PADIC_RADIX_PACK_BIG, N, s;

            for (t = 0; t < kk; t++)
            {
                slong x1 = avv[t] + bvv[t];
                slong e1 = aNN[t] + bve[t];
                slong e2 = ave[t] + bNN[t];
                slong f = aww[t] + bww[t];
                vmin = FLINT_MIN(vmin, x1);
                e1 = FLINT_MIN(e1, e2);
                Nt = FLINT_MIN(Nt, e1);
                W = FLINT_MIN(W, f);
            }

            if (vmin >= PADIC_RADIX_PACK_NONE)
            {
                /* no nonzero term: 0 + O(p^Nt), or an exact zero */
                N = (Nt >= PADIC_RADIX_PACK_NONE) ? PADIC_RADIX_EXACT : Nt;
                s = 0;
            }
            else if (Nt >= PADIC_RADIX_PACK_NONE && anyexA && anyexB)
            {
                /* only exact terms: the exact sum is kept if its a priori
                   extent is at most one limb beyond the precision; it is
                   then computed modulo B^s with a carry and a sign limb */
                slong H = -2 * PADIC_RADIX_PACK_BIG, ext, prec;
                const slong * ha = hA + i * kk, * hb = hB + j * kk;

                for (t = 0; t < kk; t++)
                    H = FLINT_MAX(H, ha[t] + hb[t]);

                N = _padic_radix_pack_prec(vmin, PADIC_RADIX_PACK_NONE, ctx);
                s = 0;
                if (N != PADIC_RADIX_EXACT && N > vmin)
                {
                    ext = (H - vmin + e - 1) / e + 1;
                    prec = (N - vmin + e - 1) / e;
                    if (ext <= prec + 1)
                    {
                        N = PADIC_RADIX_PACK_EXACT_SUM;
                        s = (H - e * W + e - 1) / e + 2;
                    }
                    else
                        s = (slong) radix_div_exp((ulong) (N - e * W) + e - 1, radix);
                }
            }
            else
            {
                N = _padic_radix_pack_prec(vmin, Nt, ctx);
                s = (N <= vmin) ? 0 : (slong) radix_div_exp((ulong) (N - e * W) + e - 1, radix);
            }

            PN[i * n + j] = N;
            PW[i * n + j] = W;
            Pslots[i * n + j] = s;
            L = FLINT_MAX(L, s);
        }
    }

    /* limb-aligned residues; at least two limbs each for the register path */
    stride = FLINT_MAX(L, 2);

    PA.U = TMP_ALLOC(sizeof(ulong) * stride * (nA + nB));
    PB.U = PA.U + stride * nA;

    for (i = 0; i < m; i++)
        for (t = 0; t < kk; t++)
        {
            slong q = i * kk + t;
            if (PA.vv[q] >= PADIC_RADIX_PACK_NONE)
            {
                PA.sz[q] = 0;
                flint_mpn_zero(PA.U + stride * q, stride);
            }
            else
            {
                PA.sz[q] = _padic_radix_pack_unit(PA.U + stride * q,
                    gr_mat_entry_srcptr(A, i, t, ctx), PA.ww[q], L, radix);
                if (stride > L)
                    PA.U[stride * q + 1] = 0;
            }
        }

    for (t = 0; t < kk; t++)
        for (j = 0; j < n; j++)
        {
            slong q = j * kk + t;
            if (PB.vv[q] >= PADIC_RADIX_PACK_NONE)
            {
                PB.sz[q] = 0;
                flint_mpn_zero(PB.U + stride * q, stride);
            }
            else
            {
                PB.sz[q] = _padic_radix_pack_unit(PB.U + stride * q,
                    gr_mat_entry_srcptr(B, t, j, ctx), PB.ww[q], L, radix);
                if (stride > L)
                    PB.U[stride * q + 1] = 0;
            }
        }

    acc = TMP_ALLOC(sizeof(ulong) * 3 * FLINT_MAX(L, 2));


    for (i = 0; i < m; i++)
    {
        const slong * aww = PA.ww + i * kk, * asz = PA.sz + i * kk;
        nn_srcptr ua = PA.U + stride * i * kk;

        for (j = 0; j < n; j++)
        {
            const slong * bww = PB.ww + j * kk, * bsz = PB.sz + j * kk;
            nn_srcptr ub = PB.U + stride * j * kk;
            padic_radix_struct * c = gr_mat_entry_ptr(C, i, j, ctx);
            slong N = PN[i * n + j], W = PW[i * n + j], nslots = Pslots[i * n + j];
            nn_ptr d;

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
                /* (no exact sums here: their windows have at least 3 limbs) */
                /* every window fits in two limbs: accumulate in registers;
                   zero entries have frame BIG and are skipped by the
                   offset test */
                ulong c0 = 0, c1 = 0, c2 = 0, d0 = 0, d1 = 0, d2 = 0, ph, pl;

                if (nslots == 1)
                {
                    for (t = 0; t < kk; t++)
                    {
                        if (aww[t] + bww[t] == W)
                        {
                            umul_ppmm(ph, pl, ua[stride * t], ub[stride * t]);
                            ADD3(c2, c1, c0, ph, pl);
                        }
                    }
                }
                else
                {
                    for (t = 0; t < kk; t++)
                    {
                        slong o = aww[t] + bww[t] - W;
                        ulong a0 = ua[stride * t], a1 = ua[stride * t + 1];
                        ulong b0 = ub[stride * t], b1 = ub[stride * t + 1];

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

                for (t = 0; t < kk; t++)
                {
                    slong an = asz[t], bn = bsz[t], o;

                    o = aww[t] + bww[t] - W;
                    if (o >= nslots || an == 0 || bn == 0)
                        continue;

                    if (an == 1 && bn == 1)
                    {
                        ulong ph, pl;
                        umul_ppmm(ph, pl, ua[stride * t], ub[stride * t]);
                        add_sssaaaaaa(acc[3 * o + 2], acc[3 * o + 1], acc[3 * o],
                                      acc[3 * o + 2], acc[3 * o + 1], acc[3 * o], 0, ph, pl);
                    }
                    else
                    {
                        _padic_radix_pack_addmul(acc, ua + stride * t, an, ub + stride * t, bn, o, nslots);
                    }
                }
            }

            d = radix_integer_fit_limbs(&c->u, nslots + 3, radix);
            c->u.size = _padic_radix_pack_normalise(d, acc, nslots, radix);

            if (N == PADIC_RADIX_PACK_EXACT_SUM)
            {
                /* the exact sum modulo B^nslots, |sum| < B^(nslots - 1) */
                slong sz = FLINT_MIN(c->u.size, nslots);
                int neg;

                flint_mpn_zero(d + sz, nslots - sz);
                neg = (d[nslots - 1] != 0);
                if (neg)
                    radix_neg(d, d, nslots, radix);
                sz = nslots;
                while (sz > 0 && d[sz - 1] == 0)
                    sz--;
                c->u.size = neg ? -sz : sz;
                c->v = e * W;
                c->N = PADIC_RADIX_EXACT;
                status |= _padic_radix_finalize(c, ctx);
                continue;
            }

            status |= _padic_radix_pack_finish(c, W, N, ctx);
        }
    }

    TMP_END;
    return status;

unable:
    TMP_END;
    return GR_UNABLE;
}

/*
    Matrix multiplication through fmpz_mat_mul (multimodular and BLAS-based
    algorithms) for large matrices. Row i of A is scaled to the base frame
    rw_i = min_t w(A_it) and column j of B to cw_j = min_t w(B_tj): the
    integer A'_it = U_it B^(w_it - rw_i) (the limb-aligned residue, truncated
    to the S_max limbs that can reach any output window), and similarly B',
    so that C_ij = (A' B')_ij B^(rw_i + cw_j). Each output is reduced modulo
    B^s (s limbs = its window) in binary before the conversion to radix B.

    The precisions are those of padic_radix_dot_strided_delayed, computed
    with blocked separable minima: with the inner index split into blocks,
    min_t (x_it + y_tj) = min_blocks (min_t x_it + min_t y_tj) whenever the
    exponents of row i of A, or those of column j of B, are constant on each
    block (e.g. uniform valuations and precisions), in O(mn * blocks).
    The remaining outputs use the exact O(k) minimum.
*/

#define MAT_FMPZ_BLOCKS 16

static int
_padic_radix_mat_mul_fmpz(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int cutoff, int alg, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong m = A->r, kk = A->c, n = B->c;
    slong i, j, t, b, nb, blk, nA, nB, Smax, Sown, binl, * hA, * hB;
    int pass;
    _pack_mat_struct PA, PB;
    padic_radix_pack_meta_struct meta;
    slong * ex, * Av, * Ave, * AN, * Bv, * Bve, * BN, * rw, * cw, * PN, * PS, * Ax, * Bx;
    char * rconst, * cconst;
    nn_ptr tmp, bin;
    fmpz_mat_t Az, Bz, Cz;
    fmpz * Bpow;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (B->r != kk || C->r != m || C->c != n)
        return GR_DOMAIN;

    if (PADIC_RADIX_CTX_PREC_REL(ctx) == PADIC_RADIX_PREC_INF
        && PADIC_RADIX_CTX_PREC_ABS(ctx) == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;

    if (m == 0 || n == 0)
        return GR_SUCCESS;

    if (kk == 0)
        return gr_mat_zero(C, ctx);

    nA = m * kk;
    nB = kk * n;
    nb = FLINT_MIN(kk, MAT_FMPZ_BLOCKS);
    blk = (kk + nb - 1) / nb;
    nb = (kk + blk - 1) / blk;

    TMP_START;

    ex = TMP_ALLOC(sizeof(slong) * (4 * (nA + nB) + 3 * nb * (m + n) + m + n + 2 * m * n));
    PA.vv = ex;           PA.ve = PA.vv + nA;   PA.NN = PA.ve + nA;   PA.ww = PA.NN + nA;
    PB.vv = PA.ww + nA;   PB.ve = PB.vv + nB;   PB.NN = PB.ve + nB;   PB.ww = PB.NN + nB;
    Av = PB.ww + nB;      Ave = Av + nb * m;    AN = Ave + nb * m;
    Bv = AN + nb * m;     Bve = Bv + nb * n;    BN = Bve + nb * n;
    rw = BN + nb * n;     cw = rw + m;
    PN = cw + n;          PS = PN + m * n;
    Ax = TMP_ALLOC(sizeof(slong) * nb * (m + n));
    Bx = Ax + nb * m;
    rconst = TMP_ALLOC(m + n);
    cconst = rconst + m;

    for (i = 0; i < m; i++)
        for (t = 0; t < kk; t++)
        {
            slong q = i * kk + t;
            if (!_padic_radix_pack_meta(&meta, gr_mat_entry_srcptr(A, i, t, ctx), radix))
                goto unable;
            PA.vv[q] = meta.vv; PA.ve[q] = meta.ve; PA.NN[q] = meta.NN; PA.ww[q] = meta.ww;
        }

    for (t = 0; t < kk; t++)
        for (j = 0; j < n; j++)
        {
            slong q = j * kk + t;
            if (!_padic_radix_pack_meta(&meta, gr_mat_entry_srcptr(B, t, j, ctx), radix))
                goto unable;
            PB.vv[q] = meta.vv; PB.ve[q] = meta.ve; PB.NN[q] = meta.NN; PB.ww[q] = meta.ww;
        }

    /* a priori extents v + e |u| of the exact nonzero entries */
    hA = TMP_ALLOC(sizeof(slong) * (nA + nB));
    hB = hA + nA;
    for (i = 0; i < m; i++)
        for (t = 0; t < kk; t++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(A, i, t, ctx);
            hA[i * kk + t] = (x->u.size != 0 && x->N == PADIC_RADIX_EXACT)
                ? x->v + e * FLINT_ABS(x->u.size) : -PADIC_RADIX_PACK_BIG;
        }
    for (t = 0; t < kk; t++)
        for (j = 0; j < n; j++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(B, t, j, ctx);
            hB[j * kk + t] = (x->u.size != 0 && x->N == PADIC_RADIX_EXACT)
                ? x->v + e * FLINT_ABS(x->u.size) : -PADIC_RADIX_PACK_BIG;
        }

    /* block minima and base frames */
    for (i = 0; i < m; i++)
    {
        rw[i] = 2 * PADIC_RADIX_PACK_BIG;
        rconst[i] = 1;
        for (b = 0; b < nb; b++)
        {
            slong x = 2 * PADIC_RADIX_PACK_BIG, y = x, z = x, xm = 0, q0 = i * kk + b * blk;
            for (t = b * blk; t < FLINT_MIN(kk, (b + 1) * blk); t++)
            {
                slong q = i * kk + t;
                x = FLINT_MIN(x, PA.vv[q]);
                xm = FLINT_MAX(xm, PA.vv[q]);
                y = FLINT_MIN(y, PA.ve[q]);
                z = FLINT_MIN(z, PA.NN[q]);
                rw[i] = FLINT_MIN(rw[i], PA.ww[q]);
                if (PA.vv[q] != PA.vv[q0] || PA.ve[q] != PA.ve[q0] || PA.NN[q] != PA.NN[q0])
                    rconst[i] = 0;
            }
Ax[i * nb + b] = xm; Av[i * nb + b] = x; Ave[i * nb + b] = y; AN[i * nb + b] = z;
        }
    }

    for (j = 0; j < n; j++)
    {
        cw[j] = 2 * PADIC_RADIX_PACK_BIG;
        cconst[j] = 1;
        for (b = 0; b < nb; b++)
        {
            slong x = 2 * PADIC_RADIX_PACK_BIG, y = x, z = x, xm = 0, q0 = j * kk + b * blk;
            for (t = b * blk; t < FLINT_MIN(kk, (b + 1) * blk); t++)
            {
                slong q = j * kk + t;
                x = FLINT_MIN(x, PB.vv[q]);
                xm = FLINT_MAX(xm, PB.vv[q]);
                y = FLINT_MIN(y, PB.ve[q]);
                z = FLINT_MIN(z, PB.NN[q]);
                cw[j] = FLINT_MIN(cw[j], PB.ww[q]);
                if (PB.vv[q] != PB.vv[q0] || PB.ve[q] != PB.ve[q0] || PB.NN[q] != PB.NN[q0])
                    cconst[j] = 0;
            }
Bx[j * nb + b] = xm; Bv[j * nb + b] = x; Bve[j * nb + b] = y; BN[j * nb + b] = z;
        }
    }

    /* in automatic mode, give up unless the dimensions are above the
       crossover with the packed kernel, estimated from an upper bound for
       the windows: v(C_ij) <= max_(t in block) v(A_it) + max_(t in block)
       v(B_tj) for every block, so the window of C_ij is at most
       rel + max_i (Ax_ib - e rw_i) + max_j (Bx_jb - e cw_j) digits */
    if (cutoff)
    {
        slong rel = PADIC_RADIX_CTX_PREC_REL(ctx), pabs = PADIC_RADIX_CTX_PREC_ABS(ctx);
        slong dim = FLINT_MIN(FLINT_MIN(m, n), kk), cut, best = WORD_MAX, nr = 0, nc = 0;
        slong rmin = WORD_MAX, cmin = WORD_MAX, Sest;

        for (i = 0; i < m; i++)
        {
            nr += !rconst[i];
            if (rw[i] < PADIC_RADIX_PACK_NONE)
                rmin = FLINT_MIN(rmin, rw[i]);
        }
        for (j = 0; j < n; j++)
        {
            nc += !cconst[j];
            if (cw[j] < PADIC_RADIX_PACK_NONE)
                cmin = FLINT_MIN(cmin, cw[j]);
        }

        if (rmin == WORD_MAX || cmin == WORD_MAX)
            goto unable;    /* a zero factor: nothing to gain */

        if (rel != PADIC_RADIX_PREC_INF)
        {
            for (b = 0; b < nb; b++)
            {
                slong X = WORD_MIN, Y = WORD_MIN;

                for (i = 0; i < m && X < PADIC_RADIX_PACK_NONE; i++)
                    if (rw[i] < PADIC_RADIX_PACK_NONE)
                        X = FLINT_MAX(X, Ax[i * nb + b] - e * rw[i]);
                for (j = 0; j < n && Y < PADIC_RADIX_PACK_NONE; j++)
                    if (cw[j] < PADIC_RADIX_PACK_NONE)
                        Y = FLINT_MAX(Y, Bx[j * nb + b] - e * cw[j]);

                if (X < PADIC_RADIX_PACK_NONE && Y < PADIC_RADIX_PACK_NONE)
                    best = FLINT_MIN(best, X + Y + rel);
            }
        }

        if (pabs != PADIC_RADIX_PREC_INF)
            best = FLINT_MIN(best, pabs - e * (rmin + cmin));

        Sest = (best >= PADIC_RADIX_PACK_NONE) ? WORD_MAX : (best + e - 1) / e;

        /* crossovers measured for uniform valuations: nmod_mat (one limb),
           mpn_mod_mat (up to 16 binary limbs), fmpz_mat */
        {
            double c, Lr;

            if (Sest == WORD_MAX)
                goto unable;

            c = (Sest <= 1) ? 8 : (Sest == 2) ? 32 : (Sest <= 8) ? 16 : (Sest <= 15) ? 28 : 48;

            /* windows wider than the precision (spread frames) make the
               products more expensive than in the packed kernel, which
               uses the frame of every output */
            Lr = (rel == PADIC_RADIX_PREC_INF) ? (double) Sest : (double) ((FLINT_MIN(rel, WORD_MAX / 4) + e - 1) / e);
            if (Sest > Lr)
                c *= pow((double) Sest / Lr, 1.5);

            /* outputs needing the exact O(k) precision computation */
            if (4 * nr * nc > m * n)
                c *= 1.5;

            cut = (slong) FLINT_MIN(c, 1e9);
        }

        if (dim < cut)
            goto unable;
    }

    /* precisions and windows: blocked minima first (exact for most
       outputs, a lower bound for the others), then the exact minimum for
       the outputs where neither the row nor the column is block-constant */
    Sown = 0;
    for (pass = 0; pass < 2; pass++)
    {
        slong nconst = 0;

        Smax = 0;
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
            {
                slong vmin = 2 * PADIC_RADIX_PACK_BIG, Nt = vmin, N, s = 0;

                if (rconst[i] || cconst[j])
                {
                    if (pass == 1)
                    {
                        Smax = FLINT_MAX(Smax, PS[i * n + j]);
                        continue;
                    }

                    for (b = 0; b < nb; b++)
                    {
                        slong x1 = Av[i * nb + b] + Bv[j * nb + b];
                        slong e1 = AN[i * nb + b] + Bve[j * nb + b];
                        slong e2 = Ave[i * nb + b] + BN[j * nb + b];
                        vmin = FLINT_MIN(vmin, x1);
                        e1 = FLINT_MIN(e1, e2);
                        Nt = FLINT_MIN(Nt, e1);
                    }
                }
                else if (pass == 0)
                {
                    nconst++;

                    for (b = 0; b < nb; b++)
                    {
                        slong x1 = Av[i * nb + b] + Bv[j * nb + b];
                        slong e1 = AN[i * nb + b] + Bve[j * nb + b];
                        slong e2 = Ave[i * nb + b] + BN[j * nb + b];
                        vmin = FLINT_MIN(vmin, x1);
                        e1 = FLINT_MIN(e1, e2);
                        Nt = FLINT_MIN(Nt, e1);
                    }
                }
                else
                {
                    const slong * avv = PA.vv + i * kk, * ave = PA.ve + i * kk, * aNN = PA.NN + i * kk;
                    const slong * bvv = PB.vv + j * kk, * bve = PB.ve + j * kk, * bNN = PB.NN + j * kk;

                    for (t = 0; t < kk; t++)
                    {
                        slong x1 = avv[t] + bvv[t];
                        slong e1 = aNN[t] + bve[t];
                        slong e2 = ave[t] + bNN[t];
                        vmin = FLINT_MIN(vmin, x1);
                        e1 = FLINT_MIN(e1, e2);
                        Nt = FLINT_MIN(Nt, e1);
                    }
                }

                if (vmin >= PADIC_RADIX_PACK_NONE || rw[i] >= PADIC_RADIX_PACK_NONE || cw[j] >= PADIC_RADIX_PACK_NONE)
                    N = (Nt >= PADIC_RADIX_PACK_NONE) ? PADIC_RADIX_EXACT : Nt;
                else if (Nt >= PADIC_RADIX_PACK_NONE)
                {
                    /* only exact terms: as in the packed version */
                    slong H = -2 * PADIC_RADIX_PACK_BIG, ext, prec;

                    for (t = 0; t < kk; t++)
                        H = FLINT_MAX(H, hA[i * kk + t] + hB[j * kk + t]);

                    N = _padic_radix_pack_prec(vmin, PADIC_RADIX_PACK_NONE, ctx);
                    if (N != PADIC_RADIX_EXACT && N > vmin)
                    {
                        ext = (H - vmin + e - 1) / e + 1;
                        prec = (N - vmin + e - 1) / e;
                        if (ext <= prec + 1)
                        {
                            N = PADIC_RADIX_PACK_EXACT_SUM;
                            s = (H - e * (rw[i] + cw[j]) + e - 1) / e + 2;
                        }
                        else
                            s = (slong) radix_div_exp((ulong) (N - e * (rw[i] + cw[j])) + e - 1, radix);
                    }
                }
                else
                {
                    N = _padic_radix_pack_prec(vmin, Nt, ctx);
                    if (N > vmin)
                        s = (slong) radix_div_exp((ulong) (N - e * (rw[i] + cw[j])) + e - 1, radix);
                }

                PN[i * n + j] = N;
                PS[i * n + j] = s;
                Smax = FLINT_MAX(Smax, s);

                /* the window of this output in its own frame */
                if (s != 0)
                    Sown = FLINT_MAX(Sown, s - (_fdiv(vmin, e) - rw[i] - cw[j]));
            }

        if (pass == 0)
        {
            if (nconst == 0)
                break;
        }
    }

    /* frames spread far beyond the precision (outputs whose windows in the
       separable frames rw_i + cw_j are much wider than in their own frames):
       the packed kernel, which uses the frame of every output, is cheaper,
       and the temporary integers could be arbitrarily large */
    if (Smax > 4 * Sown + 8)
        goto unable;

    if (Smax == 0)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * c = gr_mat_entry_ptr(C, i, j, ctx);
                radix_integer_zero(&c->u, radix);
                c->v = 0;
                c->N = PN[i * n + j];
                status |= _padic_radix_finalize(c, ctx);
            }
        goto finish;
    }

    /* single-limb windows: the products are only needed modulo B */
    if (Smax == 1)
    {
        nmod_mat_t An, Bn, Cn;
        ulong w;

        nmod_mat_init(An, m, kk, radix->B.n);
        nmod_mat_init(Bn, kk, n, radix->B.n);
        nmod_mat_init(Cn, m, n, radix->B.n);

        for (i = 0; i < m; i++)
            for (t = 0; t < kk; t++)
            {
                slong q = i * kk + t;
                if (PA.vv[q] < PADIC_RADIX_PACK_NONE && PA.ww[q] == rw[i])
                {
                    _padic_radix_pack_unit(&w, gr_mat_entry_srcptr(A, i, t, ctx), PA.ww[q], 1, radix);
                    nmod_mat_entry(An, i, t) = w;
                }
            }

        for (t = 0; t < kk; t++)
            for (j = 0; j < n; j++)
            {
                slong q = j * kk + t;
                if (PB.vv[q] < PADIC_RADIX_PACK_NONE && PB.ww[q] == cw[j])
                {
                    _padic_radix_pack_unit(&w, gr_mat_entry_srcptr(B, t, j, ctx), PB.ww[q], 1, radix);
                    nmod_mat_entry(Bn, t, j) = w;
                }
            }

        nmod_mat_mul(Cn, An, Bn);

        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * c = gr_mat_entry_ptr(C, i, j, ctx);
                slong N = PN[i * n + j];


                w = (PS[i * n + j] == 0) ? 0 : nmod_mat_entry(Cn, i, j);

                if (w == 0)
                {
                    radix_integer_zero(&c->u, radix);
                    c->v = 0;
                    c->N = N;
                    status |= _padic_radix_finalize(c, ctx);
                }
                else
                {
                    radix_integer_fit_limbs(&c->u, 1, radix)[0] = w;
                    c->u.size = 1;
                    status |= _padic_radix_pack_finish(c, rw[i] + cw[j], N, ctx);
                }
            }

        nmod_mat_clear(An);
        nmod_mat_clear(Bn);
        nmod_mat_clear(Cn);
        goto finish;
    }

    /* multi-limb windows of at most MPN_MOD_MAX_LIMBS binary limbs:
       arithmetic modulo B^Smax (classical with delayed reduction, Waksman,
       or multimodular, as chosen by mpn_mod_mat_mul) */
    if (alg == 2 || alg == 0)
    {
        fmpz_t M;
        gr_ctx_t mctx;
        gr_mat_t Am, Bm, Cm;
        slong nl;

        fmpz_init(M);
        fmpz_set_ui(M, radix->B.n);
        fmpz_pow_ui(M, M, Smax);
        nl = fmpz_size(M);

        if (nl < MPN_MOD_MIN_LIMBS || nl > MPN_MOD_MAX_LIMBS || gr_ctx_init_mpn_mod(mctx, M) != GR_SUCCESS)
        {
            fmpz_clear(M);
            if (alg == 2)
                goto unable;
            goto use_fmpz;
        }

        gr_mat_init(Am, m, kk, mctx);
        gr_mat_init(Bm, kk, n, mctx);
        gr_mat_init(Cm, m, n, mctx);

        tmp = TMP_ALLOC(sizeof(ulong) * (Smax + 1));
        bin = TMP_ALLOC(sizeof(ulong) * (radix_set_mpn_need_alloc(nl, radix) + Smax + 4));

        for (i = 0; i < m; i++)
            for (t = 0; t < kk; t++)
            {
                slong q = i * kk + t, o, bn;
                nn_ptr x = gr_mat_entry_ptr(Am, i, t, mctx);
                if (PA.vv[q] >= PADIC_RADIX_PACK_NONE)
                    continue;
                o = PA.ww[q] - rw[i];
                if (o >= Smax)
                    continue;
                flint_mpn_zero(tmp, o);
                _padic_radix_pack_unit(tmp + o, gr_mat_entry_srcptr(A, i, t, ctx), PA.ww[q], Smax - o, radix);
                bn = radix_get_mpn(bin, tmp, Smax, radix);
                flint_mpn_copyi(x, bin, FLINT_MIN(bn, nl));
            }

        for (t = 0; t < kk; t++)
            for (j = 0; j < n; j++)
            {
                slong q = j * kk + t, o, bn;
                nn_ptr x = gr_mat_entry_ptr(Bm, t, j, mctx);
                if (PB.vv[q] >= PADIC_RADIX_PACK_NONE)
                    continue;
                o = PB.ww[q] - cw[j];
                if (o >= Smax)
                    continue;
                flint_mpn_zero(tmp, o);
                _padic_radix_pack_unit(tmp + o, gr_mat_entry_srcptr(B, t, j, ctx), PB.ww[q], Smax - o, radix);
                bn = radix_get_mpn(bin, tmp, Smax, radix);
                flint_mpn_copyi(x, bin, FLINT_MIN(bn, nl));
            }

        status |= mpn_mod_mat_mul(Cm, Am, Bm, mctx);

        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * c = gr_mat_entry_ptr(C, i, j, ctx);
                nn_srcptr z = gr_mat_entry_srcptr(Cm, i, j, mctx);
                slong N = PN[i * n + j], sl = PS[i * n + j], zn = nl, dn, j0;
                nn_ptr d;

                while (zn > 0 && z[zn - 1] == 0)
                    zn--;

                if (N == PADIC_RADIX_PACK_EXACT_SUM)
                {
                    dn = (zn == 0) ? 0 : radix_set_mpn(bin, z, zn, radix);
                    status |= _exact_sum_finish(c, bin, FLINT_MIN(dn, sl), sl, rw[i] + cw[j], ctx);
                    continue;
                }

                if (sl == 0 || zn == 0)
                {
                    radix_integer_zero(&c->u, radix);
                    c->v = 0;
                    c->N = N;
                    status |= _padic_radix_finalize(c, ctx);
                    continue;
                }

                /* z < B^Smax: convert, then keep the s lowest limbs */
                dn = radix_set_mpn(bin, z, zn, radix);
                dn = FLINT_MIN(dn, sl);
                j0 = 0;
                while (j0 < dn && bin[j0] == 0)
                    j0++;
                while (dn > j0 && bin[dn - 1] == 0)
                    dn--;

                if (dn == j0)
                {
                    radix_integer_zero(&c->u, radix);
                    c->v = 0;
                    c->N = N;
                    status |= _padic_radix_finalize(c, ctx);
                    continue;
                }

                d = radix_integer_fit_limbs(&c->u, dn - j0, radix);
                flint_mpn_copyi(d, bin + j0, dn - j0);
                c->u.size = dn - j0;
                status |= _padic_radix_pack_finish(c, rw[i] + cw[j] + j0, N, ctx);
            }

        gr_mat_clear(Am, mctx);
        gr_mat_clear(Bm, mctx);
        gr_mat_clear(Cm, mctx);
        gr_ctx_clear(mctx);
        fmpz_clear(M);
        goto finish;
    }

use_fmpz:
    /* integer matrices */
    fmpz_mat_init(Az, m, kk);
    fmpz_mat_init(Bz, kk, n);
    fmpz_mat_init(Cz, m, n);

    binl = Smax + 2;
    tmp = TMP_ALLOC(sizeof(ulong) * (Smax + 1));
    bin = TMP_ALLOC(sizeof(ulong) * FLINT_MAX(binl, radix_set_mpn_need_alloc(2 * Smax + 4, radix) + 2));

    for (i = 0; i < m; i++)
        for (t = 0; t < kk; t++)
        {
            slong q = i * kk + t, o, bn;
            if (PA.vv[q] >= PADIC_RADIX_PACK_NONE)
                continue;
            o = PA.ww[q] - rw[i];
            if (o >= Smax)
                continue;
            flint_mpn_zero(tmp, o);
            _padic_radix_pack_unit(tmp + o, gr_mat_entry_srcptr(A, i, t, ctx), PA.ww[q], Smax - o, radix);
            bn = radix_get_mpn(bin, tmp, Smax, radix);
            fmpz_set_ui_array(fmpz_mat_entry(Az, i, t), bin, bn);
        }

    for (t = 0; t < kk; t++)
        for (j = 0; j < n; j++)
        {
            slong q = j * kk + t, o, bn;
            if (PB.vv[q] >= PADIC_RADIX_PACK_NONE)
                continue;
            o = PB.ww[q] - cw[j];
            if (o >= Smax)
                continue;
            flint_mpn_zero(tmp, o);
            _padic_radix_pack_unit(tmp + o, gr_mat_entry_srcptr(B, t, j, ctx), PB.ww[q], Smax - o, radix);
            bn = radix_get_mpn(bin, tmp, Smax, radix);
            fmpz_set_ui_array(fmpz_mat_entry(Bz, t, j), bin, bn);
        }

    fmpz_mat_mul(Cz, Az, Bz);

    /* outputs: reduce modulo B^s, convert, finish */
    Bpow = _fmpz_vec_init(Smax + 1);
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
            padic_radix_struct * c = gr_mat_entry_ptr(C, i, j, ctx);
            fmpz * z = fmpz_mat_entry(Cz, i, j);
            slong N = PN[i * n + j], s = PS[i * n + j], zn, dn, j0;
            nn_ptr d;

            if (s == 0)
            {
                radix_integer_zero(&c->u, radix);
                c->v = 0;
                c->N = N;
                status |= _padic_radix_finalize(c, ctx);
                continue;
            }

            if (fmpz_is_zero(Bpow + s))
            {
                fmpz_set_ui(Bpow + s, radix->B.n);
                fmpz_pow_ui(Bpow + s, Bpow + s, s);
            }

            fmpz_mod(z, z, Bpow + s);

            if (fmpz_is_zero(z) && N == PADIC_RADIX_PACK_EXACT_SUM)
            {
                status |= _exact_sum_finish(c, NULL, 0, s, rw[i] + cw[j], ctx);
                continue;
            }

            if (fmpz_is_zero(z))
            {
                radix_integer_zero(&c->u, radix);
                c->v = 0;
                c->N = N;
                status |= _padic_radix_finalize(c, ctx);
                continue;
            }

            zn = fmpz_size(z);
            fmpz_get_ui_array(bin, zn, z);

            if (N == PADIC_RADIX_PACK_EXACT_SUM)
            {
                nn_ptr dd = flint_malloc(sizeof(ulong) * (radix_set_mpn_need_alloc(zn, radix) + 1));
                dn = radix_set_mpn(dd, bin, zn, radix);
                status |= _exact_sum_finish(c, dd, FLINT_MIN(dn, s), s, rw[i] + cw[j], ctx);
                flint_free(dd);
                continue;
            }

            d = radix_integer_fit_limbs(&c->u, radix_set_mpn_need_alloc(zn, radix) + 1, radix);
            dn = radix_set_mpn(d, bin, zn, radix);

            j0 = 0;
            while (j0 < dn && d[j0] == 0)
                j0++;
            if (j0 > 0)
            {
                flint_mpn_copyi(d, d + j0, dn - j0);
                dn -= j0;
            }
            c->u.size = dn;
            status |= _padic_radix_pack_finish(c, rw[i] + cw[j] + j0, N, ctx);
        }

    _fmpz_vec_clear(Bpow, Smax + 1);
    fmpz_mat_clear(Az);
    fmpz_mat_clear(Bz);
    fmpz_mat_clear(Cz);

finish:
    TMP_END;
    return status;

unable:
    TMP_END;
    return GR_UNABLE;
}

int
padic_radix_mat_mul_fmpz(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return _padic_radix_mat_mul_fmpz(C, A, B, 0, 1, ctx);
}

int
padic_radix_mat_mul_mpn_mod(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return _padic_radix_mat_mul_fmpz(C, A, B, 0, 2, ctx);
}

int
padic_radix_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    int status;

    /* below this the dot products (or the addmul chains they use for very
       short lengths) are as fast */
    if (A->c <= 3 || A->r * B->c <= 4)
        return gr_mat_mul_classical(C, A, B, ctx);

    /* the fmpz route decides for itself after a cheap estimate of the
       windows; skip even that below the smallest crossovers, estimated from
       the relative precision of the context */
    status = GR_UNABLE;
    {
        slong dim = FLINT_MIN(FLINT_MIN(A->r, B->c), A->c);
        slong rel = PADIC_RADIX_CTX_PREC_REL(ctx), cut;
        slong e = PADIC_RADIX_CTX_RADIX(ctx)->exp;

        if (rel == PADIC_RADIX_PREC_INF)
            cut = 48;
        else
        {
            /* sample the valuations along a diagonal: a spread enlarges the
               windows and (being likely not block-constant) requires the
               exact precision bookkeeping, both raising the crossover */
            slong vlo = WORD_MAX, vhi = WORD_MIN, Lr, q;
            int any = 0;

            for (q = 0; q < FLINT_MIN(dim, 8); q++)
            {
                const padic_radix_struct * a = gr_mat_entry_srcptr(A, q, q, ctx);
                const padic_radix_struct * b = gr_mat_entry_srcptr(B, q, q, ctx);
                if (a->u.size != 0) { vlo = FLINT_MIN(vlo, a->v); vhi = FLINT_MAX(vhi, a->v); any = 1; }
                if (b->u.size != 0) { vlo = FLINT_MIN(vlo, b->v); vhi = FLINT_MAX(vhi, b->v); any = 1; }
            }

            if (!any)
                vhi = vlo = 0;
            else
                vhi = vlo + FLINT_MIN(vhi - vlo, WORD(1) << 20);

            Lr = (FLINT_MIN(rel, WORD_MAX / 4) + 2 * (vhi - vlo) + e - 1) / e;
            cut = (Lr <= 1) ? 8 : (Lr == 2) ? 32 : (Lr <= 8) ? 16 : (Lr <= 15) ? 28 : 48;
            if (vhi != vlo)
                cut = 3 * cut / 2;
        }

        if (dim >= cut)
            status = _padic_radix_mat_mul_fmpz(C, A, B, 1, 0, ctx);
    }
    if (status == GR_UNABLE)
        status = padic_radix_mat_mul_packed(C, A, B, ctx);
    if (status == GR_UNABLE)
        status = gr_mat_mul_classical(C, A, B, ctx);
    return status;
}
