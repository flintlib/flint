/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "impl.h"

/* Fixed-point sums sum_i (num_i / den_i) x_i with small integer num_i,
   den_i and (fn + 1)-limb fixed-point x_i (the terms of a q-series block,
   x_i = q^i from a shared table).

   Denominators are cleared in chunks: with D = prod den_i over a chunk,
   sum (num_i/den_i) x_i = (1/D) sum (num_i D/den_i) x_i, where the
   integer coefficients fit in one limb (D num_max < B = 2^FLINT_BITS) or
   two limbs (< B^2); the terms are accumulated with mpn_addmul_1 (one or two
   passes) and each chunk is divided once by D, the result truncated.

   The division of a long dividend by a one-limb divisor is bound by the
   latency of the remainder chain r -> (r B + x) mod D (one 2/1 step per
   limb).  The chunks give independent dividends, so the divisions of
   several chunks are interleaved: Moller-Granlund (2011) 2/1 steps for
   one-limb D, 3/2 steps for two-limb D.  Two-limb chunks hold about
   twice as many terms; they win once one-limb chunks shrink to single
   terms (2 log2 den + log2 num > FLINT_BITS), see choose_limbs.  All
   double-limb arithmetic uses the longlong.h macros. */

#define CHUNK_MAX 32
#define BATCH_K1 6      /* interleaved dividends, one-limb divisors */
#define BATCH_K2 4      /* interleaved dividends, two-limb divisors */
#define BATCH_KMAX 8

typedef struct
{
    ulong i, c;         /* terms i, ..., i + c - 1 */
    ulong d1, d0;       /* normalized divisor D 2^norm (d1 = 0: one limb) */
    ulong v;            /* inverse for the 2/1 or 3/2 division */
    unsigned int norm;
    ulong D1, D0;       /* D = D1 B + D0 (D1 = 0: one limb) */
}
chunk_t;

/* Moller-Granlund 2/1 division of u1 B + u0 (u1 < d) by normalized d */
#define DIV_2BY1(q, r, u1, u0, d, dinv)                                       \
    do {                                                                      \
        ulong _q1, _q0, _r;                                                   \
        umul_ppmm(_q1, _q0, (dinv), (u1));                                    \
        add_ssaaaa(_q1, _q0, _q1, _q0, (u1) + 1, (u0));                       \
        _r = (u0) - _q1 * (d);                                                \
        ulong _m = -(ulong) (_r > _q0);                                       \
        _q1 += _m;                                                            \
        _r += _m & (d);                                                       \
        if (FLINT_UNLIKELY(_r >= (d))) { _q1++; _r -= (d); }                  \
        (q) = _q1;                                                            \
        (r) = _r;                                                             \
    } while (0)

/* Moller-Granlund 3/2 division: (u2, u1, u0) / (d1, d0), (u2, u1) < (d1, d0),
   d1 normalized, v = floor((B^3 - 1) / d) - B */
#define DIV_3BY2(q, r1, r0, u2, u1, u0, d1, d0, v)                            \
    do {                                                                      \
        ulong _q1, _q0, _r1, _r0, _t1, _t0, _m;                               \
        umul_ppmm(_q1, _q0, (v), (u2));                                       \
        add_ssaaaa(_q1, _q0, _q1, _q0, (u2), (u1));                           \
        _r1 = (u1) - _q1 * (d1);                                              \
        umul_ppmm(_t1, _t0, (d0), _q1);                                       \
        sub_ddmmss(_r1, _r0, _r1, (u0), _t1, _t0);                            \
        sub_ddmmss(_r1, _r0, _r1, _r0, (d1), (d0));                           \
        _q1++;                                                                \
        _m = -(ulong) (_r1 >= _q0);                                           \
        _q1 += _m;                                                            \
        add_ssaaaa(_r1, _r0, _r1, _r0, (d1) & _m, (d0) & _m);                 \
        if (FLINT_UNLIKELY(_r1 >= (d1)))                                      \
        {                                                                     \
            if (_r1 > (d1) || _r0 >= (d0))                                    \
            {                                                                 \
                _q1++;                                                        \
                sub_ddmmss(_r1, _r0, _r1, _r0, (d1), (d0));                   \
            }                                                                 \
        }                                                                     \
        (q) = _q1;                                                            \
        (r1) = _r1;                                                           \
        (r0) = _r0;                                                           \
    } while (0)

static void
inverse_3by2(ulong * v, ulong d1, ulong d0)
{
    ulong ones[3], qq[2], rr[2], dd[2];
    dd[0] = d0;
    dd[1] = d1;
    ones[0] = ones[1] = ones[2] = UWORD_MAX;
    mpn_tdiv_qr(qq, rr, 0, ones, 3, dd, 2);     /* B + v */
    *v = qq[0];
}

/* (r2, r1, r0) = (a1, a0) m */
#define MUL_2X1(r2, r1, r0, a1, a0, m)                                        \
    do {                                                                      \
        ulong _h, _l;                                                         \
        umul_ppmm(_h, (r0), (a0), (m));                                       \
        umul_ppmm((r2), _l, (a1), (m));                                       \
        add_ssaaaa((r2), (r1), (r2), _l, 0, _h);                              \
    } while (0)

/* whether (e2, e1, e0) m < B^limbs */
static int
fits_limbs(ulong e2, ulong e1, ulong e0, ulong m, int limbs)
{
    ulong h0, l0, h1, l1, h2, l2, f1, f2, f3, c;

    umul_ppmm(h0, l0, e0, m);
    umul_ppmm(h1, l1, e1, m);
    umul_ppmm(h2, l2, e2, m);
    (void) l0;
    f1 = l1 + h0;
    c = (f1 < l1);
    f2 = l2 + h1;
    f3 = h2 + (f2 < l2);
    f2 += c;
    f3 += (f2 < c);
    if (limbs == 2)
        return f2 == 0 && f3 == 0;
    return f1 == 0 && f2 == 0 && f3 == 0;
}

/* the product of den_i over a chunk, as long as D tmax stays below
   B^limbs, tmax = max(1, num_i mod den_i) */
static slong
make_chunks(chunk_t * ch, const ulong * num, const ulong * den, ulong len, int limbs)
{
    ulong i, c;
    slong nch = 0;

    for (i = 0; i < len; i += c)
    {
        ulong D1 = 0, D0 = 1, tmax = 0;

        for (c = 0; i + c < len && c < CHUNK_MAX; c++)
        {
            ulong M = den[i + c];
            ulong tm = FLINT_MAX(tmax, FLINT_MAX(num[i + c] % M, 1));
            ulong E2, E1, E0;

            /* E = D M; stop before D M tmax reaches B^limbs (tm >= 1,
               so that E fits as well) */
            MUL_2X1(E2, E1, E0, D1, D0, M);
            if (c > 0 && !fits_limbs(E2, E1, E0, tm, limbs))
                break;
            D1 = E1;
            D0 = E0;
            tmax = tm;
        }

        if (ch != NULL)
        {
            ch[nch].i = i;
            ch[nch].c = c;
            ch[nch].D1 = D1;
            ch[nch].D0 = D0;
            if (limbs == 1)
            {
                ch[nch].norm = flint_clz(D0);
                ch[nch].d1 = 0;
                ch[nch].d0 = D0 << ch[nch].norm;
                ch[nch].v = n_preinvert_limb_prenorm(ch[nch].d0);
            }
            else
            {
                unsigned int norm = D1 ? flint_clz(D1) : FLINT_BITS + flint_clz(D0);
                ulong n1, n0;
                if (norm >= FLINT_BITS)
                {
                    n1 = D0 << (norm - FLINT_BITS);
                    n0 = 0;
                }
                else if (norm == 0)
                {
                    n1 = D1;
                    n0 = D0;
                }
                else
                {
                    n1 = (D1 << norm) | (D0 >> (FLINT_BITS - norm));
                    n0 = D0 << norm;
                }
                ch[nch].norm = norm;
                ch[nch].d1 = n1;
                ch[nch].d0 = n0;
                inverse_3by2(&ch[nch].v, n1, n0);
            }
        }
        nch++;
    }

    return nch;
}

slong
_acb_mzeta_fixed_ratsum_chunks_size(ulong len)
{
    return (len + 1) * sizeof(chunk_t);
}

slong
_acb_mzeta_fixed_ratsum_work_size(slong fn)
{
    return _acb_mzeta_fixed_ratsum_multi_work_size(fn, 1);
}

slong
_acb_mzeta_fixed_ratsum_multi_work_size(slong fn, slong P)
{
    /* P open accumulators, two queues of BATCH_KMAX queued dividends,
       BATCH_KMAX quotients and a zero dividend for the padding */
    return (P + 3 * BATCH_KMAX + 1) * (fn + 3);
}

/* Per-block choice between one- and two-limb divisors, from measured costs
   per limb (x86-64, GMP built for the CPU): a term costs one (two)
   mpn_addmul_1 passes at ~0.46 ns, a chunk ~1.95 ns (one-limb: batched
   2/1 division plus mpn_add_n) or ~3.1 ns (two-limb). */
int
_acb_mzeta_fixed_ratsum_choose_limbs(const ulong * num, const ulong * den, ulong len)
{
    double c1 = 0.46 * len + 1.95 * make_chunks(NULL, num, den, len, 1);
    double c2 = 0.92 * len + 3.1 * make_chunks(NULL, num, den, len, 2);
    return (c2 < 0.85 * c1) ? 2 : 1;
}

/* a queue of chunks waiting for their (interleaved) division */
typedef struct
{
    slong cnt;
    ulong * A[BATCH_KMAX];      /* dividends */
    slong p[BATCH_KMAX];        /* the sums they belong to */
    ulong d1[BATCH_KMAX], d0[BATCH_KMAX], v[BATCH_KMAX];
}
divqueue_t;

/* divides the queued dividends (padded to K), adds the truncated
   quotients to their sums and returns the buffers to the pool */
static void
queue_flush(divqueue_t * q, int limbs, ulong * bs, ulong * err, slong n,
    ulong ** Qb, ulong * zero, ulong ** pool, slong * npool)
{
    slong K = (limbs == 1) ? BATCH_K1 : BATCH_K2, i, j;

    if (q->cnt == 0)
        return;

    for (j = q->cnt; j < K; j++)
    {
        /* padding: zero dividend, any normalized divisor */
        q->A[j] = zero;
        q->d1[j] = (limbs == 2) ? UWORD(1) << (FLINT_BITS - 1) : 0;
        q->d0[j] = (limbs == 2) ? 0 : UWORD(1) << (FLINT_BITS - 1);
        if (limbs == 1)
            q->v[j] = n_preinvert_limb_prenorm(q->d0[j]);
        else
            inverse_3by2(&q->v[j], q->d1[j], q->d0[j]);
    }

    /* K interleaved long divisions, top limb down */
    if (limbs == 1)
    {
        ulong R[BATCH_KMAX] = {0};
        for (i = n - 1; i >= 0; i--)
            for (j = 0; j < BATCH_K1; j++)
                DIV_2BY1(Qb[j][i], R[j], R[j], q->A[j][i], q->d0[j], q->v[j]);
    }
    else
    {
        ulong R1[BATCH_KMAX] = {0}, R0[BATCH_KMAX] = {0};
        for (i = n - 1; i >= 0; i--)
            for (j = 0; j < BATCH_K2; j++)
                DIV_3BY2(Qb[j][i], R1[j], R0[j], R1[j], R0[j], q->A[j][i], q->d1[j], q->d0[j], q->v[j]);
    }

    /* truncated quotients: one ulp each */
    for (j = 0; j < q->cnt; j++)
    {
        mpn_add_n(bs + q->p[j] * n, bs + q->p[j] * n, Qb[j], n);
        err[q->p[j]] += 1;
        pool[(*npool)++] = q->A[j];
    }
    q->cnt = 0;
}

void
_acb_mzeta_fixed_ratsum_multi(ulong * bs, ulong * err, const ulong * xtab, slong stride,
    slong fn, const ulong * num, const ulong * den, ulong len, slong P, ulong xerr,
    const int * limbs, void * vchunks, ulong * work)
{
    chunk_t * ch = vchunks;
    slong n = fn + 3, p, j, npool = 0;
    ulong u;
    slong * cur, * nch;
    ulong ** open, ** pool, * Qb[BATCH_KMAX], * zero;
    int * used;
    divqueue_t queue[2];

    flint_mpn_zero(bs, n * P);
    for (p = 0; p < P; p++)
        err[p] = 0;

    cur = flint_malloc(sizeof(slong) * P);
    nch = flint_malloc(sizeof(slong) * P);
    used = flint_malloc(sizeof(int) * P);
    open = flint_malloc(sizeof(ulong *) * P);
    pool = flint_malloc(sizeof(ulong *) * (P + 2 * BATCH_KMAX));

    /* buffers: P open accumulators and 2 BATCH_KMAX queued ones, then the
       quotients */
    for (j = 0; j < P + 2 * BATCH_KMAX; j++)
        pool[npool++] = work + j * n;
    for (j = 0; j < BATCH_KMAX; j++)
        Qb[j] = work + (P + 2 * BATCH_KMAX + j) * n;
    zero = work + (P + 3 * BATCH_KMAX) * n;
    flint_mpn_zero(zero, n);
    queue[0].cnt = queue[1].cnt = 0;

    for (p = 0; p < P; p++)
    {
        nch[p] = make_chunks(ch + p * (len + 1), num + p * len, den + p * len, len, limbs[p]);
        cur[p] = 0;
        used[p] = 0;
        open[p] = pool[--npool];
        flint_mpn_zero(open[p], n);
    }

    /* the terms in order, all the sums for each term (the table entry is
       read from memory once) */
    for (u = 0; u < len; u++)
    {
        const ulong * tp = xtab + u * stride;

        for (p = 0; p < P; p++)
        {
            const chunk_t * C = ch + p * (len + 1) + cur[p];
            ulong M = den[p * len + u], nu = num[p * len + u], q = nu / M, t = nu % M;
            ulong * A = open[p];

            /* integer part: q x added exactly (error q xerr) */
            if (q != 0)
            {
                ulong cy = mpn_addmul_1(bs + p * n, tp, fn + 1, q);
                mpn_add_1(bs + p * n + fn + 1, bs + p * n + fn + 1, 2, cy);
                err[p] += q * xerr;
            }

            /* t < M, so t D/M < D and the shifted coefficient fits */
            if (t != 0)
            {
                if (limbs[p] == 1)
                {
                    ulong cf = (t * (C->D0 / M)) << C->norm;
                    ulong cy = mpn_addmul_1(A, tp, fn + 1, cf);
                    mpn_add_1(A + fn + 1, A + fn + 1, 2, cy);
                }
                else
                {
                    /* coef = t (D / M) 2^norm < B^2 (exact division) */
                    ulong Q[2], DD[2], c2, c1, c0, cy;
                    unsigned int nm = C->norm;
                    DD[0] = C->D0;
                    DD[1] = C->D1;
                    mpn_divexact_1(Q, DD, 2, M);
                    MUL_2X1(c2, c1, c0, Q[1], Q[0], t);
                    (void) c2;
                    if (nm >= FLINT_BITS)
                    {
                        c1 = c0 << (nm - FLINT_BITS);
                        c0 = 0;
                    }
                    else if (nm != 0)
                    {
                        c1 = (c1 << nm) | (c0 >> (FLINT_BITS - nm));
                        c0 = c0 << nm;
                    }
                    /* two mpn_addmul_1 passes (faster than mpn_addmul_2 here) */
                    cy = mpn_addmul_1(A, tp, fn + 1, c0);
                    mpn_add_1(A + fn + 1, A + fn + 1, 2, cy);
                    cy = mpn_addmul_1(A + 1, tp, fn + 1, c1);
                    A[fn + 2] += cy;
                }
                /* (t/M) times an error of xerr ulps */
                err[p] += xerr;
                used[p] = 1;
            }

            /* the end of the chunk: queue its division */
            if (u == C->i + C->c - 1)
            {
                if (used[p])
                {
                    int l = limbs[p];
                    divqueue_t * qq = queue + (l - 1);
                    slong k = qq->cnt;
                    qq->A[k] = A;
                    qq->p[k] = p;
                    qq->d1[k] = C->d1;
                    qq->d0[k] = C->d0;
                    qq->v[k] = C->v;
                    qq->cnt = k + 1;
                    if (qq->cnt == ((l == 1) ? BATCH_K1 : BATCH_K2))
                        queue_flush(qq, l, bs, err, n, Qb, zero, pool, &npool);
                    open[p] = pool[--npool];
                    flint_mpn_zero(open[p], n);
                    used[p] = 0;
                }
                cur[p]++;
            }
        }
    }

    queue_flush(queue + 0, 1, bs, err, n, Qb, zero, pool, &npool);
    queue_flush(queue + 1, 2, bs, err, n, Qb, zero, pool, &npool);

    flint_free(cur);
    flint_free(nch);
    flint_free(used);
    flint_free(open);
    flint_free(pool);
}

ulong
_acb_mzeta_fixed_ratsum(ulong * bs, const ulong * xtab, slong stride, slong fn,
    const ulong * num, const ulong * den, ulong len, ulong xerr, int limbs,
    void * vchunks, ulong * work)
{
    ulong err;
    _acb_mzeta_fixed_ratsum_multi(bs, &err, xtab, stride, fn, num, den, len, 1, xerr,
        &limbs, vchunks, work);
    return err;
}
