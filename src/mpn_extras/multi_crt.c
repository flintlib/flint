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
#include "mpn_extras.h"
#include "thread_support.h"
#include "mpn_extras/impl.h"

/* reduced residue of leaf j: (sum_t r_t coeff_t) mod l_j, which
   equals x (M_i/l_j)^-1 w_i mod l_j for the chunk containing the leaf */
FLINT_FORCE_INLINE ulong
_leaf_residue(nn_srcptr res, const flint_mpn_crt_t C, slong j)
{
    slong l0 = C->leaf_start[j];
    slong l1 = C->leaf_start[j + 1];
    ulong hi, lo, u;

    if (C->crt_use_shoup)
        return n_mulmod_shoup(C->prime_crt_coeff[l0], res[l0], C->prime_crt_shoup[l0], C->leaf[j].mod.n);

    umul_ppmm(hi, lo, res[l0], C->prime_crt_coeff[l0]);
    for (l0++; l0 < l1; l0++)
        MAC(hi, lo, res[l0], C->prime_crt_coeff[l0]);

    /* hi < leaf modulus, see crt_init.c */
    if (C->leaf[j].mod.norm != 0)
        return _n_ll_rem_l_nonfullword(hi, lo, C->leaf[j].mod.n, C->leaf[j].qhi, C->leaf[j].qlo);
    NMOD_RED2_FULLWORD(u, hi, lo, C->leaf[j].mod);
    return u;
}

/* fixed-length CRT: out = sum_j u_j (P / l_j) reduced mod P, all in
   N = prod_len limbs */
#define DEFINE_FIXED(N, M) \
static void CAT3(_crt_fixed, N, M)(nn_ptr out, nn_srcptr res, const flint_mpn_crt_t C) \
{ \
    ulong r[N], t[N]; \
    slong j, nl = C->num_leaves; \
    CAT3(_big_mul, N, M)(r, t, C->crt_leaf_mult, _leaf_residue(res, C, 0)); \
    for (j = 1; j < nl; j++) \
        CAT3(_big_addmul, N, M)(r, t, C->crt_leaf_mult + j * N, _leaf_residue(res, C, j)); \
    CAT(_reduce_big_sum, N)(r, t, C->prod); \
    for (j = 0; j < N; j++) \
        out[j] = r[j]; \
}
DEFINE_FIXED(1, 1)
DEFINE_FIXED(2, 1) DEFINE_FIXED(2, 2)
DEFINE_FIXED(3, 2) DEFINE_FIXED(3, 3)
DEFINE_FIXED(4, 3) DEFINE_FIXED(4, 4)
DEFINE_FIXED(5, 4) DEFINE_FIXED(5, 5)
DEFINE_FIXED(6, 5) DEFINE_FIXED(6, 6)
DEFINE_FIXED(7, 6) DEFINE_FIXED(7, 7)
DEFINE_FIXED(8, 7) DEFINE_FIXED(8, 8)
#undef DEFINE_FIXED

static void
_crt_fixed(nn_ptr out, nn_srcptr res, const flint_mpn_crt_t C)
{
    switch (2 * C->prod_len + (C->crt_fixed_m == C->prod_len))
    {
        case 2: case 3: _crt_fixed_1_1(out, res, C); break;
        case 4: _crt_fixed_2_1(out, res, C); break;
        case 5: _crt_fixed_2_2(out, res, C); break;
        case 6: _crt_fixed_3_2(out, res, C); break;
        case 7: _crt_fixed_3_3(out, res, C); break;
        case 8: _crt_fixed_4_3(out, res, C); break;
        case 9: _crt_fixed_4_4(out, res, C); break;
        case 10: _crt_fixed_5_4(out, res, C); break;
        case 11: _crt_fixed_5_5(out, res, C); break;
        case 12: _crt_fixed_6_5(out, res, C); break;
        case 13: _crt_fixed_6_6(out, res, C); break;
        case 14: _crt_fixed_7_6(out, res, C); break;
        case 15: _crt_fixed_7_7(out, res, C); break;
        case 16: _crt_fixed_8_7(out, res, C); break;
        case 17: _crt_fixed_8_8(out, res, C); break;
        default: FLINT_UNREACHABLE;
    }
}

/* y = x w_i mod M_i (unreduced) for chunk i; returns the length.
   Single-prime leaves: y = sum_t r_t V_t with all cofactors folded into
   the multipliers V_t, y < 2^64 b M_i. Batched leaves: the residues of each
   leaf are first combined into u_j = x (M_i/l_j)^-1 w_i mod l_j with
   single-limb arithmetic, then y = sum_j u_j (M_i/l_j) < b M_i. In both
   cases y fits in cl + 2 limbs, and the multipliers are zero-padded to cl
   limbs so that the carries can be accumulated separately. */
static slong
_crt_basecase(nn_ptr y, nn_srcptr res, const flint_mpn_crt_t C, slong i)
{
    slong cl = C->crt_chunk_limbs;
    slong yn = cl + 2;
    slong j0 = C->chunk_start[i];
    slong j1 = C->chunk_start[i + 1];
    slong j;
    ulong hi = 0, lo = 0, cy;

    flint_mpn_zero(y, cl);

    if (C->all_single)
    {
        slong t0 = C->leaf_start[j0];
        slong t1 = C->leaf_start[j1];

        for (j = t0; j < t1; j++)
        {
            FLINT_ASSERT(res[j] < C->primes[j]);
            cy = mpn_addmul_1(y, C->crt_mult + j * cl, cl, res[j]);
            add_ssaaaa(hi, lo, hi, lo, UWORD(0), cy);
        }
    }
    else
    {
        for (j = j0; j < j1; j++)
        {
            ulong u = _leaf_residue(res, C, j);
            cy = mpn_addmul_1(y, C->crt_leaf_mult + j * cl, cl, u);
            add_ssaaaa(hi, lo, hi, lo, UWORD(0), cy);
        }
    }

    y[cl] = lo;
    y[cl + 1] = hi;

    MPN_NORM(y, yn);
    return yn;
}

typedef struct
{
    const flint_mpn_crt_struct * C;
    slong k;
    nn_ptr * ptr;
    slong * len;
    nn_ptr * ptr2;
    slong * len2;
    nn_ptr buf;
    nn_ptr scratch;
    int parallel;
}
_crt_ascend_arg;

/* A = B * v_c + C * v_b for node i at level k */
static void
_crt_ascend_worker(slong i, void * varg)
{
    _crt_ascend_arg * arg = (_crt_ascend_arg *) varg;
    const flint_mpn_crt_struct * C = arg->C;
    slong k = arg->k;
    slong cl = C->level_limbs[k] + FLINT_MPN_CRT_SLOT_EXTRA;
    slong pl = C->level_limbs[k - 1];
    slong ccount = C->level_count[k - 1];
    nn_ptr A = arg->buf + i * cl;
    nn_srcptr B = arg->ptr[2 * i];
    slong bn = arg->len[2 * i];
    slong An;

    if (2 * i + 1 < ccount)
    {
        nn_srcptr Cc = arg->ptr[2 * i + 1];
        slong cn = arg->len[2 * i + 1];
        nn_srcptr vb = C->level_prod[k - 1] + (2 * i) * pl;
        slong vbn = C->level_len[k - 1][2 * i];
        nn_srcptr vc = C->level_prod[k - 1] + (2 * i + 1) * pl;
        slong vcn = C->level_len[k - 1][2 * i + 1];
        slong tn;
        nn_ptr t = arg->parallel ? FLINT_ARRAY_ALLOC(cl, ulong) : arg->scratch;

        if (bn == 0)
        {
            An = 0;
        }
        else
        {
            if (bn >= vcn)
                flint_mpn_mul(A, B, bn, vc, vcn);
            else
                flint_mpn_mul(A, vc, vcn, B, bn);
            An = bn + vcn;
        }

        if (cn != 0)
        {
            if (cn >= vbn)
                flint_mpn_mul(t, Cc, cn, vb, vbn);
            else
                flint_mpn_mul(t, vb, vbn, Cc, cn);
            tn = cn + vbn;

            if (An == 0)
            {
                flint_mpn_copyi(A, t, tn);
                An = tn;
            }
            else if (An >= tn)
            {
                ulong cy = mpn_add(A, A, An, t, tn);
                if (cy)
                    A[An++] = cy;
            }
            else
            {
                ulong cy = mpn_add(A, t, tn, A, An);
                An = tn;
                if (cy)
                    A[An++] = cy;
            }
        }

        MPN_NORM(A, An);
        FLINT_ASSERT(An <= C->level_len[k][i] + 2);

        if (arg->parallel)
            flint_free(t);
    }
    else
    {
        flint_mpn_copyi(A, B, bn);
        An = bn;
    }

    arg->ptr2[i] = A;
    arg->len2[i] = An;
}


/*
    Small-value shortcut. Reconstructions are often much smaller than P
    (e.g. matrix entries of very different sizes). We reconstruct the
    value modulo the product Q of the first few primes with the fixed-length
    code; if the result has at least a limb of slack below Q (which a random
    residue mod Q has with probability 2^-64), it is verified against all
    residues, which costs only num_primes single-limb divisions of a small
    number. Returns 1 if out has been set.
*/
static int
_crt_small_shortcut(nn_ptr out, int * negative, nn_srcptr res, const flint_mpn_crt_t C, int sign)
{
    const flint_mpn_crt_struct * S = C->crt_small;
    ulong c[FLINT_MPN_CRT_FIXED_MAX];
    slong cn = S->prod_len, t;
    int neg;

    neg = flint_mpn_multi_crt(c, res, S, sign, NULL);
    MPN_NORM(c, cn);

    if (cn >= S->prod_len)
        return 0;

    for (t = 0; t < C->num_primes; t++)
    {
        ulong r = (cn == 0) ? 0 : mpn_mod_1(c, cn, C->primes[t]);
        if (neg && r != 0)
            r = C->primes[t] - r;
        if (r != res[t])
            return 0;
    }

    flint_mpn_copyi(out, c, cn);
    flint_mpn_zero(out + cn, C->prod_len - cn);
    *negative = neg;
    return 1;
}

int
flint_mpn_multi_crt(nn_ptr out, nn_srcptr res, const flint_mpn_crt_t C, int sign, nn_ptr tmp)
{
    _crt_work_struct W;
    slong top = C->num_levels - 1;
    slong k, i;
    nn_ptr * ptr, * ptr2;
    slong * len, * len2;
    nn_ptr buf, buf2;
    slong an, pn = C->prod_len;
    nn_ptr a;
    int negative = 0;
    TMP_INIT;

    FLINT_ASSERT(C->flags & FLINT_MPN_CRT_CRT);

    if (C->crt_fixed)
    {
        _crt_fixed(out, res, C);
        goto finish;
    }

    if (C->crt_small != NULL && _crt_small_shortcut(out, &negative, res, C, sign))
        return negative;

    TMP_START;
    if (tmp == NULL)
        tmp = TMP_ALLOC(C->tmp_limbs * sizeof(ulong));

    _crt_work_init(&W, C, tmp);

    ptr = W.ptrA;
    len = W.lenA;
    ptr2 = W.ptrB;
    len2 = W.lenB;
    buf = W.bufA;
    buf2 = W.bufB;

    /* chunk values */
    {
        slong cl = C->level_limbs[0] + FLINT_MPN_CRT_SLOT_EXTRA;
        for (i = 0; i < C->num_chunks; i++)
        {
            ptr[i] = buf + i * cl;
            len[i] = _crt_basecase(buf + i * cl, res, C, i);
        }
    }

    /* ascend */
    for (k = 1; k <= top; k++)
    {
        _crt_ascend_arg arg;

        arg.C = C;
        arg.k = k;
        arg.ptr = ptr;
        arg.len = len;
        arg.ptr2 = ptr2;
        arg.len2 = len2;
        arg.buf = buf2;
        arg.scratch = W.scratch;
        arg.parallel = _crt_use_threads(C->level_count[k], C->level_limbs[k]);

        if (arg.parallel)
            flint_parallel_do(_crt_ascend_worker, &arg, C->level_count[k], -1, FLINT_PARALLEL_UNIFORM);
        else
            for (i = 0; i < C->level_count[k]; i++)
                _crt_ascend_worker(i, &arg);

        FLINT_SWAP(nn_ptr *, ptr, ptr2);
        FLINT_SWAP(slong *, len, len2);
        FLINT_SWAP(nn_ptr, buf, buf2);
    }

    /* final reduction */
    a = ptr[0];
    an = len[0];

    if (an < pn || (an == pn && mpn_cmp(a, C->prod, pn) < 0))
    {
        flint_mpn_copyi(out, a, an);
        flint_mpn_zero(out + an, pn - an);
    }
    else
    {
        FLINT_ASSERT(an - pn + 1 <= 3);
        mpn_tdiv_qr(W.scratch, out, 0, a, an, C->prod, pn);
    }

    FLINT_ASSERT(mpn_cmp(out, C->prod, pn) < 0);

    TMP_END;

finish:
    if (sign)
    {
        /* negative iff 2 x > P */
        if (mpn_cmp(out, C->prod_half, pn) > 0)
        {
            mpn_sub_n(out, C->prod, out, pn);
            negative = 1;
        }
    }

    return negative;
}

/* block of entries whose values fit comfortably in L1 */
#define VEC_BLOCK_LIMBS 1024

void
flint_mpn_multi_crt_vec(nn_ptr out, slong out_stride, int * negative, nn_srcptr res, slong res_stride,
                        slong len, const flint_mpn_crt_t C, int sign, nn_ptr tmp)
{
    slong n = C->num_primes, pn = C->prod_len;
    slong i, l, t, i0, blk;
    nn_ptr rbuf;
    TMP_INIT;

    FLINT_ASSERT(C->flags & FLINT_MPN_CRT_CRT);

    TMP_START;

    if (C->crt_fixed || C->num_levels > 1 || !C->all_single)
    {
        /* entry by entry (the fixed-length code is already table-light) */
        rbuf = TMP_ALLOC(n * sizeof(ulong));
        for (i = 0; i < len; i++)
        {
            int neg;
            for (l = 0; l < n; l++)
                rbuf[l] = res[l * res_stride + i];
            neg = flint_mpn_multi_crt(out + i * out_stride, rbuf, C, sign, tmp);
            if (negative != NULL)
                negative[i] = neg;
        }
        TMP_END;
        return;
    }

    /* single chunk: prime-outer accumulation within blocks of entries,
       then one reduction per entry (entries caught by the small-value
       shortcut are handled individually) */
    {
        slong cl = C->crt_chunk_limbs;
        slong yl = cl + 2;
        nn_ptr Y, hi, lo;
        int * done;

        blk = FLINT_MAX(1, VEC_BLOCK_LIMBS / yl);
        Y = TMP_ALLOC(blk * yl * sizeof(ulong));
        hi = TMP_ALLOC(2 * blk * sizeof(ulong));
        lo = hi + blk;
        done = TMP_ALLOC(blk * sizeof(int));
        rbuf = TMP_ALLOC(n * sizeof(ulong));

        for (i0 = 0; i0 < len; i0 += blk)
        {
            slong i1 = FLINT_MIN(i0 + blk, len);
            slong m = i1 - i0;
            slong todo = m;

            /* small-value shortcut, entry by entry */
            for (i = 0; i < m; i++)
            {
                done[i] = 0;
                if (C->crt_small != NULL)
                {
                    int neg = 0;
                    for (l = 0; l < n; l++)
                        rbuf[l] = res[l * res_stride + i0 + i];
                    if (_crt_small_shortcut(out + (i0 + i) * out_stride, &neg, rbuf, C, sign))
                    {
                        if (negative != NULL)
                            negative[i0 + i] = neg;
                        done[i] = 1;
                        todo--;
                    }
                }
            }

            if (todo == 0)
                continue;

            flint_mpn_zero(Y, m * yl);
            flint_mpn_zero(hi, 2 * blk);

            for (t = 0; t < n; t++)
            {
                nn_srcptr V = C->crt_mult + t * cl;
                nn_srcptr r = res + t * res_stride + i0;

                for (i = 0; i < m; i++)
                {
                    ulong cy;
                    if (done[i])
                        continue;
                    cy = mpn_addmul_1(Y + i * yl, V, cl, r[i]);
                    add_ssaaaa(hi[i], lo[i], hi[i], lo[i], UWORD(0), cy);
                }
            }

            for (i = 0; i < m; i++)
            {
                nn_ptr y = Y + i * yl;
                nn_ptr o = out + (i0 + i) * out_stride;
                slong yn = yl;
                int neg = 0;

                if (done[i])
                    continue;

                y[cl] = lo[i];
                y[cl + 1] = hi[i];
                MPN_NORM(y, yn);

                if (yn < pn || (yn == pn && mpn_cmp(y, C->prod, pn) < 0))
                {
                    flint_mpn_copyi(o, y, yn);
                    flint_mpn_zero(o + yn, pn - yn);
                }
                else
                {
                    ulong q[4];
                    FLINT_ASSERT(yn - pn + 1 <= 3);
                    mpn_tdiv_qr(q, o, 0, y, yn, C->prod, pn);
                }
                FLINT_ASSERT(mpn_cmp(o, C->prod, pn) < 0);

                if (sign && mpn_cmp(o, C->prod_half, pn) > 0)
                {
                    mpn_sub_n(o, C->prod, o, pn);
                    neg = 1;
                }

                if (negative != NULL)
                    negative[i0 + i] = neg;
            }
        }
    }

    TMP_END;
}
