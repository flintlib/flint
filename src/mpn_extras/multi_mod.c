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

/* r = a mod d where d is node i of level k; r has level_limbs[k] limbs */
static slong
_reduce_node(nn_ptr r, nn_srcptr a, slong an, const flint_mpn_crt_t C,
             slong k, slong i, nn_ptr scratch)
{
    slong ll = C->level_limbs[k];
    slong dn = C->level_len[k][i];
    nn_srcptr d = C->level_prod[k] + i * ll;

    FLINT_ASSERT(an >= dn);
    FLINT_ASSERT(an <= C->level_limbs[k + 1] + FLINT_MPN_CRT_SLOT_EXTRA);

    if (C->level_use_preinv[k])
    {
        ulong norm = C->level_norm[k][i];
        nn_ptr t = scratch;
        nn_ptr r2 = scratch + an + 1;
        slong tn;

        if (norm != 0)
        {
            t[an] = mpn_lshift(t, a, an, norm);
            tn = an + (t[an] != 0);
        }
        else
        {
            flint_mpn_copyi(t, a, an);
            tn = an;
        }

        flint_mpn_mod_preinvn(r2, t, tn, C->level_prod_norm[k] + i * ll, dn,
                              C->level_inv[k] + i * ll);

        if (norm != 0)
            mpn_rshift(r, r2, dn, norm);
        else
            flint_mpn_copyi(r, r2, dn);
    }
    else
    {
        mpn_tdiv_qr(scratch, r, 0, a, an, d, dn);
    }

    MPN_NORM(r, dn);
    FLINT_ASSERT(dn < C->level_len[k][i] || mpn_cmp(r, d, dn) < 0);
    return dn;
}

/* dot product of (a, n) with pw, reduced mod l; n small and fixed */
#define DEFINE_DOT2(N) \
FLINT_FORCE_INLINE ulong CAT(_dot2, N)(nn_srcptr a, nn_srcptr pw, const flint_mpn_crt_leaf_struct * L) \
{ \
    ulong h0 = 0, l0 = 0, h1 = 0, l1 = 0; \
    for (ulong k = 0; k < N; k += 2) \
    { \
        MAC(h0, l0, a[k], pw[k]); \
        if (k + 1 < N) MAC(h1, l1, a[k + 1], pw[k + 1]); \
    } \
    add_ssaaaa(h0, l0, h0, l0, h1, l1); \
    return _nmod_red2_any(h0, l0, L); \
}
DEFINE_DOT2(1) DEFINE_DOT2(2) DEFINE_DOT2(3) DEFINE_DOT2(4)
DEFINE_DOT2(5) DEFINE_DOT2(6) DEFINE_DOT2(7) DEFINE_DOT2(8)
#undef DEFINE_DOT2

#define DEFINE_DOT3(N) \
FLINT_FORCE_INLINE ulong CAT(_dot3, N)(nn_srcptr a, nn_srcptr pw, const flint_mpn_crt_leaf_struct * L) \
{ \
    ulong h = 0, m = 0, l = 0; \
    for (ulong k = 0; k < N; k++) \
    { \
        ulong p1, p0; \
        umul_ppmm(p1, p0, a[k], pw[k]); \
        add_sssaaaaaa(h, m, l, h, m, l, 0, p1, p0); \
    } \
    return _nmod_red3_any(h, m, l, L); \
}
DEFINE_DOT3(1) DEFINE_DOT3(2) DEFINE_DOT3(3) DEFINE_DOT3(4)
DEFINE_DOT3(5) DEFINE_DOT3(6) DEFINE_DOT3(7) DEFINE_DOT3(8)
#undef DEFINE_DOT3

FLINT_FORCE_INLINE ulong
_dot_mod(nn_srcptr a, slong an, nn_srcptr pw, const flint_mpn_crt_leaf_struct * L, int slack)
{
    slong k;

    if (an == 1)
        return n_mod_barrett(a[0], L->mod.n, L->npre);

    if (an == 2 && L->mod.norm == 0)
    {
        /* fullword modulus: two chained reductions beat the 3-limb accumulator */
        ulong u = n_mod_barrett(a[1], L->mod.n, L->npre);
        NMOD_RED2_FULLWORD(u, u, a[0], L->mod);
        return u;
    }

    if (slack)
    {
        ulong h0 = 0, l0 = 0, h1 = 0, l1 = 0;

        switch (an)
        {
            case 2: return _dot2_2(a, pw, L);
            case 3: return _dot2_3(a, pw, L);
            case 4: return _dot2_4(a, pw, L);
            case 5: return _dot2_5(a, pw, L);
            case 6: return _dot2_6(a, pw, L);
            case 7: return _dot2_7(a, pw, L);
            case 8: return _dot2_8(a, pw, L);
            default: break;
        }

        for (k = 0; k + 2 <= an; k += 2)
        {
            MAC(h0, l0, a[k], pw[k]);
            MAC(h1, l1, a[k + 1], pw[k + 1]);
        }
        if (k < an)
            MAC(h0, l0, a[k], pw[k]);

        add_ssaaaa(h0, l0, h0, l0, h1, l1);
        return _nmod_red2_any(h0, l0, L);
    }
    else
    {
        ulong h = 0, m = 0, l = 0;

        switch (an)
        {
            case 2: return _dot3_2(a, pw, L);
            case 3: return _dot3_3(a, pw, L);
            case 4: return _dot3_4(a, pw, L);
            case 5: return _dot3_5(a, pw, L);
            case 6: return _dot3_6(a, pw, L);
            case 7: return _dot3_7(a, pw, L);
            case 8: return _dot3_8(a, pw, L);
            default: break;
        }

        for (k = 0; k < an; k++)
        {
            ulong p1, p0;
            umul_ppmm(p1, p0, a[k], pw[k]);
            add_sssaaaaaa(h, m, l, h, m, l, 0, p1, p0);
        }

        return _nmod_red3_any(h, m, l, L);
    }
}

ulong
flint_mpn_crt_mod_leaf(nn_srcptr a, slong an, const flint_mpn_crt_t C, slong j)
{
    FLINT_ASSERT(an <= C->mod_pow_limbs);
    MPN_NORM(a, an);
    return _dot_mod(a, an, C->mod_pow + j * C->mod_pow_limbs, C->leaf + j, C->mod_pow_slack);
}

/* residues of (a, an) modulo all primes below node i at the mod base level */
static void
_mod_basecase(nn_ptr out, nn_srcptr a, slong an, const flint_mpn_crt_t C, slong i)
{
    slong L = C->mod_base_level;
    slong c0 = i << L;
    slong c1 = FLINT_MIN((i + 1) << L, C->num_chunks);
    slong j0 = C->chunk_start[c0];
    slong j1 = C->chunk_start[c1];
    slong mbl = C->mod_pow_limbs;
    slong j, k;
    int slack = C->mod_pow_slack;

    FLINT_ASSERT(an <= mbl);

    if (an == 0)
    {
        for (k = C->leaf_start[j0]; k < C->leaf_start[j1]; k++)
            out[k] = 0;
        return;
    }

    if (C->all_single)
    {
        for (j = j0; j < j1; j++)
        {
            out[j] = _dot_mod(a, an, C->mod_pow + j * mbl, C->leaf + j, slack);
            FLINT_ASSERT(out[j] < C->primes[j]);
        }
    }
    else
    {
        for (j = j0; j < j1; j++)
        {
            slong k0 = C->leaf_start[j];
            ulong u = _dot_mod(a, an, C->mod_pow + j * mbl, C->leaf + j, slack);

            FLINT_ASSERT(u < C->leaf[j].mod.n);

            if (C->leaf_start[j + 1] == k0 + 1)
                out[k0] = u;
            else
                for (k = k0; k < C->leaf_start[j + 1]; k++)
                    out[k] = n_mod_barrett(u, C->primes[k], C->prime_barrett[k]);

            for (k = k0; k < C->leaf_start[j + 1]; k++)
                FLINT_ASSERT(out[k] < C->primes[k]);
        }
    }
}


typedef struct
{
    const flint_mpn_crt_struct * C;
    slong k;
    nn_srcptr x;
    nn_ptr * ptr;
    slong * len;
    nn_ptr * ptr2;
    slong * len2;
    nn_ptr buf;
    nn_ptr scratch;
    int parallel;
}
_mod_descend_arg;

/* reduce the value of node i at level k to its children */
static void
_mod_descend_worker(slong i, void * varg)
{
    _mod_descend_arg * arg = (_mod_descend_arg *) varg;
    const flint_mpn_crt_struct * C = arg->C;
    slong k = arg->k;
    slong cl = C->level_limbs[k - 1] + FLINT_MPN_CRT_SLOT_EXTRA;
    slong ccount = C->level_count[k - 1];
    nn_srcptr a = arg->ptr[i];
    slong an = arg->len[i];
    nn_ptr buf = arg->buf;
    nn_ptr scratch;
    slong c;

    scratch = arg->parallel ? FLINT_ARRAY_ALLOC(4 * (C->prod_len + FLINT_MPN_CRT_SLOT_EXTRA) + 16, ulong) : arg->scratch;

    for (c = 2 * i; c < FLINT_MIN(2 * i + 2, ccount); c++)
    {
        slong dn = C->level_len[k - 1][c];
        nn_srcptr d = C->level_prod[k - 1] + c * C->level_limbs[k - 1];

        if (an < dn || (an == dn && mpn_cmp(a, d, dn) < 0))
        {
            /* already reduced: pass through the original input
               by pointer, otherwise copy (the ping-pong buffer
               will be overwritten two levels down) */
            if (a == arg->x)
            {
                arg->ptr2[c] = (nn_ptr) a;
            }
            else
            {
                arg->ptr2[c] = buf + c * cl;
                flint_mpn_copyi(buf + c * cl, a, an);
            }
            arg->len2[c] = an;
        }
        else
        {
            arg->ptr2[c] = buf + c * cl;
            arg->len2[c] = _reduce_node(buf + c * cl, a, an, C, k - 1, c, scratch);
        }
    }

    if (arg->parallel)
        flint_free(scratch);
}

typedef struct
{
    const flint_mpn_crt_struct * C;
    nn_ptr out;
    nn_ptr * ptr;
    slong * len;
}
_mod_base_arg;

static void
_mod_base_worker(slong i, void * varg)
{
    _mod_base_arg * arg = (_mod_base_arg *) varg;
    _mod_basecase(arg->out, arg->ptr[i], arg->len[i], arg->C, i);
}


void
flint_mpn_multi_mod(nn_ptr out, nn_srcptr x, slong xn, const flint_mpn_crt_t C, nn_ptr tmp)
{
    _crt_work_struct W;
    slong top = C->num_levels - 1;
    slong k, i, pn;
    nn_srcptr p;
    nn_ptr * ptr, * ptr2;
    slong * len, * len2;
    nn_ptr buf, buf2;
    nn_ptr tmp_alloc = NULL;
    TMP_INIT;

    FLINT_ASSERT(C->flags & FLINT_MPN_CRT_MOD);

    MPN_NORM(x, xn);

    /* fast path: no tree levels above the basecase and the input is
       small enough for the power tables */
    if (C->mod_base_level == top && xn <= C->mod_pow_limbs)
    {
        /* tiny cases: avoid the bookkeeping of the general basecase */
        if (C->num_leaves <= 4 && C->all_single)
        {
            slong j, mpl = C->mod_pow_limbs;
            for (j = 0; j < C->num_leaves; j++)
                out[j] = _dot_mod(x, xn, C->mod_pow + j * mpl, C->leaf + j, C->mod_pow_slack);
            return;
        }

        _mod_basecase(out, x, xn, C, 0);
        return;
    }

    TMP_START;
    if (tmp == NULL)
    {
        tmp = TMP_ALLOC(C->tmp_limbs * sizeof(ulong));
    }

    _crt_work_init(&W, C, tmp);

    /* reduce modulo the full product if necessary */
    if (xn > C->prod_len || (xn == C->prod_len && mpn_cmp(x, C->prod, xn) >= 0))
    {
        slong qn = xn - C->prod_len + 1;
        nn_ptr q;

        if (qn <= W.scratch_limbs)
            q = W.scratch;
        else
            q = tmp_alloc = flint_malloc(qn * sizeof(ulong));

        mpn_tdiv_qr(q, W.bufA, 0, x, xn, C->prod, C->prod_len);
        p = W.bufA;
        pn = C->prod_len;
        MPN_NORM(p, pn);
        buf = W.bufB;
        buf2 = W.bufA;
    }
    else
    {
        p = x;
        pn = xn;
        buf = W.bufA;
        buf2 = W.bufB;
    }

    ptr = W.ptrA;
    len = W.lenA;
    ptr2 = W.ptrB;
    len2 = W.lenB;
    ptr[0] = (nn_ptr) p;
    len[0] = pn;

    /* descend */
    for (k = top; k > C->mod_base_level; k--)
    {
        _mod_descend_arg arg;

        arg.C = C;
        arg.k = k;
        arg.x = x;
        arg.ptr = ptr;
        arg.len = len;
        arg.ptr2 = ptr2;
        arg.len2 = len2;
        arg.buf = buf;
        arg.scratch = W.scratch;
        arg.parallel = _crt_use_threads(C->level_count[k - 1], C->level_limbs[k - 1]);

        if (arg.parallel)
            flint_parallel_do(_mod_descend_worker, &arg, C->level_count[k], -1, FLINT_PARALLEL_UNIFORM);
        else
            for (i = 0; i < C->level_count[k]; i++)
                _mod_descend_worker(i, &arg);

        FLINT_SWAP(nn_ptr *, ptr, ptr2);
        FLINT_SWAP(slong *, len, len2);
        FLINT_SWAP(nn_ptr, buf, buf2);
    }

    /* basecase */
    {
        _mod_base_arg arg;

        arg.C = C;
        arg.out = out;
        arg.ptr = ptr;
        arg.len = len;

        if (_crt_use_threads(C->level_count[C->mod_base_level], C->mod_base_limbs))
            flint_parallel_do(_mod_base_worker, &arg, C->level_count[C->mod_base_level], -1, FLINT_PARALLEL_UNIFORM);
        else
            for (i = 0; i < C->level_count[C->mod_base_level]; i++)
                _mod_base_worker(i, &arg);
    }

    flint_free(tmp_alloc);
    TMP_END;
}

/* block of entries whose limbs fit comfortably in L1 */
#define VEC_BLOCK_LIMBS 1024

void
flint_mpn_multi_mod_vec(nn_ptr out, slong out_stride, nn_srcptr x, slong xn, slong len,
                        const flint_mpn_crt_t C, nn_ptr tmp)
{
    slong top = C->num_levels - 1;
    slong i, j, l, i0, blk;

    FLINT_ASSERT(C->flags & FLINT_MPN_CRT_MOD);

    if (!(C->mod_base_level == top && xn <= C->mod_pow_limbs && C->all_single))
    {
        /* general case: entry by entry */
        nn_ptr res;
        TMP_INIT;
        TMP_START;
        res = TMP_ALLOC(C->num_primes * sizeof(ulong));
        for (i = 0; i < len; i++)
        {
            flint_mpn_multi_mod(res, x + i * xn, xn, C, tmp);
            for (l = 0; l < C->num_primes; l++)
                out[l * out_stride + i] = res[l];
        }
        TMP_END;
        return;
    }

    /* basecase for all entries: leaf-outer loop within blocks of entries */
    blk = FLINT_MAX(1, VEC_BLOCK_LIMBS / FLINT_MAX(xn, 1));

    for (i0 = 0; i0 < len; i0 += blk)
    {
        slong i1 = FLINT_MIN(i0 + blk, len);

        for (j = 0; j < C->num_leaves; j++)
        {
            nn_srcptr pw = C->mod_pow + j * C->mod_pow_limbs;
            const flint_mpn_crt_leaf_struct * L = C->leaf + j;

            for (i = i0; i < i1; i++)
            {
                nn_srcptr a = x + i * xn;
                slong an = xn;
                MPN_NORM(a, an);
                out[j * out_stride + i] = _dot_mod(a, an, pw, L, C->mod_pow_slack);
            }
        }
    }
}
