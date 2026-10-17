/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

/* Transformed power tables for the binary splittings with a fixed
   midpoint (exp_sum_bs.c, sin_cos_sum_bs.c).

   A merge of [a, b) with step = (b - a) / 2 multiplies the table entry
   x^step by the right child, and every merge of a level uses the same
   entry, so its transform can be computed once per tree and each of
   those products saves the pack and forward transform of the entry --
   about 30% of a two-prime product.  The ranges of a level differ in
   length by at most one, so the number of merges using each entry and
   the longest right child (which bounds the partners) follow from N
   alone, and the transforms are made before the recursion, read only
   during it (the two halves may run on two threads).

   An entry is transformed only if at least two merges use it: the
   root's entry, the largest, is used once and would only add the
   storage of its transform.  The storage of all entries is one block,
   about four times the transform of the largest entry kept.

   Measured on Emerald Rapids, single-threaded: exp 8% faster at 10^6
   digits and 5% at 3 10^5, sin and cos 4% and 3% (about 1% at 10^5,
   where the bit-burst trees are a small part of the work); on two
   threads exp 8% and sin and cos 2% at 10^6.  Above the two-prime
   range the transforms are wide (16-22% per product); they add about
   3% at 10^7 digits and nothing at 10^6, where every longer table
   product is a root's (used once).  At 10^7 digits the transforms
   raise the peak memory of exp by about 17% (sin and cos 8%); at
   higher precisions MP_REAL_POWTAB_PRE_MAX_LIMBS bounds the extra
   memory, leaving the levels above it untransformed.  The
   threshold MP_REAL_POWTAB_PRE_CUTOFF made no measurable difference
   between 60 and 600 limbs. */

void
_mp_real_bs_step_uses(slong * uses, slong * rmax, const slong * xexp,
    slong len, slong N)
{
    slong n[4], c[4], k, nk, i, j;

    for (i = 0; i < len; i++)
    {
        uses[i] = 0;
        rmax[i] = 0;
    }

    /* the distinct range lengths of a level (at most two, differing
       by one) with their multiplicities; ranges of 1 or 2 terms are
       leaves */
    n[0] = N;
    c[0] = 1;
    k = 1;

    while (k > 0)
    {
        slong n2[4], c2[4];

        nk = 0;
        for (j = 0; j < k; j++)
        {
            slong s, r, t, ch[2];

            if (n[j] <= 2)
                continue;

            s = n[j] / 2;
            r = n[j] - s;

            for (i = 0; i < len && xexp[i] != s; i++)
                ;
            FLINT_ASSERT(i < len);
            if (i < len)
            {
                uses[i] += c[j];
                rmax[i] = FLINT_MAX(rmax[i], r);
            }

            ch[0] = s;
            ch[1] = r;
            for (t = 0; t < 2; t++)
            {
                slong u;
                for (u = 0; u < nk && n2[u] != ch[t]; u++)
                    ;
                if (u == nk)
                {
                    FLINT_ASSERT(nk < 4);
                    n2[nk] = ch[t];
                    c2[nk] = 0;
                    nk++;
                }
                c2[u] += c[j];
            }
        }

        for (j = 0; j < nk; j++)
        {
            n[j] = n2[j];
            c[j] = c2[j];
        }
        k = nk;
    }
}

#if FLINT_HAVE_FFT_SMALL

/* a product of l by s limbs is worth the transformed entry when it is
   in the FFT range of flint_mpn_mul with the threshold
   MP_REAL_POWTAB_PRE_CUTOFF */
static int
_pre_range(slong l, slong s)
{
    slong c = MP_REAL_POWTAB_PRE_CUTOFF;
    if (l < s)
        FLINT_SWAP(slong, l, s);
    return s >= c || (l + s > 2 * c && 3 * s > c);
}

void
_mp_real_powtab_pre_init(mp_real_powtab_pre_struct * T, nn_srcptr * xpow,
    const slong * xlen, const slong * uses, const slong * ymax, slong len)
{
    fft_small_mpn_mul_precomp_struct * pre;
    ulong * size;
    ulong total = 0;
    slong i;

    T->pre = NULL;
    T->storage = NULL;
    T->len = len;

    size = flint_malloc(len * sizeof(ulong));
    for (i = 0; i < len; i++)
    {
        size[i] = 0;
        if (uses[i] >= 2 && xlen[i] > 0 && ymax[i] > 0
            && _pre_range(xlen[i], ymax[i])
            && xlen[i] + ymax[i] <= MP_REAL_POWTAB_PRE_MAX_LIMBS)
            size[i] = _fft_small_mpn_mul_precomp_size(xlen[i], ymax[i]);
        total += size[i];
    }

    if (total != 0)
    {
        double * s;

        pre = flint_malloc(len * sizeof(fft_small_mpn_mul_precomp_struct));
        T->storage = flint_aligned_alloc(64, total * sizeof(double));
        s = T->storage;

        for (i = 0; i < len; i++)
        {
            pre[i].d = NULL;
            if (size[i] != 0)
            {
                _fft_small_mpn_mul_precomp_init(pre + i, xpow[i], xlen[i],
                    ymax[i], s);
                s += size[i];
            }
        }

        T->pre = pre;
    }

    flint_free(size);
}

void
_mp_real_powtab_pre_clear(mp_real_powtab_pre_struct * T)
{
    if (T->pre != NULL)
    {
        fft_small_mpn_mul_precomp_struct * pre = T->pre;
        slong i;
        for (i = 0; i < T->len; i++)
            if (pre[i].d != NULL)
                _fft_small_mpn_mul_precomp_clear(pre + i);
        flint_free(pre);
        flint_aligned_free(T->storage);
    }
    T->pre = NULL;
    T->storage = NULL;
}

int
_mp_real_powtab_pre_mul(nn_ptr z, const mp_real_powtab_pre_struct * T,
    slong i, nn_srcptr y, slong yn)
{
    const fft_small_mpn_mul_precomp_struct * pre;

    if (T == NULL || T->pre == NULL || yn <= 0)
        return 0;

    pre = (const fft_small_mpn_mul_precomp_struct *) T->pre + i;
    if (pre->d == NULL || (ulong) yn > pre->yn_max || !_pre_range(pre->xn, yn))
        return 0;

    _fft_small_mpn_mul_precomp(z, pre, y, yn);
    return 1;
}

#else

void
_mp_real_powtab_pre_init(mp_real_powtab_pre_struct * T, nn_srcptr * xpow,
    const slong * xlen, const slong * uses, const slong * ymax, slong len)
{
    T->pre = NULL;
    T->storage = NULL;
    T->len = len;
}

void
_mp_real_powtab_pre_clear(mp_real_powtab_pre_struct * T)
{
}

int
_mp_real_powtab_pre_mul(nn_ptr z, const mp_real_powtab_pre_struct * T,
    slong i, nn_srcptr y, slong yn)
{
    return 0;
}

#endif
