/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb.h"
#include "impl.h"

#define MAX_ENTRY (UWORD(1) << 24)

/* Indices (PARI/GP zetamultconvert): the word 0 y 1 of weight k >= 2
   corresponds to m = (1 y)_2, so the words of weight k are the
   m in [2^(k-2), 2^(k-1)). */

ulong
acb_mzeta_mzv_index(const ulong * s, slong r)
{
    unsigned char * w;
    ulong m = 1, k = 0;
    slong i;
    for (i = 0; i < r; i++)
        k += s[i];
    if (r == 0 || s[0] < 2 || k > FLINT_BITS)
        return 0;
    w = flint_malloc(k);
    _acb_mzeta_mzv_word(w, s, r);
    for (i = 1; i + 1 < (slong) k; i++)
        m = 2 * m + w[i];
    flint_free(w);
    return m;
}

slong
acb_mzeta_mzv_index_composition(ulong * s, ulong m)
{
    slong k, i, r = 0;
    ulong z = 1;      /* the leading 0 */

    if (m == 0)
        return -1;
    k = FLINT_BIT_COUNT(m) + 1;
    /* letters w[1..k-2] are the bits of m below the leading one */
    for (i = k - 3; i >= 0; i--)
    {
        if ((m >> i) & 1)
        {
            s[r++] = z + 1;
            z = 0;
        }
        else
            z++;
    }
    s[r++] = z + 1;
    return r;
}

/* -log2 of a lower bound of zeta(s): the first term prod_i (r - i)^-s_i */
static double
mzv_log2_lower(const ulong * s, slong r)
{
    double lg = 0.0;
    slong i;
    for (i = 0; i < r; i++)
        lg += s[i] * log2((double) (r - i));
    return lg;
}

/* values of nw words given by their indices m (all admissible), sharing
   the states, in batches bounded by memory; prec is the absolute
   precision, extra bits are added for the smallest values */
static void
mzv_by_index(arb_ptr res, const ulong * ms, slong nw, slong prec)
{
    slong i, j, kmax = 0, G, guard = 0;
    ulong s[FLINT_BITS + 1];
    unsigned char ** words;
    slong * lens;

    if (nw == 0)
        return;

    for (i = 0; i < nw; i++)
    {
        slong r = acb_mzeta_mzv_index_composition(s, ms[i]);
        kmax = FLINT_MAX(kmax, FLINT_BIT_COUNT(ms[i]) + 1);
        guard = FLINT_MAX(guard, (slong) mzv_log2_lower(s, r) + 2);
    }
    prec += guard;

    /* batch size: states (48 bytes) and two levels of values per state */
    {
        double per = (double) (kmax + 1) * (kmax + 2) / 2 * (48 + 16.0 * (prec / FLINT_BITS + 3));
        G = (slong) FLINT_MAX(1.0, FLINT_MIN((double) nw, 3e8 / per));
    }

    words = flint_malloc(sizeof(unsigned char *) * G);
    lens = flint_malloc(sizeof(slong) * G);
    for (i = 0; i < G; i++)
        words[i] = flint_malloc(kmax + 1);

    for (j = 0; j < nw; j += G)
    {
        slong g = FLINT_MIN(G, nw - j);
        for (i = 0; i < g; i++)
        {
            slong r = acb_mzeta_mzv_index_composition(s, ms[j + i]);
            lens[i] = _acb_mzeta_mzv_word(words[i], s, r);
        }
        _acb_mzeta_mzv_akhilesh_vec(res + j, (const unsigned char * const *) words, lens, g, prec);
    }

    for (i = 0; i < G; i++)
        flint_free(words[i]);
    flint_free(words);
    flint_free(lens);
}

static int
check_index(const ulong * s, slong r, ulong * k)
{
    slong i;
    *k = 0;
    for (i = 0; i < r; i++)
    {
        if (s[i] == 0 || s[i] > MAX_ENTRY)
            return 0;
        *k += s[i];
    }
    return r == 0 || s[0] >= 2;
}

void
acb_mzeta_mzv(arb_t res, const ulong * s, slong r, slong prec)
{
    unsigned char * w;
    ulong k;

    if (!check_index(s, r, &k))
    {
        arb_indeterminate(res);
        return;
    }
    if (r == 0)
    {
        arb_one(res);
        return;
    }
    if (r == 1)
    {
        arb_zeta_ui(res, s[0], prec);
        return;
    }

    w = flint_malloc(k);
    _acb_mzeta_mzv_word(w, s, r);
    /* quasilinear algorithm from about 25000 bits (measured crossover for
       weights 3 to 12 and various depths, one thread) */
    if (prec >= 25000 && k <= 500)
        _acb_mzeta_mzv_bsplit(res, w, k, prec + (slong) mzv_log2_lower(s, r) + 2);
    else
        _acb_mzeta_mzv_akhilesh(res, w, k, prec + (slong) mzv_log2_lower(s, r) + 2);
    arb_set_round(res, res, prec);
    flint_free(w);
}

/* c_j = sum of zeta(p) over the indices p obtained from s by replacing
   j of the r - 1 commas by plus signs, j = 0, ..., r - 1 */
static void
mzv_merge_coeffs(arb_ptr c, const ulong * s, slong r, slong prec)
{
    slong nm = WORD(1) << (r - 1), mask, i, j, wp, guard = 0;
    unsigned char ** words;
    slong * lens;
    arb_ptr v;
    ulong * p;
    ulong k = 0;

    for (i = 0; i < r; i++)
        k += s[i];

    words = flint_malloc(sizeof(unsigned char *) * nm);
    lens = flint_malloc(sizeof(slong) * nm);
    p = flint_malloc(sizeof(ulong) * r);
    v = _arb_vec_init(nm);

    for (mask = 0; mask < nm; mask++)
    {
        slong q = 0;
        p[0] = s[0];
        for (i = 1; i < r; i++)
        {
            if ((mask >> (i - 1)) & 1)
                p[q] += s[i];
            else
                p[++q] = s[i];
        }
        q++;
        guard = FLINT_MAX(guard, (slong) mzv_log2_lower(p, q) + 2);
        words[mask] = flint_malloc(k);
        lens[mask] = _acb_mzeta_mzv_word(words[mask], p, q);
    }

    wp = prec + guard + FLINT_BIT_COUNT(nm);
    _acb_mzeta_mzv_akhilesh_vec(v, (const unsigned char * const *) words, lens, nm, wp);

    _arb_vec_zero(c, r);
    for (mask = 0; mask < nm; mask++)
    {
        slong mm = mask;
        for (j = 0; mm != 0; mm &= mm - 1)
            j++;
        arb_add(c + j, c + j, v + mask, wp);
    }

    for (mask = 0; mask < nm; mask++)
        flint_free(words[mask]);
    flint_free(words);
    flint_free(lens);
    flint_free(p);
    _arb_vec_clear(v, nm);
}

void
acb_mzeta_mzv_star(arb_t res, const ulong * s, slong r, slong prec)
{
    arb_ptr c;
    ulong k;
    slong j;

    if (!check_index(s, r, &k) || r > 24)
    {
        arb_indeterminate(res);
        return;
    }
    if (r <= 1)
    {
        acb_mzeta_mzv(res, s, r, prec);
        return;
    }

    c = _arb_vec_init(r);
    mzv_merge_coeffs(c, s, r, prec);
    arb_zero(res);
    for (j = 0; j < r; j++)
        arb_add(res, res, c + j, prec + 16);
    arb_set_round(res, res, prec);
    _arb_vec_clear(c, r);
}

void
acb_mzeta_mzv_t(acb_t res, const ulong * s, slong r, const acb_t t, slong prec)
{
    arb_ptr c;
    ulong k;
    slong j;

    if (!check_index(s, r, &k) || r > 24)
    {
        acb_indeterminate(res);
        return;
    }
    if (r <= 1)
    {
        acb_mzeta_mzv(acb_realref(res), s, r, prec);
        arb_zero(acb_imagref(res));
        return;
    }

    c = _arb_vec_init(r);
    mzv_merge_coeffs(c, s, r, prec);
    /* Horner in t */
    acb_set_arb(res, c + r - 1);
    for (j = r - 2; j >= 0; j--)
    {
        acb_mul(res, res, t, prec + 16);
        arb_add(acb_realref(res), acb_realref(res), c + j, prec + 16);
    }
    acb_set_round(res, res, prec);
    _arb_vec_clear(c, r);
}

void
acb_mzeta_mzv_all(arb_ptr res, ulong k, int flags, slong prec)
{
    ulong m0, m1, m, nw, i;
    ulong * ms;

    if (k < 2 || k > FLINT_BITS - 2)
        return;

    m0 = (flags & ACB_MZETA_MZV_WEIGHT) ? (UWORD(1) << (k - 2)) : 1;
    m1 = UWORD(1) << (k - 1);
    nw = m1 - m0;

    ms = flint_malloc(sizeof(ulong) * nw);
    for (i = 0; i < nw; i++)
        ms[i] = m0 + i;

    /* the depth-1 values zeta(k) by the same algorithm (they are a small
       fraction), everything at the same absolute precision */
    mzv_by_index(res, ms, nw, prec);

    if (flags & ACB_MZETA_MZV_STAR)
    {
        /* zeta*(s) = sum over the merges p of s of zeta(p); the merges
           have the same weight, so their values are in res */
        arb_ptr st = _arb_vec_init(nw);
        ulong s[FLINT_BITS + 1], p[FLINT_BITS + 1];

        for (i = 0; i < nw; i++)
        {
            slong r = acb_mzeta_mzv_index_composition(s, ms[i]), j, q;
            ulong mask, nm = UWORD(1) << (r - 1);
            arb_zero(st + i);
            for (mask = 0; mask < nm; mask++)
            {
                ulong mp;
                p[0] = s[0];
                q = 0;
                for (j = 1; j < r; j++)
                {
                    if ((mask >> (j - 1)) & 1)
                        p[q] += s[j];
                    else
                        p[++q] = s[j];
                }
                q++;
                mp = acb_mzeta_mzv_index(p, q);
                arb_add(st + i, st + i, res + (mp - m0), prec + 16);
            }
        }
        for (i = 0; i < nw; i++)
            arb_set_round(res + i, st + i, prec);
        _arb_vec_clear(st, nw);
    }
    else
    {
        for (i = 0; i < nw; i++)
            arb_set_round(res + i, res + i, prec);
    }

    flint_free(ms);
    (void) m;
}
