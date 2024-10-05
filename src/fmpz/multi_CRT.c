/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "nmod.h"
#include "fmpz.h"

void fmpz_multi_CRT_init(fmpz_multi_CRT_t P)
{
    P->prog = NULL;
    P->moduli = NULL;
    P->fracmoduli = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
    P->temp2loc = 0;
    P->temp3loc = 0;
    P->temp4loc = 0;
    P->good = 0;
    fmpz_init(P->final_modulus);
}

void fmpz_multi_CRT_clear(fmpz_multi_CRT_t P)
{
    slong i;

    for (i = 0; i < P->alloc; i++)
    {
        fmpz_clear(P->prog[i].b_modulus);
        fmpz_clear(P->prog[i].c_modulus);
        fmpz_clear(P->moduli + i);
        fmpz_clear(P->fracmoduli + i);
    }

    flint_free(P->prog);
    flint_free(P->moduli);
    flint_free(P->fracmoduli);
    fmpz_clear(P->final_modulus);
}

int fmpz_multi_CRT(
    fmpz_t output,
    const fmpz * moduli,
    const fmpz * values,
    slong len,
    int sign)
{
    int success;
    slong i;
    fmpz_multi_CRT_t P;
    fmpz * out;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;

    fmpz_multi_CRT_init(P);
    success = fmpz_multi_CRT_precompute(P, moduli, len);

    out = TMP_ARRAY_ALLOC(P->localsize, fmpz);
    for (i = 0; i < P->localsize; i++)
        fmpz_init(out + i);

    fmpz_swap(out + 0, output);
    _fmpz_multi_CRT_precomp(out, P, values, sign);
    fmpz_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
        fmpz_clear(out + i);

    fmpz_multi_CRT_clear(P);

    TMP_END;

    return success;
}

void _fmpz_multi_CRT_precomp(
    fmpz * outputs,
    const fmpz_multi_CRT_t P,
    const fmpz * inputs,
    int sign)
{
    slong i, a, b, c;
    slong len = P->length;
    const fmpz * m = P->moduli;
    const fmpz * mf = P->fracmoduli;
    fmpz * A, * B, * C, * t1, * t2, * t3, * t4;

    t1 = outputs + P->temp1loc;
    t2 = outputs + P->temp2loc;
    t3 = outputs + P->temp3loc;
    t4 = outputs + P->temp4loc;

    FLINT_ASSERT(len < 1 || P->good);

    if (len > 0)
    {
        for (i = P->moduli_count - 1; i > 0; i--)
        {
            if (!fmpz_equal(inputs + 0, inputs + i))
                goto doit;
        }
    }

    _fmpz_smod(outputs + 0, inputs + 0, P->final_modulus, sign, t4);
    return;

doit:

    for (i = 0; i < len; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;
        c = P->prog[i].c_idx;

        A = outputs + a;
        B = outputs + b;
        C = outputs + c;

        if (b < 0)
        {
            b = -b - 1;
            B = t1;

            fmpz_mul(t3, inputs + b, mf + b);
            _fmpz_smod(B, t3, m + b, sign, t4);
        }

        if (c < 0)
        {
            c = -c - 1;
            C = t2;

            fmpz_mul(t3, inputs + c, mf + c);
            _fmpz_smod(C, t3, m + c, sign, t4);
        }

        /* A = B*c_m + C*b_m */
        fmpz_mul(A, B, P->prog[i].c_modulus);
        fmpz_mul(t3, C, P->prog[i].b_modulus);
        fmpz_add(A, A, t3);
    }

    _fmpz_smod(outputs + 0, A, P->final_modulus, sign, t4);
}

void fmpz_multi_CRT_precomp(
    fmpz_t output,
    const fmpz_multi_CRT_t P,
    const fmpz * inputs,
    int sign)
{
    slong i;
    fmpz * out;
    TMP_INIT;

    TMP_START;
    out = TMP_ARRAY_ALLOC(P->localsize, fmpz);
    for (i = 0; i < P->localsize; i++)
        fmpz_init(out + i);

    fmpz_swap(out + 0, output);
    _fmpz_multi_CRT_precomp(out, P, inputs, sign);
    fmpz_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
        fmpz_clear(out + i);

    TMP_END;
}

static void _fmpz_multi_CRT_fit_length(fmpz_multi_CRT_t P, slong k)
{
    slong i;

    k = FLINT_MAX(WORD(1), k);

    for (i = k; i < P->alloc; i++)
    {
        fmpz_clear(P->prog[i].b_modulus);
        fmpz_clear(P->prog[i].c_modulus);
        fmpz_clear(P->moduli + i);
        fmpz_clear(P->fracmoduli + i);
    }

    P->prog = FLINT_ARRAY_REALLOC(P->prog, k, _fmpz_multi_CRT_instr);
    P->moduli = FLINT_ARRAY_REALLOC(P->moduli, k, fmpz);
    P->fracmoduli = FLINT_ARRAY_REALLOC(P->fracmoduli, k, fmpz);

    for (i = P->alloc; i < k; i++)
    {
        fmpz_init(P->prog[i].b_modulus);
        fmpz_init(P->prog[i].c_modulus);
        fmpz_init(P->moduli + i);
        fmpz_init(P->fracmoduli + i);
    }

    P->alloc = k;
}

static int _fill_pfrac(
    slong * link,
    fmpz * v,
    fmpz * w,
    slong j,
    const fmpz_t A,
    fmpz_t g, /* temps */
    fmpz_t s,
    fmpz_t t)
{
    while (j >= 0)
    {
        int cmp = fmpz_cmp(v + j, v + j + 1);

        /* A/(v[j]*v[j+1]) = w[j]/v[j] + w[j+1]/v[j+1] mod 1 */

        if (fmpz_is_zero(v + j) || fmpz_is_zero(v + j + 1) ||
            fmpz_is_one(v + j) || fmpz_is_one(v + j + 1) ||
            cmp == 0)
        {
            return 0;
        }

        /* fmpz_gcdinv requires x < y AND we hit the smaller branch first */

        if (cmp > 0)
        {
            fmpz_swap(v + j, v + j + 1);
            FLINT_SWAP(slong, link[j], link[j + 1]);
        }

        fmpz_gcdinv(g, s, v + j + 0, v + j + 1);
        if (!fmpz_is_one(g))
            return 0;

        fmpz_mul(w + j + 1, A, s);
        fmpz_mod(w + j + 1, w + j + 1, v + j + 1 );

        /* w[j] = (A - v[j] w[1 + j])/v[1 + j] */

        fmpz_mul(w + j + 0, v + j + 0, w + j + 1);
        fmpz_sub(t, A, w + j + 0);
        fmpz_fdiv_qr(w + j + 0, g, t, v + j + 1);
        FLINT_ASSERT(fmpz_is_zero(g) && "division should be exact");

        fmpz_mod(w + j + 0, w + j + 0, v + j + 0);

        if (!_fill_pfrac(link, v, w, link[j + 0], w + j + 0, g, s, t))
            return 0;

        A = w + j + 1;
        j = link[j + 1];
    }

    return 1;
}

static void _fill_prog(
    fmpz_multi_CRT_t P,
    slong * link,
    fmpz * v,
    fmpz * w,
    slong j,
    slong ret_idx)
{
    slong i, b_idx, c_idx;
    slong next_ret_idx = ret_idx;

    FLINT_ASSERT(j >= 0);

    if (link[j] >= 0)
    {
        b_idx = ++next_ret_idx;
        _fill_prog(P, link, v, w, link[j], b_idx);
    }
    else
    {
        b_idx = -1 - link[j];
        FLINT_ASSERT(b_idx < P->alloc);
        fmpz_set(P->moduli + b_idx, v + j);
        fmpz_set(P->fracmoduli + b_idx, w + j);
        b_idx = -1 - b_idx;
    }

    if (link[j + 1] >= 0)
    {
        c_idx = ++next_ret_idx;
        _fill_prog(P, link, v, w, link[j + 1], c_idx);
    }
    else
    {
        c_idx = -1 - link[j + 1];
        FLINT_ASSERT(c_idx < P->alloc);
        fmpz_set(P->moduli + c_idx, v + j + 1);
        fmpz_set(P->fracmoduli + c_idx, w + j + 1);
        c_idx = -1 - c_idx;
    }

    i = P->length;
    FLINT_ASSERT(i < P->alloc);
    P->prog[i].a_idx = ret_idx;
    P->prog[i].b_idx = b_idx;
    P->prog[i].c_idx = c_idx;
    fmpz_set(P->prog[i].b_modulus, v + j);
    fmpz_set(P->prog[i].c_modulus, v + j + 1);
    P->length = i + 1;

    P->localsize = FLINT_MAX(P->localsize, 1 + next_ret_idx);
}

int fmpz_multi_CRT_precompute(
    fmpz_multi_CRT_t P,
    const fmpz * f,
    slong r)
{
    slong i, j;
    slong * link;
    fmpz * v;
    fmpz * w;
    fmpz_t one, g, s, t;

    FLINT_ASSERT(r > 0);

    _fmpz_multi_CRT_fit_length(P, r);
    P->length = 0;
    P->localsize = 1;
    P->moduli_count = r;
    P->min_modulus_bits = fmpz_bits(f + 0);

    if (r < 2)
    {
        P->good = !fmpz_is_zero(f + 0);

        if (P->good)
        {
            fmpz_abs(P->final_modulus, f + 0);
            fmpz_abs(P->moduli + 0, f + 0);
            fmpz_one(P->fracmoduli + 0);
        }

        goto done;
    }

    fmpz_init(one);
    fmpz_init(g);
    fmpz_init(s);
    fmpz_init(t);

    link = FLINT_ARRAY_ALLOC(2*r - 2, slong);
    v = FLINT_ARRAY_ALLOC(2*(2*r - 2), fmpz);
    w = v + 2*r - 2;

    for (i = 0; i < 2*(2*r - 2); i++)
        fmpz_init(v + i);

    for (i = 0; i < r; i++)
    {
        flint_bitcnt_t this_bits = fmpz_bits(f + i);
        P->min_modulus_bits = FLINT_MIN(P->min_modulus_bits, this_bits);
        fmpz_abs(v + i, f + i);
        link[i] = -i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        slong sx, minp;
        const fmpz * mind;

        minp = j;
        mind = v + j;
        for (sx = j + 1; sx < i; sx++)
        {
            if (fmpz_cmp(v + sx, mind) < 0)
            {
                mind = v + sx;
                minp = sx;
            }
        }
        fmpz_swap(v + j, v + minp);
        FLINT_SWAP(slong, link[j], link[minp]);

        minp = j + 1;
        mind = v + j + 1;
        for (sx = j + 2; sx < i; sx++)
        {
            if (fmpz_cmp(v + sx, mind) < 0)
            {
                mind = v + sx;
                minp = sx;
            }
        }
        fmpz_swap(v + j + 1, v + minp);
        FLINT_SWAP(slong, link[j + 1], link[minp]);

        fmpz_mul(v + i, v + j, v + j + 1);
        link[i] = j;
    }

    fmpz_mul(P->final_modulus, v + 2*r - 4, v + 2*r - 3);

    fmpz_one(one);
    P->good = _fill_pfrac(link, v, w, 2*r - 4, one, g, s, t);
    if (P->good)
        _fill_prog(P, link, v, w, 2*r - 4, 0);

    fmpz_clear(one);
    fmpz_clear(g);
    fmpz_clear(s);
    fmpz_clear(t);

    for (i = 0; i < 2*(2*r - 2); i++)
        fmpz_clear(v + i);

    flint_free(link);
    flint_free(v);

done:

    P->temp1loc = P->localsize++;
    P->temp2loc = P->localsize++;
    P->temp3loc = P->localsize++;
    P->temp4loc = P->localsize++;

    if (!P->good)
    {
        fmpz_one(P->final_modulus);
        P->length = 0;
    }

    return P->good;
}

#define MAC(h, l, a, b)                 \
do {                                    \
    ulong p1, p0;                   \
    umul_ppmm(p1, p0, a, b);            \
    add_ssaaaa(h, l, h, l, p1, p0);     \
} while (0)

void fmpz_multi_CRT_ui(
    fmpz_t b,
    nn_srcptr in,
    const fmpz_comb_t C,
    fmpz_comb_temp_t CT,
    int sign)
{
    slong i, j, k, l, s;
    slong klen = C->crt_klen;
    slong * step = C->step;
    crt_lut_entry * lu = C->crt_lu;
    fmpz * T = CT->T;
    fmpz * A = CT->A;
    slong * offsets = C->crt_offsets;
    const ulong * md = C->packed_multipliers;
    mpz_ptr az;
    ulong * ad;
    ulong hi, lo, t;

    for (k = 0, i = 0, l = 0; k < klen; k++)
    {
        s = step[k];
        j = offsets[k];
        az = _fmpz_promote(A + k);

        if (s < 0)
        {
            /*
                every low level combination in this chunk has 1 prime
                and md already has lu[i].i0 pre-multiplied in.
            */
            s = -s - 1;

            ad = FLINT_MPZ_REALLOC(az, s + 2);

            flint_mpn_zero(ad, s + 2);
            hi = lo = 0;

            for ( ; i < j; md += s, l++, i++)
            {
                FLINT_ASSERT(lu[i].i0 != 0);
                FLINT_ASSERT(lu[i].i1 == 0);
                FLINT_ASSERT(lu[i].i2 == 0);

                t = mpn_addmul_1(ad, md, s, in[l*1]);
                add_ssaaaa(hi, lo, hi, lo, UWORD(0), t);
            }

            ad[s] = lo;
            ad[s + 1] = hi;
        }
        else
        {
            ad = FLINT_MPZ_REALLOC(az, s + 2);

            flint_mpn_zero(ad, s + 2);

            for ( ; i < j; md += s, i++)
            {
                /* low level combination: 1, 2, or 3 small primes */
                FLINT_ASSERT(l + 1 <= C->num_primes);
                umul_ppmm(hi, lo, in[l*1], lu[i].i0); l++;

                if (lu[i].i2 != 0)
                {
                    FLINT_ASSERT(l + 2 <= C->num_primes);
                    MAC(hi, lo, in[l*1], lu[i].i1); l++;
                    MAC(hi, lo, in[l*1], lu[i].i2); l++;
                /*
                    We have lu[i].mod.n = p0*p1*p2, and each lu[i].i{0|1|2} is
                    strictly less than p1*p2*p3, and the inputs are reduced mod pi.
                    Therefore, the sum is at most (p0*p1*p2-1)*(p0-1+p1-1+p2-1).
                    Since p0*p1*p2 fits into a word, the sum fits into two words
                    and the hi word is less than p0*p1*p2.
                */
                }
                else if (lu[i].i1 != 0)
                {
                    FLINT_ASSERT(l + 1 <= C->num_primes);
                    MAC(hi, lo, in[l*1], lu[i].i1); l++;
                    /* Ditto for two */
                }

                FLINT_ASSERT(hi < lu[i].mod.n);
                NMOD_RED2(t, hi, lo, lu[i].mod);

                /* mid level combination: depends on FMPZ_CRT_UI_CUTOFF */
                hi = mpn_addmul_1(ad, md, s, t);
                add_ssaaaa(ad[s + 1], ad[s], ad[s + 1], ad[s], UWORD(0), hi);
            }
        }

        s += 2;

        MPN_NORM(ad, s);
        az->_mp_size = s;
        _fmpz_demote_val(A + k);

        _fmpz_smod(A + k, A + k, C->crt_P->moduli + k, sign, T + 0);
    }

    FLINT_ASSERT(l == C->num_primes);

    /* high level combination */
    fmpz_swap(T + 0, b);
    _fmpz_multi_CRT_precomp(T, C->crt_P, A, sign);
    fmpz_swap(T + 0, b);
}
