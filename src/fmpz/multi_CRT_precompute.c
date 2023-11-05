/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"


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
        slong s, minp;
        const fmpz * mind;

        minp = j;
        mind = v + j;
        for (s = j + 1; s < i; s++)
        {
            if (fmpz_cmp(v + s, mind) < 0)
            {
                mind = v + s;
                minp = s;
            }
        }
        fmpz_swap(v + j, v + minp);
        FLINT_SWAP(slong, link[j], link[minp]);

        minp = j + 1;
        mind = v + j + 1;
        for (s = j + 2; s < i; s++)
        {
            if (fmpz_cmp(v + s, mind) < 0)
            {
                mind = v + s;
                minp = s;
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

