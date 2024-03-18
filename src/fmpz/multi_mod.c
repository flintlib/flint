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

#include "nmod.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void fmpz_multi_mod_init(fmpz_multi_mod_t P)
{
    P->prog = NULL;
    P->moduli = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
    P->good = 0;
}

void fmpz_multi_mod_clear(fmpz_multi_mod_t P)
{
    slong i;

    for (i = 0; i < P->alloc; i++)
    {
        fmpz_clear(P->prog[i].modulus);
        fmpz_clear(P->moduli + i);
    }

    flint_free(P->prog);
    flint_free(P->moduli);
}

void _fmpz_multi_mod_precomp(
    fmpz * outputs,
    const fmpz_multi_mod_t P,
    const fmpz_t input,
    int sign,
    fmpz * T)
{
    slong i, a, b;
    slong len = P->length;
    _fmpz_multi_mod_instr * instr = P->prog;
    fmpz * t1 = T + P->temp1loc;
    unsigned char * org;
    TMP_INIT;

    TMP_START;

    /*
        Efficiently propagate small inputs without copying:
        ord[i] = 1 means T[i] should be read from input
    */
    org = TMP_ARRAY_ALLOC(P->localsize, unsigned char);

#if FLINT_WANT_ASSERT
    for (i = 0; i < P->localsize; i++)
        org[i] = 2;
#endif

    for (i = 0; i < len; i++)
    {
        a = P->prog[i].in_idx;
        b = P->prog[i].out_idx;

        FLINT_ASSERT(a < 1 || org[a] < 2);

        if (a > 0 && org[a] == 0)
        {
            /* read input from T[a] */

            if (b < 0)
            {
                _fmpz_smod(outputs - b - 1, T + a, instr[i].modulus, sign, t1);
            }
            else
            {
                org[b] = 0;
                fmpz_tdiv_qr(t1, T + b, T + a, instr[i].modulus);
            }
        }
        else
        {
            /* read input from input */

            if (b < 0)
            {
                _fmpz_smod(outputs - b - 1, input, instr[i].modulus, sign, t1);
            }
            else if (fmpz_cmpabs(instr[i].modulus, input) > 0)
            {
                org[b] = 1;
            }
            else
            {
                org[b] = 0;
                fmpz_tdiv_qr(t1, T + b, input, instr[i].modulus);
            }
        }
    }

    TMP_END;
}

void fmpz_multi_mod_precomp(
    fmpz * outputs,
    const fmpz_multi_mod_t P,
    const fmpz_t input,
    int sign)
{
    slong i;
    fmpz * tmp;
    TMP_INIT;

    TMP_START;
    tmp = TMP_ARRAY_ALLOC(P->localsize, fmpz);
    for (i = 0; i < P->localsize; i++)
        fmpz_init(tmp + i);

    _fmpz_multi_mod_precomp(outputs, P, input, sign, tmp);

    for (i = 0; i < P->localsize; i++)
        fmpz_clear(tmp + i);

    TMP_END;
}

static void _fmpz_multi_mod_fit_length(fmpz_multi_mod_t P, slong k)
{
    slong i;

    k = FLINT_MAX(WORD(1), k);

    for (i = k; i < P->alloc; i++)
    {
        fmpz_clear(P->prog[i].modulus);
        fmpz_clear(P->moduli + i);
    }

    P->prog = FLINT_ARRAY_REALLOC(P->prog, k, _fmpz_multi_mod_instr);
    P->moduli = FLINT_ARRAY_REALLOC(P->moduli, k, fmpz);

    for (i = P->alloc; i < k; i++)
    {
        fmpz_init(P->prog[i].modulus);
        fmpz_init(P->moduli + i);
    }

    P->alloc = k;
}

static int _fill_sort(slong * link, fmpz * v, slong j)
{
    while (j >= 0)
    {
        int cmp = fmpz_cmp(v + j, v + j + 1);

        if (fmpz_is_zero(v + j) || fmpz_is_zero(v + j + 1))
            return 0;

        /* hit the smaller branch first */
        if (cmp > 0)
        {
            fmpz_swap(v + j, v + j + 1);
            FLINT_SWAP(slong, link[j], link[j + 1]);
        }

        if (!_fill_sort(link, v, link[j + 0]))
            return 0;

        j = link[j + 1];
    }

    return 1;
}

static void _fill_prog(
    fmpz_multi_mod_t P,
    slong * link,
    fmpz * v,
    slong j,
    slong a_idx)
{
    slong i, b_idx, c_idx;

    FLINT_ASSERT(j >= 0);

    if (link[j] >= 0)
    {
        b_idx = a_idx + 1;
    }
    else
    {
        b_idx = -1 - link[j];
        FLINT_ASSERT(b_idx < P->alloc);
        fmpz_set(P->moduli + b_idx, v + j);
        b_idx = -1 - b_idx;
    }

    if (link[j + 1] >= 0)
    {
        c_idx = a_idx + 1;
    }
    else
    {
        c_idx = -1 - link[j + 1];
        FLINT_ASSERT(c_idx < P->alloc);
        fmpz_set(P->moduli + c_idx, v + j + 1);
        c_idx = -1 - c_idx;
    }

    i = P->length;
    FLINT_ASSERT(i < P->alloc);
    P->prog[i].in_idx = a_idx;
    P->prog[i].out_idx = b_idx;
    fmpz_set(P->prog[i].modulus, v + j);
    P->length = i + 1;

    if (link[j + 0] >= 0)
        _fill_prog(P, link, v, link[j + 0], b_idx);

    i = P->length;
    FLINT_ASSERT(i < P->alloc);
    P->prog[i].in_idx = a_idx;
    P->prog[i].out_idx = c_idx;
    fmpz_set(P->prog[i].modulus, v + j + 1);
    P->length = i + 1;

    if (link[j + 1] >= 0)
        _fill_prog(P, link, v, link[j + 1], c_idx);

    P->localsize = FLINT_MAX(P->localsize, a_idx + 1);
}

int fmpz_multi_mod_precompute(
    fmpz_multi_mod_t P,
    const fmpz * f,
    slong r)
{
    slong i, j;
    slong * link;
    fmpz * v;

    FLINT_ASSERT(r > 0);

    _fmpz_multi_mod_fit_length(P, 2*r);
    P->length = 0;
    P->localsize = 1;
    P->moduli_count = r;
    P->min_modulus_bits = fmpz_bits(f + 0);

    if (r < 2)
    {
        P->good = !fmpz_is_zero(f + 0);

        if (P->good)
        {
            fmpz_abs(P->moduli + 0, f + 0);

            P->prog[0].in_idx = 0;
            P->prog[0].out_idx = -1;
            fmpz_set(P->prog[0].modulus, P->moduli + 0);
            P->length = 1;
        }

        goto done;
    }

    link = FLINT_ARRAY_ALLOC(2*r - 2, slong);
    v = _fmpz_vec_init(2*r - 2);

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

    P->good = _fill_sort(link, v, 2*r - 4);
    if (P->good)
        _fill_prog(P, link, v, 2*r - 4, 0);

    flint_free(link);
    _fmpz_vec_clear(v, 2*r - 2);

done:

    P->temp1loc = P->localsize++;

    if (!P->good)
    {
        P->length = 0;
    }

    return P->good;
}

void fmpz_multi_mod_ui(
    mp_limb_t * out,
    const fmpz_t input,
    const fmpz_comb_t C,
    fmpz_comb_temp_t CT)
{
    slong i, k, l;
    slong stride = 1;
    fmpz * A = CT->A;
    mod_lut_entry * lu;
    slong * offsets;
    slong klen = C->mod_klen;
    fmpz_t ttt;

    /* high level split */
    if (klen == 1)
    {
        *ttt = A[0];
        A[0] = *input;
    }
    else
    {
        _fmpz_multi_mod_precomp(A, C->mod_P, input, -1, CT->T);
    }

    offsets = C->mod_offsets;
    lu = C->mod_lu;

    for (k = 0, i = 0, l = 0; k < klen; k++)
    {
        slong j = offsets[k];

        for ( ; i < j; i++)
        {
            /* mid level split: depends on FMPZ_MOD_UI_CUTOFF */
            mp_limb_t t = fmpz_get_nmod(A + k, lu[i].mod);

            /* low level split: 1, 2, or 3 small primes */
            if (lu[i].mod2.n != 0)
            {
                FLINT_ASSERT(l + 3 <= C->num_primes);
                NMOD_RED(out[l*stride], t, lu[i].mod0); l++;
                NMOD_RED(out[l*stride], t, lu[i].mod1); l++;
                NMOD_RED(out[l*stride], t, lu[i].mod2); l++;
            }
            else if (lu[i].mod1.n != 0)
            {
                FLINT_ASSERT(l + 2 <= C->num_primes);
                NMOD_RED(out[l*stride], t, lu[i].mod0); l++;
                NMOD_RED(out[l*stride], t, lu[i].mod1); l++;
            }
            else
            {
                FLINT_ASSERT(l + 1 <= C->num_primes);
                out[l*stride] = t; l++;
            }
        }
    }

    FLINT_ASSERT(l == C->num_primes);

    if (klen == 1)
        A[0] = *ttt;
}
