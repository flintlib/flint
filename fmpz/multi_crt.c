/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void fmpz_multi_crt_init(fmpz_multi_crt_t P)
{
    P->prog = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
    P->temp2loc = 0;
    P->good = 0;
}

static void _fmpz_multi_crt_fit_length(fmpz_multi_crt_t P, slong k)
{
    k = FLINT_MAX(WORD(1), k);

    if (P->alloc == 0)
    {
        FLINT_ASSERT(P->prog == NULL);
        P->prog = (_fmpz_multi_crt_prog_instr *) flint_malloc(k
                                          *sizeof(_fmpz_multi_crt_prog_instr));
        P->alloc = k;
    }
    else if (k > P->alloc)
    {
        FLINT_ASSERT(P->prog != NULL);
        P->prog = (_fmpz_multi_crt_prog_instr *) flint_realloc(P->prog, k
                                          *sizeof(_fmpz_multi_crt_prog_instr));
        P->alloc = k;
    }
}

static void _fmpz_multi_crt_set_length(fmpz_multi_crt_t P, slong k)
{
    slong i;

    FLINT_ASSERT(k <= P->length);

    for (i = k; i < P->length; i++)
    {
        fmpz_clear(P->prog[i].modulus);
        fmpz_clear(P->prog[i].idem);
    }
    P->length = k;
}

void fmpz_multi_crt_clear(fmpz_multi_crt_t P)
{
    _fmpz_multi_crt_set_length(P, 0);

    if (P->alloc > 0)
    {
        flint_free(P->prog);
    }
}

typedef struct {
    slong idx;
    flint_bitcnt_t degree;
} index_deg_pair;

/*
    combine all moduli in [start, stop)
    return index of instruction that computes the result
*/
static slong _push_prog(
    fmpz_multi_crt_t P,
    const fmpz * const * moduli,
    const index_deg_pair * perm,
    slong ret_idx,
    slong start,
    slong stop)
{
    slong i, mid;
    slong b_idx, c_idx;
    flint_bitcnt_t lefttot, righttot;
    slong leftret, rightret;
    fmpz * leftmodulus, * rightmodulus;

    /* we should have at least 2 moduli */
    FLINT_ASSERT(start + 1 < stop);

    mid = start + (stop - start)/2;

    FLINT_ASSERT(start < mid);
    FLINT_ASSERT(mid < stop);

    lefttot = 0;
    for (i = start; i < mid; i++)
    {
        lefttot += perm[i].degree;
    }

    righttot = 0;
    for (i = mid; i < stop; i++)
    {
        righttot += perm[i].degree;
    }

    /* try to balance the total degree on left and right */
    while (lefttot < righttot
            && mid + 1 < stop
            && perm[mid].degree < righttot - lefttot)
    {
        lefttot += perm[mid].degree;
        righttot -= perm[mid].degree;
        mid++;
    }

    P->localsize = FLINT_MAX(P->localsize, 1 + ret_idx);

    /* compile left [start, mid) */
    if (start + 1 < mid)
    {
        b_idx = ret_idx + 1;
        leftret = _push_prog(P, moduli, perm, b_idx, start, mid);
        if (!P->good)
        {
            return -1;
        }
        leftmodulus = P->prog[leftret].modulus;
    }
    else
    {
        b_idx = -1 - perm[start].idx;
        leftmodulus = (fmpz *) moduli[perm[start].idx];
    }

    /* compile right [mid, end) */
    if (mid + 1 < stop)
    {
        c_idx = ret_idx + 2;
        rightret = _push_prog(P, moduli, perm, c_idx, mid, stop);
        if (!P->good)
        {
            return -1;
        }
        rightmodulus = P->prog[rightret].modulus;
    }
    else
    {
        c_idx = -1 - perm[mid].idx;
        rightmodulus = (fmpz *) moduli[perm[mid].idx];
    }

    /* check if fmpz_invmod is going to throw */
    if (fmpz_is_zero(leftmodulus) || fmpz_is_zero(rightmodulus))
    {
        P->good = 0;
        return -1;
    }

    /* compile [start, end) */
    i = P->length;
    _fmpz_multi_crt_fit_length(P, i + 1);
    fmpz_init(P->prog[i].modulus);
    fmpz_init(P->prog[i].idem);
    P->good = P->good &&
                 fmpz_invmod(P->prog[i].modulus, leftmodulus, rightmodulus);
    fmpz_mul(P->prog[i].idem, leftmodulus, P->prog[i].modulus);
    fmpz_mul(P->prog[i].modulus, leftmodulus, rightmodulus);
    P->prog[i].a_idx = ret_idx;
    P->prog[i].b_idx = b_idx;
    P->prog[i].c_idx = c_idx;
    P->length = i + 1;

    return i;
}

static int index_deg_pair_cmp(
    const index_deg_pair * lhs,
    const index_deg_pair * rhs)
{
    return (lhs->degree < rhs->degree) ? -1 : (lhs->degree > rhs->degree);
}

/*
    Return 1 if moduli can be CRT'ed, 0 otherwise.
    A return of 0 means that future calls to run will leave output undefined.
*/
int fmpz_multi_crt_precompute_p(
    fmpz_multi_crt_t P,
    const fmpz * const * moduli,
    slong len)
{
    slong i;
    index_deg_pair * perm;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;
    perm = (index_deg_pair *) TMP_ALLOC(len * sizeof(index_deg_pair));

    for (i = 0; i < len; i++)
    {
        perm[i].idx = i;
        perm[i].degree = fmpz_bits(moduli[i]);
    }

    /* make perm sort the degs so that degs[perm[i-1]] <= degs[perm[i-0]] */
    qsort(perm, len, sizeof(index_deg_pair),
                        (int(*)(const void*, const void*)) index_deg_pair_cmp);
    for (i = 0; i < len; i++)
    {
        FLINT_ASSERT(perm[i].degree == fmpz_bits(moduli[perm[i].idx]));
        FLINT_ASSERT(i == 0 || perm[i - 1].degree <= perm[i].degree);
    }

    _fmpz_multi_crt_fit_length(P, FLINT_MAX(WORD(1), len - 1));
    _fmpz_multi_crt_set_length(P, 0);
    P->localsize = 1;
    P->good = 1;

    if (1 < len)
    {
        _push_prog(P, moduli, perm, 0, 0, len);
    }
    else
    {
        /*
            There is only one modulus. Lets compute as
                output[0] = input[0] + 0*(input[0] - input[0]) mod moduli[0]
        */
        i = 0;
        fmpz_init(P->prog[i].modulus);
        fmpz_init(P->prog[i].idem);
        fmpz_set(P->prog[i].modulus, moduli[0]);
        P->prog[i].a_idx = 0;
        P->prog[i].b_idx = -WORD(1);
        P->prog[i].c_idx = -WORD(1);
        P->length = i + 1;

        P->good = !fmpz_is_zero(moduli[0]);
    }

    if (!P->good)
    {
        _fmpz_multi_crt_set_length(P, 0);
    }

    /* two more spots for temporaries */
    P->temp1loc = P->localsize++;
    P->temp2loc = P->localsize++;

    TMP_END;

    return P->good;
}

int fmpz_multi_crt_precompute(
    fmpz_multi_crt_t P,
    const fmpz * moduli,
    slong len)
{
    int success;
    slong i;
    const fmpz ** m;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;

    m = (const fmpz **) TMP_ALLOC(len*sizeof(fmpz *));
    for (i = 0; i < len; i++)
    {
        m[i] = moduli + i;
    }

    success = fmpz_multi_crt_precompute_p(P, (const fmpz * const *) m, len);

    TMP_END;

    return success;
}



void fmpz_multi_crt_precomp(
    fmpz_t output,
    const fmpz_multi_crt_t P,
    const fmpz * inputs)
{
    slong i;
    fmpz * out;
    TMP_INIT;

    TMP_START;
    out = (fmpz *) TMP_ALLOC(P->localsize*sizeof(fmpz));
    for (i = 0; i < P->localsize; i++)
    {
        fmpz_init(out + i);
    }

    fmpz_swap(out + 0, output);
    _fmpz_multi_crt_run(out, P, inputs);
    fmpz_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
    {
        fmpz_clear(out + i);
    }

    TMP_END;
}

void fmpz_multi_crt_precomp_p(
    fmpz_t output,
    const fmpz_multi_crt_t P,
    const fmpz * const * inputs)
{
    slong i;
    fmpz * out;
    TMP_INIT;

    TMP_START;
    out = (fmpz *) TMP_ALLOC(P->localsize*sizeof(fmpz));
    for (i = 0; i < P->localsize; i++)
    {
        fmpz_init(out + i);
    }

    fmpz_swap(out + 0, output);
    _fmpz_multi_crt_run_p(out, P, inputs);
    fmpz_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
    {
        fmpz_clear(out + i);
    }

    TMP_END;
}


int fmpz_multi_crt(
    fmpz_t output,
    const fmpz * moduli,
    const fmpz * values,
    slong len)
{
    int success;
    slong i;
    fmpz_multi_crt_t P;
    fmpz * out;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;

    fmpz_multi_crt_init(P);
    success = fmpz_multi_crt_precompute(P, moduli, len);

    out = (fmpz *) TMP_ALLOC(P->localsize*sizeof(fmpz));
    for (i = 0; i < P->localsize; i++)
    {
        fmpz_init(out + i);
    }

    fmpz_swap(out + 0, output);
    _fmpz_multi_crt_run(out, P, values);
    fmpz_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
    {
        fmpz_clear(out + i);
    }

    fmpz_multi_crt_clear(P);

    TMP_END;

    return success;
}


/*
    If P was set with a call to fmpz_poly_crt_compile(P, m, len), return
    in outputs[0] signed integer r of smallest abs value such that
        r = inputs[0] mod m[0]
        r = inputs[1] mod m[1]
            ...
        r = inputs[len-1] mod m[len-1]
    For thread safety "outputs" is expected to have enough space for all
    temporaries, thus should be at least as long as P->localsize.
*/

void _fmpz_multi_crt_run(
    fmpz * outputs,
    const fmpz_multi_crt_t P,
    const fmpz * inputs)
{
    slong i;
    slong a, b, c;
    const fmpz * B, * C;
    fmpz * A, * t1, * t2;

    t1 = outputs + P->temp1loc;
    t2 = outputs + P->temp2loc;

    for (i = 0; i < P->length; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;
        c = P->prog[i].c_idx;
        FLINT_ASSERT(a >= 0);
        A = outputs + a;
        B = b < 0 ? inputs + (-b-1) : outputs + b;
        C = c < 0 ? inputs + (-c-1) : outputs + c;

        /* A = B + I*(C - B) mod M */
        fmpz_sub(t1, B, C);
        fmpz_mul(t2, P->prog[i].idem, t1);
        fmpz_sub(t1, B, t2);
        fmpz_smod(A, t1, P->prog[i].modulus);

        /* last calculation should write answer to outputs[0] */
        if (i + 1 >= P->length)
        {
            FLINT_ASSERT(A == outputs + 0);
        }
    }
}

void _fmpz_multi_crt_run_p(
    fmpz * outputs,
    const fmpz_multi_crt_t P,
    const fmpz * const * inputs)
{
    slong i;
    slong a, b, c;
    const fmpz * B, * C;
    fmpz * A, * t1, * t2;

    t1 = outputs + P->temp1loc;
    t2 = outputs + P->temp2loc;

    for (i = 0; i < P->length; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;
        c = P->prog[i].c_idx;
        FLINT_ASSERT(a >= 0);
        A = outputs + a;
        B = b < 0 ? inputs[-b-1] : outputs + b;
        C = c < 0 ? inputs[-c-1] : outputs + c;

        /* A = B + I*(C - B) mod M */
        fmpz_sub(t1, B, C);
        fmpz_mul(t2, P->prog[i].idem, t1);
        fmpz_sub(t1, B, t2);
        fmpz_smod(A, t1, P->prog[i].modulus);

        /* last calculation should write answer to outputs[0] */
        if (i + 1 >= P->length)
        {
            FLINT_ASSERT(A == outputs + 0);
        }
    }
}
