/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod_poly.h"

void nmod_poly_multi_crt_init(nmod_poly_multi_crt_t P)
{
    P->prog = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
    P->temp2loc = 0;
    P->good = 0;
}

static void _nmod_poly_multi_crt_fit_length(nmod_poly_multi_crt_t P, slong k)
{
    k = FLINT_MAX(WORD(1), k);

    if (P->alloc == 0)
    {
        FLINT_ASSERT(P->prog == NULL);
        P->prog = (_nmod_poly_multi_crt_prog_instr *) flint_malloc(k
                                     *sizeof(_nmod_poly_multi_crt_prog_instr));
        P->alloc = k;
    }
    else if (k > P->alloc)
    {
        FLINT_ASSERT(P->prog != NULL);
        P->prog = (_nmod_poly_multi_crt_prog_instr *) flint_realloc(P->prog, k
                                     *sizeof(_nmod_poly_multi_crt_prog_instr));
        P->alloc = k;
    }
}

static void _nmod_poly_multi_crt_set_length(nmod_poly_multi_crt_t P, slong k)
{
    slong i;

    FLINT_ASSERT(k <= P->length);

    for (i = k; i < P->length; i++)
    {
        nmod_poly_clear(P->prog[i].modulus);
        nmod_poly_clear(P->prog[i].idem);
    }
    P->length = k;
}

void nmod_poly_multi_crt_clear(nmod_poly_multi_crt_t P)
{
    _nmod_poly_multi_crt_set_length(P, 0);

    if (P->alloc > 0)
    {
        flint_free(P->prog);
    }
}

typedef struct {
    slong idx;
    slong degree;
} index_deg_pair;

/*
    combine all moduli in [start, stop)
    return index of instruction that computes the result
*/
static slong _push_prog(
    nmod_poly_multi_crt_t P,
    const nmod_poly_struct * const * moduli,
    const index_deg_pair * perm,
    slong ret_idx,
    slong start,
    slong stop)
{
    slong i, mid;
    slong b_idx, c_idx;
    slong lefttot, righttot;
    slong leftret, rightret;
    nmod_poly_struct * leftmodulus, * rightmodulus;

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
        leftmodulus = (nmod_poly_struct *) moduli[perm[start].idx];
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
        rightmodulus = (nmod_poly_struct *) moduli[perm[mid].idx];
    }

    /* check if nmod_poly_invmod is going to throw */
    if (nmod_poly_degree(leftmodulus) < 1 || nmod_poly_degree(rightmodulus) < 1)
    {
        P->good = 0;
        return -1;
    }

    /* compile [start, end) */
    i = P->length;
    _nmod_poly_multi_crt_fit_length(P, i + 1);
    nmod_poly_init_mod(P->prog[i].modulus, rightmodulus->mod);
    nmod_poly_init_mod(P->prog[i].idem, rightmodulus->mod);
    P->good = P->good && nmod_poly_invmod(P->prog[i].modulus, leftmodulus, rightmodulus);
    nmod_poly_mul(P->prog[i].idem, leftmodulus, P->prog[i].modulus);
    nmod_poly_mul(P->prog[i].modulus, leftmodulus, rightmodulus);
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
int nmod_poly_multi_crt_precompute_p(
    nmod_poly_multi_crt_t P,
    const nmod_poly_struct * const * moduli,
    slong len)
{
    slong i;
    index_deg_pair * perm;
    TMP_INIT;

    FLINT_ASSERT(len > 0);
    for (i = 1; i < len; i++)
    {
        FLINT_ASSERT(moduli[i - 1]->mod.n == moduli[i]->mod.n);
    }

    TMP_START;
    perm = (index_deg_pair *) TMP_ALLOC(len * sizeof(index_deg_pair));

    for (i = 0; i < len; i++)
    {
        perm[i].idx = i;
        perm[i].degree = nmod_poly_degree(moduli[i]);
    }

    /* make perm sort the degs so that degs[perm[i-1]] <= degs[perm[i-0]] */
    qsort(perm, len, sizeof(index_deg_pair),
                        (int(*)(const void*, const void*)) index_deg_pair_cmp);
    for (i = 0; i < len; i++)
    {
        FLINT_ASSERT(perm[i].degree == nmod_poly_degree(moduli[perm[i].idx]));
        FLINT_ASSERT(i == 0 || perm[i - 1].degree <= perm[i].degree);
    }

    _nmod_poly_multi_crt_fit_length(P, FLINT_MAX(WORD(1), len - 1));
    _nmod_poly_multi_crt_set_length(P, 0);
    P->localsize = 1;
    P->good = 1;

    if (1 < len)
    {
        _push_prog(P, moduli, perm, 0, 0, len);
    }
    else
    {
        /*
            There is only one modulus. Let's compute as
                output[0] = input[0] + 0*(input[0] - input[0]) mod moduli[0]
        */
        i = 0;
        nmod_poly_init_mod(P->prog[i].modulus, moduli[0]->mod);
        nmod_poly_init_mod(P->prog[i].idem, moduli[0]->mod);
        nmod_poly_set(P->prog[i].modulus, moduli[0]);
        P->prog[i].a_idx = 0;
        P->prog[i].b_idx = -WORD(1);
        P->prog[i].c_idx = -WORD(1);
        P->length = i + 1;

        P->good = !nmod_poly_is_zero(moduli[0]);
    }

    if (!P->good)
    {
        _nmod_poly_multi_crt_set_length(P, 0);
    }

    /* two more spots for temporaries */
    P->temp1loc = P->localsize++;
    P->temp2loc = P->localsize++;

    TMP_END;

    return P->good;
}

int nmod_poly_multi_crt_precompute(
    nmod_poly_multi_crt_t P,
    const nmod_poly_struct * moduli,
    slong len)
{
    int success;
    slong i;
    const nmod_poly_struct ** m;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;

    m = TMP_ALLOC(len * sizeof(nmod_poly_struct *));
    for (i = 0; i < len; i++)
    {
        m[i] = moduli + i;
    }

    /* GCC really wants to complain about this one */
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    success = nmod_poly_multi_crt_precompute_p(P, m, len);
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic pop
#endif
    TMP_END;

    return success;
}


void nmod_poly_multi_crt_precomp(
    nmod_poly_t output,
    const nmod_poly_multi_crt_t P,
    const nmod_poly_struct * inputs)
{
    slong i;
    nmod_poly_struct * out;
    TMP_INIT;

    TMP_START;
    out = (nmod_poly_struct *) TMP_ALLOC(P->localsize
                                                    *sizeof(nmod_poly_struct));
    for (i = 0; i < P->localsize; i++)
    {
        nmod_poly_init_mod(out + i, inputs[0].mod);
    }

    nmod_poly_swap(out + 0, output);
    _nmod_poly_multi_crt_run(out, P, inputs);
    nmod_poly_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
    {
        nmod_poly_clear(out + i);
    }

    TMP_END;
}

void nmod_poly_multi_crt_precomp_p(
    nmod_poly_t output,
    const nmod_poly_multi_crt_t P,
    const nmod_poly_struct * const * inputs)
{
    slong i;
    nmod_poly_struct * out;
    TMP_INIT;

    TMP_START;
    out = (nmod_poly_struct *) TMP_ALLOC(P->localsize
                                                    *sizeof(nmod_poly_struct));
    for (i = 0; i < P->localsize; i++)
    {
        nmod_poly_init_mod(out + i, inputs[0]->mod);
    }

    nmod_poly_swap(out + 0, output);
    _nmod_poly_multi_crt_run_p(out, P, inputs);
    nmod_poly_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
    {
        nmod_poly_clear(out + i);
    }

    TMP_END;
}

int nmod_poly_multi_crt(
    nmod_poly_t output,
    const nmod_poly_struct * moduli,
    const nmod_poly_struct * values,
    slong len)
{
    int success;
    slong i;
    nmod_poly_multi_crt_t P;
    nmod_poly_struct * out;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;

    nmod_poly_multi_crt_init(P);
    success = nmod_poly_multi_crt_precompute(P, moduli, len);

    out = (nmod_poly_struct *) TMP_ALLOC(P->localsize
                                                    *sizeof(nmod_poly_struct));
    for (i = 0; i < P->localsize; i++)
    {
        nmod_poly_init_mod(out + i, values[0].mod);
    }

    nmod_poly_swap(out + 0, output);
    _nmod_poly_multi_crt_run(out, P, values);
    nmod_poly_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
    {
        nmod_poly_clear(out + i);
    }

    nmod_poly_multi_crt_clear(P);

    TMP_END;

    return success;
}

/*
    If P was set with a call to nmod_poly_multi_crt_compile(P, m, len), return
    in outputs[0] polynomial r of smallest degree such that
        r = inputs[0] mod m[0]
        r = inputs[1] mod m[1]
            ...
        r = inputs[len-1] mod m[len-1]
    For thread safety "outputs" is expected to have enough space for all
    temporaries, thus should be at least as long as P->localsize.
*/
void _nmod_poly_multi_crt_run(
    nmod_poly_struct * outputs,
    const nmod_poly_multi_crt_t P,
    const nmod_poly_struct * inputs)
{
    slong i;
    slong a, b, c;
    const nmod_poly_struct * B, * C;
    nmod_poly_struct * A, * t1, * t2;

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

        FLINT_ASSERT(A->mod.n == P->prog[i].modulus->mod.n);
        FLINT_ASSERT(B->mod.n == P->prog[i].modulus->mod.n);
        FLINT_ASSERT(C->mod.n == P->prog[i].modulus->mod.n);

        /* A = B + I*(C - B) mod M */
        nmod_poly_sub(t1, B, C);
        nmod_poly_mul(t2, P->prog[i].idem, t1);
        nmod_poly_sub(t1, B, t2);

        if (nmod_poly_degree(t1) < nmod_poly_degree(P->prog[i].modulus))
        {
            nmod_poly_swap(A, t1);
        }
        else
        {
            nmod_poly_rem(A, t1, P->prog[i].modulus);
        }

        /* last calculation should write answer to outputs[0] */
        if (i + 1 >= P->length)
        {
            FLINT_ASSERT(A == outputs + 0);
        }
    }
}

void _nmod_poly_multi_crt_run_p(
    nmod_poly_struct * outputs,
    const nmod_poly_multi_crt_t P,
    const nmod_poly_struct * const * inputs)
{
    slong i;
    slong a, b, c;
    const nmod_poly_struct * B, * C;
    nmod_poly_struct * A, * t1, * t2;

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

        FLINT_ASSERT(A->mod.n == P->prog[i].modulus->mod.n);
        FLINT_ASSERT(B->mod.n == P->prog[i].modulus->mod.n);
        FLINT_ASSERT(C->mod.n == P->prog[i].modulus->mod.n);

        /* A = B + I*(C - B) mod M */
        nmod_poly_sub(t1, B, C);
        nmod_poly_mul(t2, P->prog[i].idem, t1);
        nmod_poly_sub(t1, B, t2);

        if (nmod_poly_degree(t1) < nmod_poly_degree(P->prog[i].modulus))
        {
            nmod_poly_swap(A, t1);
        }
        else
        {
            nmod_poly_rem(A, t1, P->prog[i].modulus);
        }

        /* last calculation should write answer to outputs[0] */
        if (i + 1 >= P->length)
        {
            FLINT_ASSERT(A == outputs + 0);
        }
    }
}
