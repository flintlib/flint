/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2017 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <string.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"
#include "gr_vec.h"
#include "gr_generic.h"
#include "gr_mpoly.h"

/*
    GR port of fmpz_mpoly_divides_monagan_pearce.

    Set Q to A/B if the division is exact, returning GR_SUCCESS, and set *lenout
    to the length of the quotient.  Return GR_DOMAIN if the division is provably
    not exact (in which case *lenout is set to 0), or GR_UNABLE if the ring
    arithmetic could not decide the question.

    The two coefficient code paths of fmpz_mpoly (machine-word vs multiprecision
    accumulation) are replaced here by the two accumulation strategies used in
    gr_mpoly_mul_heap: a fast path that gathers the divisor/quotient operands
    of each diagonal into two shallow vectors and accumulates with a single
    _gr_vec_dot (used when the ring overloads VEC_DOT), and a generic path that
    accumulates term by term as acc -= B[i]*Q[j] using a single preallocated
    product temporary (avoiding the per-call temporary of gr_submul).

    Exact division of coefficients is handled following gr_poly_divexact: if the
    leading coefficient of B is a unit we invert it once and multiply, otherwise
    we use gr_div, which over an integral domain such as Z returns GR_DOMAIN
    precisely when the division has a nonzero remainder.
*/

/* shallow copy helpers (see mul_johnson.c) */
#define SET_SHALLOW_GENERIC(a,i,b,j) memcpy(GR_ENTRY(a, i, sz), GR_ENTRY(b, j, sz), sz)
#define SET_SHALLOW1(a,i,b,j) ((uint8_t *) a)[i] = ((uint8_t *) b)[j]
#define SET_SHALLOW2(a,i,b,j) ((uint16_t *) a)[i] = ((uint16_t *) b)[j]
#define SET_SHALLOW4(a,i,b,j) ((uint32_t *) a)[i] = ((uint32_t *) b)[j]

#if FLINT_BITS == 64
#define SET_SHALLOW8(a,i,b,j) ((uint64_t *) a)[i] = ((uint64_t *) b)[j]
#define SET_SHALLOW16(a,i,b,j) ((uint64_t *) a)[2 * (i)] = ((uint64_t *) b)[2 * (j)]; ((uint64_t *) a)[2 * (i) + 1] = ((uint64_t *) b)[2 * (j) + 1]
#else
#define SET_SHALLOW8(a,i,b,j) gr_set_shallow(GR_ENTRY(a, i, sz), GR_ENTRY(b, j, sz), cctx)
#define SET_SHALLOW16(a,i,b,j) gr_set_shallow(GR_ENTRY(a, i, sz), GR_ENTRY(b, j, sz), cctx)
#endif

/* gather product node (poly3[i], p1[j]) into dot_a/dot_b; record dividend node */
#define DIVIDES_FETCH_NODE(SET_SHALLOW) \
    do { \
        store[store_len++] = x->i; \
        store[store_len++] = x->j; \
        if (x->i == -UWORD(1)) \
        { \
            have_dividend = 1; \
            dividend_j = x->j; \
        } \
        else \
        { \
            hind[x->i] |= WORD(1); \
            SET_SHALLOW(dot_a, dot_len, coeff3, x->i); \
            SET_SHALLOW(dot_b, dot_len, p1, x->j); \
            dot_len++; \
        } \
    } while (0)

/*
    Divide the accumulated value by the leading coefficient of B, storing the
    quotient coefficient in GR_ENTRY(p1, k, sz).

    On entry *qstatus is the status of forming acc.  lc_is_unit indicates that
    lc_inv = 1/coeff3[0] has been precomputed.  Returns GR_SUCCESS if the
    coefficient divided exactly, GR_DOMAIN if not, or GR_UNABLE.
*/
#define DIVIDES_COEFF(acc) \
    (lc_is_one  ? gr_set(GR_ENTRY(p1, k, sz), acc, cctx) \
     : lc_is_unit ? gr_mul(GR_ENTRY(p1, k, sz), acc, lc_inv, cctx) \
                  : gr_div(GR_ENTRY(p1, k, sz), acc, coeff3, cctx))


static int _gr_mpoly_divides_heap1(
    slong * lenout,
    gr_ptr * poly1, ulong ** exp1, slong * alloc, slong * exps_alloc,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_srcptr coeff3, const ulong * exp3, slong len3,
    flint_bitcnt_t bits, ulong maskhi,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    slong i, j, k, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    gr_ptr p1 = *poly1;
    ulong * e1 = *exp1;
    slong * hind;
    ulong mask, exp, maxexp = exp2[len2 - 1];
    int lt_divides;
    gr_ptr acc, lc_inv, pp, dot_a, dot_b;
    slong dot_len;
    int have_dividend, lc_is_unit, lc_is_one, cstatus;
    slong dividend_j = 0;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;

    GR_TMP_INIT3(acc, lc_inv, pp, cctx);
    dot_a = flint_malloc(2 * len3 * sz);
    dot_b = GR_ENTRY(dot_a, len3, sz);

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    /* precompute 1/lc(B) if it is a unit */
    lc_is_one = (gr_is_one(coeff3, cctx) == T_TRUE);
    lc_is_unit = lc_is_one || (gr_inv(lc_inv, coeff3, cctx) == GR_SUCCESS);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

        k++;
        _gr_mpoly_fit_length(&p1, alloc, &e1, exps_alloc, 1, k + 1, ctx);

        lt_divides = mpoly_monomial_divides1(e1 + k, exp, exp3[0], mask);

        /* accumulate acc = (dividend term) - sum(coeff3[i]*p1[j]) */
        dot_len = 0;
        have_dividend = 0;
        store_len = 0;

        if (have_fast_dot)
        {
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    if (sz == 1)       { DIVIDES_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { DIVIDES_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { DIVIDES_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { DIVIDES_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { DIVIDES_FETCH_NODE(SET_SHALLOW16); }
                    else               { DIVIDES_FETCH_NODE(SET_SHALLOW_GENERIC); }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);

            status |= _gr_vec_dot(acc,
                        have_dividend ? GR_ENTRY(coeff2, dividend_j, sz) : NULL,
                        1, dot_a, dot_b, dot_len, cctx);
        }
        else
        {
            status |= gr_zero(acc, cctx);
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    store[store_len++] = x->i;
                    store[store_len++] = x->j;
                    if (x->i == -UWORD(1))
                        status |= gr_add(acc, acc, GR_ENTRY(coeff2, x->j, sz), cctx);
                    else
                    {
                        hind[x->i] |= WORD(1);
                        {
                            status |= gr_mul(pp, GR_ENTRY(coeff3, x->i, sz), GR_ENTRY(p1, x->j, sz), cctx);
                            status |= gr_sub(acc, acc, pp, cctx);
                        }
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        /* process nodes taken from the heap */
        while (store_len > 0)
        {
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
            else
            {
                /* should we go right? */
                if ((i + 1 < len3) && (hind[i + 1] == 2*j + 1))
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
                /* should we go up? */
                if (j + 1 == k)
                {
                    s++;
                }
                else if (((hind[i] & 1) == 1) &&
                         ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        /* is the accumulated coefficient zero? then this term cancels */
        switch (gr_is_zero(acc, cctx))
        {
            case T_TRUE:
                k--;
                continue;
            case T_FALSE:
                break;
            default:
                status |= GR_UNABLE;
                goto unable;
        }

        /* divide the accumulated coefficient by the leading coefficient of B */
        cstatus = DIVIDES_COEFF(acc);
        if (cstatus == GR_DOMAIN)
            goto not_exact_division;
        if (cstatus != GR_SUCCESS)
        {
            status |= cstatus;
            goto unable;
        }

        if (!lt_divides || (exp^maskhi) < (maxexp^maskhi))
            goto not_exact_division;

        /* put newly generated quotient term back into the heap if necessary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = k;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
    }

    k++;

    if (status != GR_SUCCESS)
        goto unable;

    *lenout = k;
    status = GR_SUCCESS;

cleanup:

    GR_TMP_CLEAR3(acc, lc_inv, pp, cctx);
    flint_free(dot_a);

    (*poly1) = p1;
    (*exp1) = e1;

    TMP_END;

    return status;

not_exact_division:
    *lenout = 0;
    status = GR_DOMAIN;
    goto cleanup;

unable:
    *lenout = 0;
    status = GR_UNABLE;
    goto cleanup;
}


static int _gr_mpoly_divides_heap(
    slong * lenout,
    gr_ptr * poly1, ulong ** exp1, slong * alloc, slong * exps_alloc,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_srcptr coeff3, const ulong * exp3, slong len3,
    flint_bitcnt_t bits, slong N, const ulong * cmpmask,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int have_fast_dot;
    slong i, j, k, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    gr_ptr p1 = *poly1;
    ulong * e1 = *exp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    int lt_divides;
    gr_ptr acc, lc_inv, pp, dot_a, dot_b;
    slong dot_len;
    int have_dividend, lc_is_unit, lc_is_one, cstatus;
    slong dividend_j = 0;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (N == 1)
        return _gr_mpoly_divides_heap1(lenout, poly1, exp1, alloc, exps_alloc,
                   coeff2, exp2, len2, coeff3, exp3, len3, bits, cmpmask[0], ctx);

    have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);

    TMP_START;

    GR_TMP_INIT3(acc, lc_inv, pp, cctx);
    dot_a = flint_malloc(2 * len3 * sz);
    dot_b = GR_ENTRY(dot_a, len3, sz);

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    lc_is_one = (gr_is_one(coeff3, cctx) == T_TRUE);
    lc_is_unit = lc_is_one || (gr_inv(lc_inv, coeff3, cctx) == GR_SUCCESS);

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto not_exact_division;
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;
        }

        k++;
        _gr_mpoly_fit_length(&p1, alloc, &e1, exps_alloc, N, k + 1, ctx);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(e1 + k*N, exp, exp3, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(e1 + k*N, exp, exp3, N, bits);

        dot_len = 0;
        have_dividend = 0;
        store_len = 0;

        if (have_fast_dot)
        {
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    if (sz == 1)       { DIVIDES_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { DIVIDES_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { DIVIDES_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { DIVIDES_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { DIVIDES_FETCH_NODE(SET_SHALLOW16); }
                    else               { DIVIDES_FETCH_NODE(SET_SHALLOW_GENERIC); }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            status |= _gr_vec_dot(acc,
                        have_dividend ? GR_ENTRY(coeff2, dividend_j, sz) : NULL,
                        1, dot_a, dot_b, dot_len, cctx);
        }
        else
        {
            status |= gr_zero(acc, cctx);
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    store[store_len++] = x->i;
                    store[store_len++] = x->j;
                    if (x->i == -UWORD(1))
                        status |= gr_add(acc, acc, GR_ENTRY(coeff2, x->j, sz), cctx);
                    else
                    {
                        hind[x->i] |= WORD(1);
                        {
                            status |= gr_mul(pp, GR_ENTRY(coeff3, x->i, sz), GR_ENTRY(p1, x->j, sz), cctx);
                            status |= gr_sub(acc, acc, pp, cctx);
                        }
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        /* process nodes taken from the heap */
        while (store_len > 0)
        {
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
            {
                /* should we go right? */
                if ((i + 1 < len3) && (hind[i + 1] == 2*j + 1))
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                /* should we go up? */
                if (j + 1 == k)
                {
                    s++;
                }
                else if (((hind[i] & 1) == 1) &&
                         ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        switch (gr_is_zero(acc, cctx))
        {
            case T_TRUE:
                k--;
                continue;
            case T_FALSE:
                break;
            default:
                status |= GR_UNABLE;
                goto unable;
        }

        cstatus = DIVIDES_COEFF(acc);
        if (cstatus == GR_DOMAIN)
            goto not_exact_division;
        if (cstatus != GR_SUCCESS)
        {
            status |= cstatus;
            goto unable;
        }

        if (!lt_divides ||
                mpoly_monomial_gt(exp2 + (len2 - 1)*N, exp, N, cmpmask))
            goto not_exact_division;

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = k;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;
    }

    k++;

    if (status != GR_SUCCESS)
        goto unable;

    *lenout = k;
    status = GR_SUCCESS;

cleanup:

    GR_TMP_CLEAR3(acc, lc_inv, pp, cctx);
    flint_free(dot_a);

    (*poly1) = p1;
    (*exp1) = e1;

    TMP_END;

    return status;

not_exact_division:
    *lenout = 0;
    status = GR_DOMAIN;
    goto cleanup;

unable:
    *lenout = 0;
    status = GR_UNABLE;
    goto cleanup;
}


int gr_mpoly_divides_heap(
    gr_mpoly_t Q,
    const gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    slong i, N, lenq = 0;
    flint_bitcnt_t exp_bits;
    fmpz * max_fields2, * max_fields3;
    ulong * cmpmask;
    ulong * exp2 = A->exps, * exp3 = B->exps, * expq;
    int easy_exit, free2 = 0, free3 = 0;
    ulong mask;
    int status;
    TMP_INIT;

    if (B->length == 0)
        return GR_DOMAIN;     /* division by zero: not divisible */

    if (A->length == 0)
        return gr_mpoly_zero(Q, ctx);

    /* The leading monomials must be known */
    if (gr_is_zero(B->coeffs, GR_MPOLY_CCTX(ctx)) != T_FALSE ||
        gr_is_zero(A->coeffs, GR_MPOLY_CCTX(ctx)) != T_FALSE)
        return GR_UNABLE;

    TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    max_fields3 = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
    {
        fmpz_init(max_fields2 + i);
        fmpz_init(max_fields3 + i);
    }

    mpoly_max_fields_fmpz(max_fields2, A->exps, A->length, A->bits, mctx);
    mpoly_max_fields_fmpz(max_fields3, B->exps, B->length, B->bits, mctx);

    easy_exit = 0;
    for (i = 0; i < mctx->nfields; i++)
    {
        /* cannot be exact if any max field of A is less than that of B */
        if (fmpz_cmp(max_fields2 + i, max_fields3 + i) < 0)
            easy_exit = 1;
    }

    exp_bits = _fmpz_vec_max_bits(max_fields2, mctx->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, A->bits);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    for (i = 0; i < mctx->nfields; i++)
    {
        fmpz_clear(max_fields2 + i);
        fmpz_clear(max_fields3 + i);
    }

    if (easy_exit)
    {
        GR_IGNORE(gr_mpoly_zero(Q, ctx));
        TMP_END;
        return GR_DOMAIN;
    }

    N = mpoly_words_per_exp(exp_bits, mctx);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

    expq = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    /* quick check for easy case of inexact division of leading monomials */
    if (A->bits == B->bits && N == 1 && A->exps[0] < B->exps[0])
    {
        GR_IGNORE(gr_mpoly_zero(Q, ctx));
        TMP_END;
        return GR_DOMAIN;
    }

    /* ensure input exponents packed to same size as output exponents */
    if (exp_bits > A->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, A->exps, A->bits, A->length, mctx);
    }

    if (exp_bits > B->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(exp3, exp_bits, B->exps, B->bits, B->length, mctx);
    }

    /* check leading monomial divides exactly */
    mask = 0;
    if (exp_bits <= FLINT_BITS)
    {
        flint_bitcnt_t jj;
        for (jj = 0; jj < FLINT_BITS/exp_bits; jj++)
            mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));

        if (!mpoly_monomial_divides(expq, exp2, exp3, N, mask))
            easy_exit = 1;
    }
    else
    {
        if (!mpoly_monomial_divides_mp(expq, exp2, exp3, N, exp_bits))
            easy_exit = 1;
    }

    if (easy_exit)
    {
        GR_IGNORE(gr_mpoly_zero(Q, ctx));
        if (free2) flint_free(exp2);
        if (free3) flint_free(exp3);
        TMP_END;
        return GR_DOMAIN;
    }

    /* deal with aliasing and divide polynomials */
    if (Q == A || Q == B)
    {
        gr_mpoly_t T;
        gr_mpoly_init3(T, A->length/B->length + 1, exp_bits, ctx);

        status = _gr_mpoly_divides_heap(&lenq,
                            &T->coeffs, &T->exps, &T->coeffs_alloc, &T->exps_alloc,
                            A->coeffs, exp2, A->length,
                            B->coeffs, exp3, B->length, exp_bits, N, cmpmask, ctx);

        _gr_mpoly_set_length(T, lenq, ctx);
        gr_mpoly_swap(T, Q, ctx);
        gr_mpoly_clear(T, ctx);
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(Q, A->length/B->length + 1, exp_bits, ctx);

        status = _gr_mpoly_divides_heap(&lenq,
                            &Q->coeffs, &Q->exps, &Q->coeffs_alloc, &Q->exps_alloc,
                            A->coeffs, exp2, A->length,
                            B->coeffs, exp3, B->length, exp_bits, N, cmpmask, ctx);

        _gr_mpoly_set_length(Q, lenq, ctx);
    }

    if (free2)
        flint_free(exp2);

    if (free3)
        flint_free(exp3);

    TMP_END;

    return status;
}
