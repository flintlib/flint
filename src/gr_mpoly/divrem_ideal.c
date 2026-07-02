/*
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <string.h>
#include "mpoly.h"
#include "gr_vec.h"
#include "gr_generic.h"
#include "gr_mpoly.h"

/*
    GR port of fmpz_mpoly_divrem_ideal_monagan_pearce.

    Divide A by a list of divisors B[0], ..., B[len-1], producing quotients
    Q[0], ..., Q[len-1] and a remainder R with

        A = Q[0] B[0] + ... + Q[len-1] B[len-1] + R.

    Processing terms from largest to smallest, a leading term is reduced by the
    first divisor whose leading monomial divides it (and, field-like, whose
    leading coefficient divides the term coefficient exactly).  A leading term
    that no divisor can reduce is moved to R (term-to-R normal form), so the
    routine returns GR_SUCCESS rather than GR_DOMAIN and always yields the
    identity above.  As with all multivariate division by a list, R depends on
    the order of the divisors.

    Field-like (nonfield == 0): a divisor reduces a term only if its leading
    coefficient divides the coefficient exactly.

    Non-field (nonfield == 1): the coefficient is reduced with remainder by
    gr_euclidean_divrem (the analogue of fmpz_fdiv_qr), so over an integral
    domain such as Z every term whose monomial is divisible contributes a
    quotient term and its euclidean excess goes to R.
*/

/* field-like exact quotient coefficient = val / lc(B[w]) */
/* shallow copy helpers (see mul_johnson.c / divrem_monagan_pearce.c) */
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

/* gather one heap node: record (i, j, p) for later heap advancement, and for a
   product node place (divisor coeff, quotient coeff) into the shallow dot
   vectors; the single dividend node (if any) becomes the dot initial value */
#define IDEAL_FETCH_NODE(SET_SHALLOW) \
    do { \
        store[store_len++] = x->i; \
        store[store_len++] = x->j; \
        store[store_len++] = x->p; \
        if (x->i == -UWORD(1)) \
        { \
            have_dividend = 1; \
            dividend_j = x->j; \
        } \
        else \
        { \
            hinds[x->p][x->i] |= WORD(1); \
            SET_SHALLOW(dot_a, dot_len, poly3[x->p]->coeffs, x->i); \
            SET_SHALLOW(dot_b, dot_len, Q[x->p]->coeffs, x->j); \
            dot_len++; \
        } \
    } while (0)

/* field-like exact quotient coefficient = val / lc(B[w]) */
#define IDEAL_QCOEFF(dst, w, val) \
    (lc_is_one[w]  ? gr_set(dst, val, cctx) \
     : lc_is_unit[w] ? gr_mul(dst, val, GR_ENTRY(lc_inv, w, sz), cctx) \
                     : gr_div(dst, val, GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx))


/* single-word exponent version */
static int _gr_mpoly_divrem_ideal_mp1(
    gr_mpoly_struct ** Q, gr_mpoly_t R, int * overflowed, int nonfield,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_mpoly_struct * const * poly3, ulong * const * exp3, slong len,
    flint_bitcnt_t bits, ulong maskhi,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    slong i, j, p, w, r_len = 0;
    slong len3 = 0;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_nheap_t ** chains, * chains_ptr;
    slong ** hinds, * hinds_ptr;
    mpoly_nheap_t * x;
    gr_ptr r_coeff = R->coeffs;
    ulong * r_exp = R->exps;
    ulong exp, texp;
    ulong mask;
    slong * store, * store_base, store_len;
    slong * k, * s;
    gr_ptr acc, pp, rem, lc_inv, dot_a, dot_b;
    slong dot_len;
    int have_dividend, cstatus;
    slong dividend_j = 0;
    int * lc_is_one, * lc_is_unit;
    int status = GR_SUCCESS;
    TMP_INIT;

    *overflowed = 0;

    TMP_START;

    GR_TMP_INIT3(acc, pp, rem, cctx);

    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));
    for (w = 0; w < len; w++)
        len3 += poly3[w]->length;
    chains_ptr = (mpoly_nheap_t *) TMP_ALLOC(len3*sizeof(mpoly_nheap_t));
    hinds_ptr = (slong *) TMP_ALLOC(len3*sizeof(slong));

    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = chains_ptr + len3;
        hinds[w] = hinds_ptr + len3;
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    dot_a = flint_malloc(2 * len3 * sz);
    dot_b = GR_ENTRY(dot_a, len3, sz);

    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    for (w = 0; w < len; w++)
    {
        k[w] = -WORD(1);
        s[w] = poly3[w]->length;
    }

    lc_inv = flint_malloc(len * sz);
    _gr_vec_init(lc_inv, len, cctx);
    lc_is_one = (int *) TMP_ALLOC(len*sizeof(int));
    lc_is_unit = (int *) TMP_ALLOC(len*sizeof(int));
    for (w = 0; w < len; w++)
    {
        lc_is_one[w] = lc_is_unit[w] = 0;
        if (!nonfield)
        {
            lc_is_one[w] = (gr_is_one(GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx) == T_TRUE);
            lc_is_unit[w] = lc_is_one[w] ||
                (gr_inv(GR_ENTRY(lc_inv, w, sz), GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx) == GR_SUCCESS);
        }
    }

    mask = mpoly_overflow_mask_sp(bits);

    x = chains[0] + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
        {
            *overflowed = 1;
            goto cleanup;
        }

        store_len = 0;
        dot_len = 0;
        have_dividend = 0;

        if (have_fast_dot)
        {
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    if (sz == 1)       { IDEAL_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { IDEAL_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { IDEAL_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { IDEAL_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { IDEAL_FETCH_NODE(SET_SHALLOW16); }
                    else               { IDEAL_FETCH_NODE(SET_SHALLOW_GENERIC); }
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
                    store[store_len++] = x->p;
                    if (x->i == -UWORD(1))
                        status |= gr_add(acc, acc, GR_ENTRY(coeff2, x->j, sz), cctx);
                    else
                    {
                        hinds[x->p][x->i] |= WORD(1);
                        status |= gr_mul(pp, GR_ENTRY(poly3[x->p]->coeffs, x->i, sz),
                                             GR_ENTRY(Q[x->p]->coeffs, x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        /* process popped nodes */
        while (store_len > 0)
        {
            p = store_base[--store_len];
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                if (j + 1 < len2)
                {
                    x = chains[0] + 0;
                    x->i = -UWORD(1);
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, exp2[x->j], x,
                                             &next_loc, &heap_len, maskhi);
                }
            }
            else
            {
                if ((i + 1 < poly3[p]->length) && (hinds[p][i + 1] == 2*j + 1))
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[p][x->i] + Q[p]->exps[x->j], x,
                                             &next_loc, &heap_len, maskhi);
                }
                if (j == k[p])
                {
                    s[p]++;
                }
                else if (((hinds[p][i] & 1) == 1) &&
                         ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[p][x->i] + Q[p]->exps[x->j], x,
                                             &next_loc, &heap_len, maskhi);
                }
            }
        }

        switch (gr_is_zero(acc, cctx))
        {
            case T_TRUE:
                continue;
            case T_FALSE:
                break;
            default:
                status |= GR_UNABLE;
                goto cleanup;
        }

        /* try to reduce the term acc*x^exp by the divisors in order */
        {
            int div_flag = 0;

            for (w = 0; w < len; w++)
            {
                int d1, d2;

                d1 = mpoly_monomial_divides1(&texp, exp, exp3[w][0], mask);

                if (!d1)
                    continue;

                _gr_mpoly_fit_length(&Q[w]->coeffs, &Q[w]->coeffs_alloc,
                                     &Q[w]->exps, &Q[w]->exps_alloc, 1, k[w] + 2, ctx);

                if (nonfield)
                {
                    cstatus = gr_euclidean_divrem(GR_ENTRY(Q[w]->coeffs, k[w] + 1, sz),
                                     rem, acc, GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx);
                    if (cstatus != GR_SUCCESS) { status |= cstatus; goto cleanup; }
                    gr_swap(acc, rem, cctx);
                    d2 = (gr_is_zero(GR_ENTRY(Q[w]->coeffs, k[w] + 1, sz), cctx) != T_TRUE);
                }
                else
                {
                    cstatus = IDEAL_QCOEFF(GR_ENTRY(Q[w]->coeffs, k[w] + 1, sz), w, acc);
                    if (cstatus == GR_DOMAIN)
                        continue;
                    if (cstatus != GR_SUCCESS) { status |= cstatus; goto cleanup; }
                    d2 = (gr_is_zero(GR_ENTRY(Q[w]->coeffs, k[w] + 1, sz), cctx) != T_TRUE);
                }

                if (d2)
                {
                    k[w]++;
                    Q[w]->exps[k[w]] = texp;

                    if (s[w] > 1)
                    {
                        x = chains[w] + 1;
                        x->i = 1;
                        x->j = k[w];
                        x->p = w;
                        x->next = NULL;
                        hinds[w][x->i] = 2*(x->j + 1) + 0;
                        _mpoly_heap_insert1(heap, exp3[w][1] + Q[w]->exps[k[w]], x,
                                                 &next_loc, &heap_len, maskhi);
                    }
                    s[w] = 1;
                }

                if (nonfield)
                {
                    if (gr_is_zero(acc, cctx) == T_TRUE)
                    {
                        div_flag = 1;
                        break;
                    }
                }
                else
                {
                    div_flag = 1;
                    break;
                }
            }

            if (!div_flag)
            {
                _gr_mpoly_fit_length(&r_coeff, &R->coeffs_alloc, &r_exp, &R->exps_alloc, 1, r_len + 1, ctx);
                status |= gr_set(GR_ENTRY(r_coeff, r_len, sz), acc, cctx);
                r_exp[r_len] = exp;
                r_len++;
            }
        }
    }

cleanup:

    GR_TMP_CLEAR3(acc, pp, rem, cctx);
    flint_free(dot_a);
    _gr_vec_clear(lc_inv, len, cctx);
    flint_free(lc_inv);

    R->coeffs = r_coeff;
    R->exps = r_exp;

    if (*overflowed || status != GR_SUCCESS)
    {
        for (w = 0; w < len; w++)
            Q[w]->length = 0;
        R->length = 0;
    }
    else
    {
        for (w = 0; w < len; w++)
            Q[w]->length = k[w] + 1;
        R->length = r_len;
    }

    TMP_END;

    return (*overflowed) ? GR_SUCCESS : status;
}


static int _gr_mpoly_divrem_ideal_mp(
    gr_mpoly_struct ** Q, gr_mpoly_t R, int * overflowed, int nonfield,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_mpoly_struct * const * poly3, ulong * const * exp3, slong len,
    flint_bitcnt_t bits, slong N, const ulong * cmpmask,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int have_fast_dot;
    slong i, j, p, w, r_len = 0;
    slong len3 = 0;
    slong next_loc, heap_len = 2;
    mpoly_heap_s * heap;
    mpoly_nheap_t ** chains, * chains_ptr;
    slong ** hinds, * hinds_ptr;
    mpoly_nheap_t * x;
    gr_ptr r_coeff = R->coeffs;
    ulong * r_exp = R->exps;
    ulong * exp, * exps, * texp;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * store, * store_base, store_len;
    slong * k, * s;
    gr_ptr acc, pp, rem, lc_inv, dot_a, dot_b;
    slong dot_len;
    int have_dividend, cstatus;
    slong dividend_j = 0;
    int * lc_is_one, * lc_is_unit;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (N == 1)
        return _gr_mpoly_divrem_ideal_mp1(Q, R, overflowed, nonfield,
                   coeff2, exp2, len2, poly3, exp3, len, bits, cmpmask[0], ctx);

    have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);

    *overflowed = 0;

    TMP_START;

    GR_TMP_INIT3(acc, pp, rem, cctx);

    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));
    for (w = 0; w < len; w++)
        len3 += poly3[w]->length;
    chains_ptr = (mpoly_nheap_t *) TMP_ALLOC(len3*sizeof(mpoly_nheap_t));
    hinds_ptr = (slong *) TMP_ALLOC(len3*sizeof(slong));

    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = chains_ptr + len3;
        hinds[w] = hinds_ptr + len3;
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    dot_a = flint_malloc(2 * len3 * sz);
    dot_b = GR_ENTRY(dot_a, len3, sz);

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    for (w = 0; w < len; w++)
    {
        k[w] = -WORD(1);
        s[w] = poly3[w]->length;
    }

    lc_inv = flint_malloc(len * sz);
    _gr_vec_init(lc_inv, len, cctx);
    lc_is_one = (int *) TMP_ALLOC(len*sizeof(int));
    lc_is_unit = (int *) TMP_ALLOC(len*sizeof(int));
    for (w = 0; w < len; w++)
    {
        lc_is_one[w] = lc_is_unit[w] = 0;
        if (!nonfield)
        {
            lc_is_one[w] = (gr_is_one(GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx) == T_TRUE);
            lc_is_unit[w] = lc_is_one[w] ||
                (gr_inv(GR_ENTRY(lc_inv, w, sz), GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx) == GR_SUCCESS);
        }
    }

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    x = chains[0] + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
            {
                *overflowed = 1;
                goto cleanup;
            }
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
            {
                *overflowed = 1;
                goto cleanup;
            }
        }

        store_len = 0;
        dot_len = 0;
        have_dividend = 0;

        if (have_fast_dot)
        {
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    if (sz == 1)       { IDEAL_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { IDEAL_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { IDEAL_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { IDEAL_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { IDEAL_FETCH_NODE(SET_SHALLOW16); }
                    else               { IDEAL_FETCH_NODE(SET_SHALLOW_GENERIC); }
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
                    store[store_len++] = x->p;
                    if (x->i == -UWORD(1))
                        status |= gr_add(acc, acc, GR_ENTRY(coeff2, x->j, sz), cctx);
                    else
                    {
                        hinds[x->p][x->i] |= WORD(1);
                        status |= gr_mul(pp, GR_ENTRY(poly3[x->p]->coeffs, x->i, sz),
                                             GR_ENTRY(Q[x->p]->coeffs, x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        /* process popped nodes: advance in dividend and each quotient*divisor */
        while (store_len > 0)
        {
            p = store_base[--store_len];
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                if (j + 1 < len2)
                {
                    x = chains[0] + 0;
                    x->i = -UWORD(1);
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
            {
                if ((i + 1 < poly3[p]->length) && (hinds[p][i + 1] == 2*j + 1))
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3[p] + x->i*N,
                                                            Q[p]->exps + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3[p] + x->i*N,
                                                               Q[p]->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                if (j == k[p])
                {
                    s[p]++;
                }
                else if (((hinds[p][i] & 1) == 1) &&
                         ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3[p] + x->i*N,
                                                            Q[p]->exps + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3[p] + x->i*N,
                                                               Q[p]->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        /* if the accumulated coefficient is zero, the term vanishes */
        switch (gr_is_zero(acc, cctx))
        {
            case T_TRUE:
                continue;
            case T_FALSE:
                break;
            default:
                status |= GR_UNABLE;
                goto cleanup;
        }

        /* try to reduce the term acc*x^exp by the divisors in order */
        {
            int div_flag = 0;

            for (w = 0; w < len; w++)
            {
                int d1, d2;

                if (bits <= FLINT_BITS)
                    d1 = mpoly_monomial_divides(texp, exp, exp3[w], N, mask);
                else
                    d1 = mpoly_monomial_divides_mp(texp, exp, exp3[w], N, bits);

                if (!d1)
                    continue;

                _gr_mpoly_fit_length(&Q[w]->coeffs, &Q[w]->coeffs_alloc,
                                     &Q[w]->exps, &Q[w]->exps_alloc, N, k[w] + 2, ctx);

                if (nonfield)
                {
                    cstatus = gr_euclidean_divrem(GR_ENTRY(Q[w]->coeffs, k[w] + 1, sz),
                                     rem, acc, GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx);
                    if (cstatus != GR_SUCCESS) { status |= cstatus; goto cleanup; }
                    gr_swap(acc, rem, cctx);   /* acc <- euclidean remainder */
                    d2 = (gr_is_zero(GR_ENTRY(Q[w]->coeffs, k[w] + 1, sz), cctx) != T_TRUE);
                }
                else
                {
                    cstatus = IDEAL_QCOEFF(GR_ENTRY(Q[w]->coeffs, k[w] + 1, sz), w, acc);
                    if (cstatus == GR_DOMAIN)
                        continue;   /* this divisor cannot divide exactly */
                    if (cstatus != GR_SUCCESS) { status |= cstatus; goto cleanup; }
                    d2 = (gr_is_zero(GR_ENTRY(Q[w]->coeffs, k[w] + 1, sz), cctx) != T_TRUE);
                }

                if (d2)
                {
                    k[w]++;
                    mpoly_monomial_set(Q[w]->exps + k[w]*N, texp, N);

                    if (s[w] > 1)
                    {
                        x = chains[w] + 1;
                        x->i = 1;
                        x->j = k[w];
                        x->p = w;
                        x->next = NULL;
                        hinds[w][x->i] = 2*(x->j + 1) + 0;
                        if (bits <= FLINT_BITS)
                            mpoly_monomial_add(exp_list[exp_next], exp3[w] + N,
                                                                Q[w]->exps + k[w]*N, N);
                        else
                            mpoly_monomial_add_mp(exp_list[exp_next], exp3[w] + N,
                                                                Q[w]->exps + k[w]*N, N);
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                                 &next_loc, &heap_len, N, cmpmask);
                    }
                    s[w] = 1;
                }

                if (nonfield)
                {
                    /* keep reducing the euclidean remainder by later divisors */
                    if (gr_is_zero(acc, cctx) == T_TRUE)
                    {
                        div_flag = 1;
                        break;
                    }
                }
                else
                {
                    div_flag = 1;   /* exact division fully reduces the term */
                    break;
                }
            }

            if (!div_flag)
            {
                /* term is irreducible: send it to the remainder */
                _gr_mpoly_fit_length(&r_coeff, &R->coeffs_alloc, &r_exp, &R->exps_alloc, N, r_len + 1, ctx);
                status |= gr_set(GR_ENTRY(r_coeff, r_len, sz), acc, cctx);
                mpoly_monomial_set(r_exp + r_len*N, exp, N);
                r_len++;
            }
        }
    }

cleanup:

    GR_TMP_CLEAR3(acc, pp, rem, cctx);
    flint_free(dot_a);
    _gr_vec_clear(lc_inv, len, cctx);
    flint_free(lc_inv);

    R->coeffs = r_coeff;
    R->exps = r_exp;

    if (*overflowed || status != GR_SUCCESS)
    {
        for (w = 0; w < len; w++)
            Q[w]->length = 0;
        R->length = 0;
    }
    else
    {
        for (w = 0; w < len; w++)
            Q[w]->length = k[w] + 1;
        R->length = r_len;
    }

    TMP_END;

    return (*overflowed) ? GR_SUCCESS : status;
}


static int _gr_mpoly_divrem_ideal(
    gr_mpoly_struct ** Q, gr_mpoly_t R,
    const gr_mpoly_t A, gr_mpoly_struct * const * B, slong len,
    int nonfield, gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong i, N, len3 = 0;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * exp2;
    ulong ** exp3;
    int free2 = 0, * free3, overflowed;
    gr_mpoly_t TR;
    gr_mpoly_struct * r;
    int status = GR_SUCCESS;
    TMP_INIT;

    for (i = 0; i < len; i++)
    {
        if (B[i]->length == 0)
            return GR_DOMAIN;   /* division by zero */
        if (gr_is_zero(B[i]->coeffs, cctx) != T_FALSE)
            return GR_UNABLE;
        len3 = FLINT_MAX(len3, B[i]->length);
    }

    if (A->length == 0)
    {
        for (i = 0; i < len; i++)
            status |= gr_mpoly_zero(Q[i], ctx);
        status |= gr_mpoly_zero(R, ctx);
        return status;
    }

    TMP_START;

    free3 = (int *) TMP_ALLOC(len*sizeof(int));
    exp3 = (ulong **) TMP_ALLOC(len*sizeof(ulong *));

    exp_bits = A->bits;
    for (i = 0; i < len; i++)
        exp_bits = FLINT_MAX(exp_bits, B[i]->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    N = mpoly_words_per_exp(exp_bits, mctx);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

    exp2 = A->exps;
    if (exp_bits > A->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, A->exps, A->bits, A->length, mctx);
    }

    for (i = 0; i < len; i++)
    {
        exp3[i] = B[i]->exps;
        free3[i] = 0;
        if (exp_bits > B[i]->bits)
        {
            free3[i] = 1;
            exp3[i] = (ulong *) flint_malloc(N*B[i]->length*sizeof(ulong));
            mpoly_repack_monomials(exp3[i], exp_bits, B[i]->exps, B[i]->bits,
                                                            B[i]->length, mctx);
        }
        gr_mpoly_fit_length_reset_bits(Q[i], 1, exp_bits, ctx);
    }

    /* if lm(A) < lm(B[i]) for all i, the quotients are zero and R = A */
    for (i = 0; i < len; i++)
    {
        if (!mpoly_monomial_lt(exp2, exp3[i], N, cmpmask))
            break;
    }

    if (i == len)
    {
        status |= gr_mpoly_set(R, A, ctx);
        for (i = 0; i < len; i++)
            status |= gr_mpoly_zero(Q[i], ctx);
        goto cleanup;
    }

    /* handle aliasing of R with A */
    if (R == A)
    {
        gr_mpoly_init3(TR, len3, exp_bits, ctx);
        r = TR;
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(R, len3, exp_bits, ctx);
        r = R;
    }

    while (1)
    {
        r->bits = exp_bits;

        status = _gr_mpoly_divrem_ideal_mp(Q, r, &overflowed, nonfield,
                            A->coeffs, exp2, A->length,
                            B, exp3, len, exp_bits, N, cmpmask, ctx);

        if (!overflowed)
            break;

        {
            slong old_exp_bits = exp_bits;
            ulong * old_exp2 = exp2, * old_exp3;

            exp_bits = mpoly_fix_bits(exp_bits + 1, mctx);
            N = mpoly_words_per_exp(exp_bits, mctx);
            cmpmask = (ulong *) flint_realloc(cmpmask, N*sizeof(ulong));
            mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

            exp2 = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
            mpoly_repack_monomials(exp2, exp_bits, old_exp2, old_exp_bits, A->length, mctx);
            if (free2) flint_free(old_exp2);
            free2 = 1;

            gr_mpoly_fit_length_reset_bits(r, len3, exp_bits, ctx);

            for (i = 0; i < len; i++)
            {
                old_exp3 = exp3[i];
                exp3[i] = (ulong *) flint_malloc(N*B[i]->length*sizeof(ulong));
                mpoly_repack_monomials(exp3[i], exp_bits, old_exp3, old_exp_bits,
                                                            B[i]->length, mctx);
                if (free3[i]) flint_free(old_exp3);
                free3[i] = 1;
                gr_mpoly_fit_length_reset_bits(Q[i], 1, exp_bits, ctx);
            }
        }
    }

    if (r == TR)
    {
        gr_mpoly_swap(R, TR, ctx);
        gr_mpoly_clear(TR, ctx);
    }

    if (status != GR_SUCCESS)
    {
        for (i = 0; i < len; i++)
            GR_IGNORE(gr_mpoly_zero(Q[i], ctx));
        GR_IGNORE(gr_mpoly_zero(R, ctx));
    }

cleanup:

    if (free2)
        flint_free(exp2);
    for (i = 0; i < len; i++)
    {
        if (free3[i])
            flint_free(exp3[i]);
    }
    flint_free(cmpmask);

    TMP_END;

    return status;
}


/* build pointer arrays from the quotient/divisor vectors and run the kernel */
static int
_gr_mpoly_divrem_ideal_vec(gr_mpoly_vec_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_vec_t B, int nonfield, gr_mpoly_ctx_t ctx)
{
    slong w, len = B->length;
    gr_mpoly_struct ** Qptr, ** Bptr;
    int status;

    gr_mpoly_vec_set_length(Q, len, ctx);

    Qptr = (gr_mpoly_struct **) flint_malloc(len*sizeof(gr_mpoly_struct *));
    Bptr = (gr_mpoly_struct **) flint_malloc(len*sizeof(gr_mpoly_struct *));
    for (w = 0; w < len; w++)
    {
        Qptr[w] = Q->entries + w;
        Bptr[w] = B->entries + w;
    }

    status = _gr_mpoly_divrem_ideal(Qptr, R, A, Bptr, len, nonfield, ctx);

    flint_free(Qptr);
    flint_free(Bptr);

    return status;
}

int gr_mpoly_divrem_ideal(
    gr_mpoly_vec_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_vec_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_ideal_vec(Q, R, A, B, 0, ctx);
}

int gr_mpoly_divrem_ideal_weak(
    gr_mpoly_vec_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_vec_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_ideal_vec(Q, R, A, B, 1, ctx);
}
