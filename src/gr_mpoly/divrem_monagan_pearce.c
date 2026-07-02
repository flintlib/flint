/*
    Copyright (C) 2017 William Hart
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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"
#include "gr_vec.h"
#include "gr_generic.h"
#include "gr_mpoly.h"

/*
    GR port of fmpz_mpoly_divrem_monagan_pearce.

    Compute Q and R with A = Q*B + R, where no monomial of R is divisible by the
    leading monomial of B.  Following the convention of gr_poly_divrem, the
    coefficient divisions are required to be exact: the leading coefficient of B
    is inverted once if it is a unit (and the multiplication is skipped entirely
    when it is 1), otherwise gr_div is used, which over an integral domain such
    as Z returns GR_DOMAIN exactly when a quotient coefficient is not integral.
    Hence the routine returns GR_SUCCESS over any field, and over Z precisely
    when an exact quotient/remainder pair with the above property exists.

    The accumulation strategy (a single _gr_vec_dot over shallow operand vectors
    when the ring overloads VEC_DOT, otherwise gr_mul into a preallocated
    temporary followed by gr_sub) matches divides_monagan_pearce.
*/

/* shallow copy helpers (see mul_johnson.c / divides_monagan_pearce.c) */
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

#define DIVREM_FETCH_NODE(SET_SHALLOW) \
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
            SET_SHALLOW(dot_b, dot_len, q_coeff, x->j); \
            dot_len++; \
        } \
    } while (0)

/* quotient coefficient = acc / lc(B) (field-like: must be exact) */
#define DIVREM_QCOEFF(dst, acc) \
    (lc_is_one  ? gr_set(dst, acc, cctx) \
     : lc_is_unit ? gr_mul(dst, acc, lc_inv, cctx) \
                  : gr_div(dst, acc, coeff3, cctx))

/* append a remainder term (val) * x^exp; skipped entirely when R is absent
   (write_r == 0), which is how the remainder-free div variants are obtained */
#define ADD_R_TERM1(val) \
    do { \
        if (write_r) \
        { \
            _gr_mpoly_fit_length(&r_coeff, &R->coeffs_alloc, &r_exp, &R->exps_alloc, 1, r_len + 1, ctx); \
            status |= gr_set(GR_ENTRY(r_coeff, r_len, sz), val, cctx); \
            r_exp[r_len] = exp; \
            r_len++; \
        } \
    } while (0)

#define ADD_R_TERMN(val) \
    do { \
        if (write_r) \
        { \
            _gr_mpoly_fit_length(&r_coeff, &R->coeffs_alloc, &r_exp, &R->exps_alloc, N, r_len + 1, ctx); \
            status |= gr_set(GR_ENTRY(r_coeff, r_len, sz), val, cctx); \
            mpoly_monomial_set(r_exp + r_len*N, exp, N); \
            r_len++; \
        } \
    } while (0)

/*
    Decide the fate of the accumulated leading coefficient acc at monomial exp.
    Sets commit = 1 if a quotient term GR_ENTRY(q_coeff, q_len, sz) was produced.

    Field-like (nonfield == 0): the coefficient division must be exact
    (returning GR_DOMAIN otherwise); terms not reducible go to R.

    Non-field (nonfield == 1): the coefficient is divided with remainder via
    gr_euclidean_divrem (the analogue of fmpz_fdiv_qr); the excess remainder is
    added to R, so a monomial divisible by lm(B) may still contribute to R.
*/
#define DIVREM_COEFF_STEP(ADD_R_TERM) \
{ \
    if (!lt_divides) \
    { \
        ADD_R_TERM(acc); \
    } \
    else if (nonfield) \
    { \
        cstatus = gr_euclidean_divrem(GR_ENTRY(q_coeff, q_len, sz), rem, acc, coeff3, cctx); \
        if (cstatus != GR_SUCCESS) { status |= cstatus; goto cleanup; } \
        switch (gr_is_zero(rem, cctx)) \
        { \
            case T_TRUE: break; \
            case T_FALSE: ADD_R_TERM(rem); break; \
            default: status |= GR_UNABLE; goto cleanup; \
        } \
        commit = (gr_is_zero(GR_ENTRY(q_coeff, q_len, sz), cctx) != T_TRUE); \
    } \
    else \
    { \
        cstatus = DIVREM_QCOEFF(GR_ENTRY(q_coeff, q_len, sz), acc); \
        if (cstatus != GR_SUCCESS) { status |= cstatus; goto cleanup; } \
        if (gr_is_zero(GR_ENTRY(q_coeff, q_len, sz), cctx) == T_TRUE) \
            ADD_R_TERM(acc); \
        else \
            commit = 1; \
    } \
}


static int _gr_mpoly_divrem_monagan_pearce1(
    gr_mpoly_t Q, gr_mpoly_t R, int * overflowed, int nonfield,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_srcptr coeff3, const ulong * exp3, slong len3,
    flint_bitcnt_t bits, ulong maskhi,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    slong i, j, q_len, r_len, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    int write_r = (R != NULL);
    gr_ptr q_coeff = Q->coeffs;
    gr_ptr r_coeff = write_r ? R->coeffs : NULL;
    ulong * q_exp = Q->exps;
    ulong * r_exp = write_r ? R->exps : NULL;
    slong * hind;
    ulong mask, exp;
    int lt_divides, lc_is_unit = 0, lc_is_one = 0, cstatus, commit;
    gr_ptr acc, lc_inv, pp, rem, dot_a, dot_b;
    slong dot_len;
    int have_dividend;
    slong dividend_j = 0;
    int status = GR_SUCCESS;
    TMP_INIT;

    *overflowed = 0;

    TMP_START;

    GR_TMP_INIT4(acc, lc_inv, pp, rem, cctx);
    dot_a = flint_malloc(2 * len3 * sz);
    dot_b = GR_ENTRY(dot_a, len3, sz);

    next_loc = len3 + 4;
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    q_len = WORD(0);
    r_len = WORD(0);
    s = len3;

    x = chain + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    if (!nonfield)
    {
        lc_is_one = (gr_is_one(coeff3, cctx) == T_TRUE);
        lc_is_unit = lc_is_one || (gr_inv(lc_inv, coeff3, cctx) == GR_SUCCESS);
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
        {
            *overflowed = 1;
            goto cleanup;
        }

        _gr_mpoly_fit_length(&q_coeff, &Q->coeffs_alloc, &q_exp, &Q->exps_alloc, 1, q_len + 1, ctx);

        lt_divides = mpoly_monomial_divides1(q_exp + q_len, exp, exp3[0], mask);

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
                    if (sz == 1)       { DIVREM_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { DIVREM_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { DIVREM_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { DIVREM_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { DIVREM_FETCH_NODE(SET_SHALLOW16); }
                    else               { DIVREM_FETCH_NODE(SET_SHALLOW_GENERIC); }
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
                        status |= gr_mul(pp, GR_ENTRY(coeff3, x->i, sz), GR_ENTRY(q_coeff, x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
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
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = -UWORD(1);
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
            else
            {
                if ((i + 1 < len3) && (hind[i + 1] == 2*j + 1))
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[x->i] + q_exp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
                if (j + 1 == q_len)
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
                    _mpoly_heap_insert1(heap, exp3[x->i] + q_exp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        /* if accumulated coefficient is zero, this term vanishes */
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

        commit = 0;
        DIVREM_COEFF_STEP(ADD_R_TERM1);

        if (commit)
        {
            /* put the newly generated quotient term back into the heap */
            if (s > 1)
            {
                i = 1;
                x = chain + i;
                x->i = i;
                x->j = q_len;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, exp3[x->i] + q_exp[x->j], x,
                                                     &next_loc, &heap_len, maskhi);
            }
            s = 1;
            q_len++;
        }
    }

cleanup:

    GR_TMP_CLEAR4(acc, lc_inv, pp, rem, cctx);
    flint_free(dot_a);

    Q->coeffs = q_coeff;
    Q->exps = q_exp;
    if (write_r)
    {
        R->coeffs = r_coeff;
        R->exps = r_exp;
    }

    if (*overflowed || (status != GR_SUCCESS))
    {
        Q->length = 0;
        if (write_r)
            R->length = 0;
    }
    else
    {
        Q->length = q_len;
        if (write_r)
            R->length = r_len;
    }

    TMP_END;

    return (*overflowed) ? GR_SUCCESS : status;
}


static int _gr_mpoly_divrem_monagan_pearce(
    gr_mpoly_t Q, gr_mpoly_t R, int * overflowed, int nonfield,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_srcptr coeff3, const ulong * exp3, slong len3,
    flint_bitcnt_t bits, slong N, const ulong * cmpmask,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int have_fast_dot;
    slong i, j, q_len, r_len, s;
    slong next_loc, heap_len = 2;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    int write_r = (R != NULL);
    gr_ptr q_coeff = Q->coeffs;
    gr_ptr r_coeff = write_r ? R->coeffs : NULL;
    ulong * q_exp = Q->exps;
    ulong * r_exp = write_r ? R->exps : NULL;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    int lt_divides, lc_is_unit = 0, lc_is_one = 0, cstatus, commit;
    gr_ptr acc, lc_inv, pp, rem, dot_a, dot_b;
    slong dot_len;
    int have_dividend;
    slong dividend_j = 0;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (N == 1)
        return _gr_mpoly_divrem_monagan_pearce1(Q, R, overflowed, nonfield,
                   coeff2, exp2, len2, coeff3, exp3, len3, bits, cmpmask[0], ctx);

    have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);

    *overflowed = 0;

    TMP_START;

    GR_TMP_INIT4(acc, lc_inv, pp, rem, cctx);
    dot_a = flint_malloc(2 * len3 * sz);
    dot_b = GR_ENTRY(dot_a, len3, sz);

    next_loc = len3 + 4;
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

    q_len = WORD(0);
    r_len = WORD(0);
    s = len3;

    x = chain + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    if (!nonfield)
    {
        lc_is_one = (gr_is_one(coeff3, cctx) == T_TRUE);
        lc_is_unit = lc_is_one || (gr_inv(lc_inv, coeff3, cctx) == GR_SUCCESS);
    }

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

        _gr_mpoly_fit_length(&q_coeff, &Q->coeffs_alloc, &q_exp, &Q->exps_alloc, N, q_len + 1, ctx);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(q_exp + q_len*N, exp, exp3, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N, exp, exp3, N, bits);

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
                    if (sz == 1)       { DIVREM_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { DIVREM_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { DIVREM_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { DIVREM_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { DIVREM_FETCH_NODE(SET_SHALLOW16); }
                    else               { DIVREM_FETCH_NODE(SET_SHALLOW_GENERIC); }
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
                        status |= gr_mul(pp, GR_ENTRY(coeff3, x->i, sz), GR_ENTRY(q_coeff, x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        while (store_len > 0)
        {
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = -UWORD(1);
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
            {
                if ((i + 1 < len3) && (hind[i + 1] == 2*j + 1))
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               q_exp + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               q_exp + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                if (j + 1 == q_len)
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
                                                               q_exp + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               q_exp + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
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

        commit = 0;
        DIVREM_COEFF_STEP(ADD_R_TERMN);

        if (commit)
        {
            if (s > 1)
            {
                i = 1;
                x = chain + i;
                x->i = i;
                x->j = q_len;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;
                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                           q_exp + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                           q_exp + x->j*N, N);
                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                                 &next_loc, &heap_len, N, cmpmask);
            }
            s = 1;
            q_len++;
        }
    }

cleanup:

    GR_TMP_CLEAR4(acc, lc_inv, pp, rem, cctx);
    flint_free(dot_a);

    Q->coeffs = q_coeff;
    Q->exps = q_exp;
    if (write_r)
    {
        R->coeffs = r_coeff;
        R->exps = r_exp;
    }

    if (*overflowed || (status != GR_SUCCESS))
    {
        Q->length = 0;
        if (write_r)
            R->length = 0;
    }
    else
    {
        Q->length = q_len;
        if (write_r)
            R->length = r_len;
    }

    TMP_END;

    return (*overflowed) ? GR_SUCCESS : status;
}


static int _gr_mpoly_divrem_mp(
    gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B, int nonfield,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * exp2 = A->exps, * exp3 = B->exps;
    int free2 = 0, free3 = 0, overflowed;
    gr_mpoly_t TQ, TR;
    gr_mpoly_struct * q, * r;
    int status;

    if (B->length == 0)
        return GR_DOMAIN;

    if (A->length == 0)
    {
        status = gr_mpoly_zero(Q, ctx);
        if (R != NULL)
            status |= gr_mpoly_zero(R, ctx);
        return status;
    }

    if (gr_is_zero(B->coeffs, cctx) != T_FALSE)
        return GR_UNABLE;

    exp_bits = FLINT_MAX(A->bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    N = mpoly_words_per_exp(exp_bits, mctx);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

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

    /* use temporaries to handle aliasing of Q, R with A, B */
    if (Q == A || Q == B)
    {
        gr_mpoly_init3(TQ, A->length/B->length + 1, exp_bits, ctx);
        q = TQ;
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(Q, A->length/B->length + 1, exp_bits, ctx);
        q = Q;
    }

    if (R == NULL)
    {
        r = NULL;   /* remainder-free: div variants */
    }
    else if (R == A || R == B)
    {
        gr_mpoly_init3(TR, B->length, exp_bits, ctx);
        r = TR;
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(R, B->length, exp_bits, ctx);
        r = R;
    }

    while (1)
    {
        q->bits = exp_bits;
        if (r != NULL)
            r->bits = exp_bits;

        status = _gr_mpoly_divrem_monagan_pearce(q, r, &overflowed, nonfield,
                            A->coeffs, exp2, A->length,
                            B->coeffs, exp3, B->length, exp_bits, N, cmpmask, ctx);

        if (!overflowed)
            break;

        /* retry with larger exponent fields */
        {
            slong old_exp_bits = exp_bits;
            ulong * old_exp2 = exp2;
            ulong * old_exp3 = exp3;

            exp_bits = mpoly_fix_bits(exp_bits + 1, mctx);
            N = mpoly_words_per_exp(exp_bits, mctx);
            cmpmask = (ulong *) flint_realloc(cmpmask, N*sizeof(ulong));
            mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

            exp2 = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
            mpoly_repack_monomials(exp2, exp_bits, old_exp2, old_exp_bits, A->length, mctx);
            if (free2)
                flint_free(old_exp2);
            free2 = 1;

            exp3 = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
            mpoly_repack_monomials(exp3, exp_bits, old_exp3, old_exp_bits, B->length, mctx);
            if (free3)
                flint_free(old_exp3);
            free3 = 1;

            gr_mpoly_fit_length_reset_bits(q, A->length/B->length + 1, exp_bits, ctx);
            if (r != NULL)
                gr_mpoly_fit_length_reset_bits(r, B->length, exp_bits, ctx);
        }
    }

    if (q == TQ)
    {
        gr_mpoly_swap(Q, TQ, ctx);
        gr_mpoly_clear(TQ, ctx);
    }
    if (r == TR && R != NULL)
    {
        gr_mpoly_swap(R, TR, ctx);
        gr_mpoly_clear(TR, ctx);
    }

    if (status != GR_SUCCESS)
    {
        gr_mpoly_zero(Q, ctx);
        if (R != NULL)
            gr_mpoly_zero(R, ctx);
    }

    flint_free(cmpmask);
    if (free2)
        flint_free(exp2);
    if (free3)
        flint_free(exp3);

    return status;
}


int gr_mpoly_divrem_monagan_pearce(
    gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_mp(Q, R, A, B, 0, ctx);
}


int gr_mpoly_divrem(
    gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_mp(Q, R, A, B, 0, ctx);
}


int gr_mpoly_divrem_allowing_nonunit_lc(
    gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_mp(Q, R, A, B, 1, ctx);
}


int gr_mpoly_div(
    gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    /* remainder-free: the kernel skips all remainder writes when R is NULL */
    return _gr_mpoly_divrem_mp(Q, NULL, A, B, 0, ctx);
}


int gr_mpoly_div_allowing_nonunit_lc(
    gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_mp(Q, NULL, A, B, 1, ctx);
}
