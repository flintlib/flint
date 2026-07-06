/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdint.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"
#include "gr_generic.h"
#include "gr_mpoly.h"

/*
    Squaring of a sparse multivariate polynomial using a binary heap, in the
    spirit of Johnson's algorithm (see gr_mpoly_mul_heap).  A commutative
    coefficient ring is assumed, so that a_i*a_j = a_j*a_i.

    Writing f = sum_i a_i X^{e_i} (with the e_i distinct), we have

        f^2 = sum_i a_i^2 X^{2 e_i}  +  2 sum_{i<j} a_i a_j X^{e_i + e_j}.

    The heap is therefore restricted to the pairs (i,j) with i <= j: stream i
    runs over columns j = i, i+1, ..., len-1.  This generates each unordered
    pair exactly once, so only ~ len*(len+1)/2 coefficient multiplications are
    performed, about half of the len^2 done by a plain multiplication.

    For each output monomial we accumulate the off-diagonal products (i<j),
    double the result, and add the diagonal square a_d^2 whenever a diagonal
    pair (d,d) lands on the same monomial (which happens exactly when
    2 e_d = e_p + e_q for some p < q).
*/

/* Shallow-copy helpers, identical to gr_mpoly_mul_heap: place operands in
   temporary arrays so the off-diagonal part can be accumulated with a single
   _gr_vec_dot call. */

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

/* Fetch one chain of terms sharing the current exponent, sorting each term
   into either the off-diagonal dot vectors or the (unique) diagonal slot. */

#define SQR_FETCH_TERMS_INNER(SET_SHALLOW) \
    do { \
        slong _ii = x->i, _jj = x->j; \
        hind[_ii] |= WORD(1); \
        Q[Q_len++] = _ii; \
        Q[Q_len++] = _jj; \
        if (_ii == _jj) \
        { \
            diag_i = _ii; \
            have_diag = 1; \
        } \
        else \
        { \
            SET_SHALLOW(dot_a, dot_len, coeff2, _ii); \
            SET_SHALLOW(dot_b, dot_len, coeff2, _jj); \
            dot_len++; \
        } \
    } while (0); \
    while ((x = x->next) != NULL) \
    { \
        slong _ii = x->i, _jj = x->j; \
        hind[_ii] |= WORD(1); \
        Q[Q_len++] = _ii; \
        Q[Q_len++] = _jj; \
        if (_ii == _jj) \
        { \
            diag_i = _ii; \
            have_diag = 1; \
        } \
        else \
        { \
            SET_SHALLOW(dot_a, dot_len, coeff2, _ii); \
            SET_SHALLOW(dot_b, dot_len, coeff2, _jj); \
            dot_len++; \
        } \
    }

#define SQR_FETCH_TERMS1(SET_SHALLOW) \
    do { \
        while (heap_len > 1 && heap[1].exp == exp) \
        { \
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi); \
            SQR_FETCH_TERMS_INNER(SET_SHALLOW) \
        } \
    } while (0)

#define SQR_FETCH_TERMS(SET_SHALLOW) \
    do \
    { \
        exp_list[--exp_next] = heap[1].exp; \
        x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask); \
        SQR_FETCH_TERMS_INNER(SET_SHALLOW) \
    } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))

/* Combine accumulated off-diagonal dot product D with the optional diagonal
   square: result = 2*D + (a_diag)^2.  Writes into GR_ENTRY(p1, len1, sz). */

#define SQR_FINALIZE_DOT \
    do { \
        gr_ptr _out = GR_ENTRY(p1, len1, sz); \
        if (dot_len > 0) \
        { \
            status |= _gr_vec_dot(_out, NULL, 0, dot_a, dot_b, dot_len, cctx); \
            status |= gr_mul_two(_out, _out, cctx); \
            if (have_diag) \
            { \
                status |= gr_sqr(t, GR_ENTRY(coeff2, diag_i, sz), cctx); \
                status |= gr_add(_out, _out, t, cctx); \
            } \
        } \
        else \
        { \
            /* pure diagonal monomial 2*e_diag */ \
            status |= gr_sqr(_out, GR_ENTRY(coeff2, diag_i, sz), cctx); \
        } \
    } while (0)

/* Generic accumulation path (no fast dot product): off-diagonal products are
   summed with gr_mul / gr_addmul, doubled once, and the diagonal square added.
   Multiplication count is still ~ len*(len+1)/2. */

#define SQR_DOT_TERMS_GENERIC \
    do { \
        do { \
            slong _ii = x->i, _jj = x->j; \
            hind[_ii] |= WORD(1); \
            Q[Q_len++] = _ii; \
            Q[Q_len++] = _jj; \
            if (_ii == _jj) \
            { \
                status |= gr_sqr(t, GR_ENTRY(coeff2, _ii, sz), cctx); \
                have_diag = 1; \
            } \
            else if (first_off) \
            { \
                status |= gr_mul(GR_ENTRY(p1, len1, sz), GR_ENTRY(coeff2, _ii, sz), GR_ENTRY(coeff2, _jj, sz), cctx); \
                first_off = 0; \
            } \
            else \
            { \
                status |= gr_addmul(GR_ENTRY(p1, len1, sz), GR_ENTRY(coeff2, _ii, sz), GR_ENTRY(coeff2, _jj, sz), cctx); \
            } \
        } while ((x = x->next) != NULL); \
    } while (0)

#define SQR_FINALIZE_GENERIC \
    do { \
        gr_ptr _out = GR_ENTRY(p1, len1, sz); \
        if (!first_off) \
        { \
            status |= gr_mul_two(_out, _out, cctx); \
            if (have_diag) \
                status |= gr_add(_out, _out, t, cctx); \
        } \
        else \
        { \
            /* only a diagonal term contributed */ \
            status |= gr_set(_out, t, cctx); \
        } \
    } while (0)

static int _gr_mpoly_sqr_commutative_heap1(
    slong * res_len,
    gr_ptr * coeff1, ulong ** exp1, slong * alloc, slong * exps_alloc,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    ulong maskhi,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    slong i, j;
    slong next_loc;
    slong Q_len = 0, heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    slong len1;
    gr_ptr p1 = * coeff1;
    ulong * e1 = *exp1;
    ulong exp;
    slong * hind;
    gr_ptr t, dot_a, dot_b;
    slong dot_len, diag_i;
    int have_diag, first_off;
    slong sz = cctx->sizeof_elem;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;
    GR_TMP_INIT(t, cctx);

    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));

    hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
    /* triangular start: stream i's first candidate column is i */
    for (i = 0; i < len2; i++)
        hind[i] = 2*i + 1;

    dot_a = TMP_ALLOC(2 * len2 * sz);
    dot_b = GR_ENTRY(dot_a, len2, sz);

    /* put (0, 0, exp2[0] + exp2[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    HEAP_ASSIGN(heap[1], exp2[0] + exp2[0], x);
    hind[0] = 2*1 + 0;

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;
        _gr_mpoly_fit_length(&p1, alloc, &e1, exps_alloc, 1, len1 + 1, ctx);

        have_diag = 0;
        diag_i = 0;

        if (!have_fast_dot)
        {
            first_off = 1;

            while (heap_len > 1 && heap[1].exp == exp)
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                SQR_DOT_TERMS_GENERIC;
            }

            SQR_FINALIZE_GENERIC;
        }
        else
        {
            dot_len = 0;

            if (sz == 1)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW1);
            }
            else if (sz == 2)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW2);
            }
            else if (sz == 4)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW4);
            }
            else if (sz == 8)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW8);
            }
            else if (sz == 16)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW16);
            }
            else
            {
                SQR_FETCH_TERMS1(SET_SHALLOW_GENERIC);
            }

            SQR_FINALIZE_DOT;
        }

        /* set output monomial */
        e1[len1] = exp;
        len1 += (gr_is_zero(GR_ENTRY(p1, len1, sz), cctx) != T_TRUE);

        /* for each node temporarily stored */
        while (Q_len > 0)
        {
            /* take node from store */
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? spawn stream (i+1) at column j.
               The triangular initialization hind[k] = 2*k+1 makes this fire
               only for j >= i+1, so below-diagonal pairs are never created. */
            if ((i + 1 < len2) && (hind[i + 1] == 2*j + 1))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }

            /* should we go up? advance (i, j+1) within stream i */
            if ((j + 1 < len2) && ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >  2*(j + 2) + 1)
                          || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */))
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }
    }

    (*coeff1) = p1;
    (*exp1) = e1;

    TMP_END;
    GR_TMP_CLEAR(t, cctx);

    *res_len = len1;

    return status;
}

static int _gr_mpoly_sqr_commutative_heap(
    slong * res_len,
    gr_ptr * coeff1, ulong ** exp1, slong * alloc, slong * exps_alloc,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    slong i, j;
    slong next_loc;
    slong Q_len = 0, heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    slong len1;
    gr_ptr p1 = * coeff1;
    ulong * e1 = *exp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    gr_ptr t, dot_a, dot_b;
    slong dot_len, diag_i;
    int have_diag, first_off;
    slong sz = cctx->sizeof_elem;
    int status = GR_SUCCESS;
    TMP_INIT;

    /* if exponent vectors fit in single word, call special version */
    if (N == 1)
        return _gr_mpoly_sqr_commutative_heap1(res_len, coeff1, exp1, alloc, exps_alloc,
            coeff2, exp2, len2, cmpmask[0], ctx);

    TMP_START;
    GR_TMP_INIT(t, cctx);

    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));
    for (i = 0; i < len2; i++)
        exp_list[i] = exps + i*N;

    hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
    /* triangular start: stream i's first candidate column is i */
    for (i = 0; i < len2; i++)
        hind[i] = 2*i + 1;

    dot_a = TMP_ALLOC(2 * len2 * sz);
    dot_b = GR_ENTRY(dot_a, len2, sz);

    /* start with no heap nodes and no exponent vectors in use */
    exp_next = 0;

    /* put (0, 0, exp2[0] + exp2[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    mpoly_monomial_add_any_bits(heap[1].exp, exp2, exp2, N, bits);

    hind[0] = 2*1 + 0;

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _gr_mpoly_fit_length(&p1, alloc, &e1, exps_alloc, N, len1 + 1, ctx);

        mpoly_monomial_set(e1 + len1*N, exp, N);

        have_diag = 0;
        diag_i = 0;

        if (!have_fast_dot)
        {
            first_off = 1;

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                SQR_DOT_TERMS_GENERIC;
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            SQR_FINALIZE_GENERIC;
        }
        else
        {
            dot_len = 0;

            if (sz == 1)
            {
                SQR_FETCH_TERMS(SET_SHALLOW1);
            }
            else if (sz == 2)
            {
                SQR_FETCH_TERMS(SET_SHALLOW2);
            }
            else if (sz == 4)
            {
                SQR_FETCH_TERMS(SET_SHALLOW4);
            }
            else if (sz == 8)
            {
                SQR_FETCH_TERMS(SET_SHALLOW8);
            }
            else
            {
                SQR_FETCH_TERMS(SET_SHALLOW_GENERIC);
            }

            SQR_FINALIZE_DOT;
        }

        len1 += (gr_is_zero(GR_ENTRY(p1, len1, sz), cctx) != T_TRUE);

        while (Q_len > 0)
        {
            /* take node from store */
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? spawn stream (i+1) at column j.
               Triangular init keeps this on or above the diagonal. */
            if ((i + 1 < len2) && (hind[i + 1] == 2*j + 1))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                mpoly_monomial_add_any_bits(exp_list[exp_next], exp2 + x->i*N, exp2 + x->j*N, N, bits);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }

            /* should we go up? advance (i, j+1) within stream i */
            if ((j + 1 < len2) && ((hind[i] & 1) == 1) && ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1)))
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                mpoly_monomial_add_any_bits(exp_list[exp_next], exp2 + x->i*N, exp2 + x->j*N, N, bits);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }
        }
    }

    (*coeff1) = p1;
    (*exp1) = e1;

    TMP_END;
    GR_TMP_CLEAR(t, cctx);

    *res_len = len1;

    return status;
}

int gr_mpoly_sqr_commutative_heap(
    gr_mpoly_t poly1,
    const gr_mpoly_t poly2,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    slong i, N, len1 = 0;
    flint_bitcnt_t exp_bits;
    fmpz * max_fields2;
    ulong * cmpmask;
    ulong * exp2 = poly2->exps;
    int free2 = 0;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (poly2->length == 0)
        return gr_mpoly_zero(poly1, ctx);

    /* Squaring exploits a_i*a_j = a_j*a_i; a commutative (or approximately
       commutative) coefficient ring is assumed a priori by the caller. */

    TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(max_fields2 + i);

    mpoly_max_fields_fmpz(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, mctx);
    /* result exponents are bounded by 2 * (max input exponent) */
    _fmpz_vec_add(max_fields2, max_fields2, max_fields2, mctx->nfields);

    exp_bits = _fmpz_vec_max_bits(max_fields2, mctx->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(max_fields2 + i);

    N = mpoly_words_per_exp(exp_bits, mctx);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

    /* ensure input exponents are packed into same sized fields as output */
    if (exp_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                    poly2->length, mctx);
    }

    /* deal with aliasing and do the squaring */
    if (poly1 == poly2)
    {
        gr_mpoly_t temp;

        gr_mpoly_init(temp, ctx);
        gr_mpoly_fit_length_reset_bits(temp,
                                poly2->length + poly2->length, exp_bits, ctx);

        status = _gr_mpoly_sqr_commutative_heap(&len1,
                                &temp->coeffs, &temp->exps, &temp->coeffs_alloc, &temp->exps_alloc,
                                  poly2->coeffs, exp2, poly2->length,
                                      exp_bits, N, cmpmask, ctx);

        gr_mpoly_swap(temp, poly1, ctx);
        gr_mpoly_clear(temp, ctx);
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(poly1, poly2->length + poly2->length, exp_bits, ctx);

        status = _gr_mpoly_sqr_commutative_heap(&len1,
                                &poly1->coeffs, &poly1->exps, &poly1->coeffs_alloc, &poly1->exps_alloc,
                                  poly2->coeffs, exp2, poly2->length,
                                      exp_bits, N, cmpmask, ctx);
    }

    if (free2)
        flint_free(exp2);

    _gr_mpoly_set_length(poly1, len1, ctx);

    TMP_END;

    return status;
}
