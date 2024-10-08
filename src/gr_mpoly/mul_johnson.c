/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"
#include "gr_mpoly.h"

static int _gr_mpoly_mul_johnson(
    slong * res_len,
    gr_ptr * coeff1, ulong ** exp1, slong * alloc, slong * exps_alloc,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_srcptr coeff3, const ulong * exp3, slong len3,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    gr_ctx_t cctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(cctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(cctx, ADD);
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
    gr_ptr pp;
    slong sz = cctx->sizeof_elem;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;
    GR_TMP_INIT(pp, cctx);

    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));
    for (i = 0; i < len2; i++)
        exp_list[i] = exps + i*N;

    hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
    for (i = 0; i < len2; i++)
        hind[i] = 1;

    /* start with no heap nodes and no exponent vectors in use */
    exp_next = 0;

    /* put (0, 0, exp2[0] + exp3[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    if (bits <= FLINT_BITS)
        mpoly_monomial_add(heap[1].exp, exp2, exp3, N);
    else
        mpoly_monomial_add_mp(heap[1].exp, exp2, exp3, N);

    hind[0] = 2*1 + 0;

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _gr_mpoly_fit_length(&p1, alloc, &e1, exps_alloc, N, len1 + 1, cctx);

        mpoly_monomial_set(e1 + len1*N, exp, N);

        status |= gr_zero(GR_ENTRY(p1, len1, sz), cctx);

        do
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;

            /* todo: use addmul when better */
            status |= mul(pp, GR_ENTRY(coeff2, x->i, sz), GR_ENTRY(coeff3, x->j, sz), cctx);
            status |= add(GR_ENTRY(p1, len1, sz), GR_ENTRY(p1, len1, sz), pp, cctx);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                Q[Q_len++] = x->i;
                Q[Q_len++] = x->j;

                /* todo: use addmul when better */
                status |= mul(pp, GR_ENTRY(coeff2, x->i, sz), GR_ENTRY(coeff3, x->j, sz), cctx);
                status |= add(GR_ENTRY(p1, len1, sz), GR_ENTRY(p1, len1, sz), pp, cctx);
            }
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        len1 += (gr_is_zero(GR_ENTRY(p1, len1, sz), cctx) != T_TRUE);

        while (Q_len > 0)
        {
            /* take node from store */
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? */
            if (  (i + 1 < len2)
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }

            /* should we go up? */
            if (  (j + 1 < len3)
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }
        }
    }

    (* coeff1) = p1;
    (* exp1) = e1;

    TMP_END;
    GR_TMP_CLEAR(pp, cctx);

    *res_len = len1;

    return status;
}

int gr_mpoly_mul_johnson(
    gr_mpoly_t poly1,
    const gr_mpoly_t poly2,
    const gr_mpoly_t poly3,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong i, N, len1 = 0;
    flint_bitcnt_t exp_bits;
    fmpz * max_fields2, * max_fields3;
    ulong * cmpmask;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    int free2 = 0, free3 = 0;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (poly2->length == 0 || poly3->length == 0)
    {
        return gr_mpoly_zero(poly1, mctx, cctx);
    }

    if (poly3->length == 1)
    {
        return gr_mpoly_mul_monomial(poly1, poly2, poly3, mctx, cctx);
    }

    /* todo: could have a version of mul_monomial for the noncommutative case */
    if (poly2->length == 1 && gr_ctx_is_commutative_ring(cctx) == T_TRUE)
    {
        return gr_mpoly_mul_monomial(poly1, poly3, poly2, mctx, cctx);
    }

    TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    max_fields3 = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
    {
        fmpz_init(max_fields2 + i);
        fmpz_init(max_fields3 + i);
    }
    mpoly_max_fields_fmpz(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, mctx);
    mpoly_max_fields_fmpz(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, mctx);
    _fmpz_vec_add(max_fields2, max_fields2, max_fields3, mctx->nfields);

    exp_bits = _fmpz_vec_max_bits(max_fields2, mctx->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = FLINT_MAX(exp_bits, poly3->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    for (i = 0; i < mctx->nfields; i++)
    {
        fmpz_clear(max_fields2 + i);
        fmpz_clear(max_fields3 + i);
    }

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

    if (exp_bits > poly3->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
        mpoly_repack_monomials(exp3, exp_bits, poly3->exps, poly3->bits,
                                                    poly3->length, mctx);
    }

    /* deal with aliasing and do multiplication */
    if (poly1 == poly2 || poly1 == poly3)
    {
        gr_mpoly_t temp;

        gr_mpoly_init(temp, mctx, cctx);
        gr_mpoly_fit_length_reset_bits(temp,
                                poly2->length + poly3->length, exp_bits, mctx, cctx);

        if (poly2->length >= poly3->length)
        {
            status = _gr_mpoly_mul_johnson(&len1,
                                    &temp->coeffs, &temp->exps, &temp->coeffs_alloc, &temp->exps_alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                          exp_bits, N, cmpmask, cctx);
        }
        else
        {
            status = _gr_mpoly_mul_johnson(&len1,
                &temp->coeffs, &temp->exps, &temp->coeffs_alloc, &temp->exps_alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                          exp_bits, N, cmpmask, cctx);
        }

        gr_mpoly_swap(temp, poly1, mctx, cctx);
        gr_mpoly_clear(temp, mctx, cctx);
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(poly1, poly2->length + poly3->length, exp_bits, mctx, cctx);

        if (poly2->length > poly3->length)
        {
            status = _gr_mpoly_mul_johnson(&len1,
                                    &poly1->coeffs, &poly1->exps, &poly1->coeffs_alloc, &poly1->exps_alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                          exp_bits, N, cmpmask, cctx);
        }
        else
        {
            status = _gr_mpoly_mul_johnson(&len1,
                                    &poly1->coeffs, &poly1->exps, &poly1->coeffs_alloc, &poly1->exps_alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                          exp_bits, N, cmpmask, cctx);
        }
    }

    if (free2)
        flint_free(exp2);

    if (free3)
        flint_free(exp3);

    _gr_mpoly_set_length(poly1, len1, mctx, cctx);

    TMP_END;

    return status;
}
