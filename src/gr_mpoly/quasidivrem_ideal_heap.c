/*
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_vec.h"
#include "gr_mpoly.h"

/*
    GR port of fmpz_mpoly_quasidivrem_ideal_heap.

    Pseudo-division by a list of divisors: compute a coefficient-ring scalar
    "scale", quotients Q[0], ..., Q[len-1] and a remainder R with

        scale * A = Q[0] B[0] + ... + Q[len-1] B[len-1] + R,

    where every term of R has a monomial divisible by none of the leading
    monomials of the B[w].  Each leading term is reduced by the first divisor
    whose leading monomial divides it; because pseudo-division scales rather
    than requiring exact coefficient division, that first divisor always
    reduces the term (no field/non-field distinction is needed).  Over a field
    scale = 1 and this is ordinary ideal division.  Scaling uses gcd
    minimisation when the ring provides gcd, else the full lc(B[w]); quotient
    and remainder coefficients are committed at the scale in force and
    normalised at the end.  As always for division by a list, the result
    depends on the order of the divisors.
*/

static void
_gr_vec_grow(gr_ptr * vec, slong * alloc, slong len, gr_ctx_t cctx)
{
    if (len > *alloc)
    {
        slong sz = cctx->sizeof_elem;
        slong newalloc = FLINT_MAX(len, 2 * (*alloc));
        *vec = flint_realloc(*vec, newalloc * sz);
        _gr_vec_init(GR_ENTRY(*vec, *alloc, sz), newalloc - *alloc, cctx);
        *alloc = newalloc;
    }
}

static int _gr_mpoly_quasidivrem_ideal_heap(
    gr_ptr scale, gr_mpoly_struct ** Q, gr_mpoly_t R, int * overflowed,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_mpoly_struct * const * poly3, ulong * const * exp3, slong len,
    flint_bitcnt_t bits, slong N, const ulong * cmpmask,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
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
    int scaleis1, cstatus;
    gr_ptr acc, pp, tp, g, ns;
    gr_ptr * qs, rs;         /* per-term scale at commit time */
    slong * qs_alloc, rs_alloc;
    int status = GR_SUCCESS;
    TMP_INIT;

    *overflowed = 0;

    TMP_START;

    GR_TMP_INIT5(acc, pp, tp, g, ns, cctx);
    status |= gr_one(scale, cctx);
    scaleis1 = 1;

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

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));

    /* per-divisor quotient scales, plus one for the remainder */
    qs = (gr_ptr *) TMP_ALLOC(len*sizeof(gr_ptr));
    qs_alloc = (slong *) TMP_ALLOC(len*sizeof(slong));
    for (w = 0; w < len; w++)
    {
        k[w] = WORD(0);
        s[w] = poly3[w]->length;
        qs_alloc[w] = 16;
        qs[w] = flint_malloc(qs_alloc[w] * sz);
        _gr_vec_init(qs[w], qs_alloc[w], cctx);
    }
    rs_alloc = 16;
    rs = flint_malloc(rs_alloc * sz);
    _gr_vec_init(rs, rs_alloc, cctx);

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
        status |= gr_zero(acc, cctx);

        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do
            {
                store[store_len++] = x->i;
                store[store_len++] = x->j;
                store[store_len++] = x->p;
                if (x->i != -UWORD(1))
                    hinds[x->p][x->i] |= WORD(1);

                if (x->i == -UWORD(1))
                {
                    /* dividend term: + scale * poly2[j] */
                    if (scaleis1)
                        status |= gr_add(acc, acc, GR_ENTRY(coeff2, x->j, sz), cctx);
                    else
                    {
                        status |= gr_mul(pp, scale, GR_ENTRY(coeff2, x->j, sz), cctx);
                        status |= gr_add(acc, acc, pp, cctx);
                    }
                }
                else
                {
                    /* quotient term: - poly3[p][i] * (scale/qs[p][j]) * Q[p][j] */
                    if (scaleis1)
                    {
                        status |= gr_mul(pp, GR_ENTRY(poly3[x->p]->coeffs, x->i, sz),
                                             GR_ENTRY(Q[x->p]->coeffs, x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                    else
                    {
                        status |= gr_divexact(tp, scale, GR_ENTRY(qs[x->p], x->j, sz), cctx);
                        status |= gr_mul(tp, tp, GR_ENTRY(Q[x->p]->coeffs, x->j, sz), cctx);
                        status |= gr_mul(pp, GR_ENTRY(poly3[x->p]->coeffs, x->i, sz), tp, cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                }
            } while ((x = x->next) != NULL);
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
                    mpoly_monomial_add_any_bits(exp_list[exp_next], exp3[p] + x->i*N, Q[p]->exps + x->j*N, N, bits);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                if (j + 1 == k[p])
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
                    mpoly_monomial_add_any_bits(exp_list[exp_next], exp3[p] + x->i*N, Q[p]->exps + x->j*N, N, bits);
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

        /* reduce by the first divisor whose leading monomial divides exp */
        {
            int reduced = 0;

            for (w = 0; w < len; w++)
            {
                int d1;

                d1 = mpoly_monomial_divides_any_bits(texp, exp, exp3[w], N, mask, bits);

                if (!d1)
                    continue;

                _gr_mpoly_fit_length(&Q[w]->coeffs, &Q[w]->coeffs_alloc,
                                     &Q[w]->exps, &Q[w]->exps_alloc, N, k[w] + 1, ctx);
                _gr_vec_grow(&qs[w], &qs_alloc[w], k[w] + 1, cctx);

                /* pseudo-quotient coefficient, scaling up if inexact */
                cstatus = gr_div(GR_ENTRY(Q[w]->coeffs, k[w], sz), acc,
                                     GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx);
                if (cstatus == GR_SUCCESS)
                {
                    status |= gr_set(GR_ENTRY(qs[w], k[w], sz), scale, cctx);
                }
                else if (cstatus == GR_DOMAIN)
                {
                    gr_srcptr lc = GR_ENTRY(poly3[w]->coeffs, 0, sz);

                    if (gr_gcd(g, acc, lc, cctx) == GR_SUCCESS)
                    {
                        status |= gr_divexact(GR_ENTRY(Q[w]->coeffs, k[w], sz), acc, g, cctx);
                        status |= gr_divexact(ns, lc, g, cctx);
                        status |= gr_mul(scale, scale, ns, cctx);
                    }
                    else
                    {
                        status |= gr_set(GR_ENTRY(Q[w]->coeffs, k[w], sz), acc, cctx);
                        status |= gr_mul(scale, scale, lc, cctx);
                    }
                    status |= gr_set(GR_ENTRY(qs[w], k[w], sz), scale, cctx);
                    scaleis1 = 0;
                }
                else
                {
                    status |= cstatus;
                    goto cleanup;
                }

                mpoly_monomial_set(Q[w]->exps + k[w]*N, texp, N);

                if (s[w] > 1)
                {
                    x = chains[w] + 1;
                    x->i = 1;
                    x->j = k[w];
                    x->p = w;
                    x->next = NULL;
                    hinds[w][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_any_bits(exp_list[exp_next], exp3[w] + N, Q[w]->exps + k[w]*N, N, bits);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                s[w] = 1;
                k[w]++;
                reduced = 1;
                break;
            }

            if (!reduced)
            {
                /* no divisor reduces this term: send it to the remainder */
                _gr_mpoly_fit_length(&r_coeff, &R->coeffs_alloc, &r_exp, &R->exps_alloc, N, r_len + 1, ctx);
                _gr_vec_grow(&rs, &rs_alloc, r_len + 1, cctx);
                status |= gr_set(GR_ENTRY(r_coeff, r_len, sz), acc, cctx);
                status |= gr_set(GR_ENTRY(rs, r_len, sz), scale, cctx);
                mpoly_monomial_set(r_exp + r_len*N, exp, N);
                r_len++;
            }
        }
    }

    /* normalise all committed terms to the final scale */
    for (w = 0; w < len && status == GR_SUCCESS; w++)
    {
        for (i = 0; i < k[w] && status == GR_SUCCESS; i++)
        {
            status |= gr_divexact(tp, scale, GR_ENTRY(qs[w], i, sz), cctx);
            status |= gr_mul(GR_ENTRY(Q[w]->coeffs, i, sz), GR_ENTRY(Q[w]->coeffs, i, sz), tp, cctx);
        }
    }
    for (i = 0; i < r_len && status == GR_SUCCESS; i++)
    {
        status |= gr_divexact(tp, scale, GR_ENTRY(rs, i, sz), cctx);
        status |= gr_mul(GR_ENTRY(r_coeff, i, sz), GR_ENTRY(r_coeff, i, sz), tp, cctx);
    }

cleanup:

    GR_TMP_CLEAR5(acc, pp, tp, g, ns, cctx);
    for (w = 0; w < len; w++)
    {
        _gr_vec_clear(qs[w], qs_alloc[w], cctx);
        flint_free(qs[w]);
    }
    _gr_vec_clear(rs, rs_alloc, cctx);
    flint_free(rs);

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
            Q[w]->length = k[w];
        R->length = r_len;
    }

    TMP_END;

    return (*overflowed) ? GR_SUCCESS : status;
}


static int _gr_mpoly_quasidivrem_ideal_arr(
    gr_ptr scale, gr_mpoly_struct ** Q, gr_mpoly_t R,
    const gr_mpoly_t A, gr_mpoly_struct * const * B, slong len,
    gr_mpoly_ctx_t ctx)
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
            return GR_DOMAIN;
        if (gr_is_zero(B[i]->coeffs, cctx) != T_FALSE)
            return GR_UNABLE;
        len3 = FLINT_MAX(len3, B[i]->length);
    }

    if (A->length == 0)
    {
        status = gr_one(scale, cctx);
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

    MPOLY_GET_CMPMASK_FLINT_MALLOC(cmpmask, N, exp_bits, mctx);

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

    /* if lm(A) < lm(B[i]) for all i, quotients are zero and R = A */
    for (i = 0; i < len; i++)
    {
        if (!mpoly_monomial_lt(exp2, exp3[i], N, cmpmask))
            break;
    }

    if (i == len)
    {
        status |= gr_one(scale, cctx);
        status |= gr_mpoly_set(R, A, ctx);
        for (i = 0; i < len; i++)
            status |= gr_mpoly_zero(Q[i], ctx);
        goto cleanup;
    }

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

        status = _gr_mpoly_quasidivrem_ideal_heap(scale, Q, r, &overflowed,
                            A->coeffs, exp2, A->length,
                            B, exp3, len, exp_bits, N, cmpmask, ctx);

        if (!overflowed)
            break;

        {
            slong old_exp_bits = exp_bits;
            ulong * old_exp2 = exp2, * old_exp3;

            exp2 = mpoly_monomials_repack_wider_cmpmask(&exp_bits,
                                                        &N, &cmpmask, old_exp2,
                                                        old_exp_bits,
                                                        A->length, A->length,
                                                        mctx);
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
        GR_IGNORE(gr_one(scale, cctx));
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


int gr_mpoly_quasidivrem_ideal(
    gr_ptr scale, gr_mpoly_vec_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_vec_t B,
    gr_mpoly_ctx_t ctx)
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

    status = _gr_mpoly_quasidivrem_ideal_arr(scale, Qptr, R, A, Bptr, len, ctx);

    flint_free(Qptr);
    flint_free(Bptr);

    return status;
}
