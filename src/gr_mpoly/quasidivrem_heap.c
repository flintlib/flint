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
    GR port of fmpz_mpoly_quasidivrem_heap.

    Pseudo-division: compute a scalar "scale" (a coefficient-ring element), a
    quotient Q and a remainder R with

        scale * A = Q * B + R,

    where no monomial of R is divisible by the leading monomial of B.  Over a
    field every nonzero leading coefficient is a unit, so scale = 1 and this
    reduces to ordinary division.  Over an integral domain such as Z, whenever a
    quotient coefficient would not be exact we scale the running computation and
    proceed, mirroring classical pseudo-division.  If the coefficient ring
    provides gcd we scale by the minimal factor lc(B)/gcd(acc, lc(B)) (as the
    fmpz version does), keeping the scale small; otherwise we fall back to
    scaling by the full lc(B).  Quotient and remainder coefficients are committed
    at the scale in force at the time and normalised to the final scale at the
    end.
*/

/* grow a scratch gr vector, initialising the new entries */
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

static int _gr_mpoly_quasidivrem_heap(
    gr_ptr scale, gr_mpoly_t Q, gr_mpoly_t R, int * overflowed,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    gr_srcptr coeff3, const ulong * exp3, slong len3,
    flint_bitcnt_t bits, slong N, const ulong * cmpmask,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    slong i, j, s = len3;
    slong q_len = 0, r_len = 0;
    slong next_loc, heap_len = 2;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    slong * hind;
    int write_r = (R != NULL);
    gr_ptr q_coeff = Q->coeffs;
    gr_ptr r_coeff = write_r ? R->coeffs : NULL;
    ulong * q_exp = Q->exps;
    ulong * r_exp = write_r ? R->exps : NULL;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    int lt_divides, scaleis1, cstatus;
    gr_srcptr lc = coeff3;   /* leading coefficient of B */
    gr_ptr acc, pp, tp, g, ns;
    gr_ptr qs, rs;           /* per-term scale at commit time */
    slong qs_alloc, rs_alloc;
    int status = GR_SUCCESS;
    TMP_INIT;

    *overflowed = 0;

    TMP_START;

    GR_TMP_INIT5(acc, pp, tp, g, ns, cctx);
    status |= gr_one(scale, cctx);
    scaleis1 = 1;

    qs_alloc = 16;
    qs = flint_malloc(qs_alloc * sz);
    _gr_vec_init(qs, qs_alloc, cctx);
    rs_alloc = 16;
    rs = flint_malloc(rs_alloc * sz);
    _gr_vec_init(rs, rs_alloc, cctx);

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

    x = chain + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    while (heap_len > 1)
    {
        _gr_mpoly_fit_length(&q_coeff, &Q->coeffs_alloc, &q_exp, &Q->exps_alloc, N, q_len + 1, ctx);
        _gr_vec_grow(&qs, &qs_alloc, q_len + 1, cctx);
        if (write_r)
        {
            _gr_mpoly_fit_length(&r_coeff, &R->coeffs_alloc, &r_exp, &R->exps_alloc, N, r_len + 1, ctx);
            _gr_vec_grow(&rs, &rs_alloc, r_len + 1, cctx);
        }

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
            {
                *overflowed = 1;
                goto cleanup;
            }
            lt_divides = mpoly_monomial_divides(q_exp + q_len*N, exp, exp3, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
            {
                *overflowed = 1;
                goto cleanup;
            }
            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N, exp, exp3, N, bits);
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
                    hind[x->i] |= WORD(1);
                    /* quotient term: - poly3[i] * (scale/qs[j]) * q_coeff[j] */
                    if (scaleis1)
                    {
                        status |= gr_mul(pp, GR_ENTRY(coeff3, x->i, sz), GR_ENTRY(q_coeff, x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                    else
                    {
                        status |= gr_divexact(tp, scale, GR_ENTRY(qs, x->j, sz), cctx);
                        status |= gr_mul(tp, tp, GR_ENTRY(q_coeff, x->j, sz), cctx);
                        status |= gr_mul(pp, GR_ENTRY(coeff3, x->i, sz), tp, cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                }
            } while ((x = x->next) != NULL);
        }

        /* process popped nodes */
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
                    mpoly_monomial_add_any_bits(exp_list[exp_next], exp3 + x->i*N, q_exp + x->j*N, N, bits);
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
                    mpoly_monomial_add_any_bits(exp_list[exp_next], exp3 + x->i*N, q_exp + x->j*N, N, bits);
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

        if (!lt_divides)
        {
            /* remainder term, committed at the current scale (skipped for the
               remainder-free quasidiv variant) */
            if (write_r)
            {
                status |= gr_set(GR_ENTRY(r_coeff, r_len, sz), acc, cctx);
                mpoly_monomial_set(r_exp + r_len*N, exp, N);
                status |= gr_set(GR_ENTRY(rs, r_len, sz), scale, cctx);
                r_len++;
            }
            continue;
        }

        /* quotient coefficient = acc / lc(B); scale up if inexact */
        cstatus = gr_div(GR_ENTRY(q_coeff, q_len, sz), acc, lc, cctx);
        if (cstatus == GR_SUCCESS)
        {
            /* exact: no scaling needed */
            status |= gr_set(GR_ENTRY(qs, q_len, sz), scale, cctx);
        }
        else if (cstatus == GR_DOMAIN)
        {
            /*
                Inexact.  If the ring provides a gcd, scale by the minimal
                factor ns = lc/gcd(acc,lc): then q_coeff = acc/gcd is exact and
                the scale grows as little as possible.  Otherwise fall back to
                scaling by the full leading coefficient (q_coeff = acc).
            */
            if (gr_gcd(g, acc, lc, cctx) == GR_SUCCESS)
            {
                /* g divides both acc and lc exactly */
                status |= gr_divexact(GR_ENTRY(q_coeff, q_len, sz), acc, g, cctx);
                status |= gr_divexact(ns, lc, g, cctx);
                status |= gr_mul(scale, scale, ns, cctx);
            }
            else
            {
                status |= gr_set(GR_ENTRY(q_coeff, q_len, sz), acc, cctx);
                status |= gr_mul(scale, scale, lc, cctx);
            }
            status |= gr_set(GR_ENTRY(qs, q_len, sz), scale, cctx);
            scaleis1 = 0;
        }
        else
        {
            status |= cstatus;
            goto cleanup;
        }

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            mpoly_monomial_add_any_bits(exp_list[exp_next], exp3 + x->i*N, q_exp + x->j*N, N, bits);
            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                     &next_loc, &heap_len, N, cmpmask);
        }
        q_len++;
        s = 1;
    }

    /* normalise every committed term to the final scale */
    for (i = 0; i < q_len && status == GR_SUCCESS; i++)
    {
        status |= gr_divexact(tp, scale, GR_ENTRY(qs, i, sz), cctx);
        status |= gr_mul(GR_ENTRY(q_coeff, i, sz), GR_ENTRY(q_coeff, i, sz), tp, cctx);
    }
    for (i = 0; i < r_len && status == GR_SUCCESS; i++)
    {
        status |= gr_divexact(tp, scale, GR_ENTRY(rs, i, sz), cctx);
        status |= gr_mul(GR_ENTRY(r_coeff, i, sz), GR_ENTRY(r_coeff, i, sz), tp, cctx);
    }

cleanup:

    GR_TMP_CLEAR5(acc, pp, tp, g, ns, cctx);
    _gr_vec_clear(qs, qs_alloc, cctx);
    flint_free(qs);
    _gr_vec_clear(rs, rs_alloc, cctx);
    flint_free(rs);

    Q->coeffs = q_coeff;
    Q->exps = q_exp;
    if (write_r)
    {
        R->coeffs = r_coeff;
        R->exps = r_exp;
    }

    if (*overflowed || status != GR_SUCCESS)
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


int gr_mpoly_quasidivrem_heap(
    gr_ptr scale, gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B,
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
        status = gr_one(scale, cctx);
        status |= gr_mpoly_zero(Q, ctx);
        if (R != NULL)
            status |= gr_mpoly_zero(R, ctx);
        return status;
    }

    if (gr_is_zero(B->coeffs, cctx) != T_FALSE)
        return GR_UNABLE;

    exp_bits = FLINT_MAX(A->bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    MPOLY_GET_CMPMASK_FLINT_MALLOC(cmpmask, N, exp_bits, mctx);

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

    /* if lm(A) < lm(B) the quotient is zero and the remainder is A */
    if (mpoly_monomial_lt(exp2, exp3, N, cmpmask))
    {
        status = gr_one(scale, cctx);
        if (R != NULL)
            status |= gr_mpoly_set(R, A, ctx);
        status |= gr_mpoly_zero(Q, ctx);
        flint_free(cmpmask);
        if (free2) flint_free(exp2);
        if (free3) flint_free(exp3);
        return status;
    }

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
        r = NULL;   /* remainder-free: quasidiv */
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

        status = _gr_mpoly_quasidivrem_heap(scale, q, r, &overflowed,
                            A->coeffs, exp2, A->length,
                            B->coeffs, exp3, B->length, exp_bits, N, cmpmask, ctx);

        if (!overflowed)
            break;

        {
            slong old_exp_bits = exp_bits;
            ulong * old_exp2 = exp2;
            ulong * old_exp3 = exp3;

            exp2 = mpoly_monomials_repack_wider_cmpmask(&exp_bits,
                                                        &N, &cmpmask, old_exp2,
                                                        old_exp_bits,
                                                        A->length, A->length,
                                                        mctx);
            if (free2) flint_free(old_exp2);
            free2 = 1;

            exp3 = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
            mpoly_repack_monomials(exp3, exp_bits, old_exp3, old_exp_bits, B->length, mctx);
            if (free3) flint_free(old_exp3);
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
        GR_IGNORE(gr_mpoly_zero(Q, ctx));
        if (R != NULL)
            GR_IGNORE(gr_mpoly_zero(R, ctx));
        GR_IGNORE(gr_one(scale, cctx));
    }

    flint_free(cmpmask);
    if (free2) flint_free(exp2);
    if (free3) flint_free(exp3);

    return status;
}


int gr_mpoly_quasidivrem(
    gr_ptr scale, gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return gr_mpoly_quasidivrem_heap(scale, Q, R, A, B, ctx);
}


int gr_mpoly_quasidiv(
    gr_ptr scale, gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    /* remainder-free: the kernel skips all remainder writes when R is NULL */
    return gr_mpoly_quasidivrem_heap(scale, Q, NULL, A, B, ctx);
}
