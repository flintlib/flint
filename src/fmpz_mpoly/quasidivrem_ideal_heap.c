/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


slong _fmpz_mpoly_quasidivrem_ideal_heap1(fmpz_t scale, fmpz_mpoly_struct ** polyq,
  fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2,
     const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3,
                     ulong * const * exp3, slong len, slong bits, 
                                      const fmpz_mpoly_ctx_t ctx, ulong maskhi)
{
    slong i, j, p, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    ulong exp, texp;
    ulong mask;
    fmpz ** qs, * rs;
    slong * qs_alloc, rs_alloc;
    slong * q_len, * s;
    fmpz * r_coeff = *polyr;
    ulong * r_exp = *expr;
    slong r_len;
    fmpz_t ns, gcd, acc_lg, tp;
    TMP_INIT;

    TMP_START;

    fmpz_init(ns);
    fmpz_init(gcd);
    fmpz_init(acc_lg);
    fmpz_init(tp);

    fmpz_one(scale);

    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    q_len = (slong *) TMP_ALLOC(len*sizeof(slong));
    qs_alloc = (slong *) TMP_ALLOC(len*sizeof(slong));
    qs = (fmpz **) TMP_ALLOC(len*sizeof(fmpz *));
    for (w = 0; w < len; w++)
    {
        q_len[w] = WORD(0);
        qs_alloc[w] = 16;
        qs[w] = (fmpz *) flint_calloc(qs_alloc[w], sizeof(fmpz));
        s[w] = poly3[w]->length;
    }
    r_len = WORD(0);
    rs_alloc = 64;
    rs = (fmpz *) flint_calloc(rs_alloc, sizeof(fmpz));

    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));

    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((poly3[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((poly3[w]->length)*sizeof(slong));
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    mask = mpoly_overflow_mask_sp(bits);

    x = chains[0] + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    while (heap_len > 1)
    {
        exp = heap[1].exp;
        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        fmpz_zero(acc_lg);
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;
                if (x->i != -WORD(1))
                    hinds[x->p][x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    fmpz_addmul(acc_lg, scale, poly2 + x->j);
                } else
                {
                    fmpz_divexact(tp, scale, qs[x->p] + x->j);
                    fmpz_mul(tp, tp, polyq[x->p]->coeffs + x->j);
                    fmpz_submul(acc_lg, poly3[x->p]->coeffs + x->i, tp);
                }

            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        while (store > store_base)
        {
            p = *--store;
            j = *--store;
            i = *--store;
            if (i == -WORD(1))
            {
                if (j + 1 < len2)
                {
                    x = chains[0] + 0;
                    x->i = -WORD(1);
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, exp2[x->j], x, &next_loc, &heap_len, maskhi);
                }
            } else
            {
                if ( (i + 1 < poly3[p]->length)
                   && (hinds[p][i + 1] == 2*j + 1)
                   )
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[x->p][x->i] +
                                polyq[x->p]->exps[x->j], x, &next_loc, &heap_len, maskhi);
                }
                if (j + 1 == q_len[p])
                {
                    s[p]++;
                } else if (  ((hinds[p][i] & 1) == 1)
                          && ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[x->p][x->i] +
                                polyq[x->p]->exps[x->j], x, &next_loc, &heap_len, maskhi);
                }
            }
        }

        if (fmpz_is_zero(acc_lg))
            continue;

        for (w = 0; w < len; w++)
        {
            if (mpoly_monomial_divides1(&texp, exp, exp3[w][0], mask))
            {
                fmpz_mpoly_fit_length(polyq[w], q_len[w] + 1, ctx);
                if (q_len[w] + 1 > qs_alloc[w])
                {
                    slong len = FLINT_MAX(q_len[w] + 1, 2*qs_alloc[w]);
                    qs[w] = (fmpz *) flint_realloc(qs[w], len*sizeof(fmpz));
                    flint_mpn_zero((mp_ptr) (qs[w] + qs_alloc[w]), len - qs_alloc[w]);
                    qs_alloc[w] = len;
                }

                fmpz_gcd(gcd, acc_lg, poly3[w]->coeffs + 0);
                fmpz_divexact(ns, poly3[w]->coeffs + 0, gcd);
                fmpz_abs(ns, ns);
                fmpz_mul(polyq[w]->coeffs + q_len[w], ns, acc_lg);
                fmpz_divexact(polyq[w]->coeffs + q_len[w],
                            polyq[w]->coeffs + q_len[w], poly3[w]->coeffs + 0);
                fmpz_mul(qs[w] + q_len[w], scale, ns);
                fmpz_set(scale, qs[w] + q_len[w]);
                polyq[w]->exps[q_len[w]] = texp;

                if (s[w] > 1)
                {
                    i = 1;
                    x = chains[w] + i;
                    x->i = i;
                    x->j = q_len[w];
                    x->p = w;
                    x->next = NULL;
                    hinds[w][x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[w][x->i] +
                                          polyq[w]->exps[x->j], x,
                                         &next_loc, &heap_len, maskhi);
                }
                s[w] = 1;

                q_len[w]++;

                goto break_continue; /* break out of w for loop and continue in heap loop */
            }
        }

        /* if get here, no leading terms divided */

        _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, 1);
        if (r_len + 1 > rs_alloc)
        {
            slong len = FLINT_MAX(r_len + 1, 2*rs_alloc);
            rs = (fmpz *) flint_realloc(rs, len*sizeof(fmpz));
            flint_mpn_zero((mp_ptr) (rs + rs_alloc), len - rs_alloc);
            rs_alloc = len;
        }
        fmpz_set(r_coeff + r_len, acc_lg);
        fmpz_set(rs + r_len, scale);
        r_exp[r_len] = exp;

        r_len++;

break_continue:
        (void)(NULL);
    }


cleanup2:

    for (w = 0; w < len; w++)
    {
        _fmpz_mpoly_set_length(polyq[w], q_len[w], ctx); 
        for (i = 0; i < q_len[w]; i++)
        {
            fmpz_divexact(tp, scale, qs[w] + i);
            fmpz_mul(polyq[w]->coeffs + i, polyq[w]->coeffs + i, tp);
        }
        for (i = 0; i < qs_alloc[w]; i++)
            fmpz_clear(qs[w] + i);
        flint_free(qs[w]);
    }

    (*polyr) = r_coeff;
    (*expr) = r_exp;
    for (i = 0; i < r_len; i++)
    {
        fmpz_divexact(tp, scale, rs + i);
        fmpz_mul(r_coeff + i, r_coeff + i, tp);
    }
    for (i = 0; i < rs_alloc; i++)
        fmpz_clear(rs + i);
    flint_free(rs);

    fmpz_clear(ns);
    fmpz_clear(gcd);
    fmpz_clear(acc_lg);
    fmpz_clear(tp);

    TMP_END;

    return r_len;

exp_overflow:
    for (i = 0; i < r_len; i++)
       _fmpz_demote(r_coeff + i);
    r_len = -WORD(1);

    for (w = 0; w < len; w++)
    {
        for (i = 0; i < q_len[w]; i++)
           _fmpz_demote(polyq[w]->coeffs + i);
        q_len[w] = WORD(0);
    }

    goto cleanup2;
}



slong _fmpz_mpoly_quasidivrem_ideal_heap(fmpz_t scale, fmpz_mpoly_struct ** polyq,
  fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2,
     const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3,
                     ulong * const * exp3, slong len, slong N, slong bits, 
                        const fmpz_mpoly_ctx_t ctx, const ulong * cmpmask)
{
    slong i, j, p, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    ulong * exp, * exps, * texp;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    fmpz ** qs, * rs;
    slong * qs_alloc, rs_alloc;
    slong * q_len, * s;
    fmpz * r_coeff = *polyr;
    ulong * r_exp = *expr;
    slong r_len;
    fmpz_t ns, gcd, acc_lg, tp;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_quasidivrem_ideal_heap1(scale, polyq, polyr, expr,
               allocr, poly2, exp2, len2, poly3, exp3, len, bits, ctx, cmpmask[0]);

    TMP_START;

    fmpz_init(ns);
    fmpz_init(gcd);
    fmpz_init(acc_lg);
    fmpz_init(tp);

    fmpz_one(scale);

    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    q_len = (slong *) TMP_ALLOC(len*sizeof(slong));
    qs_alloc = (slong *) TMP_ALLOC(len*sizeof(slong));
    qs = (fmpz **) TMP_ALLOC(len*sizeof(fmpz *));
    for (w = 0; w < len; w++)
    {
        q_len[w] = WORD(0);
        qs_alloc[w] = 16;
        qs[w] = (fmpz *) flint_calloc(qs_alloc[w], sizeof(fmpz));
        s[w] = poly3[w]->length;
    }
    r_len = WORD(0);
    rs_alloc = 64;
    rs = (fmpz *) flint_calloc(rs_alloc, sizeof(fmpz));

    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));
    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((poly3[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((poly3[w]->length)*sizeof(slong));
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;
   
    x = chains[0] + 0;
    x->i = -WORD(1);
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
                goto exp_overflow;
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto exp_overflow;
        }

        fmpz_zero(acc_lg);
        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;
                if (x->i != -WORD(1))
                    hinds[x->p][x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    fmpz_addmul(acc_lg, scale, poly2 + x->j);
                } else
                {
                    fmpz_divexact(tp, scale, qs[x->p] + x->j);
                    fmpz_mul(tp, tp, polyq[x->p]->coeffs + x->j);
                    fmpz_submul(acc_lg, poly3[x->p]->coeffs + x->i, tp);
                }

            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        while (store > store_base)
        {
            p = *--store;
            j = *--store;
            i = *--store;
            if (i == -WORD(1))
            {
                if (j + 1 < len2)
                {
                    x = chains[0] + 0;
                    x->i = -WORD(1);
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            } else
            {
                /* should we go right */
                if ( (i + 1 < poly3[p]->length)
                   && (hinds[p][i + 1] == 2*j + 1)
                   )
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3[x->p] + x->i*N,
                                                polyq[x->p]->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }

                /* should we go up? */
                if (j + 1 == q_len[p])
                {
                    s[p]++;
                } else if (  ((hinds[p][i] & 1) == 1)
                          && ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3[x->p] + x->i*N,
                                                polyq[x->p]->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if (fmpz_is_zero(acc_lg))
            continue;

        for (w = 0; w < len; w++)
        {
            int divides;

            if (bits <= FLINT_BITS)
                divides = mpoly_monomial_divides(texp, exp, exp3[w] + N*0, N, mask);
            else
                divides = mpoly_monomial_divides_mp(texp, exp, exp3[w] + N*0, N, bits);

            if (divides)
            {
                fmpz_mpoly_fit_length(polyq[w], q_len[w] + 1, ctx);
                if (q_len[w] + 1 > qs_alloc[w])
                {
                    slong len = FLINT_MAX(q_len[w] + 1, 2*qs_alloc[w]);
                    qs[w] = (fmpz *) flint_realloc(qs[w], len*sizeof(fmpz));
                    flint_mpn_zero((mp_ptr) (qs[w] + qs_alloc[w]), len - qs_alloc[w]);
                    qs_alloc[w] = len;
                }

                fmpz_gcd(gcd, acc_lg, poly3[w]->coeffs + 0);
                fmpz_divexact(ns, poly3[w]->coeffs + 0, gcd);
                fmpz_abs(ns, ns);
                fmpz_mul(polyq[w]->coeffs + q_len[w], ns, acc_lg);
                fmpz_divexact(polyq[w]->coeffs + q_len[w],
                            polyq[w]->coeffs + q_len[w], poly3[w]->coeffs + 0);
                fmpz_mul(qs[w] + q_len[w], scale, ns);
                fmpz_set(scale, qs[w] + q_len[w]);
                mpoly_monomial_set(polyq[w]->exps + N*q_len[w], texp, N);

                if (s[w] > 1)
                {
                    i = 1;
                    x = chains[w] + i;
                    x->i = i;
                    x->j = q_len[w];
                    x->p = w;
                    x->next = NULL;
                    hinds[w][x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], exp3[w] + N*x->i, 
                                                    polyq[w]->exps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                s[w] = 1;

                q_len[w]++;

                goto break_continue; /* break out of w for loop and continue in heap loop */
            }
        }

        /* if get here, no leading terms divided */

        _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, N);
        if (r_len + 1 > rs_alloc)
        {
            slong len = FLINT_MAX(r_len + 1, 2*rs_alloc);
            rs = (fmpz *) flint_realloc(rs, len*sizeof(fmpz));
            flint_mpn_zero((mp_ptr) (rs + rs_alloc), len - rs_alloc);
            rs_alloc = len;
        }
        fmpz_set(r_coeff + r_len, acc_lg);
        fmpz_set(rs + r_len, scale);
        mpoly_monomial_set(r_exp + r_len*N, exp, N);

        r_len++;

break_continue:
        (void)(NULL);
    }


cleanup2:

    for (w = 0; w < len; w++)
    {
        _fmpz_mpoly_set_length(polyq[w], q_len[w], ctx); 
        for (i = 0; i < q_len[w]; i++)
        {
            fmpz_divexact(tp, scale, qs[w] + i);
            fmpz_mul(polyq[w]->coeffs + i, polyq[w]->coeffs + i, tp);
        }
        for (i = 0; i < qs_alloc[w]; i++)
            fmpz_clear(qs[w] + i);
        flint_free(qs[w]);
    }

    (*polyr) = r_coeff;
    (*expr) = r_exp;
    for (i = 0; i < r_len; i++)
    {
        fmpz_divexact(tp, scale, rs + i);
        fmpz_mul(r_coeff + i, r_coeff + i, tp);
    }
    for (i = 0; i < rs_alloc; i++)
        fmpz_clear(rs + i);
    flint_free(rs);

    fmpz_clear(ns);
    fmpz_clear(gcd);
    fmpz_clear(acc_lg);
    fmpz_clear(tp);

    TMP_END;

    return r_len;

exp_overflow:
    for (i = 0; i < r_len; i++)
       _fmpz_demote(r_coeff + i);
    r_len = -WORD(1);

    for (w = 0; w < len; w++)
    {
        for (i = 0; i < q_len[w]; i++)
           _fmpz_demote(polyq[w]->coeffs + i);
        q_len[w] = WORD(0);
    }

    goto cleanup2;
}

/* Assumes divisor polys don't alias any output polys */
void fmpz_mpoly_quasidivrem_ideal_heap(fmpz_t scale,
                                 fmpz_mpoly_struct ** q, fmpz_mpoly_t r,
                const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3,
                                         slong len, const fmpz_mpoly_ctx_t ctx)
{
    slong i, exp_bits, N, lenr = 0;
    slong len3 = 0;
    ulong * cmpmask;
    ulong * exp2;
    ulong ** exp3;
    int free2 = 0;
    int * free3;
    fmpz_mpoly_t temp2;
    fmpz_mpoly_struct * tr;
    TMP_INIT;

    /* check none of the divisor polynomials is zero */
    for (i = 0; i < len; i++)
    {  
        if (poly3[i]->length == 0)
            flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divrem_ideal_monagan_pearce");

        len3 = FLINT_MAX(len3, poly3[i]->length);
    }

    fmpz_one(scale);

    /* dividend is zero, write out quotients and remainder */
    if (poly2->length == 0)
    {
        for (i = 0; i < len; i++)
        {
            fmpz_mpoly_zero(q[i], ctx);
        }
        fmpz_mpoly_zero(r, ctx);
        return;
    }


    TMP_START;

    free3 = (int *) TMP_ALLOC(len*sizeof(int));
    exp3 = (ulong **) TMP_ALLOC(len*sizeof(ulong *));

    /* compute maximum degrees that can occur in any input or output polys */
    exp_bits = poly2->bits;
    for (i = 0; i < len; i++)
        exp_bits = FLINT_MAX(exp_bits, poly3[i]->bits);

    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    exp2 = poly2->exps;
    free2 = 0;
    if (exp_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
    }

    for (i = 0; i < len; i++)
    {
        exp3[i] = poly3[i]->exps;
        free3[i] = 0;
        if (exp_bits > poly3[i]->bits)
        {
            free3[i] = 1;
            exp3[i] = (ulong *) flint_malloc(N*poly3[i]->length*sizeof(ulong));
            mpoly_repack_monomials(exp3[i], exp_bits, poly3[i]->exps,
                                 poly3[i]->bits, poly3[i]->length, ctx->minfo);
        }
        fmpz_mpoly_fit_length(q[i], 1, ctx);
        fmpz_mpoly_fit_bits(q[i], exp_bits, ctx);
        q[i]->bits = exp_bits;
    }

    /* check leading mon. of at least one divisor is at most that of dividend */
    for (i = 0; i < len; i++)
    {
        if (!mpoly_monomial_lt(exp2, exp3[i], N, cmpmask))
            break;
    }

    if (i == len)
    {
        fmpz_mpoly_set(r, poly2, ctx);
        for (i = 0; i < len; i++)
            fmpz_mpoly_zero(q[i], ctx);

        goto cleanup3;
    }

    /* take care of aliasing */
    if (r == poly2)
    {
        fmpz_mpoly_init2(temp2, len3, ctx);
        fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);
        temp2->bits = exp_bits;
        tr = temp2;
    } else
    {
        fmpz_mpoly_fit_length(r, len3, ctx);
        fmpz_mpoly_fit_bits(r, exp_bits, ctx);
        r->bits = exp_bits;
        tr = r;
    }

    /* do division with remainder */
    while (1)
    {
        slong old_exp_bits = exp_bits;
        ulong * old_exp2 = exp2, * old_exp3;

        lenr = _fmpz_mpoly_quasidivrem_ideal_heap(scale,
                          q, &tr->coeffs, &tr->exps, &tr->alloc,
                             poly2->coeffs, exp2, poly2->length,
                           poly3, exp3, len, N, exp_bits, ctx, cmpmask);

        if (lenr >= 0) /* check if division was successful */
            break;

        exp_bits = mpoly_fix_bits(exp_bits + 1, ctx->minfo);
        N = mpoly_words_per_exp(exp_bits, ctx->minfo);
        cmpmask = (ulong *) flint_realloc(cmpmask, N*sizeof(ulong));
        mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, old_exp2, old_exp_bits,
                                                    poly2->length, ctx->minfo);

        if (free2)
            flint_free(old_exp2);

        free2 = 1;

        fmpz_mpoly_fit_bits(tr, exp_bits, ctx);
        tr->bits = exp_bits;

        for (i = 0; i < len; i++)
        {
            old_exp3 = exp3[i];

            exp3[i] = (ulong *) flint_malloc(N*poly3[i]->length*sizeof(ulong));
            mpoly_repack_monomials(exp3[i], exp_bits, old_exp3, old_exp_bits,
                                                 poly3[i]->length, ctx->minfo);
   
            if (free3[i])
                flint_free(old_exp3);

            free3[i] = 1; 

            fmpz_mpoly_fit_bits(q[i], exp_bits, ctx);
            q[i]->bits = exp_bits;
        }
    }

    /* take care of aliasing */
    if (r == poly2)
    {
        fmpz_mpoly_swap(temp2, r, ctx);
        fmpz_mpoly_clear(temp2, ctx);
    } 

    _fmpz_mpoly_set_length(r, lenr, ctx);

cleanup3:

    if (free2)
        flint_free(exp2);

    for (i = 0; i < len; i++)
    {
        if (free3[i])
            flint_free(exp3[i]);
    }

    flint_free(cmpmask);

    TMP_END;
}
