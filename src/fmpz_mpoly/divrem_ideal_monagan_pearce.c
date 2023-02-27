/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "longlong.h"

/* ensure rounding is towards -\infty like fmpz_fdiv_qr */
#define fdiv_qrnnd(qxx, rxx, nhixx, nloxx, dxx)        \
   do {                                                \
      slong __dxx = (dxx);                             \
      sdiv_qrnnd(qxx, rxx, nhixx, nloxx, dxx);         \
      if (((slong) (rxx) < 0 && (slong) (dxx) > 0) ||  \
          ((slong) (rxx) > 0 && (slong) (dxx) < 0))    \
      {                                                \
         (rxx) += __dxx;                               \
         (qxx)--;                                      \
      }                                                \
   } while (0)

/*
   As for divrem_monagan_pearce1 except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
slong _fmpz_mpoly_divrem_ideal_monagan_pearce1(fmpz_mpoly_struct ** polyq,
  fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2,
     const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3,
       ulong * const * exp3, slong len, slong bits, const fmpz_mpoly_ctx_t ctx,
                                                                  ulong maskhi)
{
    slong i, j, p, l, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    fmpz * p2 = *polyr;
    ulong * e2 = *expr;
    ulong exp, texp;
    ulong c[3]; /* for accumulating coefficients */
    ulong mask;
    ulong * ub;
    slong * k, * s;
    fmpz_t qc, q;
    fmpz * mb;
    int small;
    slong bits2, bits3;
    int d1, d2, div_flag;
    TMP_INIT;

    TMP_START;

    fmpz_init(q);
    fmpz_init(qc);

    bits2 = _fmpz_vec_max_bits(poly2, len2);

    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));

    bits3 = 0;
    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((poly3[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((poly3[w]->length)*sizeof(slong));
        bits3 = FLINT_MAX(bits3, FLINT_ABS(fmpz_mpoly_max_bits(poly3[w])));
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }
      
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX)
          && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    ub = (ulong *) TMP_ALLOC(len*sizeof(ulong));
    mb = (fmpz * ) TMP_ALLOC(len*sizeof(fmpz));

    mask = mpoly_overflow_mask_sp(bits);

    for (w = 0; w < len; w++)
    {
        k[w] = -WORD(1);
        s[w] = poly3[w]->length;
    }
    l = -WORD(1);
   
    x = chains[0] + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    for (i = 0; i < len; i++)
    {
        fmpz_init(mb + i);
        fmpz_neg(mb + i, poly3[i]->coeffs);
        ub[i] = ((ulong) FLINT_ABS(mb[i])) >> 1; /* abs(poly3[0])/2 */
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;
        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        c[0] = c[1] = c[2] = 0;
        fmpz_zero(qc);
        while (heap_len > 1 && heap[1].exp == exp)
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;
                if (x->i != -WORD(1))
                    hinds[x->p][x->i] |= WORD(1);

                if (small)
                {
                    if (x->i == -WORD(1))
                      _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
                    else
                      _fmpz_mpoly_addmul_uiuiui_fmpz(c,
                             poly3[x->p]->coeffs[x->i], polyq[x->p]->coeffs[x->j]);
                } else
                {
                    if (x->i == -WORD(1))
                        fmpz_sub(qc, qc, poly2 + x->j);
                    else
                        fmpz_addmul(qc, poly3[x->p]->coeffs + x->i,
                                                   polyq[x->p]->coeffs + x->j);
                }
            } while ((x = x->next) != NULL);
        }

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
                if (j == k[p])
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

        if ((small && (c[2] != 0 || c[1] != 0 || c[0] != 0)) ||
                                                 (!small && !fmpz_is_zero(qc)))
        {
            div_flag = 0;
            for (w = 0; w < len; w++)
            {
                d1 = mpoly_monomial_divides1(&texp, exp, exp3[w][0], mask);
                if (d1)
                {
                    d2 = 0;
                    if (small)
                    {
                        ulong d[3];
                        if (0 > (slong) c[2])
                            mpn_neg(d, c, 3);
                        else
                            flint_mpn_copyi(d, c, 3);

                        if (d[2] != 0 || ub[w] <= d[1] ||
                          (ub[w] == 0 && 0 > (slong) d[0])) /* quotient not a small */
                        {
                            fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                            small = 0;
                        } else /* quotient fits a small */
                        {
                            slong r1;
                            slong tq;
                            fdiv_qrnnd(tq, r1, c[1], c[0], mb[w]);
                            if (tq > COEFF_MAX || tq < COEFF_MIN)
                            {
                                small = 0;
                                fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                            } else
                            {
                                div_flag = (r1 == 0);
                                d2 = tq != 0;
                                if (d2)
                                {
                                    k[w]++;
                                    fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);
                                    fmpz_set_si(polyq[w]->coeffs + k[w], tq);
                                    polyq[w]->exps[k[w]] = texp;
                                }
                                c[0] = r1;
                                c[2] = c[1] = -(slong)(r1 < 0);
                            }
                        }
                    } 
                    /* quotient non-small case */
                    if (!small)
                    {
                        fmpz_fdiv_qr(q, qc, qc, mb + w);
                        div_flag = fmpz_is_zero(qc);
                        d2 = !fmpz_is_zero(q);
                        if (d2)
                        {
                            k[w]++;
                            fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);
                            fmpz_set(polyq[w]->coeffs + k[w], q);                     
                            polyq[w]->exps[k[w]] = texp;
                        }
                    }
                    if (d2)
                    {
                        if (s[w] > 1)
                        {
                            i = 1;
                            x = chains[w] + i;
                            x->i = i;
                            x->j = k[w];
                            x->p = w;
                            x->next = NULL;
                            hinds[w][x->i] = 2*(x->j + 1) + 0;
                            _mpoly_heap_insert1(heap, exp3[w][i] +
                                                  polyq[w]->exps[k[w]], x,
                                                 &next_loc, &heap_len, maskhi);
                        }
                        s[w] = 1;
                    }
                }
            }
            if (!div_flag)
            {
                l++;
                _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

                if (small)
                {
                    fmpz_set_signed_uiuiui(p2 + l, c[2], c[1], c[0]);
                    fmpz_neg(p2 + l, p2 + l);
                } else
                {
                    fmpz_neg(p2 + l, qc);
                }
                e2[l] = exp;
            }
        }
    }

cleanup:

   for (i = 0; i < len; i++)
      _fmpz_mpoly_set_length(polyq[i], k[i] + 1, ctx); 
   for (i = 0; i < len; i++)
      fmpz_clear(mb + i);
   fmpz_clear(qc);
   fmpz_clear(q);

   (*polyr) = p2;
   (*expr) = e2;
   
   TMP_END;
   return l + 1;

exp_overflow:
    for (i = 0; i < l; i++)
        _fmpz_demote(p2 + i);
    for (w = 0; w < len; w++)
    {
        for (i = 0; i < k[w]; i++)
            _fmpz_demote(polyq[w]->coeffs + i);
        k[w] = -WORD(1); /* we add 1 to this before exit */
    }
    l = -WORD(2); /* we add 1 to this upon return */
    goto cleanup;
}

/*
   As for divrem_monagan_pearce except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
slong _fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** polyq,
  fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2,
     const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3,
                     ulong * const * exp3, slong len, slong N, slong bits, 
                        const fmpz_mpoly_ctx_t ctx, const ulong * cmpmask)
{
    slong i, j, p, l, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    fmpz * p2 = *polyr;
    ulong * e2 = *expr;
    ulong * exp, * exps, * texp;
    ulong ** exp_list;
    ulong c[3]; /* for accumulating coefficients */
    slong exp_next;
    ulong mask = 0;
    ulong * ub;
    slong * k, * s;
    fmpz_t qc, q;
    fmpz * mb;
    int small;
    slong bits2, bits3;
    int d1, d2, div_flag;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_divrem_ideal_monagan_pearce1(polyq, polyr, expr,
           allocr, poly2, exp2, len2, poly3, exp3, len, bits, ctx, cmpmask[0]);

    TMP_START;

    fmpz_init(q);
    fmpz_init(qc);

    bits2 = _fmpz_vec_max_bits(poly2, len2);
   
    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));
    bits3 = 0;
    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((poly3[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((poly3[w]->length)*sizeof(slong));
        bits3 = FLINT_MAX(bits3, FLINT_ABS(fmpz_mpoly_max_bits(poly3[w])));
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }
      
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) +
           FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX) &&
           FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    ub = (ulong *) TMP_ALLOC(len*sizeof(ulong));
    mb = (fmpz * ) TMP_ALLOC(len*sizeof(fmpz));

    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    for (w = 0; w < len; w++)
    {
        k[w] = -WORD(1);
        s[w] = poly3[w]->length;
    }
    l = -WORD(1);
   
    x = chains[0] + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    for (i = 0; i < len; i++)
    {
        fmpz_init(mb + i);
        fmpz_neg(mb + i, poly3[i]->coeffs);
        ub[i] = ((ulong) FLINT_ABS(mb[i])) >> 1; /* abs(poly3[0])/2 */
    }

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        c[0] = c[1] = c[2] = 0;
        fmpz_zero(qc);

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

        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
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

                if (small)
                {
                    if (x->i == -WORD(1))
                        _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
                    else
                        _fmpz_mpoly_addmul_uiuiui_fmpz(c,
                             poly3[x->p]->coeffs[x->i], polyq[x->p]->coeffs[x->j]);
                } else
                {
                    if (x->i == -WORD(1))
                        fmpz_sub(qc, qc, poly2 + x->j);
                    else
                        fmpz_addmul(qc, poly3[x->p]->coeffs + x->i,
                                                       polyq[x->p]->coeffs + x->j);
                }
            } while ((x = x->next) != NULL);
        }

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
                if (j == k[p])
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

        if ((small && (c[2] != 0 || c[1] != 0 || c[0] != 0)) ||
                                                 (!small && !fmpz_is_zero(qc)))
        {
            div_flag = 0;
            for (w = 0; w < len; w++)
            {
                if (bits <= FLINT_BITS)
                    d1 = mpoly_monomial_divides(texp, exp, exp3[w], N, mask);
                else
                    d1 = mpoly_monomial_divides_mp(texp, exp, exp3[w], N, bits);

                if (d1)
                {
                    d2 = 0;
                    if (small)
                    {
                        ulong d[3];
                        if (0 > (slong) c[2])
                            mpn_neg(d, c, 3);
                        else
                            flint_mpn_copyi(d, c, 3);

                        if (d[2] != 0 || ub[w] <= d[1] ||
                          (ub[w] == 0 && 0 > (slong) d[0])) /* quotient not a small */
                        {
                            small = 0;
                            fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                        } else /* quotient fits a small */
                        {
                            slong r1;
                            slong tq;
                            fdiv_qrnnd(tq, r1, c[1], c[0], mb[w]);
                            if (tq > COEFF_MAX || tq < COEFF_MIN)
                            {
                                small = 0;
                                fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                            } else
                            {
                                d2 = tq != 0;
                                div_flag = r1 == 0;
                                if (d2)
                                {
                                    k[w]++;
                                    fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);
                                    fmpz_set_si(polyq[w]->coeffs + k[w], tq);
                                    mpoly_monomial_set(polyq[w]->exps + k[w]*N, texp, N);
                                }
                                c[0] = r1;
                                c[2] = c[1] = -(slong)(r1 < 0);
                            }
                        }
                    } 
                    /* quotient non-small case */
                    if (!small)
                    {
                        fmpz_fdiv_qr(q, qc, qc, mb + w);
                        d2 = !fmpz_is_zero(q);
                        div_flag = fmpz_is_zero(qc);
                        if (d2)
                        {
                            k[w]++;
                            fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);
                            fmpz_set(polyq[w]->coeffs + k[w], q);
                            mpoly_monomial_set(polyq[w]->exps + k[w]*N, texp, N);
                        }
                    }
                    if (d2)
                    {
                        if (s[w] > 1)
                        {
                            i = 1;
                            x = chains[w] + i;
                            x->i = i;
                            x->j = k[w];
                            x->p = w;
                            x->next = NULL;
                            hinds[w][x->i] = 2*(x->j + 1) + 0;
                            mpoly_monomial_add_mp(exp_list[exp_next], exp3[w] + i*N, 
                                                   polyq[w]->exps + k[w]*N, N);
                            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                        }
                        s[w] = 1;
                    }
                }
            }
            if (!div_flag)
            {
                l++;
                _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, N);
                if (small)
                {
                    fmpz_set_signed_uiuiui(p2 + l, c[2], c[1], c[0]);
                    fmpz_neg(p2 + l, p2 + l);
                } else
                {
                    fmpz_neg(p2 + l, qc);
                }
                mpoly_monomial_set(e2 + l*N, exp, N);
            }
        } 
    }

cleanup2:

   for (i = 0; i < len; i++)
      _fmpz_mpoly_set_length(polyq[i], k[i] + 1, ctx); 

   for (i = 0; i < len; i++)
      fmpz_clear(mb + i);
   fmpz_clear(qc);
   fmpz_clear(q);

   (*polyr) = p2;
   (*expr) = e2;
   
   TMP_END;

   return l + 1;

exp_overflow:
    for (i = 0; i < l; i++)
       _fmpz_demote(p2 + i);
    for (w = 0; w < len; w++)
    {
        for (i = 0; i < k[w]; i++)
           _fmpz_demote(polyq[w]->coeffs + i);
        k[w] = -WORD(1);
    }
    l = -WORD(2);
    goto cleanup2;

}

/* Assumes divisor polys don't alias any output polys */
void fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** q, fmpz_mpoly_t r,
    const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len,
                                                    const fmpz_mpoly_ctx_t ctx)
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
         flint_throw(FLINT_DIVZERO,
                   "Divide by zero in fmpz_mpoly_divrem_ideal_monagan_pearce");

      len3 = FLINT_MAX(len3, poly3[i]->length);
   }

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

      lenr = _fmpz_mpoly_divrem_ideal_monagan_pearce(q, &tr->coeffs, &tr->exps,
                         &tr->alloc, poly2->coeffs, exp2, poly2->length,
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
