/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include <assert.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "longlong.h"


slong _fmpz_mpoly_quasidivrem_heapV2(fmpz_t scale, slong * lenr,
  fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr,
                  ulong ** expr, slong * allocr, const fmpz * poly2,
   const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, 
                   slong len3, slong bits, slong N, ulong maskhi, ulong masklo)
{
    slong a, b, i, j, k, l, s;
    slong next_loc;
    slong next_free, Q_len = 0;
    slong reuse_len = 0, heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    mpoly_heap_t ** Q, ** reuse;
    mpoly_heap_t * x, * x2;
    fmpz * p1 = *polyq;
    fmpz * p2 = *polyr;
    ulong * e1 = *expq;
    ulong * e2 = *expr;
    ulong * exp, * exps;
    ulong ** exp_list;
    ulong c[3]; /* for accumulating coefficients */
    slong exp_next;
    ulong mask = 0, ub;
    fmpz_t mb, sf, gcd, qc, r, tp;
    int small;
    slong bits2, bits3;
    int divides, scaleis1;
    fmpz * rs, * qs, * qsf[32], * qsu, * qsr;
    slong qs_alloc, rs_alloc, qsf_alloc[32], ta;
    

slong maxiter = 0;
slong itercnt = 1;
slong avgiter = 0;

    TMP_INIT;


    /* if exponent vectors fit in one word, call specialised version */

    TMP_START;

    fmpz_init(tp);
    fmpz_init(sf);
    fmpz_init(gcd);
    fmpz_set_ui(scale, 1);
    scaleis1 = 1;

    fmpz_init(mb);
    fmpz_init(qc);
    fmpz_init(r);

    qs_alloc = 64;
    qs  = (fmpz *) flint_calloc(qs_alloc, sizeof(fmpz));
    ta = qs_alloc;
    for (i=0; i<32; i++) {
        qsf[i] = (fmpz *) flint_calloc(ta, sizeof(fmpz));
        qsf_alloc[i] = ta;
        ta = FLINT_MAX(1, ta/2);
    }
    qsr = (fmpz *) flint_calloc(qs_alloc, sizeof(fmpz));
    qsu = (long *) flint_calloc(qs_alloc, sizeof(long));

    rs_alloc = 64;
    rs = (fmpz *) flint_calloc(rs_alloc, sizeof(fmpz));
    

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) +
           FLINT_BITS - 2) && FLINT_ABS(bits3) <= FLINT_BITS - 2;

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    /* alloc array of heap nodes which can be chained together */
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    /* space for temporary storage of pointers to heap nodes */
    Q = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));
    /* space for pointers to heap nodes which can be reused */
    reuse = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));
    /* array of exponents of N words each */
    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    /* array of pointers to unused exponent vectors */
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    /* space to save copy of current exponent vector */
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    /* set up list of available exponent vectors */
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    /* start with no heap nodes and no exponent vectors in use */
    next_free = 0;
    exp_next = 0;

    /* mask with high bit set in each field of exponent vector */
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    /* quotient and remainder poly indices start at -1 */
    k = -WORD(1);
    l = -WORD(1);

    /* see description of divisor heap division in paper */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + next_free++;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    mpoly_monomial_set(heap[1].exp, exp2, N);

    /* precompute -c_n where c_n is the leading coeff of poly3 */
    fmpz_neg(mb, poly3);

/*
flint_printf("numberator length: %wd  denominator length: %wd  ",len2,len3);
printf("leading coeff: "); fmpz_print(mb); printf("\n");
*/
    ub = ((ulong) FLINT_ABS(*mb)) >> 1; /* abs(lc(poly3))/2 */

    /* while heap is nonempty */
    while (heap_len > 1)
    {
        /* make temporary copy of exponent at top of heap */
        mpoly_monomial_set(exp, heap[1].exp, N);

        /* check there has been no overflow */
        if (mpoly_monomial_overflows(exp, N, mask))
        {
            for (i = 0; i < k; i++)
                _fmpz_demote(p1 + i);
            for (i = 0; i < l; i++)
                _fmpz_demote(p2 + i);

            k = 0;
            l = 0;

            goto cleanup2;
        }
      
        /* realloc quotient poly, for next quotient term */
        k++;
        _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, N);
        if (k + 1 > qs_alloc)
        {
            slong len;
            len = 2*qs_alloc;
            qs  = (fmpz *) flint_realloc(qs, len*sizeof(fmpz));
            qsr = (fmpz *) flint_realloc(qsr, len*sizeof(fmpz));
            qsu = (long *) flint_realloc(qsu, len*sizeof(long));
            flint_mpn_zero((mp_ptr) (qs  + qs_alloc), len - qs_alloc);
            flint_mpn_zero((mp_ptr) (qsr + qs_alloc), len - qs_alloc);
            ta = len;
            for (i=0; i<32; i++) {
                qsf[i] = (fmpz *) flint_realloc(qsf[i], ta*sizeof(fmpz));
                flint_mpn_zero((mp_ptr) (qsf[i] + qsf_alloc[i]), ta - qsf_alloc[i]);
                qsf_alloc[i] = ta;
                ta = FLINT_MAX(1, ta/2);
            }
            qs_alloc = len;
        }


        /* set temporary coeff to zero */
        c[0] = c[1] = c[2] = 0;
        fmpz_zero(qc);  

        /* while heap nonempty and contains chain with current output exponent */
        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
        {
            /* put pointer to exponent at heap top into list of available exps */
            exp_list[--exp_next] = heap[1].exp;

            /* pop chain from heap */
            x = _mpoly_heap_pop(heap, &heap_len, N, maskhi, masklo);

            if (small)
            {
                if (x->i == -WORD(1)) {
                   _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
                } else
                {
                    assert(k>0);
                    qsu[x->j] = k - 1;
                    _fmpz_mpoly_addmul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
                }
            }
            else
            {
                if (scaleis1)
                {
                    if (x->i == -WORD(1))
                    {
                        fmpz_sub(qc, qc, poly2 + x->j);
                    }
                    else 
                    {
                        assert(k>0);
                        qsu[x->j] = k - 1;
                        fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
                    }
                }
                else
                {
                    if (x->i == -WORD(1))
                        fmpz_submul(qc, scale, poly2 + x->j);
                    else
                    {

fmpz_set(sf, qsr + x->j);

if (k - qsu[x->j] > maxiter)
    maxiter = k - qsu[x->j];
itercnt++;
avgiter += k - qsu[x->j];


/*
for (a = qsu[x->j] + 1; a < k; a++)
{
    fmpz_mul(sf, sf, qsf[0] + a);
}
*/


a = qsu[x->j] + 1;
b = k;
i=0;
while (a<b) {
    if (a&1) {
        fmpz_mul(sf,sf, qsf[i]+a);
        a++;
    }
    if (b&1) {
        b--;
        fmpz_mul(sf,sf, qsf[i]+b);
    }
    a=a/2;
    b=b/2;
    i++;
}


qsu[x->j] = k - 1;
fmpz_set(qsr + x->j, sf);


/*
fmpz_divexact(tp, scale, qs + x->j);
assert(fmpz_equal(sf, tp));
*/

                        fmpz_mul(sf, sf, p1 + x->j);
                        fmpz_addmul(qc, poly3 + x->i, sf);

                    }
                }
            }

            /* temporarily store pointer to this node, or designate for reuse */
            if (x->i != -WORD(1) || x->j < len2 - 1)
                Q[Q_len++] = x;
            else
                reuse[reuse_len++] = x;

            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
                if (small)
                {
                    if (x->i == -WORD(1))
                       _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
                    else
                       _fmpz_mpoly_addmul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
                }
                else
                {
                    if (scaleis1)
                    {

                        if (x->i == -WORD(1))
                            fmpz_sub(qc, qc, poly2 + x->j);
                        else
                            fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
                    }
                    else
                    {
                        if (x->i == -WORD(1))
                            fmpz_submul(qc, scale, poly2 + x->j);
                        else
                        {
                            fmpz_divexact(tp, scale, qs + x->j);
                            fmpz_mul(tp, tp, p1 + x->j);
                            fmpz_addmul(qc, poly3 + x->i, tp);
                        }
                    }
                }
                /* temporarily store pointer to node, or designate for reuse */
                if (x->i != -WORD(1) || x->j < len2 - 1)
                    Q[Q_len++] = x;
                else
                    reuse[reuse_len++] = x;
            }
        }

        /* for each node temporarily stored */
        while (Q_len > 0)
        {
            x = Q[--Q_len];
            if (x->i == -WORD(1))
            {
                x->j++;
                x->next = NULL;
                mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, maskhi, masklo))
                    exp_next--;
            } else if (x->j < k - 1)
            {
                x->j++;
                x->next = NULL;
                mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N, e1 + x->j*N, N);
                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, maskhi, masklo))
                    exp_next--;
            } else if (x->j == k - 1)
            {
                s++;
                reuse[reuse_len++] = x;
            }
        }

        /* check current exp divisible by leading exp of poly3... */
        divides = mpoly_monomial_divides(e1 + k*N, exp, exp3, N, mask);


        /* if accumulated coeff is zero, no output coeffs to be written */
        if (small)
        {
            if (c[2] == 0 && c[1] == 0 && c[0] == 0)
            {
                k--;
                continue;
            }
            else if (!divides)
            {
                l++;
                _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, N);

                if (l + 1 > rs_alloc)
                {
                    slong len;
                    len = FLINT_MAX(l + 1, 2*rs_alloc);
                    rs = (fmpz *) flint_realloc(rs, len*sizeof(fmpz));
                    if (len > rs_alloc)
                        flint_mpn_zero((mp_ptr) (rs + rs_alloc), len - rs_alloc);
                    rs_alloc = len;
                }


                fmpz_set_signed_uiuiui(p2 + l, c[2], c[1], c[0]);
                fmpz_neg(p2 + l, p2 + l);
                mpoly_monomial_set(e2 + l*N, exp, N);
                fmpz_set_ui(rs + l, 1);
                k--;
                continue;
            }
            else
            {
                ulong d[3];
                if (0 > (slong) c[2])
                    mpn_neg(d, c, 3);
                else
                    flint_mpn_copyi(d, c, 3);

                if (d[2] != 0 || ub <= d[1] || (ub == 0 && 0 > (slong) d[0]))
                {
                    /* upgrade to multiprecision accumulated coeffs */
                    fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                    small = 0;
                } else /* quotient maybe fits a small */
                {
                    ulong q, r1;

                    /* compute quotient and remainder coeff */
                    sdiv_qrnnd(q, r1, c[1], c[0], *mb);

                    if (COEFF_IS_MPZ(FLINT_ABS(q))) /* quotient too large */
                    {
                        /* upgrade to multiprecision accumulated coeffs */
                        fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                        small = 0;
                    }
                    else if (r1 != 0)
                    {
                        /* upgrade to multiprecision accumulated coeffs */
                        fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                        small = 0;
                    }
                    else
                    {
                        fmpz_set_si(p1 + k, q);            
                        fmpz_set_ui(qs + k, 1);
                        fmpz_set_ui(qsr + k, 1);
                        qsu[k] = k;

                        fmpz_set_ui(qsf[0] + k, 1);
                        fmpz_set_ui(sf, 1);

        a = k;
        i = 0;
looper1:
        fmpz_set_ui(qsf[i] + a, 1);
        if (a & 1) {
            a = a/2;
            i++;
            goto looper1;
        }


                    }
                }
            }
        }
        else /* not small */
        {
            if (fmpz_is_zero(qc))
            {
                k--;
                continue;
            }
            else if (!divides)
            {
                l++;
                _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, N);

                if (l + 1 > rs_alloc)
                {
                    slong len;
                    len = FLINT_MAX(l + 1, 2*rs_alloc);
                    rs = (fmpz *) flint_realloc(rs, len*sizeof(fmpz));
                    if (len > rs_alloc)
                        flint_mpn_zero((mp_ptr) (rs + rs_alloc), len - rs_alloc);
                    rs_alloc = len;
                }

                fmpz_neg(p2 + l, qc); 
                mpoly_monomial_set(e2 + l*N, exp, N);
                fmpz_set(rs + l, scale);                
                k--;
                continue;
            }
        }

        if (!small)
        {
            fmpz_gcd(gcd, qc, mb);
            fmpz_abs(sf, mb);
            fmpz_divexact(sf, sf, gcd);
            if (!fmpz_is_one(sf))
                scaleis1 = 0;
            fmpz_mul(p1 + k, sf, qc);
            fmpz_divexact(p1 + k, p1 + k, mb);
            fmpz_mul(scale, scale, sf);

            fmpz_set(qs + k, scale);
            fmpz_set_ui(qsr + k, 1);
            qsu[k] = k;
/*
            fmpz_set(qsf[0] + k, sf);
*/
        a = k;
        i = 0;
looper:
        fmpz_set(qsf[i] + a, sf);
        if (a & 1) {
            fmpz_mul(sf, sf, qsf[i]+a-1);
            a = a/2;
            i++;
            goto looper;
        }

        }




        /* see paper */
        for (i = 1; i < s; i++)
        {
            /* get an empty node, from reuse array if possible */
            if (reuse_len != 0)
                x2 = reuse[--reuse_len];
            else
                x2 = chain + next_free++;
            x2->i = i;
            x2->j = k;
            x2->next = NULL;
            mpoly_monomial_add(exp_list[exp_next], exp3 + i*N, e1 + k*N, N);
            /* insert (i, k, exp3[i] + e1[k]) */
            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x2,
                                      &next_loc, &heap_len, N, maskhi, masklo))
                exp_next--;
        }
        s = 1;   
    }

    k++;
    l++;

cleanup2:

    fmpz_clear(mb);
    fmpz_clear(qc);
    fmpz_clear(r);

    (*polyq) = p1;
    (*expq) = e1;
    (*polyr) = p2;
    (*expr) = e2;

    /* set remainder poly length */
    (*lenr) = l;

    TMP_END;


    fmpz_clear(tp);
    fmpz_clear(gcd);
    fmpz_clear(sf);

    for (i = 0; i < k; i++)
    {
        fmpz_divexact(tp, scale, qs + i);
        fmpz_mul(p1 + i, p1 + i, tp);
    }
    for (i = 0; i < qs_alloc; i++) {
        fmpz_clear(qs + i);
        fmpz_clear(qsr + i);
    }
    for (j=0; j<32; j++) {
        for (i=0; i<qsf_alloc[j]; i++) {
            fmpz_clear(qsf[j] + i);
        }
        flint_free(qsf[j]);
    }
    flint_free(qs);
    flint_free(qsr);
    flint_free(qsu);


    for (i = 0; i < l; i++)
    {
        fmpz_divexact(tp, scale, rs + i);
        fmpz_mul(p2 + i, p2 + i, tp);
    }
    for (i = 0; i < rs_alloc; i++)
        fmpz_clear(rs + i);
    flint_free(rs);



/*
flint_printf("maxiter: %wd  avg iter: %f\n",maxiter,(double)(avgiter)/(double)(itercnt));
printf("final scale: "); fmpz_print(scale); printf("\n");
*/

    /* return quotient poly length */
    return k;
}

void fmpz_mpoly_quasidivrem_heapV2(fmpz_t scale, fmpz_mpoly_t q, fmpz_mpoly_t r,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong exp_bits, N, lenq = 0, lenr = 0;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    ulong maskhi, masklo;
    int free2 = 0, free3 = 0;
    fmpz_mpoly_t temp1, temp2;
    fmpz_mpoly_struct * tq, * tr;


    /* check divisor is nonzero */
    if (poly3->length == 0)
        flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_quasidivrem_heap");

    fmpz_set_ui(scale, 1);

    /* dividend zero, write out quotient and remainder */
    if (poly2->length == 0)
    {
        fmpz_mpoly_zero(q, ctx);
        fmpz_mpoly_zero(r, ctx);
        return;
    }

    /* compute maximum degree appearing in inputs */

    /* maximum bits in quotient and remainder exps is max for poly2 and poly3 */
    exp_bits = FLINT_MAX(poly2->bits, poly3->bits);

    masks_from_bits_ord(maskhi, masklo, exp_bits, ctx->ord);
    /* number of words required for exponent vectors */
    N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

    /* ensure input exponents packed to same size as output exponents */
    if (exp_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_unpack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                        poly2->length, ctx->n);
    }

    if (exp_bits > poly3->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
        mpoly_unpack_monomials(exp3, exp_bits, poly3->exps, poly3->bits,
                                                        poly3->length, ctx->n);
    }

    /* check divisor leading monomial is at most that of the dividend */
    if (mpoly_monomial_lt(exp3, exp2, N, maskhi, masklo))
    {
        fmpz_mpoly_set(r, poly2, ctx);
        fmpz_mpoly_zero(q, ctx);
        goto cleanup3;
    }

   /* take care of aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_init2(temp1, FLINT_MAX(poly2->length/poly3->length + 1, 1),
                                                                          ctx);
      fmpz_mpoly_fit_bits(temp1, exp_bits, ctx);
      temp1->bits = exp_bits;

      tq = temp1;
   } else
   {
      fmpz_mpoly_fit_length(q, FLINT_MAX(poly2->length/poly3->length + 1, 1),
                                                                          ctx);
      fmpz_mpoly_fit_bits(q, exp_bits, ctx);
      q->bits = exp_bits;

      tq = q;
   }

   if (r == poly2 || r == poly3)
   {
      fmpz_mpoly_init2(temp2, poly3->length, ctx);
      fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);
      temp2->bits = exp_bits;

      tr = temp2;
   } else
   {
      fmpz_mpoly_fit_length(r, poly3->length, ctx);
      fmpz_mpoly_fit_bits(r, exp_bits, ctx);
      r->bits = exp_bits;

      tr = r;
   }

   /* do division with remainder */
   while ((lenq = _fmpz_mpoly_quasidivrem_heapV2(scale, &lenr, &tq->coeffs, &tq->exps,
         &tq->alloc, &tr->coeffs, &tr->exps, &tr->alloc, poly2->coeffs, exp2, 
         poly2->length, poly3->coeffs, exp3, poly3->length, exp_bits,
                                                       N, maskhi, masklo)) == 0
         && lenr == 0 && exp_bits < FLINT_BITS)
   {
      ulong * old_exp2 = exp2, * old_exp3 = exp3;

      masks_from_bits_ord(maskhi, masklo, 2*exp_bits, ctx->ord);
      N = (2*exp_bits*ctx->n - 1)/FLINT_BITS + 1;

      exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
      mpoly_unpack_monomials(exp2, 2*exp_bits, old_exp2, exp_bits,
                                                        poly2->length, ctx->n);

      exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
      mpoly_unpack_monomials(exp3, 2*exp_bits, old_exp3, exp_bits,
                                                        poly3->length, ctx->n);

      exp_bits *= 2;

      if (free2)
         flint_free(old_exp2);

      if (free3)
         flint_free(old_exp3);

      free2 = free3 = 1; 

      fmpz_mpoly_fit_bits(tq, exp_bits, ctx);
      tq->bits = exp_bits;

      fmpz_mpoly_fit_bits(tr, exp_bits, ctx);
      tr->bits = exp_bits;
   }
   
   if (lenq == 0 && lenr == 0)
      flint_throw(FLINT_EXPOF,
                      "Exponent overflow in fmpz_mpoly_quasidivrem_heap");

   /* deal with aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_swap(temp1, q, ctx);

      fmpz_mpoly_clear(temp1, ctx);
   } 

   if (r == poly2 || r == poly3)
   {
      fmpz_mpoly_swap(temp2, r, ctx);

      fmpz_mpoly_clear(temp2, ctx);
   } 

   _fmpz_mpoly_set_length(q, lenq, ctx);
   _fmpz_mpoly_set_length(r, lenr, ctx);

cleanup3:

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);
}
