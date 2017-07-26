/*
    Copyright (C) 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "longlong.h"

/*
   Set poly1 to poly2/poly3 if the division is exact, and return the length
   of the quotient. Otherwise return 0. This version of the function assumes
   the exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input polys
   are nonzero. Implements "Polynomial division using dynamic arrays, heaps
   and packed exponents" by Michael Monagan and Roman Pearce [1], except that
   we divide from right to left and use a heap with smallest exponent at head.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf 
*/
slong _fmpz_mpoly_divides_monagan_pearce1(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                const fmpz * poly3, const ulong * exp3, slong len3, slong bits,
                                                                  ulong maskhi)
{
   slong i, k, s;
   slong next_free, Q_len = 0, reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x, * x2;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong exp, maxexp = exp2[len2 - 1];
   ulong c[3]; /* for accumulating coefficients */
   int first, d1, d2;
   ulong mask = 0, ub;
   fmpz_t mb, qc, r;
   int small;
   slong bits2, bits3;
   TMP_INIT;

   TMP_START;

   fmpz_init(mb);
   fmpz_init(qc);
   fmpz_init(r);

   /* whether intermediate computations q - a*b will fit in three words */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* allow one bit for sign, one bit for subtraction */
   small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) +
           FLINT_BITS - 2) && FLINT_ABS(bits3) <= FLINT_BITS - 2;

   heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));
   /* space for pointers to heap nodes which can be reused */
   reuse = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));

   /* start with no heap nodes in use */
   next_free = 0;

   /* mask with high bit set in each field of exponent vector */
   for (i = 0; i < FLINT_BITS/bits; i++)
      mask = (mask << bits) + (UWORD(1) << (bits - 1));

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);
   
   /* see description of divisor heap division in paper */
   s = len3;
   
   /* insert (-1, 0, exp2[0]) into heap */
   x = chain + next_free++;
   x->i = -WORD(1);
   x->j = 0;
   x->next = NULL;

   HEAP_ASSIGN(heap[1], exp2[0], x);

   /* precompute -c_0 where c_i is the i-th coeff of poly3 */
   fmpz_neg(mb, poly3 + 0);

   ub = ((ulong) FLINT_ABS(*mb)) >> 1; /* abs(poly3[0])/2 */
 
   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get exponent field of heap top */
      exp = heap[1].exp;
      
      /* check for overflow in exponent: not an exact division */
      if (mpoly_monomial_overflows1(exp, mask))
      {
            for (i = 0; i < k; i++)
               _fmpz_demote(p1 + i);

            k = 0;

            goto cleanup;
      }

      /* realloc output poly ready for next quotient term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

      /* whether we are on first heap node for this exponent */
      first = 1;
      
      /* whether current exponent is divisible by exp3[0] */
      d1 = 0;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && heap[1].exp == exp)
      {
         /* pop chain from heap */
         x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
         
         /* if first heap node for this exp, check it's divisible by exp3[0] */
         if (first)
         {
            d1 = mpoly_monomial_divides1(e1 + k, exp, exp3[0], mask);

            first = 0; 
         }

         /* if accumulated coeffs will fit in three words */
         if (small)
         {
            if (x->i == -WORD(1))
            {
               /* subtract poly2[j] from accumulated three word coeff */
               _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
            } else
            {
               /* subtract poly3[i]*q[j] from accumulated three word coeff */
               _fmpz_mpoly_submul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
            }
         } else /* accumulated coeffs are multiprecision */
         {
            if (x->i == -WORD(1))
            {
               /* subtract poly2[j] from accumulated coeff */
               fmpz_sub(qc, qc, poly2 + x->j);
            } else
            {
               /* subtract poly3[i]*q[j] from accumulated coeff */
               fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
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
            /* if accumulated coeffs will fit in three words */
            if (small)
            {
               if (x->i == -WORD(1))
               {
                  /* subtract poly2[j] from accumulated three word coeff */
                  _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
               } else
               {
                  /* subtract poly3[i]*q[j] from accum. three word coeff */
                  _fmpz_mpoly_submul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
               }
            } else /* accumulated coeffs are multiprecision */
            {
               if (x->i == -WORD(1))
               {
                  /* subtract poly2[j] from accumulated coeff */
                     fmpz_sub(qc, qc, poly2 + x->j);
               } else
               {
                  /* subtract poly3[i]*q[j] from accum. coeff */
                  fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
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
         /* take node from store */
         x = Q[--Q_len];
     
         if (x->i == -WORD(1))
         {
            x->j++;
            x->next = NULL;

            /* insert (x->i, x->j + 1, exp2[x->j]) in heap */
            _mpoly_heap_insert1(heap, exp2[x->j], x, &heap_len, maskhi);
         } else if (x->j < k - 1)
         {
            x->j++;
            x->next = NULL;

            /* insert (x->i, x->j + 1, exp3[x->j] + e1[x->j]) in heap */
            _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x, &heap_len,
                                                                       maskhi);
         } else if (x->j == k - 1)
         {
            s++;

            /* node x no longer needed, designate for reuse */
            reuse[reuse_len++] = x;
         }
      }

      /* if accumulated coeff is zero, no output coeff to be written */
      if ((small && (c[2] == 0 && c[1] == 0 && c[0] == 0)) ||
          (!small && fmpz_is_zero(qc)))
         k--;
      else /* accumulated coeff is nonzero */
      {
         /* if accumulated coeff fit in three words */
         if (small)
         {
            ulong d[3];

            /* set d = abs(c) */
            if (0 > (slong) c[2])
               mpn_neg(d, c, 3);
            else
               flint_mpn_copyi(d, c, 3);

            /* check quotient of accumulated coeff by -poly3[0] is small */
            if (d[2] != 0 || ub <= d[1] || 
               (ub == 0 && 0 > (slong) d[0])) /* quotient not a small */
            {
               /* convert three words to multiprecision value */
               fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

               /* continue in non-small case from now on */
               fmpz_fdiv_qr(p1 + k, r, qc, mb); /* quotient and remainder */

               /* coeff division exact if remainder zero */
               d2 = fmpz_is_zero(r);

               small = 0;
            } else /* quotient may fit a small */
            {
               ulong q, r1;

               /* write out quotient and compute remainder */
               sdiv_qrnnd(q, r1, c[1], c[0], *mb);
               
               /* check quotient really fit a small */
               if (!COEFF_IS_MPZ(FLINT_ABS((slong) q)))
               {
                  _fmpz_demote(p1 + k);
                  p1[k] = q;
               }
               else
               {
                  fmpz_set_si(p1 + k, q);
                  small = 0;
               }

               /* coeff division exact if remainder zero */
               d2 = r1 == 0;
            }
         } else /* multiprecision case */
         {
            /* write out quotient and compute remainder */
            fmpz_fdiv_qr(p1 + k, r, qc, mb);

            /* coeff division exact if remainder zero */
            d2 = fmpz_is_zero(r);
         }

         /* if coeffs or monomials don't divide or exponent too large */
         if (!d1 || !d2 || (exp^maskhi) < (maxexp^maskhi)) /* not an exact division */
         {
            for (i = 0; i <= k; i++)
               _fmpz_demote(p1 + i);

            k = 0;

            goto cleanup;
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

            /* insert (i, k, exp3[i] + e1[k]) in heap */
            _mpoly_heap_insert1(heap, exp3[i] + e1[k], x2, &heap_len, maskhi);
         }

         s = 1;
      } 
      
      /* zero temporary accumulator */
      fmpz_zero(qc);  
   }

   k++;

cleanup:

   fmpz_clear(mb);
   fmpz_clear(qc);
   fmpz_clear(r);

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   /* return length of quotient, or zero if division not exact */
   return k;
}

/*
   Set poly1 to poly2/poly3 if the division is exact, and return the length
   of the quotient. Otherwise return 0. This version allows exponent vectors
   that each fit in "N" word. The exponent vectors are assumed to have fields
   with the given number of bits. Assumes input polys are nonzero. Implements
   "Polynomial division using dynamic arrays, heaps and packed exponents" by
   Michael Monagan and Roman Pearce [1], except that we divide from right to
   left and use a heap with smallest exponent at head.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf 
*/
slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
       const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N,
                                                    ulong maskhi, ulong masklo)
{
   slong i, k, s;
   slong next_free, Q_len = 0;
   slong reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x, * x2;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong * exp, * exps;
   ulong ** exp_list;
   ulong c[3]; /* for accumulating coefficients */
   slong exp_next;
   int first, d1, d2;
   ulong mask = 0, ub;
   fmpz_t mb, qc, r;
   int small;
   slong bits2, bits3;
   TMP_INIT;

   /* if exponent vectors are all one word, call specialised version */
   if (N == 1)
      return _fmpz_mpoly_divides_monagan_pearce1(poly1, exp1, alloc,
                           poly2, exp2, len2, poly3, exp3, len3, bits, maskhi);

   TMP_START;

   fmpz_init(mb);
   fmpz_init(qc);
   fmpz_init(r);

   /* whether intermediate computations q - a*b will fit in three words */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* allow one bit for sign, one bit for subtraction */
   small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + FLINT_BITS - 2) &&
           FLINT_ABS(bits3) <= FLINT_BITS - 2;

   heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));
   /* space for pointers to heap nodes which can be reused */
   reuse = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));
   /* array of exponent vectors, each of "N" words */
   exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
   /* list of pointers to available exponent vectors */
   exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));

   /* set up list of available exponent vectors */
   for (i = 0; i < len3; i++)
      exp_list[i] = exps + i*N;

   /* start with no heap nodes or exponents in use */
   next_free = 0;
   exp_next = 0;

   /* mask with high bit set in each word of each field of exponent vector */
   for (i = 0; i < FLINT_BITS/bits; i++)
      mask = (mask << bits) + (UWORD(1) << (bits - 1));

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);
   
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

   /* precompute -c_0 where c_i is the i-th coeff of poly3 */
   fmpz_neg(mb, poly3 + 0);

   ub = ((ulong) FLINT_ABS(*mb)) >> 1; /* abs(poly3[0])/2 */
   
   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get pointer to exponent field of heap top */
      exp = heap[1].exp;

      /* check for overflow in exponent: not an exact division */
      if (mpoly_monomial_overflows(exp, N, mask))
      {
            for (i = 0; i < k; i++)
               _fmpz_demote(p1 + i);

            k = 0;

            goto cleanup;
      }
      
      /* realloc output poly ready for next quotient term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

      /* whether we are on first heap node for this exponent */
      first = 1;
      
      /* whether current exponent is divisible by exp3[0] */
      d1 = 0;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         /* put pointer to exponent on heap top into list of available exps */
         exp_list[--exp_next] = heap[1].exp;

         /* pop chain from heap */
         x = _mpoly_heap_pop(heap, &heap_len, N, maskhi, masklo);
         
         /* if first heap node for this exp, check it's divisible by exp3[0] */
         if (first)
         {
            d1 = mpoly_monomial_divides(e1 + k*N, exp, exp3, N, mask);

            first = 0; 
         }

         /* if accumulated coeffs will fit in three words */
         if (small)
         {
            if (x->i == -WORD(1))
            {
               /* subtract poly2[j] from accumulated three word coeff */
               _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
            } else
            {
               /* subtract poly3[i]*q[j] from accumulated three word coeff */
               _fmpz_mpoly_submul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
            }
         } else
         {
            if (x->i == -WORD(1))
            {
               /* subtract poly2[j] from accumulated coeff */
               fmpz_sub(qc, qc, poly2 + x->j);
            } else
            {
               /* subtract poly3[i]*q[j] from accumulated coeff */
               fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
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
            /* if accumulated coeffs will fit in three words */
            if (small)
            {
               if (x->i == -WORD(1))
               {
                  /* subtract poly2[j] from accumulated three word coeff */
                  _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
               } else
               {
                  /* subtract poly3[i]*q[j] from accum. three word coeff */
                  _fmpz_mpoly_submul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
               }
            } else /* accumulated coeffs are multiprecision */
            {
               if (x->i == -WORD(1))
               {
                  /* subtract poly2[j] from accumulated coeff */
                  fmpz_sub(qc, qc, poly2 + x->j);
               } else
               {
                  /* subtract poly3[i]*q[j] from accum. coeff */
                  fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
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
         /* take node from store */
         x = Q[--Q_len];
     
         if (x->i == -WORD(1))
         {
            x->j++;
            x->next = NULL;

            mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);

            /* insert (x->i, x->j + 1, exp2[x->j]) in heap */
            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len,
                                                            N, maskhi, masklo))
               exp_next--;
         } else if (x->j < k - 1)
         {
            x->j++;
            x->next = NULL;

            mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N, e1 + x->j*N, N);

            /* insert (x->i, x->j + 1, exp3[x->j] + e1[x->j]) in heap */
            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len,
                                                            N, maskhi, masklo))
               exp_next--;
         } else if (x->j == k - 1)
         {
            s++;
            
            /* node x no longer needed, designate for reuse */
            reuse[reuse_len++] = x;
         }
      }

      /* if accumulated coeff is zero, no output coeff to be written */
      if ((small && (c[2] == 0 && c[1] == 0 && c[0] == 0))
            || (!small && fmpz_is_zero(qc)))
         k--;
      else
      {
         /* if accumulated coeff fit in three words */
         if (small)
         {
            ulong d[3];

            /* set d = abs(c) */
            if (0 > (slong) c[2])
               mpn_neg(d, c, 3);
            else
               flint_mpn_copyi(d, c, 3);

            /* check quotient of accumulated coeff by -poly3[0] is small */
            if (d[2] != 0 || ub <= d[1] || (ub == 0 && 0 > (slong) d[0])) /* quotient not a small */
            {
               /* convert three words to multiprecision value */
               fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

               /* continue in non-small case from now on */
               fmpz_fdiv_qr(p1 + k, r, qc, mb); /* quotient and remainder */

               /* coeff division exact if remainder zero */
               d2 = fmpz_is_zero(r);

               small = 0;
            } else /* quotient may fit a small */
            {
               ulong q, r1;

               sdiv_qrnnd(q, r1, c[1], c[0], *mb);
               
               /* check quotient really fit a small */
               if (!COEFF_IS_MPZ(FLINT_ABS((slong) q)))
               {
                  _fmpz_demote(p1 + k);
                  p1[k] = q;
               }
               else
               {
                  fmpz_set_si(p1 + k, q);
                  small = 0;
               }

               /* coeff division exact if remainder zero */
               d2 = r1 == 0;
            }
         } else /* multiprecision case */
         {
            /* write out quotient and compute remainder */
            fmpz_fdiv_qr(p1 + k, r, qc, mb);

            /* coeff division exact if remainder zero */
            d2 = fmpz_is_zero(r);
         }

         /* if coeffs or monomials don't divide, or exponent too large */
         if (!d1 || !d2 ||
          mpoly_monomial_gt(exp, exp2 + (len2 - 1)*N, N, maskhi, masklo)) /* inexact division */
         {
            for (i = 0; i <= k; i++)
               _fmpz_demote(p1 + i);

            k = 0;

            goto cleanup;
         }

         /* see paper */
         for (i = 1; i < s; i++)
         {
            if (reuse_len != 0)
               x2 = reuse[--reuse_len];
            else
               x2 = chain + next_free++;
            
            x2->i = i;
            x2->j = k;
            x2->next = NULL;

            mpoly_monomial_add(exp_list[exp_next], exp3 + i*N, e1 + k*N, N);

            /* insert (i, k, exp3[i] + e1[k]) in heap */
            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x2, &heap_len,
                                                            N, maskhi, masklo))
               exp_next--;
         }

         s = 1;
      } 
      
      /* zero temporary accumulator */
      fmpz_zero(qc);  
   }

   k++;

cleanup:

   fmpz_clear(mb);
   fmpz_clear(qc);
   fmpz_clear(r);

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   /* return length of quotient, or zero if division not exact */
   return k;
}

/* return 1 if quotient is exact */
int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0;
   ulong * max_degs2, * max_degs3;
   ulong max = 0;
   ulong maskhi, masklo;
   ulong * exp2 = poly2->exps, * exp3 = poly3->exps, * expq;
   int free2 = 0, free3 = 0;
   ulong mask = 0;
   TMP_INIT;

   /* check divisor is nonzero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divides_monagan_pearce");

   /* dividend zero, write out quotient */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return 1;
   }

   TMP_START;

   max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
   max_degs3 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   /* compute maximum degree appearing in inputs and outputs */
   
   fmpz_mpoly_max_degrees(max_degs2, poly2, ctx);
   fmpz_mpoly_max_degrees(max_degs3, poly3, ctx);

   for (i = 0; i < ctx->n; i++)
   {
      if (max_degs2[i] > max)
         max = max_degs2[i];

      /* cannot be exact division if poly2 degrees less than those of poly3 */
      if (max_degs2[i] < max_degs3[i])
      {
         len = 0;

         goto cleanup;
      }
   }

   /* compute number of bits required for exponent fields */
   bits = FLINT_BIT_COUNT(max);

   exp_bits = 8;
   while (bits >= exp_bits)
      exp_bits *= 2;

   exp_bits = FLINT_MAX(exp_bits, poly2->bits);
   exp_bits = FLINT_MAX(exp_bits, poly3->bits);

   masks_from_bits_ord(maskhi, masklo, exp_bits, ctx->ord);
   /* number of words required for exponent vectors */
   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   /* temporary space to check leading monomials divide */
   expq = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   /* quick check for easy case of inexact division of leading monomials */
   if (poly2->bits == poly3->bits && N == 1 && 
       poly2->exps[0] < poly3->exps[0])
   {
      goto cleanup;
   }

   /* ensure input exponents packed to same size as output exponents */
   if (exp_bits > poly2->bits)
   {
      free2 = 1;
      exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
      mpoly_unpack_monomials_noalloc(exp2, exp_bits, poly2->exps, poly2->bits,
                                                        poly2->length, ctx->n);
   }

   if (exp_bits > poly3->bits)
   {
      free3 = 1;
      exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
      mpoly_unpack_monomials_noalloc(exp3, exp_bits, poly3->exps, poly3->bits,
                                                        poly3->length, ctx->n);
   }

   /* mask with high bit of each exponent vector field set */
   for (i = 0; i < FLINT_BITS/exp_bits; i++)
      mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));

   /* check leading monomial divides exactly */
   if (!mpoly_monomial_divides(expq, exp2, exp3, N, mask))
   {
      len = 0;

      goto cleanup;
   }

   /* deal with aliasing and divide polynomials */
   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      len = _fmpz_mpoly_divides_monagan_pearce(&temp->coeffs, &temp->exps,
                            &temp->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                               maskhi, masklo);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      len = _fmpz_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                               maskhi, masklo);
   }

cleanup:

   _fmpz_mpoly_set_length(poly1, len, ctx);

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   TMP_END;

   /* division is exact if len is nonzero */
   return (len != 0);
}
