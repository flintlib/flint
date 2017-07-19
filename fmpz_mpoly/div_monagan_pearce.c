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
   Set polyq to the quotient poly2 by poly3 discarding remainder (with notional
   remainder coeffs reduced modulo the leading coeff of poly3), and return
   the length of the quotient. This version of the function assumes the
   exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input polys
   are nonzero. Implements "Polynomial division using dynamic arrays, heaps
   and packed exponents" by Michael Monagan and Roman Pearce [1], except that
   we use a heap with smallest exponent at head. Note that if a < b then
   (n - b) < (n - b) where n is the maximum value a and b can take. The word
   "maxn" is set to an exponent vector whose fields are all set to such a
   value n. This allows division from left to right with a heap with smallest
   exponent at the head. Quotient poly is written in reverse order.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf 
*/
slong _fmpz_mpoly_div_monagan_pearce1(fmpz ** polyq, ulong ** expq,
                slong * allocq, const fmpz * poly2, const ulong * exp2,
            slong len2, const fmpz * poly3, const ulong * exp3, slong len3,
                                                      slong bits, ulong maskhi)
{
   slong i, k, s;
   slong next_free, Q_len = 0;
   slong reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x, * x2;
   fmpz * p1 = *polyq;
   ulong * e1 = *expq;
   ulong exp, exp_temp, exp3_lead, tmp = 0;
   ulong c[3]; /* for accumulating coefficients */
   ulong mask = 0, ub;
   fmpz_t mb, qc;
   int small;
   slong bits2, bits3;
   int d1, divides;
   TMP_INIT;

   TMP_START;

   fmpz_init(mb);
   fmpz_init(qc);

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

   /* quotient poly indices start at -1 */
   k = -WORD(1);
   
   /* see description of divisor heap division in paper */
   s = len3;
   
   x = chain + next_free++;
   x->i = -WORD(1);
   x->j = 0;
   x->next = NULL;

   /* insert (-1, 0, exp2[0]) into heap */
   HEAP_ASSIGN(heap[1], exp2[0], x);

   /* precompute -c_n where c_n is the leading coeff of poly3 */
   fmpz_neg(mb, poly3);

   ub = ((ulong) FLINT_ABS(*mb)) >> 1; /* abs(lc(poly3))/2 */
   
   /* precompute leading exponent of poly3 */
   exp3_lead = exp3[0];

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get exponent field of heap top */
      exp = heap[1].exp;
      
      /* check there has been no overflow */
      if ((exp & mask) != 0)
      {
         for (i = 0; i < k; i++)
            _fmpz_demote(p1 + i);

         k = -WORD(1);

         goto cleanup;
      }

      /* only do arithmetic for monomials that divide */
      divides = mpoly_monomial_divides1(&tmp, exp, exp3_lead, mask);

      /* realloc quotient poly ready for next quotient term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, 1);

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && heap[1].exp == exp)
      {
         /* pop chain from heap */
         x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
         
         if (divides)
         {
            /* if accumulated coeffs will fit in three words */
            if (small)
            {
               if (x->i == -WORD(1))
               {
                  /* subtract poly2 coeff from accumulated three word coeff */
                  _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
               } else
               {
                  /* subtract q[j]*poly3 coeff from accum. three word coeff */
                  _fmpz_mpoly_submul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
               }
            } else /* accumulated coeffs are multiprecision */
            {
               if (x->i == -WORD(1))
               {
                  /* subtract poly2 coeff from accum. coeff */
                  fmpz_sub(qc, qc, poly2 + x->j);
               } else
               {
                  /* subtract q[j]*poly3 coeff from accum. coeff */
                  fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
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
            if (divides)
            {
               /* if accumulated coeffs will fit in three words */
               if (small)
               {
                  if (x->i == -WORD(1))
                  {
                     /* subtract poly2 coeff from accum. three word coeff */
                     _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
                  } else
                  {
                     /* subtract q[j]*poly3 coeff from accum. three word coeff */
                     _fmpz_mpoly_submul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
                  }
               } else /* accumulated coeffs are multiprecision */
               {
                  if (x->i == -WORD(1))
                  {
                     /* subtract poly2 coeff from accum. coeff */
                     fmpz_sub(qc, qc, poly2 + x->j);
                  } else
                  {
                     /* subtract q[j]*poly3 coeff from accum. coeff */
                     fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
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
         /* take node from store */
         x = Q[--Q_len];
     
         if (x->i == -WORD(1))
         {
            x->j++;
            x->next = NULL;

            /* insert (x->i, x->j + 1, exp2[x->j]) in heap */
            exp_temp = exp2[x->j];

            /* but only if monomial is big enough */
            if ((exp_temp^maskhi) >= (exp3_lead^maskhi))
               _mpoly_heap_insert1(heap, exp_temp, x, &heap_len, maskhi);
         } else if (x->j < k - 1)
         {
            x->j++;
            x->next = NULL;

            /* insert (x->i, x->j + 1, exp3[x->i] + e1[x->j]) in heap */
            exp_temp = exp3[x->i] + e1[x->j];

            /* but only if monomial is big enough */
            if ((exp_temp^maskhi) >= (exp3_lead^maskhi))
               _mpoly_heap_insert1(heap, exp_temp, x, &heap_len, maskhi);
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
         /* check current exp divisible by leading exp of poly3... */
         d1 = mpoly_monomial_divides1(e1 + k, exp, exp3_lead, mask);

         /* ... if not, no quotient term */
         if (!d1)
            k--;
         else /* monomial exact division, so quotient term */
         {
            /* if accumulated coeff in three words */
            if (small)
            {
               ulong d[3];

               /* compute d = abs(c) */
               if (0 > (slong) c[2])
                  mpn_neg(d, c, 3);
               else
                  flint_mpn_copyi(d, c, 3);

               /* check quotient of accumulated coeff by lc(poly3) is small */
               if (d[2] != 0 || ub <= d[1] ||
                  (ub == 0 && 0 > (slong) d[0])) /* quotient not a small */
               {
                  /* upgrade to multiprecision accumulated coeffs */
                  fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

                  small = 0;
               } else /* quotient fits a small */
               {
                  ulong q, r1;

                  /* compute quotient and discard remainder coeff */
                  sdiv_qrnnd(q, r1, c[1], c[0], *mb);

                  fmpz_set_si(p1 + k, q);

                  if (COEFF_IS_MPZ(FLINT_ABS(q))) /* quotient too large */
                  {
                     /* upgrade to multiprecision accumulated coeffs */
                     fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

                     small = 0;
                  }
               }
            } 

            /* if accumulated coeff is multiprecision, incl. from upgrade */
            if (!small)
            {
               /* write out quotient coeff */
               fmpz_fdiv_q(p1 + k, qc, mb);
            }

            /* if the quotient coeff was nonzero */
            if (!fmpz_is_zero(p1 + k))
            {
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

                  /* insert (i, k, exp3[i] + e1[k]) */
                  exp_temp = exp3[i] + e1[k];

                  /* but only if monomial is big enough */
                  if ((exp_temp^maskhi) >= (exp3_lead^maskhi))
                     _mpoly_heap_insert1(heap, exp_temp, x2, &heap_len, maskhi);
               }
               s = 1;
            } else /* quotient coeff was zero, nothing to write out */
               k--;
         }
      } 
      
      /* zero temporary accumulator */
      fmpz_zero(qc);  
   }

   k++;

cleanup:

   fmpz_clear(mb);
   fmpz_clear(qc);

   (*polyq) = p1;
   (*expq) = e1;

   TMP_END;

   /* return quotient poly length */
   return k;
}

/*
   Set polyq to the quotient poly2 by poly3 discarding remainder (with notional
   remainder coeffs reduced modulo the leading coeff of poly3), and return
   the length of the quotient. This version of the function assumes the
   exponent vectors all fit in "N" words. The exponent vectors are assumed to
   have fields with the given number of bits. Assumes input polys are nonzero.
   Implements "Polynomial division using dynamic arrays, heaps and packed
   exponents" by Michael Monagan and Roman Pearce [1], except that we use a
   heap with smallest exponent at head. Note that if a < b then
   (n - b) < (n - b) where n is the maximum value a and b can take. The word
   "maxn" is set to an exponent vector whose fields are all set to such a
   value n. This allows division from left to right with a heap with smallest
   exponent at the head. Quotient poly is written in reverse order.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf 
*/
slong _fmpz_mpoly_div_monagan_pearce(fmpz ** polyq,
           ulong ** expq, slong * allocq, const fmpz * poly2,
   const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, 
                   slong len3, slong bits, slong N, ulong maskhi, ulong masklo)
{
   slong i, k, s;
   slong next_free, Q_len = 0;
   slong reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x, * x2;
   fmpz * p1 = *polyq;
   ulong * e1 = *expq;
   ulong * exp, * exps, * texp, * tmp;
   ulong ** exp_list;
   ulong c[3]; /* for accumulating coefficients */
   slong exp_next;
   ulong mask = 0, ub;
   fmpz_t mb, qc;
   int small;
   slong bits2, bits3;
   int d1, divides;
   TMP_INIT;

   /* if exponent vectors fit in one word, call specialised version */
   if (N == 1)
      return _fmpz_mpoly_div_monagan_pearce1(polyq, expq, allocq,
                           poly2, exp2, len2, poly3, exp3, len3, bits, maskhi);

   TMP_START;

   fmpz_init(mb);
   fmpz_init(qc);

   /* whether intermediate computations q - a*b will fit in three words */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* allow one bit for sign, one bit for subtraction */
   small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + 
           FLINT_BIT_COUNT(len3) + FLINT_BITS - 2) &&
           FLINT_ABS(bits3) <= FLINT_BITS - 2;

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
   /* space for temporary exponent vectors, for exponent computations */
   texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
   texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
   tmp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
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

   /* quotient and poly index starts at -1 */
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

   /* precompute -c_n where c_n is the leading coeff of poly3 */
   fmpz_neg(mb, poly3);

   ub = ((ulong) FLINT_ABS(*mb)) >> 1; /* abs(lc(poly3))/2 */
   
   /* precompute leading exponent of poly3 */
   mpoly_monomial_set(texp, exp3, N);

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

         k = -WORD(1);

         goto cleanup2;
      }

      divides = mpoly_monomial_divides(tmp, exp, texp, N, mask);

      /* realloc quotient poly, for next quotient term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, N);

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         /* put pointer to exponent at heap top into list of available exps */
         exp_list[--exp_next] = heap[1].exp;

         /* pop chain from heap */
         x = _mpoly_heap_pop(heap, &heap_len, N, maskhi, masklo);

         if (divides)
         {
            /* if accumulated coeffs will fit in three words */
            if (small)
            {
               if (x->i == -WORD(1))
               {
                  /* subtract poly2 coeff from accumulated three word coeff */
                  _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
               } else
               {
                  /* subtract q[j]*poly3 coeff from accum. three word coeff */
                  _fmpz_mpoly_submul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
               }
            } else /* accumulated coeffs are multiprecision */
            {
              if (x->i == -WORD(1))
               {
                  /* subtract poly2 coeff from accum. coeff */
                  fmpz_sub(qc, qc, poly2 + x->j);
               } else
               {
                  /* subtract q[j]*poly3 coeff from accum. coeff */
                  fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
               }
            }
         }

         /* temporarily store pointer to this node, or designate for reuse */
         if (x->i != -WORD(1) || x->j < len2 - 1)
            Q[Q_len++] = x;
         else
            reuse[reuse_len++] = x;

         while ((x = x->next) != NULL)
         {
            if (divides)
            {
               /* if accumulated coeffs will fit in three words */
               if (small)
               {
                  if (x->i == -WORD(1))
                  {
                     /* subtract poly2 coeff from accumulated three word coeff */
                     _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
                  } else
                  {
                     /* subtract q[j]*poly3 coeff from accum. three word coeff */
                     _fmpz_mpoly_submul_uiuiui_fmpz(c, poly3[x->i], p1[x->j]);
                  }
               } else /* accumulated coeffs are multiprecision */
               {
                  if (x->i == -WORD(1))
                  {
                     /* subtract poly2 coeff from accum. coeff */
                     fmpz_sub(qc, qc, poly2 + x->j);
                  } else
                  {
                     /* subtract q[j]*poly3 coeff from accum. coeff */
                     fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
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
         /* take node from store */
         x = Q[--Q_len];
     
         if (x->i == -WORD(1))
         {
            x->j++;
            x->next = NULL;

            /* insert (x->i, x->j + 1, exp2[x->j]) in heap */
            mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);

            /* but only if the exponent is big enough */
            if (!mpoly_monomial_lt(texp, exp_list[exp_next], N, maskhi, masklo))
            {
               if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                                 &heap_len, N, maskhi, masklo))
                  exp_next--;
            }
         } else if (x->j < k - 1)
         {
            x->j++;
            x->next = NULL;

            /* insert (x->i, x->j + 1, exp3[x->i] + e1[x->j]) in heap */
            mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N, e1 + x->j*N, N);

            /* but only if exponent is big enough */
            if (!mpoly_monomial_lt(texp, exp_list[exp_next], N, maskhi, masklo))
            {
               if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len,
                                                             N, maskhi, masklo))
               exp_next--;
            }
         } else if (x->j == k - 1)
         {
            s++;
            
            /* node x no longer needed, designate for reuse */
            reuse[reuse_len++] = x;
         }
      }

      /* if accumulated coeff is zero, no quotient coeff to be written */
      if ((small && (c[2] == 0 && c[1] == 0 && c[0] == 0)) ||
         (!small && fmpz_is_zero(qc)))
         k--;
      else /* accumulated coeff is nonzero */
      {
         /* check current exp divisible by leading exp of poly3... */
         d1 = mpoly_monomial_divides(e1 + k*N, exp, texp, N, mask);

         /* ... if not, no quotient term to be written */
         if (!d1)
            k--;
         else
         {
            /* if accumulated coeff in three words */
            if (small)
            {
               ulong d[3];

               /* compute d = abs(c) */
               if (0 > (slong) c[2])
                  mpn_neg(d, c, 3);
               else
                  flint_mpn_copyi(d, c, 3);

               /* check quotient of accumulated coeff by lc(poly3) is small */
               if (d[2] != 0 || ub <= d[1] ||
                  (ub == 0 && 0 > (slong) d[0])) /* quotient not a small */
               {
                  /* upgrade to multiprecision accumulated coeffs */
                  fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

                  small = 0;
               } else /* quotient fits a small */
               {
                  ulong q, r1;

                  /* write out quotient coefficient and discard remainder */
                  sdiv_qrnnd(q, r1, c[1], c[0], *mb);

                  fmpz_set_si(p1 + k, q);

                  if (COEFF_IS_MPZ(FLINT_ABS(q))) /* quotient too large */
                  {
                     /* upgrade to multiprecision accumulated coeffs */
                     fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

                     small = 0;
                  }
               }
            } 

            /* if accumulated coeff is multiprecision, incl. from upgrade */
            if (!small)
            {
               /* write out quotient coefficient */
               fmpz_fdiv_q(p1 + k, qc, mb);
            }

            /* if quotient coefficient nonzero */
            if (!fmpz_is_zero(p1 + k))
            {
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

                  /* insert (i, k, exp3[i] + e1[k]) */
                  mpoly_monomial_add(exp_list[exp_next], exp3 + i*N, e1 + k*N, N);

                  /* but only if exponent is big enough */
                  if (!mpoly_monomial_lt(texp, exp_list[exp_next], N,
                                                               maskhi, masklo))
                  {
                     if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x2,
                                                 &heap_len, N, maskhi, masklo))
                        exp_next--;
                  }
               }
               s = 1;
            } else /* quotient coeff was zero, nothing to write out */
               k--;
         }
      } 
      
      /* zero temporary accumulator */
      fmpz_zero(qc);  
   }

   k++;

cleanup2:

   fmpz_clear(mb);
   fmpz_clear(qc);

   (*polyq) = p1;
   (*expq) = e1;
   
   TMP_END;

   /* return quotient poly length */
   return k;
}

void fmpz_mpoly_div_monagan_pearce(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, 
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
   slong exp_bits, N, lenq = 0;
   ulong * exp2, * exp3;
   ulong maskhi, masklo;
   int free2 = 0, free3 = 0;
   fmpz_mpoly_t temp1;
   fmpz_mpoly_struct * tq;

   /* check divisor is nonzero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_div_monagan_pearce");

   /* dividend zero, write out quotient */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(q, ctx);
 
      return;
   }

   /* maximum bits in quotient exps and inputs is max for poly2 and poly3 */
   exp_bits = FLINT_MAX(poly2->bits, poly3->bits);

   masks_from_bits_ord(maskhi, masklo, exp_bits, ctx->ord);
   /* number of words required for exponent vectors */
   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   /* ensure input exponents packed to same size as output exponents */
   exp2 = mpoly_unpack_monomials(exp_bits, poly2->exps, 
                                           poly2->length, ctx->n, poly2->bits);

   free2 = exp2 != poly2->exps;

   exp3 = mpoly_unpack_monomials(exp_bits, poly3->exps, 
                                           poly3->length, ctx->n, poly3->bits);
   
   free3 = exp3 != poly3->exps;

   /* check divisor leading monomial is at most that of the dividend */
   if (mpoly_monomial_lt(exp3, exp2, N, maskhi, masklo))
   {
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

   /* do division with remainder */
   while ((lenq = _fmpz_mpoly_div_monagan_pearce(&tq->coeffs, &tq->exps,
                         &tq->alloc, poly2->coeffs, exp2, poly2->length, 
                                     poly3->coeffs, exp3, poly3->length,
                                      exp_bits, N, maskhi, masklo)) == -WORD(1)
            && exp_bits < FLINT_BITS)
   {
      ulong * old_exp2 = exp2, * old_exp3 = exp3;

      exp2 = mpoly_unpack_monomials(2*exp_bits, old_exp2, 
                                              poly2->length, ctx->n, exp_bits);

      exp3 = mpoly_unpack_monomials(2*exp_bits, old_exp3, 
                                              poly3->length, ctx->n, exp_bits);
   
      exp_bits *= 2;

      if (free2)
         flint_free(old_exp2);

      if (free3)
         flint_free(old_exp3);

      free2 = free3 = 1; 

      masks_from_bits_ord(maskhi, masklo, exp_bits, ctx->ord);
      N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

      fmpz_mpoly_fit_bits(tq, exp_bits, ctx);
      tq->bits = exp_bits;
   }

   if (lenq == -WORD(1))
      flint_throw(FLINT_EXPOF,
                      "Exponent overflow in fmpz_mpoly_div_monagan_pearce");

   /* take care of aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_swap(temp1, q, ctx);

      fmpz_mpoly_clear(temp1, ctx);
   } 

   _fmpz_mpoly_set_length(q, lenq, ctx);

cleanup3:

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);
}
