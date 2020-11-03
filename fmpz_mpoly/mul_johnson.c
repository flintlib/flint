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

/*
   Set poly1 to poly2*poly3 using Johnson's heap method. The function
   realocates its output and returns the length of the product. This
   version of the function assumes the exponent vectors all fit in a
   single word. Assumes input polys are nonzero.
*/
slong _fmpz_mpoly_mul_johnson1(fmpz ** poly1, ulong ** exp1, slong * alloc,
              const fmpz * poly2, const ulong * exp2, slong len2,
              const fmpz * poly3, const ulong * exp3, slong len3, ulong maskhi)
{
   slong i, j, k;
   slong next_loc;
   slong Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_heap_t * chain;
   slong * Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong * hind;
   ulong exp, cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   int first, small;
   TMP_INIT;

   TMP_START;

   /* whether input coeffs are small, thus output coeffs fit in three words */
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   next_loc = len2 + 4;   /* something bigger than heap can ever be */
   heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
    for (i = 0; i < len2; i++)
        hind[i] = 1;

   /* put (0, 0, exp2[0] + exp3[0]) on heap */
   x = chain + 0;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   HEAP_ASSIGN(heap[1], exp2[0] + exp3[0], x);
   hind[0] = 2*1 + 0;

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get exponent field of heap top */
      exp = heap[1].exp;
      
      /* realloc output poly ready for next product term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && heap[1].exp == exp)
      {
         /* pop chain from heap */
         x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
         
         /* take node out of heap and put into store */
         hind[x->i] |= WORD(1);
         Q[Q_len++] = x->i;
         Q[Q_len++] = x->j;

         /* if output coeffs will fit in three words */
         if (small)
         {
            /* compute product of input poly coeffs */
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               /* set output monomial */
               e1[k] = exp;
               first = 0; 
            } else /* addmul product of input poly coeffs */
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }

            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         } else /* output coeffs require multiprecision */
         {
            if (first) /* compute product of input poly coeffs */
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               e1[k] = exp;
               first = 0; 
            } else
            {  /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
            }

            /* for each node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         }
      }
      
      /* for each node temporarily stored */
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
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
         }

         /* should we go up? */
         if (  (j + 1 < len3)
            && ((hind[i] & 1) == 1)
            && (  (i == 0)
               || (hind[i - 1] >  2*(j + 2) + 1)
               || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */
               )
            )
         {
            x = chain + i;
            x->i = i;
            x->j = j + 1;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
         }
      }

      /* set output poly coeff from temporary accumulation, if not multiprec */
      if (small)
         fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);

      if (fmpz_is_zero(p1 + k))
         k--;
   }

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

/*
   Set poly1 to poly2*poly3 using Johnson's heap method. The function
   realocates its output and returns the length of the product. This
   version of the function assumes the exponent vectors take N words.
*/
slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
                              flint_bitcnt_t bits, slong N, const ulong * cmpmask)
{
   slong i, j, k;
   slong next_loc;
   slong Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   slong * Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   ulong * exp, * exps;
   ulong ** exp_list;
   slong exp_next;
   slong * hind;
   int first, small;
   TMP_INIT;

   /* if exponent vectors fit in single word, call special version */
   if (N == 1)
      return _fmpz_mpoly_mul_johnson1(poly1, exp1, alloc,
                             poly2, exp2, len2, poly3, exp3, len3, cmpmask[0]);

   TMP_START;

   /* whether input coeffs are small, thus output coeffs fit in three words */
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   next_loc = len2 + 4;   /* something bigger than heap can ever be */
   heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
   /* allocate space for exponent vectors of N words */
   exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
   /* list of pointers to allocated exponent vectors */
   exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));
   for (i = 0; i < len2; i++)
      exp_list[i] = exps + i*N;

   /* space for heap indices */
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

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get pointer to exponent field of heap top */
      exp = heap[1].exp;

      /* realloc output poly ready for next product term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         /* pop chain from heap and set exponent field to be reused */
         exp_list[--exp_next] = heap[1].exp;

         x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

         /* take node out of heap and put into store */
         hind[x->i] |= WORD(1);
         Q[Q_len++] = x->i;
         Q[Q_len++] = x->j;

         /* if output coeffs will fit in three words */
         if (small)
         {
            /* compute product of input poly coeffs */
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               /* set output monomial */
               mpoly_monomial_set(e1 + k*N, exp, N);

               first = 0; 
            } else /* addmul product of input poly coeffs */
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }
      
            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         } else /* output coeffs require multiprecision */
         {
            if (first) /* compute product of input poly coeffs */
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               /* set output monomial */
               mpoly_monomial_set(e1 + k*N, exp, N);

               first = 0; 
            } else
            {  /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
            }

            /* for each node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         }
      }

      /* for each node temporarily stored */
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
               || (hind[i - 1] >  2*(j + 2) + 1)
               || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */
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

      /* set output poly coeff from temporary accumulation, if not multiprec */
      if (small)
         fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);

      if (fmpz_is_zero(p1 + k))
         k--;
   }

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

/* maxBfields gets clobbered */
void _fmpz_mpoly_mul_johnson_maxfields(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B, fmpz * maxBfields,
    const fmpz_mpoly_t C, fmpz * maxCfields,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N, Alen;
    flint_bitcnt_t Abits;
    ulong * cmpmask;
    ulong * Bexp, * Cexp;
    int freeBexp, freeCexp;
    TMP_INIT;

    TMP_START;

    _fmpz_vec_add(maxBfields, maxBfields, maxCfields, ctx->minfo->nfields);

    Abits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(MPOLY_MIN_BITS, Abits + 1);
    Abits = FLINT_MAX(Abits, B->bits);
    Abits = FLINT_MAX(Abits, C->bits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
    freeBexp = 0;
    Bexp = B->exps;
    if (Abits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, Abits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    freeCexp = 0;
    Cexp = C->exps;
    if (Abits > C->bits)
    {
        freeCexp = 1;
        Cexp = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexp, Abits, C->exps, C->bits,
                                                        C->length, ctx->minfo);
    }

    /* deal with aliasing and do multiplication */
    if (A == B || A == C)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, B->length + C->length, Abits, ctx);

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
        {
            Alen = _fmpz_mpoly_mul_johnson(&T->coeffs, &T->exps, &T->alloc,
                                                  C->coeffs, Cexp, C->length,
                                                  B->coeffs, Bexp, B->length,
                                                         Abits, N, cmpmask);
        }
        else
        {
            Alen = _fmpz_mpoly_mul_johnson(&T->coeffs, &T->exps, &T->alloc,
                                                  B->coeffs, Bexp, B->length,
                                                  C->coeffs, Cexp, C->length,
                                                         Abits, N, cmpmask);
        }

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, B->length + C->length, Abits, ctx);

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
        {
            Alen = _fmpz_mpoly_mul_johnson(&A->coeffs, &A->exps, &A->alloc,
                                                  C->coeffs, Cexp, C->length,
                                                  B->coeffs, Bexp, B->length,
                                                         Abits, N, cmpmask);
        }
        else
        {
            Alen = _fmpz_mpoly_mul_johnson(&A->coeffs, &A->exps, &A->alloc,
                                                  B->coeffs, Bexp, B->length,
                                                  C->coeffs, Cexp, C->length,
                                                         Abits, N, cmpmask);
        }
    }

    if (freeBexp)
        flint_free(Bexp);

    if (freeCexp)
        flint_free(Cexp);

    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}



void fmpz_mpoly_mul_johnson(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    _fmpz_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
}
