/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "long_extras.h"
#include "fmpz_mod_mpoly.h"

/*
f = f0 + f1*x + .... + fd*x^d
g = f^k
  = g0 + g1*x +      + g_kd*x^d

i*g_i*f_0 = sum_{0 <= j <= min(i,d)} ((k+1)*j-i)*f_j*g_{i-j}


g_kd = f_d^k
for i from k*d-1 to 0
    c = 0
    for 1 <= j <= min(d,k*d-i)
        c += (i+j-(d-j)*k)*g_{i+j}*f_{d-j} [*x^(i+d)]
    end
    g_i = c/((k*d-i)*f_d [*x^d])
*/


/*
output g = (f1 + f2 + ... + ft)^k:

g = g1^k
H = {2,1,f2*g1}
while #H > 0 && H1 >= f1
    M = H1
    C = 0
    Q = {};
    while #H > 0 && H1 == M
        (i,j) = popmax(H)
        C += (exp(gj)-k*exp(fi))*coeff(fi)*coeff(gj)
        Q += (i,j)
    end
    for (i,j) in Q
        ...
    end
    if C != 0
        g += C/((exp(g1)-M)*coeff(f1))*x^(M-exp(f1))
        ...
    end
end
*/




/*
g_{(k-1)d} = f_d^(k-1)
for i from k*d-1 to 0
    c = 0
    for 1 <= j <= min(d,k*d-i)
        c += (i+j-(d-j)*k)*g_{i+j}*f_{d-j} [*x^(i+d)]
    end
    g_i = c/((k*d-i)*f_d [*x^d])
*/
/*  
output h = (f1 + f2 + ... + ft)^k:

h = 0;
g = g1^(k-1)
H = {2,1,f2*g1}
while #H > 0
    M = H1
    C = 0
    S = 0
    Q = {};
    while #H > 0 && H1 == M
        (i,j) = popmax(H)
        S += coeff(fi)*coeff(gj)
        if M >= exp(f1)
            C += (exp(gj)-k*exp(fi))*coeff(fi)*coeff(gj)
        end
        Q += (i,j)
    end
    for (i,j) in Q
        ...
    end
    if C != 0
        C /= (exp(g1)-M)*coeff(f1)
        S += C*coeff(f1)
        g += C*x^(M-exp(f1))
        ...
    end
    h += S*x^M
end

*/




slong _fmpz_mpoly_pow_fps1(
    fmpz ** poly1, ulong ** exp1, slong * Aalloc,
    const fmpz * fcoeffs, const ulong * fexps, slong flen,
    ulong k,
    flint_bitcnt_t bits,
    ulong cmpmask)
{
    slong i, j, Alen, galloc, glen;
    slong next_loc;
    slong Qlen = 0, heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    fmpz * Acoeffs = *poly1, * gc = NULL;
    ulong * Aexps = *exp1, * ge;
    ulong exp;
    slong * hind;
    fmpz_t t1, temp1;
    fmpz * S, * C;
    slong gdemote = 0;
    ulong ofmask = mpoly_overflow_mask_sp(bits);
    TMP_INIT;

    TMP_START;

    next_loc = flen + 4;   /* something bigger than heap can ever be */
    heap = TMP_ARRAY_ALLOC(flen + 1, mpoly_heap1_s);
    chain = TMP_ARRAY_ALLOC(flen, mpoly_heap_t);
    Q = TMP_ARRAY_ALLOC(3*flen, slong);
    hind = Q + 2*flen; /* space for heap indices */    
    for (i = 0; i < flen; i++)
        hind[i] = 1;

    fmpz_init(t1);
    fmpz_init(temp1);

    galloc = k*(flen - 1) + 1;
    ge = FLINT_ARRAY_ALLOC(galloc, ulong);
    glen = 1;
    Alen = 1;

    ge[0] = fexps[0]*(k - 1);

    Aexps[0] = fexps[0]*k;

    gc = (fmpz *) flint_calloc(galloc, sizeof(fmpz));
    fmpz_pow_ui(gc + 0, fcoeffs + 0, k - 1);
    fmpz_mul(Acoeffs + 0, gc + 0, fcoeffs + 0);

    x = chain + 1;
    x->i = 1;
    x->j = 0;
    x->next = NULL;

    HEAP_ASSIGN(heap[1], fexps[1] + ge[0], x);

    hind[1] = 2*1 + 0;

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (Alen >= *Aalloc)
        {
            Aexps = (ulong *) flint_realloc(Aexps, 2*sizeof(ulong)*(*Aalloc));
            Acoeffs = (fmpz *) flint_realloc(Acoeffs, 2*sizeof(fmpz)*(*Aalloc));
            flint_mpn_zero(Acoeffs + *Aalloc, *Aalloc);
            (*Aalloc) *= 2;
        }

        if (glen >= galloc)
        {
            ge = FLINT_ARRAY_REALLOC(ge, 2*galloc, ulong);
            gc = FLINT_ARRAY_REALLOC(gc, 2*galloc, fmpz);
            flint_mpn_zero(gc + galloc, galloc);
            galloc *= 2;
        }

        Aexps[Alen] = exp;
        S = Acoeffs + Alen;
        C = gc + glen;

        fmpz_zero(C);
        fmpz_zero(S);

        Qlen = 0;

        ge[glen] = exp - fexps[0];

        if ((ge[glen] & ofmask) == 0)
        {
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do {
                    fmpz fi;
                    slong dd, ddd;
                    hind[x->i] |= 1;
                    Q[Qlen++] = i = x->i;
                    Q[Qlen++] = j = x->j;
                    fi = fcoeffs[i];

                    dd = (k - 1)*fexps[i] - ge[j];

                    if (!COEFF_IS_MPZ(fi) && !z_mul_checked(&ddd, dd, fi))
                    {
                        fmpz_addmul_si(S, gc + j, fi);
                        fmpz_addmul_si(C, gc + j, ddd);
                    }
                    else
                    {
                        fmpz_mul(t1, fcoeffs + i, gc + j);
                        fmpz_add(S, S, t1);
                        fmpz_addmul_si(C, t1, dd);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }
        else
        {
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do {
                    hind[x->i] |= 1;
                    Q[Qlen++] = x->i;
                    Q[Qlen++] = x->j;
                    /*FLINT_ASSERT(x->j >= gdemote);*/
                    fmpz_addmul(S, fcoeffs + x->i, gc + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        FLINT_ASSERT(Qlen <= 2*flen);
      
        while (Qlen > 0)
        {
            /* take node from store */
            j = Q[--Qlen];
            i = Q[--Qlen];

            /* should we go right? */
            if (i + 1 < flen && hind[i + 1] == 2*j + 1)
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*j + 2;
                _mpoly_heap_insert1(heap, fexps[i + 1] + ge[j], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            /* should we go up */
            if (j + 1 < glen && hind[i] < 2*j + 4)
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*j + 4;
                _mpoly_heap_insert1(heap, fexps[i] + ge[j + 1], x,
                                                &next_loc, &heap_len, cmpmask);
            }
        }

        j = hind[flen - 1]/2 - 1;
        for ( ; gdemote < j; gdemote++)
            fmpz_clear(gc + gdemote);

        if (!fmpz_is_zero(C))
        {
            if (fmpz_is_one(fcoeffs + 0))
            {
                fmpz_divexact_si(C, C, exp - k*fexps[0]);
                fmpz_add(S, S, C);
            }
            else
            {
                fmpz_divexact_si(temp1, C, exp - k*fexps[0]);
                fmpz_add(S, S, temp1);
                fmpz_divexact(C, temp1, fcoeffs + 0);
            }

            if ((hind[1] & 1) != 0)
            {
                x = chain + 1;

                x->i = 1;
                x->j = glen;
                x->next = NULL;
                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, fexps[1] + ge[glen], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            glen += 1;
        }

        Alen += !fmpz_is_zero(S);
    }

    (*poly1) = Acoeffs;
    (*exp1) = Aexps;

    fmpz_clear(t1);
    fmpz_clear(temp1);

    flint_free(ge);
    for ( ; gdemote < galloc; gdemote++)
        fmpz_clear(gc + gdemote);
    flint_free(gc);

    TMP_END;

    return Alen;
}


slong _fmpz_mpoly_pow_fps(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2, ulong k,
                              flint_bitcnt_t bits, slong N, const ulong * cmpmask)
{
   const slong topbit = (WORD(1) << (FLINT_BITS - 1));
   const slong mask = ~topbit;
   slong i, rnext, galloc, gnext, exp_next;
   slong next_loc;
   slong next_free, Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x;
   fmpz * Acoeffs = *poly1, * gc = NULL;
   ulong * Aexps = *exp1, * ge, * fik, * exp, * exps, * exp_copy;
   ulong ** exp_list;
   ulong * finalexp, * temp2;
   slong * largest;
   fmpz_t t1, t2, C, S, temp1;
   int first;
   TMP_INIT;

   if (N == 1)
      return _fmpz_mpoly_pow_fps1(poly1, exp1, alloc, poly2, exp2, len2, k,
                                                             bits, cmpmask[0]);

   TMP_START;

   next_loc = len2 + 4;   /* something bigger than heap can ever be */
   heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
   /* 2x as we pull from heap and insert more before processing pulled ones */
   chain = (mpoly_heap_t *) TMP_ALLOC(2*len2*sizeof(mpoly_heap_t));
   reuse = (mpoly_heap_t **) TMP_ALLOC(2*len2*sizeof(mpoly_heap_t *));
   Q = (mpoly_heap_t **) TMP_ALLOC(len2*sizeof(mpoly_heap_t *));
   /* we add 1 to all entries of largest to free up the value 0 */
   largest = (slong *) TMP_ALLOC(len2*sizeof(slong));
   exps = (ulong *) TMP_ALLOC((len2 + 1)*N*sizeof(ulong));
   exp_list = (ulong **) TMP_ALLOC((len2 + 1)*sizeof(ulong *));
   finalexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
   temp2 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
   exp_copy = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   fmpz_init(t1);
   fmpz_init(t2);
   fmpz_init(C);
   fmpz_init(S);
   fmpz_init(temp1);

   for (i = 0; i < len2 + 1; i++)
      exp_list[i] = exps + i*N;

   exp_next = 0;

   for (i = 0; i < 2*len2; i++)
      reuse[i] = chain + i;

   galloc = k*(len2 - 1) + 1;
   ge = (ulong *) flint_malloc(galloc*sizeof(ulong)*N);
   gnext = 0;
   rnext = 0;


    mpoly_monomial_mul_ui_mp(ge + 0, exp2 + 0, N, k - 1);
    mpoly_monomial_mul_ui_mp(Aexps + 0, exp2 + 0, N, k);

   gc = (fmpz *) flint_calloc(galloc, sizeof(fmpz));
   fmpz_pow_ui(gc + 0, poly2 + 0, k - 1);
   fmpz_mul(Acoeffs + 0, gc + 0, poly2 + 0);

   next_free = 0;

   x = reuse[next_free++];
   x->i = 1;
   x->j = 0;
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];

    mpoly_monomial_add_mp(heap[1].exp, exp2 + N, ge + 0, N);

    for (i = 0; i < len2; i++)
        largest[i] = topbit;
    largest[1] = 1;

    fik = (ulong *) TMP_ALLOC(N*len2*sizeof(ulong));

    for (i = 0; i < len2; i++)
        mpoly_monomial_mul_ui_mp(fik + i*N, exp2 + i*N, N, k - 1);

    mpoly_monomial_set(finalexp, exp2, N);

   while (heap_len > 1)
   {
      exp = heap[1].exp;
      mpoly_monomial_set(exp_copy, exp, N);

      rnext++;
      gnext++;

      if (rnext >= *alloc)
      {
         Acoeffs = (fmpz *) flint_realloc(Acoeffs, 2*sizeof(fmpz)*(*alloc));
         Aexps = (ulong *) flint_realloc(Aexps, 2*N*sizeof(ulong)*(*alloc));
         flint_mpn_zero(Acoeffs + *alloc, *alloc);
         (*alloc) *= 2;
      }

      if (gnext >= galloc)
      {
         ge = (ulong *) flint_realloc(ge, 2*N*sizeof(ulong)*galloc);
         gc = (fmpz *) flint_realloc(gc, 2*sizeof(fmpz)*galloc);
         flint_mpn_zero(gc + galloc, galloc);
         galloc *= 2;
      }

      first = 1;

      fmpz_zero(C);
      fmpz_zero(S);

      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         exp_list[--exp_next] = heap[1].exp;
         x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

         largest[x->i] |= topbit;

         fmpz_mul(t1, poly2 + x->i, gc + x->j);
         fmpz_add(S, S, t1);

         if (!mpoly_monomial_gt(finalexp, exp, N, cmpmask))
         {
            mpn_sub_n(temp2, fik + x->i*N, ge + x->j*N, N);
            fmpz_set_signed_ui_array(t2, temp2, N);
            fmpz_addmul(C, t1, t2);
         }

         if (first)
         {
            mpoly_monomial_sub_mp(ge + gnext*N, exp, exp2 + 0, N);
            first = 0; 
         }
      
         Q[Q_len++] = x;

         while ((x = x->next) != NULL)
         {
            largest[x->i] |= topbit;

            fmpz_mul(t1, poly2 + x->i, gc + x->j);
            fmpz_add(S, S, t1);

            if (!mpoly_monomial_gt(finalexp, exp, N, cmpmask))
            {
               mpn_sub_n(temp2, fik + x->i*N, ge + x->j*N, N);
               fmpz_set_signed_ui_array(t2, temp2, N);
               fmpz_addmul(C, t1, t2);
            }

            Q[Q_len++] = x;
         }
      }
      
      while (Q_len > 0)
      {
         slong i, j;

         x = Q[--Q_len];
         i = x->i;
         j = x->j;

         if (i < len2 - 1 && largest[i + 1] == (j | topbit))
         {
            x->i++;
            x->next = NULL;

            mpoly_monomial_add_mp(exp_list[exp_next], exp2 + (i + 1)*N, ge + j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            largest[i + 1] = j + 1;
         }
         else
         {
            reuse[--next_free] = x;
         }

         if (j < gnext - 1 && (largest[i] & mask) < j + 2)
         {
            x = reuse[next_free++];

            x->i = i;
            x->j = j + 1;
            x->next = NULL;

            mpoly_monomial_add_mp(exp_list[exp_next], exp2 + i*N, ge + (j + 1)*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            largest[i] = j + 2;
         }
      }

      if (!fmpz_is_zero(C))
      {
         mpoly_monomial_mul_ui_mp(temp2, exp2 + 0, N, k);

         mpn_sub_n(temp2, exp_copy, temp2, N);
         fmpz_set_signed_ui_array(t2, temp2, N);
         fmpz_divexact(temp1, C, t2);
         fmpz_add(S, S, temp1);
         fmpz_divexact(gc + gnext, temp1, poly2 + 0);

         if ((largest[1] & topbit) != 0)
         {
            x = reuse[next_free++];

            x->i = 1;
            x->j = gnext;
            x->next = NULL;

            mpoly_monomial_add_mp(exp_list[exp_next], exp2 + N, ge + gnext*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
 
            largest[1] = gnext + 1;
         }
      }

      if (!fmpz_is_zero(S))
      {
         fmpz_set(Acoeffs + rnext, S);
         mpoly_monomial_add_mp(Aexps + rnext*N, ge + gnext*N, exp2 + 0, N);

      } else
         rnext--;

      if (fmpz_is_zero(C))
         gnext--;
   }

   rnext++;

   (*poly1) = Acoeffs;
   (*exp1) = Aexps;
   
   fmpz_clear(t1);
   fmpz_clear(t2);
   fmpz_clear(C);
   fmpz_clear(S);
   fmpz_clear(temp1);

   flint_free(ge);
   for (i = 0; i < galloc; i++)
      fmpz_clear(gc + i);
   flint_free(gc);

   TMP_END;

   return rnext;
}

void fmpz_mpoly_pow_fps(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           ulong k, const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, len = 0;
    fmpz * maxBfields;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Bexp = B->exps;
    int freeBexp = 0;
    TMP_INIT;

    FLINT_ASSERT(k >= 2);
    FLINT_ASSERT(B->length > 0);

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    _fmpz_vec_scalar_mul_ui(maxBfields, maxBfields, ctx->minfo->nfields, k);

    exp_bits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    if (exp_bits > B->bits)
    {
       freeBexp = 1;
       Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
       mpoly_repack_monomials(Bexp, exp_bits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    if (B->length == 1)
    {
        fmpz_mpoly_fit_length(A, 1, ctx);
        fmpz_mpoly_fit_bits(A, exp_bits, ctx);
        A->bits = exp_bits;

        fmpz_pow_ui(A->coeffs + 0, B->coeffs + 0, k);

        if (exp_bits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps, Bexp, N, k);
        else
            mpoly_monomial_mul_ui_mp(A->exps, Bexp, N, k);

        len = 1;
        goto cleanup;
    }

    if (A == B)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, k*(B->length - 1) + 1, exp_bits, ctx);
        len = _fmpz_mpoly_pow_fps(&T->coeffs, &T->exps, &T->alloc,
                          B->coeffs, Bexp, B->length, k, exp_bits, N, cmpmask);

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, k*(B->length - 1) + 1, exp_bits, ctx);
        len = _fmpz_mpoly_pow_fps(&A->coeffs, &A->exps, &A->alloc,
                          B->coeffs, Bexp, B->length, k, exp_bits, N, cmpmask);
    }

cleanup:

    if (freeBexp)
        flint_free(Bexp);

    _fmpz_mpoly_set_length(A, len, ctx);

    TMP_END;
}
