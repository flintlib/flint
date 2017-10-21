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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "assert.h"


void fmpz_mpoly_univar_init(fmpz_mpoly_univar_t poly, const fmpz_mpoly_ctx_t ctx)
{
    poly->coeffs = NULL;
    poly->exps = NULL;
    poly->alloc = 0;
    poly->length = 0;
    poly->var = -WORD(1);
}


void fmpz_mpoly_univar_clear(fmpz_mpoly_univar_t poly, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->alloc; i++)
        fmpz_mpoly_clear(poly->coeffs + i, ctx);
    flint_free(poly->coeffs);
    flint_free(poly->exps);
}


void fmpz_mpoly_univar_swap(fmpz_mpoly_univar_t poly1, fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpoly_univar_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}



void fmpz_mpoly_univar_fit_length(fmpz_mpoly_univar_t poly, slong length, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = poly->alloc;
    slong new_alloc = FLINT_MAX(length, 2*poly->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            poly->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            poly->coeffs = (fmpz_mpoly_struct *) flint_malloc(new_alloc*sizeof(fmpz_mpoly_struct));
        } else
        {
            poly->exps = (ulong *) flint_realloc(poly->exps, new_alloc*sizeof(ulong));
            poly->coeffs = (fmpz_mpoly_struct *) flint_realloc(poly->coeffs, new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpoly_init(poly->coeffs + i, ctx);
        }
        poly->alloc = new_alloc;
    }

}

/* if the coefficient doesn't exist, a new one is created */
fmpz_mpoly_struct * fmpz_mpoly_univar_get_coeff(fmpz_mpoly_univar_t poly, ulong pow, slong bits, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_mpoly_struct * xk;

    for (i = 0; i < poly->length && poly->exps[i] >= pow; i++)
    {
        if (poly->exps[i] == pow) 
        {
            return poly->coeffs + i;
        }
    }

    fmpz_mpoly_univar_fit_length(poly, poly->length + 1, ctx);

    for (j = poly->length; j > i; j--)
    {
        poly->exps[j] = poly->exps[j - 1];
        fmpz_mpoly_swap(poly->coeffs + j, poly->coeffs + j - 1, ctx);
    }
    
    poly->length++;
    poly->exps[i] = pow;
    xk = poly->coeffs + i;
    xk->length = 0;
    fmpz_mpoly_fit_bits(xk, bits, ctx);
    xk->bits = bits;
    return xk;
}



void mpoly_off_shift(slong * off, slong * shift, slong idx, int deg, int rev, slong fpw, slong nfields, slong bits)
{
    if (rev)
    {
        *off   = (nfields - 1 - idx)/fpw;
        *shift = (nfields - 1 - idx)%fpw;
    } else
    {
        *off   = (deg + idx)/fpw;
        *shift = (deg + idx)%fpw;
    }
    *shift = (fpw - 1 - *shift) * bits;
}

void mpoly_univar_exp(ulong * one, slong var, int deg, slong N, slong off, slong shift, slong fpw, slong bits)
{
    slong i;
    for (i = 0; i < N; i++)
        one[i] = 0;
    one[off] = WORD(1) << shift;
    if (deg)
        one[0] |= WORD(1) << ((fpw - 1)*bits);
}


void fmpz_mpoly_univar_print(fmpz_mpoly_univar_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf("+");
        flint_printf("(");
        fmpz_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        if (x==NULL)
            flint_printf(")*x%wd^%wd", poly->var+1,poly->exps[i]);
        else
            flint_printf(")*%s^%wd", x[poly->var],poly->exps[i]);
    }
}


void mpoly_monomial_madd(ulong * exp1, const ulong * exp2, ulong scalar,
                                                   const ulong * exp3, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp1[i] = exp2[i] + scalar*exp3[i];
}

void mpoly_monomial_msub(ulong * exp1, const ulong * exp2, ulong scalar,
                                                   const ulong * exp3, slong N)
{
   slong i;
   for (i = 0; i < N; i++)
      exp1[i] = exp2[i] - scalar*exp3[i];
}


void fmpz_mpoly_to_univar(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, shift, off, bits, fpw, N;
    ulong k;
    int deg, rev;
    ulong mask;
    slong poly1_old_length = poly1->length;
    slong len = poly2->length;
    fmpz * coeff = poly2->coeffs;
    ulong * exp = poly2->exps;
    ulong * one;
    fmpz_mpoly_struct * xk;
    slong xk_len;
    TMP_INIT;

    TMP_START;
    
    bits = poly2->bits;
    fpw = FLINT_BITS/bits;
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    N = words_per_exp(ctx->n, poly2->bits);
    degrev_from_ord(deg, rev, ctx->ord);
    mpoly_off_shift(&off, &shift, var, deg, rev, fpw, ctx->n, bits);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_univar_exp(one, var, deg, N, off, shift, fpw, bits);

    poly1->length = 0;
    poly1->var = var;

    for (i = 0; i < len; i++)
    {
        k = (exp[N*i + off] >> shift) & mask;
        xk = fmpz_mpoly_univar_get_coeff(poly1, k, bits, ctx);
        xk_len = xk->length;
        fmpz_mpoly_fit_length(xk, xk_len + 1, ctx);
        fmpz_set(xk->coeffs + xk_len, coeff + i);
        mpoly_monomial_msub(xk->exps + N*xk_len, exp + N*i, k, one, N);
        xk->length = xk_len + 1;
    }

    /* demote remaining coefficients of coefficients */
    for (i = 0; i < poly1->length; i++)
    {
        xk = poly1->coeffs + i;
        for (j = xk->length; j < xk->alloc; j++)
        {
            _fmpz_demote(xk->coeffs + j);
        }
    }

    /* demote remaining coefficients */
    for (i = poly1->length; i < poly1_old_length; i++)
    {
        fmpz_mpoly_clear(poly1->coeffs + i, ctx);
        fmpz_mpoly_init(poly1->coeffs + i, ctx);
    }    

    TMP_END;
}


/*
    currently this function doesn't work if the cofficients depend on the main variable
    the assertion x->next == NULL would need to be removed and a loop put in place
    other asserts would need to be removed as well
*/
void fmpz_mpoly_from_univar(fmpz_mpoly_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong i, shift, off, bits, fpw, N;
    ulong k;
    slong next_loc, heap_len = 1;
    ulong maskhi, masklo;
    int deg, rev;
    slong total_len, p_len;
    fmpz * p_coeff;
    ulong * p_exp;
    slong p_alloc;
    slong var = poly2->var;
    mpoly_heap_s * heap;
    ulong ** poly2_exps;
    ulong * exp;
    ulong * one;
    mpoly_heap_t * chain, * x;
    TMP_INIT;

    TMP_START;

    if (poly2->length == 0)
    {
        fmpz_mpoly_zero(poly1, ctx);
        return;        
    }

    bits = 0;
    for (i = 0; i < poly2->length; i++)
        bits = FLINT_MAX(bits, (poly2->coeffs + i)->bits);

    fpw = FLINT_BITS/bits;
    N = words_per_exp(ctx->n, bits);
    masks_from_bits_ord(maskhi, masklo, bits, ctx->ord);
    degrev_from_ord(deg, rev, ctx->ord);
    mpoly_off_shift(&off, &shift, var, deg, rev, fpw, ctx->n, bits);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_univar_exp(one, var, deg, N, off, shift, fpw, bits);

    poly2_exps = (ulong **) TMP_ALLOC(poly2->length*sizeof(ulong));
    total_len = 0;
    for (i = 0; i < poly2->length; i++)
    {
        total_len += (poly2->coeffs + i)->length;
        poly2_exps[i] = (poly2->coeffs + i)->exps;
        if (bits != (poly2->coeffs + i)->bits)
        {
            poly2_exps[i] = (ulong *) flint_malloc(N*(poly2->coeffs + i)->length*sizeof(ulong));
            mpoly_unpack_monomials(poly2_exps[i], bits, (poly2->coeffs + i)->exps, (poly2->coeffs + i)->bits, (poly2->coeffs + i)->length, ctx->n);
        }
    }

    fmpz_mpoly_fit_length(poly1, total_len, ctx);
    fmpz_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    p_alloc = poly1->alloc;

    next_loc = poly2->length + 2;
    heap = (mpoly_heap_s *) TMP_ALLOC((poly2->length + 1)*sizeof(mpoly_heap_s));
    exp = (ulong *) TMP_ALLOC(poly2->length*N*sizeof(ulong));
    chain = (mpoly_heap_t *) TMP_ALLOC(poly2->length*sizeof(mpoly_heap_t));

    for (i = 0; i < poly2->length; i++)
    {
        k = poly2->exps[i];
        x = chain + i;
        x->i = i;
        x->j = 0;
        x->next = NULL;
        mpoly_monomial_madd(exp + N*i, poly2_exps[x->i] + N*x->j, k, one, N);
        _mpoly_heap_insert(heap, exp + N*i, x, &next_loc, &heap_len, N, maskhi, masklo);
    }

    p_len = 0;
    while (heap_len > 1)
    {
        _fmpz_mpoly_fit_length(&p_coeff, &p_exp, &p_alloc, p_len + 1, N);
        mpoly_monomial_set(p_exp + N*p_len, heap[1].exp, N);
        x = _mpoly_heap_pop(heap, &heap_len, N, maskhi, masklo);
        fmpz_set(p_coeff + p_len, (poly2->coeffs + x->i)->coeffs + x->j);
        p_len++;

        assert(x->next == NULL);
        
        if (x->j + 1 < (poly2->coeffs + x->i)->length) {
            k = poly2->exps[x->i];
            x->j = x->j + 1;
            x->next = NULL;
            mpoly_monomial_madd(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);
            _mpoly_heap_insert(heap, exp + N*x->i, x, &next_loc, &heap_len, N, maskhi, masklo);
        }
    }

    assert(total_len == p_len);
    poly1->coeffs = p_coeff;
    poly1->exps = p_exp;
    poly1->alloc = p_alloc;
    _fmpz_mpoly_set_length(poly1, p_len, ctx);

    for (i = 0; i < poly2->length; i++)
    {
        if (poly2_exps[i] != (poly2->coeffs + i)->exps)
            flint_free(poly2_exps[i]);
    }

    TMP_END;
}

void fmpz_mpoly_univar_test(fmpz_mpoly_univar_t poly, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (!mpoly_monomials_test(poly->exps, poly->length, WORD(1), WORD(0), WORD(0)))
        flint_throw(FLINT_ERROR, "Univariate polynomial exponents invalid");
    for (i = 0; i < poly->length; i++)
        fmpz_mpoly_test(poly->coeffs + i, ctx);
}

int fmpz_mpoly_univar_equal(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (poly1->length != poly2->length)
        return 0;
    for (i = 0; i < poly1->length; i++)
    {
        if (poly1->exps[i] != poly2->exps[i])
            return 0;
    }
    for (i = 0; i < poly1->length; i++)
    {
        if (!fmpz_mpoly_equal(poly1->coeffs + i, poly2->coeffs + i, ctx))
            return 0;
    }
    return 1;
}


void fmpz_mpoly_univar_add(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_univar_t poly3, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_mpoly_struct * coeff1, * coeff2, * coeff3;
    ulong * exp1, * exp2, * exp3;
    slong len1, len2, len3;

    assert(poly2->var == poly3->var);

    if (poly1 == poly2 || poly1 == poly3)
    {
        fmpz_mpoly_univar_t temp;
        fmpz_mpoly_univar_init(temp, ctx);
        fmpz_mpoly_univar_add(temp, poly2, poly3, ctx);
        fmpz_mpoly_univar_swap(poly1, temp, ctx);
        fmpz_mpoly_univar_clear(temp, ctx);
        return;
    }

    poly1->var = poly2->var;

    len1 = 0;
    len2 = poly2->length;
    len3 = poly3->length;
    fmpz_mpoly_univar_fit_length(poly1, len2 + len3, ctx);
    coeff1 = poly1->coeffs;
    coeff2 = poly2->coeffs;
    coeff3 = poly3->coeffs;
    exp1 = poly1->exps;
    exp2 = poly2->exps;
    exp3 = poly3->exps;

    i = j = 0;
    while (i < len2 && j < len3)
    {
        if (exp2[i] > exp3[j])
        {
            fmpz_mpoly_set(coeff1 + len1, coeff2 + i, ctx);
            exp1[len1++] = exp2[i++];
        } else if (exp2[i] == exp3[j])
        {
            fmpz_mpoly_add(coeff1 + len1, coeff2 + i, coeff3 + j, ctx);
            exp1[len1] = exp2[i];
            len1 += !fmpz_mpoly_is_zero(coeff1 + len1, ctx);
            i++;
            j++;
        } else
        {
            fmpz_mpoly_set(coeff1 + len1, coeff3 + j, ctx);
            exp1[len1++] = exp3[j++];
        }
    }

    while (i < len2)
    {
        fmpz_mpoly_set(coeff1 + len1, coeff2 + i, ctx);
        exp1[len1++] = exp2[i++];
    }

    while (j < len3)
    {
        fmpz_mpoly_set(coeff1 + len1, coeff3 + j, ctx);
        exp1[len1++] = exp3[j++];
    }

    /* demote remaining coefficients */
    for (i = len1; i < poly1->length; i++)
    {
        fmpz_mpoly_clear(poly1->coeffs + i, ctx);
        fmpz_mpoly_init(poly1->coeffs + i, ctx);
    }
    poly1->length = len1;
}


void fmpz_mpoly_univar_mul(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_univar_t poly3, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong next_loc, heap_len;
    fmpz_mpoly_struct * coeff1, * coeff2, * coeff3;
    ulong exp, * exp1, * exp2, * exp3;
    slong len1, len2, len3;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain, * x;
    slong * store, * store_base;
    int first;
    fmpz_mpoly_t tcoeff;
    TMP_INIT;

    assert(poly2->var == poly3->var);

    len1 = 0;
    if (poly2->length == 0 || poly3->length == 0)
        goto done;
    if (poly1 == poly2 || poly1 == poly3)
    {
        fmpz_mpoly_univar_t temp;
        fmpz_mpoly_univar_init(temp, ctx);
        fmpz_mpoly_univar_mul(temp, poly2, poly3, ctx);
        fmpz_mpoly_univar_swap(poly1, temp, ctx);
        fmpz_mpoly_univar_clear(temp, ctx);
        return;
    }

    poly1->var = poly2->var;

    len2 = poly2->length;
    len3 = poly3->length;
    fmpz_mpoly_univar_fit_length(poly1, len2*len3, ctx);
    coeff1 = poly1->coeffs;
    coeff2 = poly2->coeffs;
    coeff3 = poly3->coeffs;
    exp1 = poly1->exps;
    exp2 = poly2->exps;
    exp3 = poly3->exps;

    TMP_START;

    next_loc = len2 + 4;
    heap_len = 1;
    heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len2*sizeof(slong));

    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x, &next_loc, &heap_len, WORD(0));

    fmpz_mpoly_init(tcoeff, ctx);

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;
        exp1[len1] = exp;
        first = 1;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, WORD(0));
            *store++ = x->i;
            *store++ = x->j;

            if (first)
            {
                first = 0; 
                fmpz_mpoly_mul_johnson(coeff1 + len1, coeff2 + x->i, coeff3 + x->j, ctx);
            } else
            {
                fmpz_mpoly_mul_johnson(tcoeff, coeff2 + x->i, coeff3 + x->j, ctx);
                fmpz_mpoly_add(coeff1 + len1, coeff1 + len1, tcoeff, ctx);
            }
            while ((x = x->next) != NULL)
            {
                *store++ = x->i;
                *store++ = x->j;
                fmpz_mpoly_mul_johnson(tcoeff, coeff2 + x->i, coeff3 + x->j, ctx);
                fmpz_mpoly_add(coeff1 + len1, coeff1 + len1, tcoeff, ctx);
            }
        } while (heap_len > 1 && heap[1].exp == exp);

        len1 += !fmpz_mpoly_is_zero(coeff1 + len1, ctx);

        while (store > store_base)
        {
            j = *--store;
            i = *--store;
            if ((j == 0) && (i + 1 < len2))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x, &next_loc, &heap_len, WORD(0));
            }
            if (j + 1 < len3)
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x, &next_loc, &heap_len, WORD(0));
            }
        }
    }

    fmpz_mpoly_clear(tcoeff, ctx);

    TMP_END;

done:

    /* demote remaining coefficients */
    for (i = len1; i < poly1->length; i++)
    {
        fmpz_mpoly_clear(poly1->coeffs + i, ctx);
        fmpz_mpoly_init(poly1->coeffs + i, ctx);
    }
    poly1->length = len1;

}

