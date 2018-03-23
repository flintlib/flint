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
#include "fmpz_poly.h"
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


void fmpz_mpoly_univar_swap(fmpz_mpoly_univar_t poly1,
                         fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpoly_univar_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}

void fmpz_mpoly_univar_fit_length(fmpz_mpoly_univar_t poly,
                                      slong length, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = poly->alloc;
    slong new_alloc = FLINT_MAX(length, 2*poly->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            poly->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            poly->coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        } else
        {
            poly->exps = (ulong *) flint_realloc(poly->exps,
                                                      new_alloc*sizeof(ulong));
            poly->coeffs = (fmpz_mpoly_struct *) flint_realloc(poly->coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpoly_init(poly->coeffs + i, ctx);
        }
        poly->alloc = new_alloc;
    }

}

void fmpz_mpoly_univar_set(fmpz_mpoly_univar_t poly1,
                   const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_struct * coeff1, * coeff2;
    ulong * exp1, * exp2;
    slong len1, len2;

    poly1->var = poly2->var;

    len1 = 0;
    len2 = poly2->length;
    fmpz_mpoly_univar_fit_length(poly1, len2, ctx);
    coeff1 = poly1->coeffs;
    coeff2 = poly2->coeffs;
    exp1 = poly1->exps;
    exp2 = poly2->exps;

    i = 0;
    while (i < len2)
    {
        fmpz_mpoly_set(coeff1 + len1, coeff2 + i, ctx);
        exp1[len1++] = exp2[i++];
    }
    len1 = len2;

    /* demote remaining coefficients */
    for (i = len1; i < poly1->length; i++)
    {
        fmpz_mpoly_clear(poly1->coeffs + i, ctx);
        fmpz_mpoly_init(poly1->coeffs + i, ctx);
    }
    poly1->length = len1;
}


/* if the coefficient doesn't exist, a new one is created */
fmpz_mpoly_struct * _fmpz_mpoly_univar_get_coeff(fmpz_mpoly_univar_t poly,
                             ulong pow, slong bits, const fmpz_mpoly_ctx_t ctx)
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

void fmpz_mpoly_univar_print_pretty(const fmpz_mpoly_univar_t poly,
                                   const char ** x, const fmpz_mpoly_ctx_t ctx)
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


int fmpz_mpoly_to_univar(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_t poly2,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, shift, off, bits, N;
    ulong k;
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

    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));

    if (bits <= FLINT_BITS)
    {
        mpoly_gen_oneexp_offset_shift(one, &off, &shift, var, N, bits, ctx->minfo);

        poly1->length = 0;
        poly1->var = var;
        for (i = 0; i < len; i++)
        {
            ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
            k = (exp[N*i + off] >> shift) & mask;
            xk = _fmpz_mpoly_univar_get_coeff(poly1, k, bits, ctx);
            xk_len = xk->length;
            fmpz_mpoly_fit_length(xk, xk_len + 1, ctx);
            fmpz_set(xk->coeffs + xk_len, coeff + i);
            xk->length = xk_len + 1;
            mpoly_monomial_msub(xk->exps + N*xk_len, exp + N*i, k, one, N);
        }
    } else
    {
        fmpz_t c;
        fmpz_init(c);

        mpoly_gen_oneexp_offset_mp(one, &off, var, N, bits, ctx->minfo);

        poly1->length = 0;
        poly1->var = var;
        for (i = 0; i < len; i++)
        {
            fmpz_set_ui_array(c, exp + N*i + off, bits/FLINT_BITS);
            if (!fmpz_fits_si(c))
                goto failed;
            k = fmpz_get_si(c);
            xk = _fmpz_mpoly_univar_get_coeff(poly1, k, bits, ctx);
            xk_len = xk->length;
            fmpz_mpoly_fit_length(xk, xk_len + 1, ctx);
            fmpz_set(xk->coeffs + xk_len, coeff + i);
            xk->length = xk_len + 1;
            mpoly_monomial_msub_mp(xk->exps + N*xk_len, exp + N*i, k, one, N);
        }        

        fmpz_clear(c);
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
    return 1;

failed:
    fmpz_mpoly_univar_clear(poly1, ctx);
    fmpz_mpoly_univar_init(poly1, ctx);
    TMP_END;
    return 0;
}


/*
    currently this function doesn't work if the coefficients depend on the main variable
    the assertion x->next == NULL would need to be removed and a loop put in place
    other asserts would need to be removed as well
*/
void fmpz_mpoly_from_univar(fmpz_mpoly_t poly1, const fmpz_mpoly_univar_t poly2,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, bits, N;
    ulong k;
    slong next_loc, heap_len = 1;
    ulong * cmpmask;
    slong total_len, p_len;
    fmpz * p_coeff;
    fmpz * gen_fields, * tmp_fields, * max_fields;
    ulong * p_exp;
    slong p_alloc;
    slong var = poly2->var;
    mpoly_heap_s * heap;
    ulong ** poly2_exps;
    ulong * exp;
    ulong * one;
    mpoly_heap_t * chain, * x;
    TMP_INIT;

    if (poly2->length == 0)
    {
        fmpz_mpoly_zero(poly1, ctx);
        return;        
    }

    TMP_START;

    /* find bits required to represent result */
    gen_fields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    tmp_fields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    max_fields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(gen_fields + i);
        fmpz_init(tmp_fields + i);
        fmpz_init_set_ui(max_fields + i, 0);
    }

    mpoly_gen_fields_fmpz(gen_fields, var, ctx->minfo);

    for (i = 0; i < poly2->length; i++)
    {
        fmpz_mpoly_struct * pi = poly2->coeffs + i;
        mpoly_max_fields_fmpz(tmp_fields, pi->exps, pi->length, pi->bits, ctx->minfo);
        _fmpz_vec_scalar_addmul_si(tmp_fields, gen_fields, ctx->minfo->nfields, poly2->exps[i]);
        _fmpz_vec_max_inplace(max_fields, tmp_fields, ctx->minfo->nfields);
    }
    bits = _fmpz_vec_max_bits(max_fields, ctx->minfo->nfields);
    bits = FLINT_MAX(MPOLY_MIN_BITS, bits + 1);
    bits = mpoly_fix_bits(bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(gen_fields + i);
        fmpz_clear(tmp_fields + i);
        fmpz_clear(max_fields + i);
    }

    /* pack everything into bits */
    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_pack_vec_fmpz(one, gen_fields, bits, ctx->minfo->nfields, 1);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    poly2_exps = (ulong **) TMP_ALLOC(poly2->length*sizeof(ulong));
    total_len = 0;
    for (i = 0; i < poly2->length; i++)
    {
        total_len += (poly2->coeffs + i)->length;
        poly2_exps[i] = (poly2->coeffs + i)->exps;
        if (bits != (poly2->coeffs + i)->bits)
        {
            poly2_exps[i] = (ulong *) flint_malloc(
                                  N*(poly2->coeffs + i)->length*sizeof(ulong));
            mpoly_repack_monomials(poly2_exps[i], bits,
                    (poly2->coeffs + i)->exps, (poly2->coeffs + i)->bits,
                                      (poly2->coeffs + i)->length, ctx->minfo);
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

        if (bits <= FLINT_BITS)
            mpoly_monomial_madd(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);
        else
            mpoly_monomial_madd_mp(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);

        _mpoly_heap_insert(heap, exp + N*i, x, &next_loc, &heap_len, N,
                                                               cmpmask);
    }

    p_len = 0;
    while (heap_len > 1)
    {
        _fmpz_mpoly_fit_length(&p_coeff, &p_exp, &p_alloc, p_len + 1, N);
        mpoly_monomial_set(p_exp + N*p_len, heap[1].exp, N);
        x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
        fmpz_set(p_coeff + p_len, (poly2->coeffs + x->i)->coeffs + x->j);
        p_len++;

        assert(x->next == NULL);
        
        if (x->j + 1 < (poly2->coeffs + x->i)->length)
        {
            k = poly2->exps[x->i];
            x->j = x->j + 1;
            x->next = NULL;

            if (bits <= FLINT_BITS)
                mpoly_monomial_madd(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);
            else
                mpoly_monomial_madd_mp(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);

            _mpoly_heap_insert(heap, exp + N*x->i, x, &next_loc, &heap_len, N,
                                                               cmpmask);
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

void fmpz_mpoly_univar_assert_canonical(fmpz_mpoly_univar_t poly, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (poly->length > 0 && (slong)(poly->exps[0]) < WORD(0))
        flint_throw(FLINT_ERROR, "Univariate polynomial exponents overflow");

    for (i = 0; i + 1 < poly->length; i++)
    {
        if (poly->exps[i] <= poly->exps[i + 1])
            flint_throw(FLINT_ERROR, "Univariate polynomial exponents out of order");
    }

    for (i = 0; i < poly->length; i++)
        fmpz_mpoly_assert_canonical(poly->coeffs + i, ctx);
}

int fmpz_mpoly_univar_equal(fmpz_mpoly_univar_t poly1,
                   const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)
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


void fmpz_mpoly_univar_add(fmpz_mpoly_univar_t poly1,
            const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_univar_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
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
    while (i < len2 || j < len3)
    {
        if (i < len2 && (j >= len3 || exp2[i] > exp3[j]))
        {
            fmpz_mpoly_set(coeff1 + len1, coeff2 + i, ctx);
            exp1[len1++] = exp2[i++];
        }
        else if (j < len3 && (i >= len2 || exp3[j] > exp2[i]))
        {
            fmpz_mpoly_set(coeff1 + len1, coeff3 + j, ctx);
            exp1[len1++] = exp3[j++];
        }
        else if (i < len2 && j < len3 && (exp2[i] == exp3[j]))
        {
            fmpz_mpoly_add(coeff1 + len1, coeff2 + i, coeff3 + j, ctx);
            exp1[len1] = exp2[i];
            len1 += !fmpz_mpoly_is_zero(coeff1 + len1, ctx);
            i++;
            j++;
        } else 
        {
            assert(0);
        }
    }

    /* demote remaining coefficients */
    for (i = len1; i < poly1->length; i++)
    {
        fmpz_mpoly_clear(poly1->coeffs + i, ctx);
        fmpz_mpoly_init(poly1->coeffs + i, ctx);
    }
    poly1->length = len1;
}


int fmpz_mpoly_univar_mul(fmpz_mpoly_univar_t poly1,
        const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_univar_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
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
    {
        goto done;
    }

    if ((slong)(poly2->exps[0] + poly3->exps[0]) < WORD(0))
    {
        success = 0;
        goto done;
    }

    if (poly1 == poly2 || poly1 == poly3)
    {
        fmpz_mpoly_univar_t temp;
        fmpz_mpoly_univar_init(temp, ctx);
        success = fmpz_mpoly_univar_mul(temp, poly2, poly3, ctx);
        fmpz_mpoly_univar_swap(poly1, temp, ctx);
        fmpz_mpoly_univar_clear(temp, ctx);
        return success;
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

    _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x, &next_loc, &heap_len,
                                                                      WORD(0));

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
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                &next_loc, &heap_len, WORD(0));
            }
            if (j + 1 < len3)
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                &next_loc, &heap_len, WORD(0));
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

    return success;
}

/*
    A = prem(A, -B)
    C is used for working space
*/
void _fmpz_mpoly_univar_prem(fmpz_mpoly_univar_t polyA,
            const fmpz_mpoly_univar_t polyB, fmpz_mpoly_univar_t polyC,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong a_len, b_len, c_len;
    slong a_deg, b_deg;
    ulong * a_exp, * b_exp, * c_exp;
    fmpz_mpoly_struct * a_coeff, * b_coeff, * c_coeff;
    slong i, j, delta, delta_org;
    fmpz_mpoly_t u, v;

    assert(polyA != polyB);
    assert(polyB != polyC);
    assert(polyC != polyA);
    assert(polyA->length != 0);
    assert(polyB->length != 0);
    assert(polyA->exps[0] >= polyB->exps[0]);

    delta_org = polyA->exps[0] - polyB->exps[0] + 1;

    fmpz_mpoly_init(u, ctx);
    fmpz_mpoly_init(v, ctx);

    fmpz_mpoly_univar_fit_length(polyA, polyA->length + polyB->exps[0], ctx);
    fmpz_mpoly_univar_fit_length(polyC, polyA->length + polyB->exps[0], ctx);

    b_len = polyB->length;
    b_deg = polyB->exps[0];
    b_exp = polyB->exps;
    b_coeff = polyB->coeffs;

looper:
    a_len = polyA->length;
    a_deg = polyA->exps[0];
    a_exp = polyA->exps;
    a_coeff = polyA->coeffs;

    c_exp = polyC->exps;
    c_coeff = polyC->coeffs;

    delta = a_deg - b_deg;

    if (a_len == 0 || delta < 0)
        goto done;

    c_len = 0;
    i = 1;
    j = 1;
    while (i < a_len || j < b_len)
    {
        if (i < a_len && (j >= b_len || a_exp[i] > b_exp[j] + delta))
        {
            fmpz_mpoly_mul_johnson(c_coeff + c_len, a_coeff + i, b_coeff + 0, ctx);
            fmpz_mpoly_neg(c_coeff + c_len, c_coeff + c_len, ctx);
            c_exp[c_len++] = a_exp[i++];
        }
        else if (j < b_len && (i >= a_len || b_exp[j] + delta > a_exp[i]))
        {
            fmpz_mpoly_mul_johnson(c_coeff + c_len, a_coeff + 0, b_coeff + j, ctx);
            c_exp[c_len++] = b_exp[j++] + delta;
        }
        else if (i < a_len && j < b_len && a_exp[i] == b_exp[j] + delta)
        {
            fmpz_mpoly_mul_johnson(u, a_coeff + i, b_coeff + 0, ctx);
            fmpz_mpoly_mul_johnson(v, a_coeff + 0, b_coeff + j, ctx);
            fmpz_mpoly_sub(c_coeff + c_len, v, u, ctx);
            c_exp[c_len] = a_exp[i];
            c_len += !fmpz_mpoly_is_zero(c_coeff + c_len, ctx);
            i++;
            j++;
        } else
        {
            assert(0);
        }
    }

    polyC->length = c_len;
    fmpz_mpoly_univar_swap(polyA, polyC, ctx);
    delta_org--;
    goto looper;

done:
    if (delta_org != 0)
    {
        assert(delta_org > 0);
        fmpz_mpoly_neg(v, b_coeff + 0, ctx);
        fmpz_mpoly_pow_fps(u, v, delta_org, ctx);
        for (i = 0; i < polyA->length; i++)
            fmpz_mpoly_mul_johnson(polyA->coeffs + i, polyA->coeffs + i, u, ctx);
    }

    fmpz_mpoly_clear(u, ctx);
    fmpz_mpoly_clear(v, ctx);
}




void _fmpz_mpoly_univar_pgcd(fmpz_mpoly_univar_t poly1,
            const fmpz_mpoly_univar_t polyP, const fmpz_mpoly_univar_t polyQ,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, d, e;
    fmpz_mpoly_t u, v, w, s;
    fmpz_mpoly_univar_t A, B, C, D;
    fmpz_mpoly_univar_struct * last;

    assert(polyP->length != 0);
    assert(polyQ->length != 0);
    assert(polyP->exps[0] >= polyQ->exps[0]);
    assert(polyQ->exps[0] >= 1);


    fmpz_mpoly_init(u,ctx);
    fmpz_mpoly_init(v,ctx);
    fmpz_mpoly_init(w,ctx);
    fmpz_mpoly_init(s,ctx);
    fmpz_mpoly_univar_init(A,ctx);
    fmpz_mpoly_univar_init(B,ctx);
    fmpz_mpoly_univar_init(C,ctx);
    fmpz_mpoly_univar_init(D,ctx);
    i = FLINT_MAX(polyP->exps[0], polyQ->exps[0]);
    fmpz_mpoly_univar_fit_length(A, 1 + i, ctx);
    fmpz_mpoly_univar_fit_length(B, 1 + i, ctx);
    fmpz_mpoly_univar_fit_length(C, 1 + i, ctx);
    fmpz_mpoly_univar_fit_length(D, 1 + i, ctx);

    fmpz_mpoly_univar_set(B, polyP, ctx);
    fmpz_mpoly_univar_set(A, polyQ, ctx);

    last = A;

    C->var = polyP->var;
    D->var = polyP->var;

    fmpz_mpoly_pow_fps(s, polyQ->coeffs + 0, polyP->exps[0] - polyQ->exps[0], ctx);

    _fmpz_mpoly_univar_prem(B, A, D, ctx);

looper:

    d = A->exps[0]; e = B->exps[0];
    if (B->length == 0)
        goto done;

    if (d - e > 1)
    {
        fmpz_mpoly_pow_fps(u, B->coeffs + 0, d - e - 1, ctx);
        fmpz_mpoly_pow_fps(v, s, d - e - 1, ctx);
        for (i = 0; i < B->length; i++)
        {
            fmpz_mpoly_mul_johnson(w, u, B->coeffs + i, ctx);
            fmpz_mpoly_divides_monagan_pearce(C->coeffs + i, w, v, ctx);
            C->exps[i] = B->exps[i];
        }
        C->length = B->length;
        last = C;
        fmpz_mpoly_mul_johnson(w, s, A->coeffs + 0, ctx);
        fmpz_mpoly_mul_johnson(u, v, w, ctx);
    } else {
        for (i = 0; i < B->length; i++)
        {
            fmpz_mpoly_set(C->coeffs + i, B->coeffs + i, ctx);
            C->exps[i] = B->exps[i];
        }
        C->length = B->length;
        fmpz_mpoly_mul_johnson(u, s, A->coeffs + 0, ctx);
    }

    last = C;

    if (e == 0)
        goto done;

    _fmpz_mpoly_univar_prem(A, B, D, ctx);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_divides_monagan_pearce(B->coeffs + i, A->coeffs + i, u, ctx);
        B->exps[i] = A->exps[i];
    }
    B->length = A->length;

    fmpz_mpoly_univar_swap(A, C, ctx);
    fmpz_mpoly_set(s, A->coeffs + 0, ctx);

    last = A;

    goto looper;

done:

    fmpz_mpoly_univar_swap(poly1, last, ctx);

    fmpz_mpoly_clear(u,ctx);
    fmpz_mpoly_clear(v,ctx);
    fmpz_mpoly_clear(w,ctx);
    fmpz_mpoly_clear(s,ctx);
    fmpz_mpoly_univar_clear(A,ctx);
    fmpz_mpoly_univar_clear(B,ctx);
    fmpz_mpoly_univar_clear(C,ctx);
    fmpz_mpoly_univar_clear(D,ctx);
    return;
}



void _fmpz_mpoly_univar_pgcd_ducos(fmpz_mpoly_univar_t poly1,
        const fmpz_mpoly_univar_t polyP, const fmpz_mpoly_univar_t polyQ,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    ulong exp;
    slong i, j, k, d, e;
    slong alpha, n, J, aJ, ae;
    slong a_len, b_len, c_len, d_len, h_len, t_len;
    ulong * a_exp, * b_exp, * c_exp, * d_exp, * h_exp, * t_exp;
    fmpz_mpoly_struct * a_coeff, * b_coeff, * c_coeff, * d_coeff, * h_coeff,
                                                                     * t_coeff;
    int iexists, jexists, kexists;
    fmpz_mpoly_t u, v, w, s;
    fmpz_mpoly_univar_t A, B, C, D, H, T;
    fmpz_mpoly_univar_struct * last;

    assert(polyP->length != 0);
    assert(polyQ->length != 0);
    assert(polyP->exps[0] >= polyQ->exps[0]);
    assert(polyQ->exps[0] >= 1);

    fmpz_mpoly_init(u,ctx);
    fmpz_mpoly_init(v,ctx);
    fmpz_mpoly_init(w,ctx);
    fmpz_mpoly_init(s,ctx);
    fmpz_mpoly_univar_init(A,ctx);
    fmpz_mpoly_univar_init(B,ctx);
    fmpz_mpoly_univar_init(C,ctx);
    fmpz_mpoly_univar_init(D,ctx);
    fmpz_mpoly_univar_init(H,ctx);
    fmpz_mpoly_univar_init(T,ctx);
    i = FLINT_MAX(polyP->exps[0], polyQ->exps[0]);
    fmpz_mpoly_univar_fit_length(A, 1 + i, ctx);
    fmpz_mpoly_univar_fit_length(B, 1 + i, ctx);
    fmpz_mpoly_univar_fit_length(C, 1 + i, ctx);
    fmpz_mpoly_univar_fit_length(D, 1 + i, ctx);
    fmpz_mpoly_univar_fit_length(H, 1 + i, ctx);
    fmpz_mpoly_univar_fit_length(T, 1 + i, ctx);

    fmpz_mpoly_univar_set(B, polyP, ctx);
    fmpz_mpoly_univar_set(A, polyQ, ctx);
    C->var = polyP->var;
    D->var = polyP->var;
    H->var = polyP->var;
    T->var = polyP->var;

    last = A;

    fmpz_mpoly_pow_fps(s, polyQ->coeffs + 0, polyP->exps[0] - polyQ->exps[0], ctx);

    _fmpz_mpoly_univar_prem(B, A, D, ctx);

looper:

    d = A->exps[0]; e = B->exps[0];
    if (B->length == 0)
        goto done;

    last = B;

    if (d - e == 1)
    {
        a_len = A->length;
        a_exp = A->exps;
        a_coeff = A->coeffs;

        b_len = B->length;
        b_exp = B->exps;
        b_coeff = B->coeffs;

        d_len = D->length;
        d_exp = D->exps;
        d_coeff = D->coeffs;

        if (e == 0)
            goto done;

        /* D = (B[e]*A - A[e]*B)/A[d] */
        /*           i        j       */
        i = 1;
        j = 1;
        if (a_len > 1 && a_exp[1] == e)
            i++;
        else
            j = b_len;
        d_len = 0;
        while (i < a_len || j < b_len)
        {
            if (i < a_len && j < b_len && a_exp[i] == b_exp[j])
            {
                fmpz_mpoly_mul_johnson(u, a_coeff + i, b_coeff + 0, ctx);
                fmpz_mpoly_mul_johnson(v, a_coeff + 1, b_coeff + j, ctx);
                fmpz_mpoly_sub(w, u, v, ctx);
                fmpz_mpoly_divides_monagan_pearce(d_coeff + d_len, w, a_coeff + 0, ctx);
                d_exp[d_len] = a_exp[i];
                d_len += !fmpz_mpoly_is_zero(d_coeff + d_len, ctx);
                i++;
                j++;                
            } else if (i < a_len && (j >= b_len || a_exp[i] > b_exp[j]))
            {
                fmpz_mpoly_mul_johnson(u, a_coeff + i, b_coeff + 0, ctx);
                fmpz_mpoly_divides_monagan_pearce(d_coeff + d_len, u, a_coeff + 0, ctx);
                d_exp[d_len++] = a_exp[i];
                i++;
            } else if (j < b_len && (i >= a_len || b_exp[j] > a_exp[i]))
            {
                fmpz_mpoly_mul_johnson(v, a_coeff + 1, b_coeff + j, ctx);
                fmpz_mpoly_divides_monagan_pearce(d_coeff + d_len, v, a_coeff + 0, ctx);
                fmpz_mpoly_neg(d_coeff + d_len, d_coeff + d_len, ctx);
                d_exp[d_len++] = b_exp[j];
                j++;
            } else
            {
                assert(0);
            }
        }
        D->length = d_len;

        /* A = (B[e]*(D - B*x) + B[e-1]*B)/s */
        /*            i    j            k    */
        i = 0;
        if (b_len > 1 && b_exp[1] == e - 1) {
            j = 2;            
            k = 1;
        } else {
            j = 1;
            k = b_len;
        }
        a_len = 0;
        while (i < d_len || j < b_len || k < b_len)
        {
            exp = 0;
            if (i < d_len)
                exp = FLINT_MAX(exp, d_exp[i]);
            if (j < b_len)
                exp = FLINT_MAX(exp, b_exp[j] + 1);
            if (k < b_len)
                exp = FLINT_MAX(exp, b_exp[k]);

            a_exp[a_len] = exp;

            iexists = (i < d_len) && (exp == d_exp[i]);
            jexists = (j < b_len) && (exp == b_exp[j] + 1);
            kexists = (k < b_len) && (exp == b_exp[k]);

            assert(iexists || jexists || kexists);

            if (iexists) {
                if (jexists) {
                    fmpz_mpoly_sub(w, d_coeff + i, b_coeff + j, ctx);
                    fmpz_mpoly_mul_johnson(u, b_coeff + 0, w, ctx);
                } else {
                    fmpz_mpoly_mul_johnson(u, b_coeff + 0, d_coeff + i, ctx);
                }
                if (kexists) {
                    fmpz_mpoly_mul_johnson(v, b_coeff + 1, b_coeff + k, ctx);
                    fmpz_mpoly_add(w, u, v, ctx);
                    fmpz_mpoly_divides_monagan_pearce(a_coeff + a_len, w, s, ctx);
                } else {
                    fmpz_mpoly_divides_monagan_pearce(a_coeff + a_len, u, s, ctx);
                }
            } else {
                if (kexists) {
                    fmpz_mpoly_mul_johnson(u, b_coeff + 1, b_coeff + k, ctx);
                    if (jexists) {
                        fmpz_mpoly_mul_johnson(v, b_coeff + 0, b_coeff + j, ctx);
                        fmpz_mpoly_sub(w, u, v, ctx);
                        fmpz_mpoly_divides_monagan_pearce(a_coeff + a_len, w, s, ctx);
                    } else {
                        fmpz_mpoly_divides_monagan_pearce(a_coeff + a_len, u, s, ctx);
                    }
                } else {
                    fmpz_mpoly_mul_johnson(u, b_coeff + 0, b_coeff + j, ctx);
                    fmpz_mpoly_divides_monagan_pearce(a_coeff + a_len, u, s, ctx);
                    fmpz_mpoly_neg(a_coeff + a_len, a_coeff + a_len, ctx);
                }
            }

            a_len += !fmpz_mpoly_is_zero(a_coeff + a_len, ctx);

            i += iexists;
            j += jexists;
            k += kexists;
        }
        A->length = a_len;

        /* A <-> B */
        fmpz_mpoly_univar_swap(A, B, ctx);

        fmpz_mpoly_set(s, A->coeffs + 0, ctx);
        last = A;

    } else {

        a_len = A->length;
        a_exp = A->exps;
        a_coeff = A->coeffs;
        b_len = B->length;
        b_exp = B->exps;
        b_coeff = B->coeffs;
        c_len = C->length;
        c_exp = C->exps;
        c_coeff = C->coeffs;
        d_len = D->length;
        d_exp = D->exps;
        d_coeff = D->coeffs;
        h_len = H->length;
        h_exp = H->exps;
        h_coeff = H->coeffs;
        t_len = T->length;
        t_exp = T->exps;
        t_coeff = T->coeffs;

        n = d - e - 1;
        assert(n > 0);
        alpha = 1;
        while (2*alpha <= n)
            alpha = 2*alpha;

        fmpz_mpoly_set(u, b_coeff + 0, ctx);
        n = n - alpha;
        while (alpha > 1)
        {
            alpha = alpha/2;
            fmpz_mpoly_mul_johnson(v, u, u, ctx);
            fmpz_mpoly_divides_monagan_pearce(u, v, s, ctx);
            if (n >= alpha)
            {
                fmpz_mpoly_mul_johnson(v, u, b_coeff + 0, ctx);
                fmpz_mpoly_divides_monagan_pearce(u, v, s, ctx);                
                n = n - alpha;
            }
        }
        for (i = 0; i < b_len; i++)
        {
            fmpz_mpoly_mul_johnson(v, u, b_coeff + i, ctx);
            fmpz_mpoly_divides_monagan_pearce(c_coeff + i, v, s, ctx);
            c_exp[i] = b_exp[i];
        }
        c_len = b_len;
        C ->length = c_len;

        last = C;

        if (e == 0)
            goto done;

        /* H = C - C[e]*x^e */
        for (i = 1; i < c_len; i++)
        {
            fmpz_mpoly_set(h_coeff + i - 1, c_coeff + i, ctx);
            h_exp[i - 1] = c_exp[i];
        }
        h_len = c_len - 1;
        H->length = h_len;

        /* D = C[e]*A - A[e]*H  (truncated to powers of x < e) */
        i = 0;
        j = h_len;
        ae = a_len;
        while (i < a_len && a_exp[i] >= e)
        {
            if (a_exp[i] == e)
            {
                j = 0;
                ae = i;
            }
            i++;
        }
        d_len = 0;
        while (i < a_len || j < h_len)
        {
            if (i < a_len && j < h_len && a_exp[i] == h_exp[j])
            {
                fmpz_mpoly_mul_johnson(u, a_coeff + i, c_coeff + 0, ctx);
                fmpz_mpoly_mul_johnson(v, a_coeff + ae, h_coeff + j, ctx);
                fmpz_mpoly_sub(d_coeff + d_len, u, v, ctx);
                d_exp[d_len] = a_exp[i];
                d_len += !fmpz_mpoly_is_zero(d_coeff + d_len, ctx);
                i++;
                j++;                
            } else if (i < a_len && (j >= h_len || a_exp[i] > h_exp[j]))
            {
                fmpz_mpoly_mul_johnson(d_coeff + d_len, a_coeff + i, c_coeff + 0, ctx);
                d_exp[d_len++] = a_exp[i];
                i++;
            } else if (j < h_len && (i >= a_len || h_exp[j] > a_exp[i]))
            {
                fmpz_mpoly_mul_johnson(d_coeff + d_len, a_coeff + ae, h_coeff + j, ctx);
                fmpz_mpoly_neg(d_coeff + d_len, d_coeff + d_len, ctx);
                d_exp[d_len++] = h_exp[j];
                j++;
            } else
            {
                assert(0);
            }
        }
        D->length = d_len;


        for (J = e + 1; J < d; J++)
        {

            if (h_len == 0)
                break;

            /* H = H*x - H[e-1]*B/B[e] */
            if (h_exp[0] == e - 1)
            {
                i = 1;
                j = 1;
                t_len = 0;
                while (i < h_len || j < b_len)
                {
                    if (i < h_len && j < b_len && h_exp[i] + 1 == b_exp[j])
                    {
                        fmpz_mpoly_mul_johnson(u, h_coeff + 0, b_coeff + j, ctx);
                        fmpz_mpoly_divides_monagan_pearce(v, u, b_coeff + 0, ctx);
                        fmpz_mpoly_sub(t_coeff + t_len, h_coeff + i, v, ctx);
                        t_exp[t_len] = b_exp[j];
                        t_len += !fmpz_mpoly_is_zero(t_coeff + t_len, ctx);
                        i++;
                        j++;                
                    } else if (i < h_len && (j >= b_len || h_exp[i] + 1 > b_exp[j]))
                    {
                        fmpz_mpoly_swap(t_coeff + t_len, h_coeff + i, ctx);
                        t_exp[t_len++] = h_exp[i] + 1;
                        i++;
                    } else if (j < b_len && (i >= h_len || b_exp[j] > h_exp[i] + 1))
                    {
                        fmpz_mpoly_mul_johnson(u, h_coeff + 0, b_coeff + j, ctx);
                        fmpz_mpoly_divides_monagan_pearce(t_coeff + t_len, u, b_coeff + 0, ctx);
                        fmpz_mpoly_neg(t_coeff + t_len, t_coeff + t_len, ctx);
                        t_exp[t_len++] = b_exp[j];
                        j++;                
                    } else
                    {
                        assert(0);
                    }
                }
                T->length = t_len;

                fmpz_mpoly_univar_swap(H, T, ctx);

                h_len = H->length;
                h_exp = H->exps;
                h_coeff = H->coeffs;

                t_len = T->length;
                t_exp = T->exps;
                t_coeff = T->coeffs;
                                
            } else
            {
                assert(h_exp[0] < e - 1);
                for (i = 0; i < h_len; i++)
                    h_exp[i]++;
            }

            /* find coefficient of x^J in A */
            aJ = 0;
            while (aJ < a_len && a_exp[aJ] != J)
                aJ++;
            if (aJ >= a_len)
                continue;

            /* D = D - A[J]*H */
            i = 0;
            j = 0;
            t_len = 0;
            while (i < d_len || j < h_len)
            {
                if (i < d_len && j < h_len && d_exp[i] == h_exp[j])
                {
                    fmpz_mpoly_mul_johnson(u, h_coeff + j, a_coeff + aJ, ctx);
                    fmpz_mpoly_sub(t_coeff + t_len, d_coeff + i, u, ctx);
                    t_exp[t_len] = d_exp[i];
                    t_len += !fmpz_mpoly_is_zero(t_coeff + t_len, ctx);
                    i++;
                    j++;                
                } else if (i < d_len && (j >= h_len || d_exp[i] > h_exp[j]))
                {
                    fmpz_mpoly_swap(t_coeff + t_len, d_coeff + i, ctx);
                    t_exp[t_len++] = d_exp[i];
                    i++;
                } else if (j < h_len && (i >= d_len || h_exp[j] > d_exp[i]))
                {
                    fmpz_mpoly_mul_johnson(t_coeff + t_len, h_coeff + j, a_coeff + aJ, ctx);
                    fmpz_mpoly_neg(t_coeff + t_len, t_coeff + t_len, ctx);
                    t_exp[t_len++] = h_exp[j];
                    j++;
                } else
                {
                    assert(0);
                }
            }
            T->length = t_len;
            fmpz_mpoly_univar_swap(D, T, ctx);
            d_len = D->length;
            d_exp = D->exps;
            d_coeff = D->coeffs;
            t_len = T->length;
            t_exp = T->exps;
            t_coeff = T->coeffs;
        }

        /* B = (-1)^(d-e+1) * (B[e]*(D/A[d] - H*x) +  H[e-1]*B)/s */
        i = 0;
        if (h_len > 0 && h_exp[0] == e - 1) {
            j = 1;
            k = 1;
        } else {
            j = 0;
            k = b_len;
        }
        t_len = 0;
        while (i < d_len || j < h_len || k < b_len)
        {
            exp = 0;
            if (i < d_len)
                exp = FLINT_MAX(exp, d_exp[i]);
            if (j < h_len)
                exp = FLINT_MAX(exp, h_exp[j] + 1);
            if (k < b_len)
                exp = FLINT_MAX(exp, b_exp[k]);

            t_exp[t_len] = exp;

            iexists = (i < d_len && exp == d_exp[i]);
            jexists = (j < h_len && exp == h_exp[j] + 1);
            kexists = (k < b_len && exp == b_exp[k]);

            assert(iexists || jexists || kexists);

            if (iexists) {
                if (jexists) {
                    fmpz_mpoly_divides_monagan_pearce(u, d_coeff + i, a_coeff + 0, ctx);
                    fmpz_mpoly_sub(w, u, h_coeff + j, ctx);
                    fmpz_mpoly_mul_johnson(u, b_coeff + 0, w, ctx);
                } else {
                    fmpz_mpoly_divides_monagan_pearce(u, d_coeff + i, a_coeff + 0, ctx);
                    fmpz_mpoly_mul_johnson(u, b_coeff + 0, u, ctx);
                }
                if (kexists) {
                    fmpz_mpoly_mul_johnson(v, h_coeff + 0, b_coeff + k, ctx);
                    fmpz_mpoly_add(w, u, v, ctx);
                    fmpz_mpoly_divides_monagan_pearce(t_coeff + t_len, w, s, ctx);
                } else {
                    fmpz_mpoly_divides_monagan_pearce(t_coeff + t_len, u, s, ctx);
                }
            } else {
                if (kexists) {
                    fmpz_mpoly_mul_johnson(u, h_coeff + 0, b_coeff + k, ctx);
                    if (jexists) {
                        fmpz_mpoly_mul_johnson(v, b_coeff + 0, h_coeff + j, ctx);
                        fmpz_mpoly_sub(w, u, v, ctx);
                        fmpz_mpoly_divides_monagan_pearce(t_coeff + t_len, w, s, ctx);
                    } else {
                        fmpz_mpoly_divides_monagan_pearce(t_coeff + t_len, u, s, ctx);
                    }
                } else {
                    fmpz_mpoly_mul_johnson(u, b_coeff + 0, h_coeff + j, ctx);
                    fmpz_mpoly_divides_monagan_pearce(t_coeff + t_len, u, s, ctx);
                    fmpz_mpoly_neg(t_coeff + t_len, t_coeff + t_len, ctx);
                }
            }

            if (((d - e) & 1) == 0)
                fmpz_mpoly_neg(t_coeff + t_len, t_coeff + t_len, ctx);            

            t_len += !fmpz_mpoly_is_zero(t_coeff + t_len, ctx);

            i += iexists;
            j += jexists;
            k += kexists;
        }
        T->length = t_len;
        fmpz_mpoly_univar_swap(B, T, ctx);
        b_len = T->length;
        b_exp = T->exps;
        b_coeff = T->coeffs;
        t_len = T->length;
        t_exp = T->exps;
        t_coeff = T->coeffs;

        /* A <-> C */
        fmpz_mpoly_univar_swap(A, C, ctx);
        a_len = A->length;
        a_exp = A->exps;
        a_coeff = A->coeffs;
        c_len = C->length;
        c_exp = C->exps;
        c_coeff = C->coeffs;

        fmpz_mpoly_set(s, A->coeffs + 0, ctx);
        last = A;
    }

    goto looper;

done:

    fmpz_mpoly_univar_swap(poly1, last, ctx);

    fmpz_mpoly_clear(u,ctx);
    fmpz_mpoly_clear(v,ctx);
    fmpz_mpoly_clear(w,ctx);
    fmpz_mpoly_clear(s,ctx);
    fmpz_mpoly_univar_clear(A,ctx);
    fmpz_mpoly_univar_clear(B,ctx);
    fmpz_mpoly_univar_clear(C,ctx);
    fmpz_mpoly_univar_clear(D,ctx);
    fmpz_mpoly_univar_clear(H,ctx);
    fmpz_mpoly_univar_clear(T,ctx);
    return;
}



void fmpz_mpoly_univar_derivative(fmpz_mpoly_univar_t poly1,
                   const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_struct * coeff1, * coeff2;
    ulong * exp1, * exp2;
    slong len1, len2;

    poly1->var = poly2->var;

    len2 = poly2->length;
    coeff2 = poly2->coeffs;
    exp2 = poly2->exps;

    len1 = len2 - ((len2 > 0) && (exp2[len2 - 1] == 0));
    fmpz_mpoly_univar_fit_length(poly1, len1, ctx);
    coeff1 = poly1->coeffs;
    exp1 = poly1->exps;

    for (i = 0; i < len1; i++)
    {
        assert(exp2[i] > 0);
        fmpz_mpoly_scalar_mul_ui(coeff1 + i, coeff2 + i, exp2[i], ctx);
        exp1[i] = exp2[i] - 1;
    }

    /* demote remaining coefficients */
    for (i = len1; i < poly1->length; i++)
    {
        fmpz_mpoly_clear(poly1->coeffs + i, ctx);
        fmpz_mpoly_init(poly1->coeffs + i, ctx);
    }
    poly1->length = len1;
}



void fmpz_mpoly_to_fmpz_poly(fmpz_poly_t poly1, slong * poly1_shift,
               const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong i, shift, off, bits, N;
    ulong k;
    slong _shift = 0, len = poly2->length;
    fmpz * coeff = poly2->coeffs;
    ulong * exp = poly2->exps;

    bits = poly2->bits;
    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift(&off, &shift, var, N, bits, ctx->minfo);

    fmpz_poly_zero(poly1);
    if (len > 0)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        _shift = (exp[N*(len - 1)] >> shift) & mask;
        for (i = 0; i < len; i++)
        {
            k = (exp[N*i + off] >> shift) & mask;
            k -= _shift;
            assert(((slong)k) >= 0);
            fmpz_poly_set_coeff_fmpz(poly1, k, coeff + i);
        }
    }

    *poly1_shift = _shift;
}

void fmpz_mpoly_from_fmpz_poly(fmpz_mpoly_t poly1, const fmpz_poly_t poly2,
                           slong shift2, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong shift, off, bits, N;
    slong k;
    slong p_len;
    fmpz * p_coeff;
    ulong * p_exp;
    slong p_alloc;
    ulong * one;
    TMP_INIT;

    TMP_START;

    bits = fmpz_poly_degree(poly2);
    bits = 1 + FLINT_BIT_COUNT(FLINT_MAX(WORD(1), shift2 + bits));
    if (bits > FLINT_BITS)
        flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_from_fmpz_poly");
    bits = mpoly_fix_bits(bits, ctx->minfo);
    
    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(one, &off, &shift, var, N, bits, ctx->minfo);

    fmpz_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    p_alloc = poly1->alloc;
    p_len = 0;
    for (k = fmpz_poly_degree(poly2); k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&p_coeff, &p_exp, &p_alloc, p_len + 1, N);
        mpoly_monomial_mul_si(p_exp + N*p_len, one, N, k + shift2);
        fmpz_poly_get_coeff_fmpz(p_coeff + p_len, poly2, k);
        p_len += !fmpz_is_zero(p_coeff + p_len);
    }

    poly1->coeffs = p_coeff;
    poly1->exps = p_exp;
    poly1->alloc = p_alloc;
    _fmpz_mpoly_set_length(poly1, p_len, ctx);

    TMP_END;
}


