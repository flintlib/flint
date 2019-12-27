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

/*
    the coefficients of poly1 should be constructed with the same bitcounts
    as that of poly2
*/
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
        mpoly_gen_monomial_offset_shift_sp(one, &off, &shift,
                                                        var, bits, ctx->minfo);

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

        off = mpoly_gen_monomial_offset_mp(one, var, bits, ctx->minfo);

        poly1->length = 0;
        poly1->var = var;
        for (i = 0; i < len; i++)
        {
            fmpz_set_ui_array(c, exp + N*i + off, bits/FLINT_BITS);
            if (!fmpz_fits_si(c))
            {
                fmpz_clear(c);
                goto failed;
            }
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
void fmpz_mpoly_from_univar_bits(fmpz_mpoly_t poly1, flint_bitcnt_t bits1,
                   const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong i, bits, N;
    ulong k;
    slong next_loc, heap_len = 1;
    ulong * cmpmask;
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

    bits = bits1;
    FLINT_ASSERT(bits <= FLINT_BITS);

    if (poly2->length == 0)
    {
        fmpz_mpoly_fit_bits(poly1, bits, ctx);
        poly1->bits = bits;
        _fmpz_mpoly_set_length(poly1, 0, ctx);
        return;
    }

    TMP_START;

    /* pack everything into bits */
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(one, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    poly2_exps = (ulong **) TMP_ALLOC(poly2->length*sizeof(ulong*));
    total_len = 0;
    for (i = 0; i < poly2->length; i++)
    {
        total_len += (poly2->coeffs + i)->length;
        poly2_exps[i] = (poly2->coeffs + i)->exps;
        if (bits != (poly2->coeffs + i)->bits)
        {
            poly2_exps[i] = (ulong *) flint_malloc(
                                  N*(poly2->coeffs + i)->length*sizeof(ulong));
            if (!mpoly_repack_monomials(poly2_exps[i], bits,
                    (poly2->coeffs + i)->exps, (poly2->coeffs + i)->bits,
                                      (poly2->coeffs + i)->length, ctx->minfo))
            {
                FLINT_ASSERT(0 && "repack does not fit");
            }
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
        mpoly_monomial_madd(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);
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

        FLINT_ASSERT(x->next == NULL);

        if (x->j + 1 < (poly2->coeffs + x->i)->length)
        {
            k = poly2->exps[x->i];
            x->j = x->j + 1;
            x->next = NULL;
            mpoly_monomial_madd(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);
            _mpoly_heap_insert(heap, exp + N*x->i, x, &next_loc, &heap_len, N,
                                                               cmpmask);
        }
    }

    FLINT_ASSERT(total_len == p_len);
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

    /* pack everything into bits */
    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_pack_vec_fmpz(one, gen_fields, bits, ctx->minfo->nfields, 1);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(gen_fields + i);
        fmpz_clear(tmp_fields + i);
        fmpz_clear(max_fields + i);
    }

    poly2_exps = (ulong **) TMP_ALLOC(poly2->length*sizeof(ulong*));
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

        FLINT_ASSERT(x->next == NULL);
        
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

    FLINT_ASSERT(total_len == p_len);
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

