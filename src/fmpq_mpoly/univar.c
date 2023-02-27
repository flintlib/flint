/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_univar_init(fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
}


void fmpq_mpoly_univar_clear(fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpq_mpoly_clear(A->coeffs + i, ctx);
        fmpz_clear(A->exps + i);
    }

    if (A->coeffs)
        flint_free(A->coeffs);

    if (A->exps)
        flint_free(A->exps);
}


void fmpq_mpoly_univar_fit_length(fmpq_mpoly_univar_t A,
                                      slong length, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (fmpz *) flint_malloc(new_alloc*sizeof(fmpz));
            A->coeffs = (fmpq_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpq_mpoly_struct));
        }
        else
        {
            A->exps = (fmpz *) flint_realloc(A->exps, new_alloc*sizeof(fmpz));
            A->coeffs = (fmpq_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpq_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_init(A->exps + i);
            fmpq_mpoly_init(A->coeffs + i, ctx);
        }
        A->alloc = new_alloc;
    }
}


void fmpq_mpoly_univar_assert_canonical(fmpq_mpoly_univar_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i + 1 < A->length; i++)
    {
        if (fmpz_cmp(A->exps + i, A->exps + i + 1) <= 0
            || fmpz_sgn(A->exps + i) < 0
            || fmpz_sgn(A->exps + i + 1) < 0)
        {
            flint_throw(FLINT_ERROR, "Univariate polynomial exponents out of order");
        }
    }

    for (i = 0; i < A->length; i++)
        fmpq_mpoly_assert_canonical(A->coeffs + i, ctx);
}


void fmpq_mpoly_univar_print_pretty(const fmpq_mpoly_univar_t A,
                                   const char ** x, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    if (A->length == 0)
        flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        if (i != 0)
            flint_printf("+");
        flint_printf("(");
        fmpq_mpoly_print_pretty(A->coeffs + i,x,ctx);
        flint_printf(")*X^");
        fmpz_print(A->exps + i);
    }
}


void fmpq_mpoly_to_univar(fmpq_mpoly_univar_t A, const fmpq_mpoly_t B,
                                         slong var, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_univar_t Z;
    fmpz_mpoly_univar_init(Z, ctx->zctx);
    fmpz_mpoly_to_univar(Z, B->zpoly, var, ctx->zctx);
    fmpq_mpoly_univar_fit_length(A, Z->length, ctx);
    A->length = Z->length;
    for (i = Z->length - 1; i >= 0; i--)
    {
        fmpz_swap(A->exps + i, Z->exps + i);
        fmpz_mpoly_swap(A->coeffs[i].zpoly, Z->coeffs + i, ctx->zctx);
        fmpq_set(A->coeffs[i].content, B->content);
        fmpq_mpoly_reduce(A->coeffs + i, ctx);
    }
    fmpz_mpoly_univar_clear(Z, ctx->zctx);
}

/*
    Currently this function does not work if the coefficients depend on "var".
    The assertion x->next == NULL would need to be replaced by a loop.
    Other asserts would need to be removed as well.
*/
void fmpq_mpoly_from_univar_bits(fmpq_mpoly_t A, flint_bitcnt_t Abits,
            const fmpq_mpoly_univar_t B, slong var, const fmpq_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->zctx->minfo);
    slong i;
    slong next_loc, heap_len = 1;
    ulong * cmpmask;
    slong total_len;
    mpoly_heap_s * heap;
    slong Alen;
    ulong ** Btexp;
    fmpz * Bscales;
    fmpz_t t;
    ulong * exp;
    ulong * one;
    mpoly_heap_t * chain, * x;
    TMP_INIT;

    if (B->length == 0)
    {
        fmpz_mpoly_fit_bits(A->zpoly, Abits, ctx->zctx);
        A->zpoly->bits = Abits;
        _fmpz_mpoly_set_length(A->zpoly, 0, ctx->zctx);
        fmpq_zero(A->content);
        return;
    }

    TMP_START;

    /* pack everything into Abits */

    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    Btexp = (ulong **) TMP_ALLOC(B->length*sizeof(ulong*));
    Bscales = (fmpz *) TMP_ALLOC(B->length*sizeof(ulong*));

    total_len = 0;
    fmpq_zero(A->content);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bi;

        fmpz_init(Bscales + i);

        fmpq_gcd(A->content, A->content, B->coeffs[i].content);

        Bi = B->coeffs[i].zpoly;
        total_len += Bi->length;
        Btexp[i] = Bi->exps;
        if (Abits != Bi->bits)
        {
            Btexp[i] = (ulong *) flint_malloc(N*Bi->length*sizeof(ulong));
            if (!mpoly_repack_monomials(Btexp[i], Abits,
                             Bi->exps, Bi->bits, Bi->length, ctx->zctx->minfo))
            {
                FLINT_ASSERT(0 && "repack does not fit");
            }
        }
    }

    fmpz_init(t);
    if (!fmpq_is_zero(A->content))
    {
        for (i = 0; i < B->length; i++)
        {
            _fmpq_div(Bscales + i, t, fmpq_numref(B->coeffs[i].content),
                                      fmpq_denref(B->coeffs[i].content),
                                      fmpq_numref(A->content),
                                      fmpq_denref(A->content));
            FLINT_ASSERT(fmpz_is_one(t));
        }
    }
    fmpz_clear(t);

    fmpz_mpoly_fit_length(A->zpoly, total_len, ctx->zctx);
    fmpz_mpoly_fit_bits(A->zpoly, Abits, ctx->zctx);
    A->zpoly->bits = Abits;

    next_loc = B->length + 2;
    heap = (mpoly_heap_s *) TMP_ALLOC((B->length + 1)*sizeof(mpoly_heap_s));
    exp = (ulong *) TMP_ALLOC(B->length*N*sizeof(ulong));
    chain = (mpoly_heap_t *) TMP_ALLOC(B->length*sizeof(mpoly_heap_t));

    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->zctx->minfo);

    Alen = 0;
    if (Abits <= FLINT_BITS)
    {
        mpoly_gen_monomial_sp(one, var, Abits, ctx->zctx->minfo);

        for (i = 0; i < B->length; i++)
        {
            FLINT_ASSERT(fmpz_fits_si(B->exps + i));
            x = chain + i;
            x->i = i;
            x->j = 0;
            x->next = NULL;
            mpoly_monomial_madd(exp + N*x->i, Btexp[x->i] + N*x->j,
                                             fmpz_get_si(B->exps + i), one, N);
            _mpoly_heap_insert(heap, exp + N*i, x, &next_loc, &heap_len,
                                                                   N, cmpmask);
        }

        while (heap_len > 1)
        {
            FLINT_ASSERT(Alen < A->zpoly->alloc);
            mpoly_monomial_set(A->zpoly->exps + N*Alen, heap[1].exp, N);
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            fmpz_mul(A->zpoly->coeffs + Alen, Bscales + x->i,
                                     (B->coeffs + x->i)->zpoly->coeffs + x->j);
            Alen++;

            FLINT_ASSERT(x->next == NULL);

            if (x->j + 1 < (B->coeffs + x->i)->zpoly->length)
            {
                FLINT_ASSERT(fmpz_fits_si(B->exps + x->i));
                x->j = x->j + 1;
                x->next = NULL;
                mpoly_monomial_madd(exp + N*x->i, Btexp[x->i] + N*x->j,
                                          fmpz_get_si(B->exps + x->i), one, N);
                _mpoly_heap_insert(heap, exp + N*x->i, x, &next_loc, &heap_len,
                                                                   N, cmpmask);
            }
        }
    }
    else
    {
        mpoly_gen_monomial_offset_mp(one, var, Abits, ctx->zctx->minfo);

        for (i = 0; i < B->length; i++)
        {
            x = chain + i;
            x->i = i;
            x->j = 0;
            x->next = NULL;
            mpoly_monomial_madd_fmpz(exp + N*x->i, Btexp[x->i] + N*x->j,
                                                          B->exps + i, one, N);
            _mpoly_heap_insert(heap, exp + N*i, x, &next_loc, &heap_len,
                                                                   N, cmpmask);
        }

        while (heap_len > 1)
        {
            FLINT_ASSERT(Alen < A->zpoly->alloc);
            mpoly_monomial_set(A->zpoly->exps + N*Alen, heap[1].exp, N);
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            fmpz_mul(A->zpoly->coeffs + Alen, Bscales + x->i,
                                     (B->coeffs + x->i)->zpoly->coeffs + x->j);
            Alen++;

            FLINT_ASSERT(x->next == NULL);

            if (x->j + 1 < (B->coeffs + x->i)->zpoly->length)
            {
                x->j = x->j + 1;
                x->next = NULL;
                mpoly_monomial_madd_fmpz(exp + N*x->i, Btexp[x->i] + N*x->j,
                                                       B->exps + x->i, one, N);
                _mpoly_heap_insert(heap, exp + N*x->i, x, &next_loc, &heap_len,
                                                                   N, cmpmask);
            }
        }
    }

    for (i = 0; i < B->length; i++)
    {
        fmpz_clear(Bscales + i);
        if (Btexp[i] != (B->coeffs + i)->zpoly->exps)
            flint_free(Btexp[i]);
    }
    TMP_END;

    _fmpz_mpoly_set_length(A->zpoly, Alen, ctx->zctx);
    fmpq_mpoly_reduce(A, ctx);
}


void fmpq_mpoly_from_univar(fmpq_mpoly_t A, const fmpq_mpoly_univar_t B,
                                         slong var, const fmpq_mpoly_ctx_t ctx)
{
    slong n = ctx->zctx->minfo->nfields;
    flint_bitcnt_t bits;
    slong i;
    fmpz * gen_fields, * tmp_fields, * max_fields;
    TMP_INIT;

    if (B->length == 0)
    {
        fmpq_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    /* find bits required to represent result */
    gen_fields = (fmpz *) TMP_ALLOC(ctx->zctx->minfo->nfields*sizeof(fmpz));
    tmp_fields = (fmpz *) TMP_ALLOC(ctx->zctx->minfo->nfields*sizeof(fmpz));
    max_fields = (fmpz *) TMP_ALLOC(ctx->zctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->zctx->minfo->nfields; i++)
    {
        fmpz_init(gen_fields + i);
        fmpz_init(tmp_fields + i);
        fmpz_init(max_fields + i);
    }

    mpoly_gen_fields_fmpz(gen_fields, var, ctx->zctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bi = B->coeffs[i].zpoly;
        mpoly_max_fields_fmpz(tmp_fields, Bi->exps, Bi->length, Bi->bits, ctx->zctx->minfo);
        _fmpz_vec_scalar_addmul_fmpz(tmp_fields, gen_fields, n, B->exps + i);
        _fmpz_vec_max_inplace(max_fields, tmp_fields, n);
    }
    bits = _fmpz_vec_max_bits(max_fields, n);
    bits = FLINT_MAX(MPOLY_MIN_BITS, bits + 1);
    bits = mpoly_fix_bits(bits, ctx->zctx->minfo);

    for (i = 0; i < ctx->zctx->minfo->nfields; i++)
    {
        fmpz_clear(gen_fields + i);
        fmpz_clear(tmp_fields + i);
        fmpz_clear(max_fields + i);
    }
    TMP_END;

    fmpq_mpoly_from_univar_bits(A, bits, B, var, ctx);
}

