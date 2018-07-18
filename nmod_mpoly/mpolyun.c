/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpolyun_init(nmod_mpolyun_t A, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void nmod_mpolyun_clear(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpolyn_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void nmod_mpolyun_swap(nmod_mpolyun_t A, nmod_mpolyun_t B)
{
   nmod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyun_zero(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        nmod_mpolyn_clear(A->coeffs + i, ctx);
        nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
    }
    A->length = 0;
}

void nmod_mpolyun_print_pretty(const nmod_mpolyun_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        FLINT_ASSERT((poly->coeffs + i)->bits == poly->bits);
        nmod_mpolyn_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void nmod_mpolyun_fit_length(nmod_mpolyun_t A, slong length, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpolyn_struct *) flint_malloc(
                                          new_alloc*sizeof(nmod_mpolyn_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpolyn_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
        }
        A->alloc = new_alloc;
    }
}

void nmod_mpolyun_shift_right(nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void nmod_mpolyun_shift_left(nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

slong nmod_mpolyun_lastdeg(nmod_mpolyun_t A, nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            deg = FLINT_MAX(deg, nmod_poly_degree((A->coeffs + i)->coeffs + j));
        }
    }
    FLINT_ASSERT(deg >= 0);
    return deg;
}

void nmod_mpolyun_set(nmod_mpolyun_t A, const nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpolyn_set(Acoeff + i, Bcoeff + i, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_mpolyn_clear(Acoeff + i, ctx);
        nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


void nmod_mpolyun_mul_poly(
    nmod_mpolyun_t A,
    const nmod_mpolyun_t B,
    const nmod_poly_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpolyn_mul_poly(Acoeff + i, Bcoeff + i, c, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_mpolyn_clear(Acoeff + i, ctx);
        nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


void nmod_mpolyun_content_last(nmod_poly_t a, nmod_mpolyun_t B,
                                                          nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    nmod_poly_zero(a);
    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            nmod_poly_gcd(a, a, (B->coeffs + i)->coeffs + j);
        }
    }
}

void nmod_mpolyun_divexact_last(nmod_mpolyun_t A, nmod_poly_t b,
                                                          nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    nmod_poly_t r;
    nmod_poly_init(r, ctx->ffinfo->mod.n);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            nmod_poly_divrem((A->coeffs + i)->coeffs + j, r, (A->coeffs + i)->coeffs + j, b);
            FLINT_ASSERT(nmod_poly_is_zero(r));
        }
    }
    nmod_poly_clear(r);
}
