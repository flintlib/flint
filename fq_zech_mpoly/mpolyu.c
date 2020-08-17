/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"


int fq_zech_mpolyu_is_canonical(const fq_zech_mpolyu_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < A->length; i++)
    {
        if ((slong)(A->exps[i]) < 0)
        {
            return 0;
        }

        if (i > 0 && A->exps[i - 1] <= A->exps[i])
        {
            return 0;
        }

        if (!fq_zech_mpoly_is_canonical(A->coeffs + i, ctx))
        {
            return 0;
        }

        if (fq_zech_mpoly_is_zero(A->coeffs + i, ctx))
        {
            return 0;
        }
    }
    return 1;
}

void fq_zech_mpolyu_init(fq_zech_mpolyu_t A, flint_bitcnt_t bits,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void fq_zech_mpolyu_clear(fq_zech_mpolyu_t A, const fq_zech_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_zech_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fq_zech_mpolyu_swap(fq_zech_mpolyu_t A, fq_zech_mpolyu_t B)
{
   fq_zech_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_zech_mpolyu_zero(fq_zech_mpolyu_t A, const fq_zech_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fq_zech_mpoly_clear(A->coeffs + i, uctx);
        fq_zech_mpoly_init(A->coeffs + i, uctx);
    }
    A->length = 0;
}

int fq_zech_mpolyu_is_one(fq_zech_mpolyu_t A, const fq_zech_mpoly_ctx_t uctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
        return 0;

    return fq_zech_mpoly_is_one(A->coeffs + 0, uctx);
}

void fq_zech_mpolyu_print_pretty(const fq_zech_mpolyu_t poly,
                                const char ** x, const fq_zech_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fq_zech_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fq_zech_mpolyu_fit_length(fq_zech_mpolyu_t A, slong length,
                                                const fq_zech_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fq_zech_mpoly_struct *) flint_malloc(
                                       new_alloc*sizeof(fq_zech_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fq_zech_mpoly_struct *) flint_realloc(A->coeffs,
                                       new_alloc*sizeof(fq_zech_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fq_zech_mpoly_init(A->coeffs + i, uctx);
            fq_zech_mpoly_fit_bits(A->coeffs + i, A->bits, uctx);
            (A->coeffs + i)->bits = A->bits;
        }
        A->alloc = new_alloc;
    }
}

void fq_zech_mpolyu_one(fq_zech_mpolyu_t A, const fq_zech_mpoly_ctx_t uctx)
{
    fq_zech_mpolyu_fit_length(A, 1, uctx);
    A->exps[0] = 0;
    fq_zech_mpoly_one(A->coeffs + 0, uctx);
    A->length = 1;
}

/* if the coefficient doesn't exist, a new one is created (and set to zero) */
fq_zech_mpoly_struct * _fq_zech_mpolyu_get_coeff(fq_zech_mpolyu_t A,
                                     ulong pow, const fq_zech_mpoly_ctx_t uctx)
{
    slong i, j;
    fq_zech_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
        if (A->exps[i] == pow) 
        {
            return A->coeffs + i;
        }
    }

    fq_zech_mpolyu_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        fq_zech_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
    }
    
    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}

