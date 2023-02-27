/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_poly.h"

void
_fq_nmod_poly_mul_classical(fq_nmod_struct * rop,
                                  const fq_nmod_struct * op1, slong len1,
                                  const fq_nmod_struct * op2, slong len2,
                                  const fq_nmod_ctx_t ctx)
{
    if (len1 == 1 && len2 == 1)
    {
        fq_nmod_mul(rop, op1, op2, ctx);
    }
    else
    {
        slong i, j;
        nmod_poly_t t;

	nmod_poly_init_mod(t, (rop + 0)->mod);

        /* Set res[i] = poly1[i]*poly2[0] */
	for (j = 0; j < len1; j++)
	   nmod_poly_mul(rop + j, op1 + j, op2 + 0);

        /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
        for (j = 0; j < len2 - 1; j++)
	   nmod_poly_mul(rop + len1 + j, op2 + j + 1, op1 + len1 - 1);

        /* out[i+j] += in1[i]*in2[j] */
        for (i = 0; i < len1 - 1; i++)
	{
	    for (j = 0; j < len2 - 1; j++)
	    {
		 nmod_poly_mul(t, op2 + j + 1, op1 + i);
		 nmod_poly_add(rop + i + j + 1, rop + i + j + 1, t);
	    }
	}

	for (i = 0; i < len1 + len2 - 1; i++)
	   fq_nmod_reduce(rop + i, ctx);

	nmod_poly_clear(t);
    }
}

void
fq_nmod_poly_mul_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op1,
                             const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
{
    const slong len = op1->length + op2->length - 1;

    if (op1->length == 0 || op2->length == 0)
    {
        fq_nmod_poly_zero(rop, ctx);
        return;
    }

    if (rop == op1 || rop == op2)
    {
        fq_nmod_poly_t t;

        fq_nmod_poly_init2(t, len, ctx);
        _fq_nmod_poly_mul_classical(t->coeffs, op1->coeffs, op1->length,
                                          op2->coeffs, op2->length, ctx);
        fq_nmod_poly_swap(rop, t, ctx);
        fq_nmod_poly_clear(t, ctx);
    }
    else
    {
        fq_nmod_poly_fit_length(rop, len, ctx);
        _fq_nmod_poly_mul_classical(rop->coeffs, op1->coeffs,
                                          op1->length, op2->coeffs,
                                          op2->length, ctx);
    }

    _fq_nmod_poly_set_length(rop, len, ctx);
}

