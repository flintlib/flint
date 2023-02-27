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

void _fq_nmod_poly_mul(fq_nmod_struct * rop,
		       const fq_nmod_struct * op1, slong len1,
		       const fq_nmod_struct * op2, slong len2,
		                                  const fq_nmod_ctx_t ctx)
{
   if (len1 <= 1 || len2 <= 1 ||
      (fq_nmod_ctx_degree(ctx) == 2 && FLINT_MAX(len1, len2) <= 2) ||
      FLINT_BIT_COUNT(fmpz_get_ui(fq_nmod_ctx_prime(ctx)))*
         fq_nmod_ctx_degree(ctx)* FLINT_MAX(len1, len2) <= 8)
      _fq_nmod_poly_mul_classical(rop, op1, len1, op2, len2, ctx);
   else
      _fq_nmod_poly_mul_univariate(rop, op1, len1, op2, len2, ctx);
}

void fq_nmod_poly_mul(fq_nmod_poly_t rop, const fq_nmod_poly_t op1,
		const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
{
   slong len1 = fq_nmod_poly_length(op1, ctx);
   slong len2 = fq_nmod_poly_length(op2, ctx);

   if (len1 <= 1 || len2 <= 1 ||
      (fq_nmod_ctx_degree(ctx) == 2 && FLINT_MAX(len1, len2) <= 2) ||
      FLINT_BIT_COUNT(fmpz_get_ui(fq_nmod_ctx_prime(ctx)))*
         fq_nmod_ctx_degree(ctx)* FLINT_MAX(len1, len2) <= 8)
      fq_nmod_poly_mul_classical(rop, op1, op2, ctx);
   else
      fq_nmod_poly_mul_univariate(rop, op1, op2, ctx);
}
