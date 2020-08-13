/*
    Copyright (C) 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_vec.h"

void _fq_nmod_vec_dot(fq_nmod_t res, const fq_nmod_struct * vec1,
         const fq_nmod_struct * vec2, slong len2, const fq_nmod_ctx_t ctx)
{
   slong i;
   nmod_poly_t t;

   if (len2 == 0)
   {
      fq_nmod_zero(res, ctx);

      return;
   }

   nmod_poly_init(t, ctx->p);

   nmod_poly_mul(res, vec1 + 0, vec2 + 0);

   for (i = 1; i < len2; i++)
   {
      nmod_poly_mul(t, vec1 + i, vec2 + i);

      nmod_poly_add(res, res, t);
   }

   fq_nmod_reduce(res, ctx);

   nmod_poly_clear(t);
}

