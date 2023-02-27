/*
    Copyright (C) 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_vec.h"

void _fq_vec_dot(fq_t res, const fq_struct * vec1,
         const fq_struct * vec2, slong len2, const fq_ctx_t ctx)
{
   slong i;
   fmpz_poly_t t;

   if (len2 == 0)
   {
      fq_zero(res, ctx);

      return;
   }

   fmpz_poly_init(t);

   fmpz_poly_mul(res, vec1 + 0, vec2 + 0);

   for (i = 1; i < len2; i++)
   {
      fmpz_poly_mul(t, vec1 + i, vec2 + i);

      fmpz_poly_add(res, res, t);
   }

   fq_reduce(res, ctx);

   fmpz_poly_clear(t);
}

