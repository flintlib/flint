/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default_poly.h"

void fq_default_poly_set_fmpz_poly(fq_default_poly_t rop,
                              const fmpz_poly_t op, const fq_default_ctx_t ctx)
{
   fmpz const * p;
   fmpz_mod_ctx_t mod;
   fmpz_mod_poly_t mod_poly;
   if (ctx->type == 1)
   {
      p = fq_zech_ctx_prime(ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      p = fq_nmod_ctx_prime(ctx->ctx.fq_nmod);
   } else
   {
      p = fq_ctx_prime(ctx->ctx.fq);
   }
   fmpz_mod_ctx_init(mod, p);
   fmpz_mod_poly_init(mod_poly, mod);
   fmpz_mod_poly_set_fmpz_poly(mod_poly, op, mod);
   fq_default_poly_set_fmpz_mod_poly(rop, mod_poly, ctx);
   fmpz_mod_poly_clear(mod_poly, mod);
   fmpz_mod_ctx_clear(mod);
}

