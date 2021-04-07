/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default.h"

void fq_default_get_fmpz_mod_poly(fmpz_mod_poly_t poly,
                             const fq_default_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      nmod_poly_t p;
      ulong mod = fmpz_get_ui(fq_zech_ctx_prime(ctx->ctx.fq_zech));
      nmod_poly_init(p, mod);
      fq_zech_get_nmod_poly(p, op->fq_zech, ctx->ctx.fq_zech);
      fmpz_mod_poly_set_nmod_poly(poly, p);
      nmod_poly_clear(p);
   } else if (ctx->type == 2)
   {
      nmod_poly_t p;
      ulong mod = fmpz_get_ui(fq_nmod_ctx_prime(ctx->ctx.fq_nmod));
      nmod_poly_init(p, mod);
      fq_nmod_get_nmod_poly(p, op->fq_nmod, ctx->ctx.fq_nmod);
      fmpz_mod_poly_set_nmod_poly(poly, p);
      nmod_poly_clear(p);
   } else
   {
      fq_get_fmpz_mod_poly(poly, op->fq, ctx->ctx.fq);
   }
}

