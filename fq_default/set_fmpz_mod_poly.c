/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default.h"

void fq_default_set_fmpz_mod_poly(fq_default_t op,
                        const fmpz_mod_poly_t poly, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      nmod_poly_t p;
      ulong mod = fmpz_get_ui(fq_zech_ctx_prime(ctx->ctx.fq_zech));
      nmod_poly_init(p, mod);
      fmpz_mod_poly_get_nmod_poly(p, poly);
      fq_zech_set_nmod_poly(op->fq_zech, p, ctx->ctx.fq_zech);
      nmod_poly_clear(p);
   } else if (ctx->type == 2)
   {
      nmod_poly_t p;
      ulong mod = fmpz_get_ui(fq_nmod_ctx_prime(ctx->ctx.fq_nmod));
      nmod_poly_init(p, mod);
      fmpz_mod_poly_get_nmod_poly(p, poly);
      fq_nmod_set_nmod_poly(op->fq_nmod, p, ctx->ctx.fq_nmod);
      nmod_poly_clear(p);
   } else
   {
      fq_set_fmpz_mod_poly(op->fq, poly, ctx->ctx.fq);
   }
}

