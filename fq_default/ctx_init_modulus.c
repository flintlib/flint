/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default.h"

void fq_default_ctx_init_modulus_type(fq_default_ctx_t ctx,
                const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx,
                                                    const char * var, int type)
{
   fmpz const * p = fmpz_mod_ctx_modulus(mod_ctx);
   int bits = fmpz_bits(p);
   int d = fmpz_mod_poly_degree(modulus, mod_ctx);

   if (type == 1 ||
       (type == 0 && bits*d <= 16 &&
        n_pow(fmpz_get_ui(p), d) < (UWORD(1) << 16)))
   {
      nmod_poly_t nmodulus;
      fq_nmod_ctx_struct * fq_nmod_ctx;
      ctx->type = 1;
      nmod_poly_init(nmodulus, fmpz_get_ui(p));
      fmpz_mod_poly_get_nmod_poly(nmodulus, modulus);
      fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));
      fq_nmod_ctx_init_modulus(fq_nmod_ctx, nmodulus, var);
      if (fq_zech_ctx_init_fq_nmod_ctx_check(ctx->ctx.fq_zech, fq_nmod_ctx))
         ctx->ctx.fq_zech->owns_fq_nmod_ctx = 1;
      else
      {
         *ctx->ctx.fq_nmod = *fq_nmod_ctx;
         flint_free(fq_nmod_ctx);
         ctx->type = 2;
      }
      nmod_poly_clear(nmodulus);
   } else if (type == 2 || (type == 0 && fmpz_abs_fits_ui(p)))
   {
      nmod_poly_t nmodulus;
      ctx->type = 2;
      nmod_poly_init(nmodulus, fmpz_get_ui(p));
      fmpz_mod_poly_get_nmod_poly(nmodulus, modulus);
      fq_nmod_ctx_init_modulus(ctx->ctx.fq_nmod, nmodulus, var);
      nmod_poly_clear(nmodulus);
   } else
   {
      ctx->type = 3;
      fq_ctx_init_modulus(ctx->ctx.fq, modulus, mod_ctx, var);
   }
}

void fq_default_ctx_init_modulus(fq_default_ctx_t ctx,
       const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var)
{
   fq_default_ctx_init_modulus_type(ctx, modulus, mod_ctx, var, 0);
}

