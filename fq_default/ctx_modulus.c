/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default.h"

void fq_default_ctx_modulus(fmpz_mod_poly_t p, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      nmod_poly_t mod;
      nmod_poly_init(mod, fmpz_get_ui(fq_zech_ctx_prime(ctx->ctx.fq_zech)));
      nmod_poly_set(mod, fq_zech_ctx_modulus(ctx->ctx.fq_zech));
      fmpz_mod_poly_set_nmod_poly(p, mod);
      nmod_poly_clear(mod);
   } else if (ctx->type == 2)
   {
      nmod_poly_t mod;
      nmod_poly_init(mod, fmpz_get_ui(fq_nmod_ctx_prime(ctx->ctx.fq_nmod)));
      nmod_poly_set(mod, fq_nmod_ctx_modulus(ctx->ctx.fq_nmod));
      fmpz_mod_poly_set_nmod_poly(p, mod);
      nmod_poly_clear(mod);
   } else
   {
      fmpz_mod_poly_struct const * modulus = fq_ctx_modulus(ctx->ctx.fq);
      slong i, len = modulus->length;

      /* fit_length */
      if (p->alloc < len)
      {
         slong alloc = FLINT_MAX(2*p->alloc, len);
         p->coeffs = (fmpz *) flint_realloc(p->coeffs, alloc*sizeof(fmpz));
         flint_mpn_zero((mp_ptr) (p->coeffs + p->alloc), alloc - p->alloc);
      }

      _fmpz_vec_set(p->coeffs, modulus_coeffs, len);

      _fmpz_mod_poly_set_length(p, len);
   }
}

