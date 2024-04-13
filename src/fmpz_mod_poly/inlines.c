/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FMPZ_MOD_POLY_INLINES_C

#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif

#include "fmpz_mod_poly.h"

#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

#include "fmpz.h"
#include "fmpz_mod.h"

/* FIXME: These should be moved to their own files */

void fmpz_mod_poly_add_si(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                             slong c, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   fmpz_set_si(d, c);

   if (fmpz_size(fmpz_mod_ctx_modulus(ctx)) > 1)
   {
      if (c < 0)
         fmpz_add(d, d, fmpz_mod_ctx_modulus(ctx));
   } else
      fmpz_mod(d, d, fmpz_mod_ctx_modulus(ctx));

   if (poly->length == 0)
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   else
   {
      fmpz_mod_poly_set(res, poly, ctx);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_cmp(res->coeffs + 0, fmpz_mod_ctx_modulus(ctx)) >= 0)
         fmpz_sub(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_sub_si(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                             slong c, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   fmpz_set_si(d, c);

   if (fmpz_size(fmpz_mod_ctx_modulus(ctx)) > 1)
   {
      if (c < 0)
         fmpz_add(d, d, fmpz_mod_ctx_modulus(ctx));
   } else
      fmpz_mod(d, d, fmpz_mod_ctx_modulus(ctx));

   if (poly->length == 0)
   {
      fmpz_sub(d, fmpz_mod_ctx_modulus(ctx), d);
      if (fmpz_cmp(d, fmpz_mod_ctx_modulus(ctx)) == 0)
         fmpz_zero(d);
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   }
   else
   {
      fmpz_mod_poly_set(res, poly, ctx);

      fmpz_sub(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_sgn(res->coeffs + 0) < 0)
         fmpz_add(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_si_sub(fmpz_mod_poly_t res, slong c,
                          const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   fmpz_set_si(d, c);

   if (fmpz_size(fmpz_mod_ctx_modulus(ctx)) > 1)
   {
      if (c < 0)
         fmpz_add(d, d, fmpz_mod_ctx_modulus(ctx));
   } else
      fmpz_mod(d, d, fmpz_mod_ctx_modulus(ctx));


   if (poly->length == 0)
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   else
   {
      fmpz_mod_poly_neg(res, poly, ctx);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_cmp(res->coeffs + 0, fmpz_mod_ctx_modulus(ctx)) >= 0)
         fmpz_sub(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_add_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                     const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   if (fmpz_sgn(c) < 0 || fmpz_cmp(c, fmpz_mod_ctx_modulus(ctx)) >= 0)
      fmpz_mod(d, c, fmpz_mod_ctx_modulus(ctx));
   else
      fmpz_set(d, c);

   if (poly->length == 0)
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   else
   {
      fmpz_mod_poly_set(res, poly, ctx);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_cmp(res->coeffs + 0, fmpz_mod_ctx_modulus(ctx)) >= 0)
         fmpz_sub(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_sub_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                      const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   if (fmpz_sgn(c) < 0 || fmpz_cmp(c, fmpz_mod_ctx_modulus(ctx)) >= 0)
      fmpz_mod(d, c, fmpz_mod_ctx_modulus(ctx));
   else
      fmpz_set(d, c);

   if (poly->length == 0)
   {
      fmpz_sub(d, fmpz_mod_ctx_modulus(ctx), d);
      if (fmpz_cmp(d, fmpz_mod_ctx_modulus(ctx)) == 0)
         fmpz_zero(d);
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   } else
   {
      fmpz_mod_poly_set(res, poly, ctx);

      fmpz_sub(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_sgn(res->coeffs + 0) < 0)
         fmpz_add(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_fmpz_sub(fmpz_mod_poly_t res, const fmpz_t c,
                          const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   if (fmpz_sgn(c) < 0 || fmpz_cmp(c, fmpz_mod_ctx_modulus(ctx)) >= 0)
      fmpz_mod(d, c, fmpz_mod_ctx_modulus(ctx));
   else
      fmpz_set(d, c);


   if (poly->length == 0)
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   else
   {
      fmpz_mod_poly_neg(res, poly, ctx);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_cmp(res->coeffs + 0, fmpz_mod_ctx_modulus(ctx)) >= 0)
         fmpz_sub(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}
