/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#define FMPZ_POLY_INLINES_C

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_add_si(fmpz_poly_t res, const fmpz_poly_t poly, slong c)
{
   if (poly->length == 0)
      fmpz_poly_set_si(res, c);
   else
   {
      fmpz_poly_set(res, poly);

      if (c < 0)
         fmpz_sub_ui(res->coeffs + 0, res->coeffs + 0, -c);
      else
         fmpz_add_ui(res->coeffs + 0, res->coeffs + 0, c);

      _fmpz_poly_normalise(res);
   }
}

void fmpz_poly_sub_si(fmpz_poly_t res, const fmpz_poly_t poly, slong c)
{
   if (poly->length == 0)
      fmpz_poly_set_si(res, -c);
   else
   {
      fmpz_poly_set(res, poly);

      if (c < 0)
         fmpz_add_ui(res->coeffs + 0, res->coeffs + 0, -c);
      else
         fmpz_sub_ui(res->coeffs + 0, res->coeffs + 0, c);

      _fmpz_poly_normalise(res);
   }
}

void fmpz_poly_si_sub(fmpz_poly_t res, slong c, const fmpz_poly_t poly)
{
   if (poly->length == 0)
      fmpz_poly_set_si(res, c);
   else
   {
      fmpz_poly_neg(res, poly);

      if (c < 0)
         fmpz_sub_ui(res->coeffs + 0, res->coeffs + 0, -c);
      else
         fmpz_add_ui(res->coeffs + 0, res->coeffs + 0, c);

      _fmpz_poly_normalise(res);
   }
}

void fmpz_poly_add_fmpz(fmpz_poly_t res, const fmpz_poly_t poly, fmpz_t c)
{
   if (poly->length == 0)
      fmpz_poly_set_fmpz(res, c);
   else
   {
      fmpz_poly_set(res, poly);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, c);

      _fmpz_poly_normalise(res);
   }
}

void fmpz_poly_sub_fmpz(fmpz_poly_t res, const fmpz_poly_t poly, fmpz_t c)
{
   if (poly->length == 0)
   {
      fmpz_poly_set_fmpz(res, c);
      fmpz_neg(res->coeffs + 0, res->coeffs + 0);
   } else
   {
      fmpz_poly_set(res, poly);

      fmpz_sub(res->coeffs + 0, res->coeffs + 0, c);

      _fmpz_poly_normalise(res);
   }
}

void fmpz_poly_fmpz_sub(fmpz_poly_t res, fmpz_t c, const fmpz_poly_t poly)
{
   if (poly->length == 0)
      fmpz_poly_set_fmpz(res, c);
   else
   {
      fmpz_poly_neg(res, poly);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, c);

      _fmpz_poly_normalise(res);
   }
}
