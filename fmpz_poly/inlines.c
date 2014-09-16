/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009 William Hart

******************************************************************************/

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
