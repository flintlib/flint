/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

/* Setters ********************************************************************/

void fmpz_mod_poly_set_coeff_si(fmpz_mod_poly_t poly, slong n, slong x,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, n + 1, ctx);

    if (n + 1 > poly->length)
    {
        flint_mpn_zero((mp_ptr) (poly->coeffs + poly->length), n - poly->length);
        poly->length = n + 1;
    }

    fmpz_set_si(poly->coeffs + n, x);
    fmpz_mod(poly->coeffs + n, poly->coeffs + n, fmpz_mod_ctx_modulus(ctx));
    _fmpz_mod_poly_normalise(poly);
}

void fmpz_mod_poly_set_coeff_ui(fmpz_mod_poly_t poly, slong n, ulong x,
                                                      const fmpz_mod_ctx_t ctx)
{
    if (x == 0)
    {
       if (n >= poly->length)
          return;

       fmpz_zero(poly->coeffs + n);
    }
    else
    {
        fmpz_mod_poly_fit_length(poly, n + 1, ctx);

        if (n + 1 > poly->length)
        {
            flint_mpn_zero((mp_ptr) (poly->coeffs + poly->length), n - poly->length);
            poly->length = n + 1;
        }

        fmpz_set_ui(poly->coeffs + n, x);
        fmpz_mod(poly->coeffs + n, poly->coeffs + n, fmpz_mod_ctx_modulus(ctx));
    }

    if (n == poly->length - 1)
        _fmpz_mod_poly_normalise(poly); /* we may have set the leading coefficient to 0 */
}

void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, slong n, const fmpz_t x,
                                                      const fmpz_mod_ctx_t ctx)
{
    if (fmpz_is_zero(x))
    {
       if (n >= poly->length)
          return;

       fmpz_zero(poly->coeffs + n);
    }
    else
    {
        fmpz_mod_poly_fit_length(poly, n + 1, ctx);

        if (n + 1 > poly->length)
        {
            flint_mpn_zero((mp_ptr) (poly->coeffs + poly->length), n - poly->length);
            poly->length = n + 1;
        }

        fmpz_mod_set_fmpz(poly->coeffs + n, x, ctx);
    }

    if (n == poly->length - 1)
        _fmpz_mod_poly_normalise(poly); /* we may have set the leading coefficient to 0 */
}

/* Getters ********************************************************************/

void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly,
                                             slong n, const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
{
    if (n < poly->length)
        fmpz_set(x, poly->coeffs + n);
    else
        fmpz_zero(x);
}
