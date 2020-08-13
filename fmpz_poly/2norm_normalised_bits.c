/*
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

flint_bitcnt_t _fmpz_poly_2norm_normalised_bits(const fmpz * poly, slong len)
{
   fmpz_t norm;
   flint_bitcnt_t bits;
   fmpz_init(norm);

	_fmpz_poly_2norm(norm, poly, len);

	bits = fmpz_bits(norm);
	fmpz_clear(norm);
   
   return bits - fmpz_bits(poly + len - 1) + 1;
}
