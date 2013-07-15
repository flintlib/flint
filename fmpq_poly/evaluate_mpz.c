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

    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpq.h"

void fmpq_poly_evaluate_mpz(mpq_t res, const fmpq_poly_t poly, const mpz_t a)
{
    fmpq_t r;
    fmpz_t b;

    fmpq_init(r);
    fmpz_init(b);
    fmpz_set_mpz(b, a);
    fmpq_poly_evaluate_fmpz(r, poly, b);
    fmpq_get_mpq(res, r);
    fmpq_clear(r);
    fmpz_clear(b);
}
