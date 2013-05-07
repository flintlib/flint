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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"

/* TODO: combine helper functions with partition function code */

static __inline__ void
mpfr_set_fmpz(mpfr_t c, const fmpz_t b)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_set_z(c, COEFF_TO_PTR(*b), MPFR_RNDN);
    else
        mpfr_set_si(c, *b, MPFR_RNDN);
}

static __inline__ void
mpfr_mul_fmpz(mpfr_t c, mpfr_srcptr a, const fmpz_t b)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_mul_z(c, a, COEFF_TO_PTR(*b), MPFR_RNDN);
    else
        mpfr_mul_si(c, a, *b, MPFR_RNDN);
}

static __inline__ void
mpfr_div_fmpz(mpfr_t c, mpfr_srcptr a, const fmpz_t b)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_div_z(c, a, COEFF_TO_PTR(*b), MPFR_RNDN);
    else
        mpfr_div_si(c, a, *b, MPFR_RNDN);
}

void
fmpq_bsplit_get_mpfr(mpfr_t x, const fmpq_bsplit_t s)
{
    if (fmpz_is_zero(s->Q))
    {
        mpfr_set_ui(x, 0, MPFR_RNDN);
    }
    else if (fmpz_is_zero(s->B))
    {
        mpfr_set_fmpz(x, s->T);
        mpfr_div_fmpz(x, x, s->Q);
    }
    else if (fmpz_is_zero(s->D))
    {
        mpfr_set_fmpz(x, s->T);
        mpfr_div_fmpz(x, x, s->Q);
        if (!fmpz_is_one(s->B))
            mpfr_div_fmpz(x, x, s->B);
    }
    else
    {
        mpfr_set_fmpz(x, s->V);
        mpfr_div_fmpz(x, x, s->Q);
        mpfr_div_fmpz(x, x, s->D);
        if (!fmpz_is_one(s->B))
            mpfr_div_fmpz(x, x, s->B);
    }
}
