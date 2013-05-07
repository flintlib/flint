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
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"

void
fmpq_bsplit_get_fmpq(fmpq_t x, const fmpq_bsplit_t s)
{
    if (fmpz_is_zero(s->Q))
    {
        fmpq_zero(x);
    }
    else if (fmpz_is_zero(s->B))
    {
        fmpq_set_fmpz_frac(x, s->T, s->Q);
    }
    else if (fmpz_is_zero(s->D))
    {
        /* T / (B * Q) */
        fmpq_set_fmpz_frac(x, s->T, s->Q);
        fmpq_div_fmpz(x, x, s->B);
    }
    else
    {
        /* V / (D * B * Q) */
        fmpq_set_fmpz_frac(x, s->V, s->Q);
        fmpq_div_fmpz(x, x, s->D);
        fmpq_div_fmpz(x, x, s->B);
    }
}
