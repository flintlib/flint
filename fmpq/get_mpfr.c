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


int
fmpq_get_mpfr(mpfr_t r, const fmpq_t x, mpfr_rnd_t rnd)
{
    __mpq_struct mpq;
    fmpz p, q;
    mp_limb_t pp, qq;

    p = *fmpq_numref(x);
    q = *fmpq_denref(x);

    if (p == 0)
        return mpfr_set_ui(r, 0, rnd);

    if (COEFF_IS_MPZ(p))
        mpq._mp_num = *COEFF_TO_PTR(p);
    else
    {
        pp = FLINT_ABS(p);
        mpq._mp_num._mp_alloc = 1;
        mpq._mp_num._mp_size = (p < 0) ? -1 : 1;
        mpq._mp_num._mp_d = &pp;
    }

    if (COEFF_IS_MPZ(q))
        mpq._mp_den = *COEFF_TO_PTR(q);
    else
    {
        qq = q;
        mpq._mp_den._mp_alloc = 1;
        mpq._mp_den._mp_size = 1;
        mpq._mp_den._mp_d = &qq;
    }

    return mpfr_set_q(r, &mpq, rnd);
}
