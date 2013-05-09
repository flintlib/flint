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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"

void _fmpq_pow_si(fmpz_t rnum, fmpz_t rden, 
                  const fmpz_t opnum, const fmpz_t opden, len_t e)
{
    if (e >= 0)
    {
        fmpz_pow_ui(rnum, opnum, e);
        fmpz_pow_ui(rden, opden, e);
    }
    else
    {
        if (rnum == opnum)
        {
            fmpz t;

            fmpz_pow_ui(rnum, opnum, -e);
            fmpz_pow_ui(rden, opden, -e);

            t     = *rnum;
            *rnum = *rden;
            *rden = t;
        }
        else
        {
            fmpz_pow_ui(rden, opnum, -e);
            fmpz_pow_ui(rnum, opden, -e);
        }

        if (fmpz_sgn(rden) < 0)
        {
            fmpz_neg(rnum, rnum);
            fmpz_neg(rden, rden);
        }
    }
}

void fmpq_pow_si(fmpq_t rop, const fmpq_t op, len_t e)
{
    _fmpq_pow_si(fmpq_numref(rop), fmpq_denref(rop), 
                 fmpq_numref(op), fmpq_denref(op), e);
}

