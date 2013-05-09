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
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include "flint.h"
#include "fmpz_mod_poly_factor.h"

void
fmpz_mod_poly_factor_set(fmpz_mod_poly_factor_t res,
                         const fmpz_mod_poly_factor_t fac)
{
    if (res != fac)
    {
        if (fac->num == 0)
        {
            fmpz_mod_poly_factor_clear(res);
            fmpz_mod_poly_factor_init(res);
        }
        else
        {
            len_t i;

            fmpz_mod_poly_factor_fit_length(res, fac->num);
            for (i = 0; i < fac->num; i++)
            {
                fmpz_mod_poly_set(res->poly + i, fac->poly + i);
                res->exp[i] = fac->exp[i];
            }
            for (; i < res->num; i++)
            {
                fmpz_mod_poly_zero(res->poly + i);
                res->exp[i] = 0;
            }
            res->num = fac->num;
        }
    }
}
