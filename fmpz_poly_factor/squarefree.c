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

    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz
   
******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, const fmpz_poly_t F)
{
    fmpz_poly_content(&(fac->c), F);
    if (fmpz_sgn(fmpz_poly_lead(F)) < 0)
        fmpz_neg(&(fac->c), &(fac->c));

    if (F->length > 1)
    {
        fmpz_poly_t f, d, t1;

        fmpz_poly_init(f);
        fmpz_poly_init(d);
        fmpz_poly_init(t1);

        fmpz_poly_scalar_divexact_fmpz(f, F, &(fac->c));

        fmpz_poly_derivative(t1, f);
        fmpz_poly_gcd(d, f, t1);

        if (d->length == 1)
        {
            fmpz_poly_factor_insert(fac, f, 1);
        }
        else
        {
            fmpz_poly_t v, w, s;
            long i;

            fmpz_poly_init(v);
            fmpz_poly_init(w);
            fmpz_poly_init(s);

            fmpz_poly_div(v, f, d);
            fmpz_poly_div(w, t1, d);

            for (i = 1; ; i++)
            {
                fmpz_poly_derivative(t1, v);
                fmpz_poly_sub(s, w, t1);

                if (s->length == 0)
                {
                    if (v->length > 1)
                        fmpz_poly_factor_insert(fac, v, i);
                    break;
                }

                fmpz_poly_gcd(d, v, s);
                fmpz_poly_div(v, v, d);
                fmpz_poly_div(w, s, d);

                if (d->length > 1)
                    fmpz_poly_factor_insert(fac, d, i);
            }

            fmpz_poly_clear(v);
            fmpz_poly_clear(w);
            fmpz_poly_clear(s);
        }
        fmpz_poly_clear(f);
        fmpz_poly_clear(d);
        fmpz_poly_clear(t1);
    }
}
