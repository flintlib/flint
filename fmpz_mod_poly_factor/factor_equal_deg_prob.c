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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
fmpz_mod_poly_factor_equal_deg_prob(fmpz_mod_poly_t factor,
          flint_rand_t state, const fmpz_mod_poly_t pol, long d)
{
    fmpz_mod_poly_t a, b;
    fmpz_t exp, t;
    int res = 1;

    if (pol->length <= 1)
    {
        printf("Exception (fmpz_mod_poly_factor_equal_deg_prob): \n");
        printf("Input polynomial is linear.\n");
        abort();
    }

    fmpz_mod_poly_init(a, &pol->p);

    do {
        fmpz_mod_poly_randtest(a, state, pol->length - 1);
    } while (a->length <= 1);

    fmpz_mod_poly_gcd(factor, a, pol);

    if (factor->length != 1)
    {
        fmpz_mod_poly_clear(a);
        return 1;
    }

    fmpz_mod_poly_init(b, &pol->p);

    fmpz_init(exp);
    fmpz_pow_ui(exp, &pol->p, d);
    fmpz_sub_ui(exp, exp, 1);
    fmpz_fdiv_q_2exp(exp, exp, 1);

    fmpz_mod_poly_powmod_fmpz_binexp(b, a, exp, pol);
    fmpz_clear(exp);

    fmpz_init(t);
    fmpz_sub_ui(t, &(b->coeffs[0]), 1);
    fmpz_mod(t, t, &pol->p);
    fmpz_mod_poly_set_coeff_fmpz(b, 0, t);
    fmpz_clear(t);

    fmpz_mod_poly_gcd(factor, b, pol);

    if ((factor->length <= 1) || (factor->length == pol->length)) res = 0;

    fmpz_mod_poly_clear(a);
    fmpz_mod_poly_clear(b);

    return res;
}
