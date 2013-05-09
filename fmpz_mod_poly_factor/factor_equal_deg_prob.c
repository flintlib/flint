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
                                    flint_rand_t state,
                                    const fmpz_mod_poly_t pol, len_t d)
{
    fmpz_mod_poly_t a, b, c;
    fmpz_t exp, t, p;
    int res = 1;
    len_t i;

    if (pol->length <= 1)
    {
        printf("Exception (fmpz_mod_poly_factor_equal_deg_prob): \n");
        printf("Input polynomial is linear.\n");
        abort();
    }

    fmpz_init_set(p, &pol->p);

    fmpz_mod_poly_init(a, p);

    do
    {
        fmpz_mod_poly_randtest(a, state, pol->length - 1);
    } while (a->length <= 1);

    fmpz_mod_poly_gcd(factor, a, pol);

    if (factor->length != 1)
    {
        fmpz_mod_poly_clear(a);
        return 1;
    }

    fmpz_mod_poly_init(b, p);

    fmpz_init(exp);
    if (fmpz_cmp_ui(p, 2) > 0)
    {
        /* compute a^{(p^d-1)/2} rem pol */
        fmpz_pow_ui(exp, p, d);
        fmpz_sub_ui(exp, exp, 1);
        fmpz_fdiv_q_2exp(exp, exp, 1);

        fmpz_mod_poly_powmod_fmpz_binexp(b, a, exp, pol);
    }
    else
    {
        /* compute b = (a^{2^{d-1}}+a^{2^{d-2}}+...+a^4+a^2+a) rem pol */
        fmpz_mod_poly_rem(b, a, pol);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_set(c, b);
        for (i = 1; i < d; i++)
        {
            /* c = a^{2^i} = (a^{2^{i-1}})^2 */
            fmpz_mod_poly_powmod_ui_binexp(c, c, 2, pol);
            fmpz_mod_poly_add(b, b, c);
        }
        fmpz_mod_poly_rem(b, b, pol);
        fmpz_mod_poly_clear(c);
    }
    fmpz_clear(exp);

    fmpz_init(t);
    fmpz_sub_ui(t, &(b->coeffs[0]), 1);
    fmpz_mod(t, t, p);
    fmpz_mod_poly_set_coeff_fmpz(b, 0, t);
    fmpz_clear(t);

    fmpz_mod_poly_gcd(factor, b, pol);

    if ((factor->length <= 1) || (factor->length == pol->length))
        res = 0;

    fmpz_mod_poly_clear(a);
    fmpz_mod_poly_clear(b);
    fmpz_clear(p);

    return res;
}
