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

#include "nmod_poly.h"
#include "ulong_extras.h"

int
nmod_poly_factor_equal_deg_prob(nmod_poly_t factor,
    flint_rand_t state, const nmod_poly_t pol, len_t d)
{
    nmod_poly_t a, b, c;
    mpz_t exp;
    int res = 1;
    len_t i;

    if (pol->length <= 1)
    {
        printf("Exception (nmod_poly_factor_equal_deg_prob). \n");
        printf("Input polynomial is linear.\n");
        abort();
    }

    nmod_poly_init_preinv(a, pol->mod.n, pol->mod.ninv);

    do {
        nmod_poly_randtest(a, state, pol->length - 1);
    } while (a->length <= 1);

    nmod_poly_gcd(factor, a, pol);

    if (factor->length != 1)
    {
        nmod_poly_clear(a);
        return 1;
    }

    nmod_poly_init_preinv(b, pol->mod.n, pol->mod.ninv);

    mpz_init(exp);
    if (pol->mod.n > 2)
    {
        /* compute a^{(p^d-1)/2} rem pol */
        mpz_ui_pow_ui(exp, pol->mod.n, d);
        mpz_sub_ui(exp, exp, 1);
        mpz_tdiv_q_2exp(exp, exp, 1);

        nmod_poly_powmod_mpz_binexp(b, a, exp, pol);
    }
    else
    {
        /* compute b = (a^{2^{d-1}}+a^{2^{d-2}}+...+a^4+a^2+a) rem pol */
        nmod_poly_rem(b, a, pol);
        nmod_poly_init_preinv(c, pol->mod.n, pol->mod.ninv);
        nmod_poly_set(c, b);
        for (i = 1; i < d; i++)
        {
            /* c = a^{2^i} = (a^{2^{i-1}})^2 */
            nmod_poly_powmod_ui_binexp(c, c, 2, pol);
            nmod_poly_add(b, b, c);
        }
        nmod_poly_rem(b, b, pol);
        nmod_poly_clear(c);
    }
    mpz_clear(exp);

    b->coeffs[0] = n_submod(b->coeffs[0], 1, pol->mod.n);

    nmod_poly_gcd(factor, b, pol);

    if ((factor->length <= 1) || (factor->length == pol->length)) res = 0;

    nmod_poly_clear(a);
    nmod_poly_clear(b);

    return res;
}
