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

******************************************************************************/

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
nmod_poly_factor_equal_deg_prob(nmod_poly_t factor,
    flint_rand_t state, const nmod_poly_t pol, ulong d)
{
    nmod_poly_t a, b;
    mpz_t exp;
    int res = 1;

    if (pol->length <= 1)
    {
        printf("Attempt to factor a linear polynomial "
                "in nmod_poly_factor_equal_deg_prob\n");
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
    mpz_ui_pow_ui(exp, pol->mod.n, d);
    mpz_sub_ui(exp, exp, 1);
    mpz_tdiv_q_2exp(exp, exp, 1);

    nmod_poly_powmod_mpz_binexp(b, a, exp, pol);
    mpz_clear(exp);

    b->coeffs[0] = n_submod(b->coeffs[0], 1, pol->mod.n);

    nmod_poly_gcd(factor, b, pol);

    if ((factor->length <= 1) || (factor->length == pol->length)) res = 0;

    nmod_poly_clear(a);
    nmod_poly_clear(b);

    return res;
}
