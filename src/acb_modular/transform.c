/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_modular.h"

void acb_modular_transform(acb_t w, const psl2z_t g, const acb_t z, slong prec)
{
#define a (&g->a)
#define b (&g->b)
#define c (&g->c)
#define d (&g->d)
#define x acb_realref(z)
#define y acb_imagref(z)

    if (fmpz_is_zero(c))
    {
        /* (az+b)/d, where we must have a = d = 1 */
        acb_add_fmpz(w, z, b, prec);
    }
    else if (fmpz_is_zero(a))
    {
        /* b/(cz+d), where -bc = 1, c = 1 => -1/(z+d) */
        acb_add_fmpz(w, z, d, prec);
        acb_inv(w, w, prec);
        acb_neg(w, w);
    }
    else if (0)
    {
        acb_t t, u;

        acb_init(t);
        acb_init(u);

        acb_set_fmpz(t, b);
        acb_addmul_fmpz(t, z, a, prec);

        acb_set_fmpz(u, d);
        acb_addmul_fmpz(u, z, c, prec);

        acb_div(w, t, u, prec);

        acb_clear(t);
        acb_clear(u);
    }
    else
    {
        /* (az+b)/(cz+d) = (re+im*i)/den where

            re = bd + (bc+ad)x + ac(x^2+y^2)
            im = (ad-bc)y
            den = c^2(x^2+y^2) + 2cdx + d^2
        */

        fmpz_t t;
        arb_t re, im, den;

        arb_init(re);
        arb_init(im);
        arb_init(den);
        fmpz_init(t);

        arb_mul(im, x, x, prec);
        arb_addmul(im, y, y, prec);

        fmpz_mul(t, b, d);
        arb_set_fmpz(re, t);
        fmpz_mul(t, b, c);
        fmpz_addmul(t, a, d);
        arb_addmul_fmpz(re, x, t, prec);
        fmpz_mul(t, a, c);
        arb_addmul_fmpz(re, im, t, prec);

        fmpz_mul(t, d, d);
        arb_set_fmpz(den, t);
        fmpz_mul(t, c, d);
        fmpz_mul_2exp(t, t, 1);
        arb_addmul_fmpz(den, x, t, prec);
        fmpz_mul(t, c, c);
        arb_addmul_fmpz(den, im, t, prec);

        fmpz_mul(t, a, d);
        fmpz_submul(t, b, c);
        arb_mul_fmpz(im, y, t, prec);

        arb_div(acb_realref(w), re, den, prec);
        arb_div(acb_imagref(w), im, den, prec);

        arb_clear(re);
        arb_clear(im);
        arb_clear(den);
        fmpz_clear(t);
    }

#undef a
#undef b
#undef c
#undef d
#undef x
#undef y
}

