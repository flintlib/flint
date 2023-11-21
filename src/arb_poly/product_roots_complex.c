/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_poly.h"

void
_arb_poly_product_roots_complex(arb_ptr poly, arb_srcptr r, slong rn,
    acb_srcptr c, slong cn, slong prec)
{
    if (rn == 0 && cn == 0)
    {
        arb_one(poly);
    }
    else if (rn == 1 && cn == 0)
    {
        arb_neg(poly, r);
        arb_one(poly + 1);
    }
    else if (rn == 2 && cn == 0)
    {
        arb_mul(poly, r + 0, r + 1, prec);
        arb_add(poly + 1, r + 0, r + 1, prec);
        arb_neg(poly + 1, poly + 1);
        arb_one(poly + 2);
    }
    else if (rn == 3 && cn == 0)
    {
        arb_mul(poly + 1, r, r + 1, prec);
        arb_mul(poly, poly + 1, r + 2, prec);
        arb_neg(poly, poly);
        arb_add(poly + 2, r, r + 1, prec);
        arb_addmul(poly + 1, poly + 2, r + 2, prec);
        arb_add(poly + 2, poly + 2, r + 2, prec);
        arb_neg(poly + 2, poly + 2);
        arb_one(poly + 3);
    }
    else if (rn == 0 && cn == 1)
    {
        arb_mul(poly, acb_realref(c), acb_realref(c), prec);
        arb_addmul(poly, acb_imagref(c), acb_imagref(c), prec);
        arb_mul_2exp_si(poly + 1, acb_realref(c), 1);
        arb_neg(poly + 1, poly + 1);
        arb_one(poly + 2);
    }
    else if (rn == 1 && cn == 1)
    {
        arb_mul(poly + 1, acb_realref(c), acb_realref(c), prec);
        arb_addmul(poly + 1, acb_imagref(c), acb_imagref(c), prec);
        arb_mul(poly, poly + 1, r, prec);
        arb_neg(poly, poly);
        arb_mul_2exp_si(poly + 2, acb_realref(c), 1);
        arb_addmul(poly + 1, poly + 2, r, prec);
        arb_add(poly + 2, poly + 2, r, prec);
        arb_neg(poly + 2, poly + 2);
        arb_one(poly + 3);
    }
    else
    {
        slong rm, cm, rm2, cm2;
        arb_ptr tmp, tmp2;

        rm = (rn + 1) / 2;
        cm = cn / 2;
        rm2 = rn - rm;
        cm2 = cn - cm;

        tmp = _arb_vec_init(rn + 2 * cn + 2);
        tmp2 = tmp + rm + (2 * cm) + 1;

        _arb_poly_product_roots_complex(tmp, r, rm, c, cm, prec);
        _arb_poly_product_roots_complex(tmp2, r + rm, rm2, c + cm, cm2, prec);
        _arb_poly_mul_monic(poly, tmp, rm + 2 * cm + 1, tmp2, rm2 + 2 * cm2 + 1, prec);

        _arb_vec_clear(tmp, rn + 2 * cn + 2);
    }
}

void
arb_poly_product_roots_complex(arb_poly_t poly,
        arb_srcptr r, slong rn, acb_srcptr c, slong cn, slong prec)
{
    arb_poly_fit_length(poly, rn + 2 * cn + 1);
    _arb_poly_product_roots_complex(poly->coeffs, r, rn, c, cn, prec);
    _arb_poly_set_length(poly, rn + 2 * cn + 1);
}

