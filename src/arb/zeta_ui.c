/*
    Copyright (C) 2012, 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arith.h"
#include "arb.h"

/* The constant factor is nearly optimal up to at least 300000 bits. */
static double
euler_product_cutoff(double prec)
{
    if (prec > 200 && prec < 15000)  /* This range has a slight "bulge". */
        return 0.39 * pow(prec, 0.8);
    else
        return 7 + 0.535 * prec / log(prec);
}

void
arb_zeta_ui(arb_t x, ulong n, slong prec)
{
    if (n == 0)
    {
        arb_set_si(x, -1);
        arb_mul_2exp_si(x, x, -1);
    }
    else if (n == 1)
    {
        arb_indeterminate(x);
    }
    /* fast detection of asymptotic case */
    else if (n > 0.7 * prec)
    {
        arb_zeta_ui_asymp(x, n, prec);
    }
    else
    {
        /* even */
        if (n % 2 == 0)
        {
            if (((prec < 10000) && (n < 40 + 0.11*prec)) ||
                ((prec >= 10000) && (arith_bernoulli_number_size(n) * 0.9 < prec)))
            {
                arb_zeta_ui_bernoulli(x, n, prec);
            }
            else
            {
                arb_zeta_ui_euler_product(x, n, prec);
            }
        }
        else
        {
            if (n == 3)
            {
                arb_const_apery(x, prec);
            }
            else if (n < prec * 0.0006)
            {
                /* small odd n, extremely high precision */
                arb_zeta_ui_borwein_bsplit(x, n, prec);
            }
            else if (n > euler_product_cutoff(prec))
            {
                /* large n */
                arb_zeta_ui_euler_product(x, n, prec);
            }
            else
            {
                /* fallback */
                arb_zeta_ui_vec_borwein(x, n, 1, 0, prec);
            }
        }
    }
}

