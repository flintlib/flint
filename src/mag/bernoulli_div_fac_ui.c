/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

const int mag_bernoulli_div_fac_ui_tab[32] =
{
    1, 536870912,
    -3, 715827883,
    -9, 763549742,
    -14, 581752185,
    -20, 930803495,
    -25, 752164440,
    -30, 609225645,
    -36, 987456906,
    -41, 800365631,
    -46, 648744507,
    -52, 1051701970,
    -57, 852476898,
    -62, 690991623,
    -67, 560096688,
    -73, 907994540,
    -78, 735992650,
};

/* computes upper bound for B_n / n! */
void
mag_bernoulli_div_fac_ui(mag_t z, ulong n)
{
    if (n % 2 == 1)
    {
        if (n == 1)
        {
            mag_one(z);
            mag_mul_2exp_si(z, z, -1);
        }
        else
        {
            mag_zero(z);
        }
    }
    else if (n < 32)
    {
        _fmpz_demote(MAG_EXPREF(z));
        MAG_EXP(z) = mag_bernoulli_div_fac_ui_tab[n];
        MAG_MAN(z) = mag_bernoulli_div_fac_ui_tab[n + 1];
    }
    else
    {
        /* upper bound for 1/(2pi) */
        mag_set_ui_2exp_si(z, 683565276, -32);

        /* 2 (1/(2pi))^n */
        mag_pow_ui(z, z, n);
        mag_mul_2exp_si(z, z, 1);

        /* we should multiply by an upper bound for zeta(n),
           but the error in the upper bound for 1/(2pi)
           used above is already larger than this factor
           when n >= 32, so there is no need */
    }
}

