/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
sp2gz_fundamental_g2(fmpz_mat_t mat, slong j)
{
    slong g = 2;
    fmpz_mat_t a, b, c, d;

    fmpz_mat_init(a, g, g);
    fmpz_mat_init(b, g, g);
    fmpz_mat_init(c, g, g);
    fmpz_mat_init(d, g, g);

    if (j < 15)
    {
        fmpz_mat_zero(a);
        fmpz_mat_one(c);
        fmpz_mat_neg(b, c);
    }
    if (15 <= j && j < 17)
    {
        fmpz_mat_one(a);
        fmpz_mat_neg(b, a);
    }
    if (17 <= j && j < 19)
    {
        fmpz_mat_zero(b);
        fmpz_one(fmpz_mat_entry(c, 0, 0));
        fmpz_set_si(fmpz_mat_entry(c, 0, 1), -1);
        fmpz_set_si(fmpz_mat_entry(c, 1, 0), -1);
        fmpz_one(fmpz_mat_entry(c, 1, 1));
    }

    switch (j)
    {
        case 0:
            fmpz_mat_zero(d);
            break;
        case 1:
            fmpz_one(fmpz_mat_entry(d, 0, 0));
            break;
        case 2:
            fmpz_set_si(fmpz_mat_entry(d, 0, 0), -1);
            break;
        case 3:
            fmpz_one(fmpz_mat_entry(d, 1, 1));
            break;
        case 4:
            fmpz_set_si(fmpz_mat_entry(d, 1, 1), -1);
            break;
        case 5:
            fmpz_mat_one(d);
            break;
        case 6:
            fmpz_mat_one(d);
            fmpz_mat_neg(d, d);
            break;
        case 7:
            fmpz_set_si(fmpz_mat_entry(d, 0, 0), -1);
            fmpz_one(fmpz_mat_entry(d, 1, 1));
            break;
        case 8:
            fmpz_one(fmpz_mat_entry(d, 0, 0));
            fmpz_set_si(fmpz_mat_entry(d, 1, 1), -1);
            break;
        case 9:
            fmpz_one(fmpz_mat_entry(d, 0, 1));
            fmpz_one(fmpz_mat_entry(d, 1, 0));
            break;
        case 10:
            fmpz_set_si(fmpz_mat_entry(d, 0, 1), -1);
            fmpz_set_si(fmpz_mat_entry(d, 1, 0), -1);
            break;
        case 11:
            fmpz_one(fmpz_mat_entry(d, 0, 0));
            fmpz_one(fmpz_mat_entry(d, 0, 1));
            fmpz_one(fmpz_mat_entry(d, 1, 0));
            break;
        case 12:
            fmpz_set_si(fmpz_mat_entry(d, 0, 0), -1);
            fmpz_set_si(fmpz_mat_entry(d, 0, 1), -1);
            fmpz_set_si(fmpz_mat_entry(d, 1, 0), -1);
            break;
        case 13:
            fmpz_one(fmpz_mat_entry(d, 0, 1));
            fmpz_one(fmpz_mat_entry(d, 1, 0));
            fmpz_one(fmpz_mat_entry(d, 1, 1));
            break;
        case 14:
            fmpz_set_si(fmpz_mat_entry(d, 0, 1), -1);
            fmpz_set_si(fmpz_mat_entry(d, 1, 0), -1);
            fmpz_set_si(fmpz_mat_entry(d, 1, 1), -1);
            break;
        case 15:
            fmpz_one(fmpz_mat_entry(c, 0, 0));
            fmpz_one(fmpz_mat_entry(d, 1, 1));
            break;
        case 16:
            fmpz_one(fmpz_mat_entry(c, 1, 1));
            fmpz_one(fmpz_mat_entry(d, 0, 0));
            break;
        case 17:
            fmpz_mat_one(a);
            fmpz_mat_one(d);
            break;
        default: //case 18:
            fmpz_mat_one(a);
            fmpz_mat_neg(a, a);
            fmpz_mat_set(d, a);
    }

    sp2gz_set_abcd(mat, a, b, c, d);

    fmpz_mat_clear(a);
    fmpz_mat_clear(b);
    fmpz_mat_clear(c);
    fmpz_mat_clear(d);
}

void
sp2gz_fundamental(fmpz_mat_t mat, slong j)
{
    slong g = sp2gz_dim(mat);

    if (g == 1)
    {
        sp2gz_j(mat);
    }
    else if (g == 2)
    {
        sp2gz_fundamental_g2(mat, j);
    }
    else
    {
        fmpz_mat_one(mat);
        fmpz_zero(fmpz_mat_entry(mat, 0, 0));
        fmpz_zero(fmpz_mat_entry(mat, g, g));
        fmpz_one(fmpz_mat_entry(mat, g, 0));
        fmpz_set_si(fmpz_mat_entry(mat, 0, g), -1);
    }
}
