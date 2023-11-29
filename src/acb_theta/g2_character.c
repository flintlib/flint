/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "acb_theta.h"

/* See Cléry, Faber, van der Geer, "Covariants of binary sextics and modular
   forms of degree 2 with character", §12 */

static void
g2_block_coeffs_mod_2(slong * coeffs, const fmpz_mat_t w)
{
    fmpz_t x;

    fmpz_init(x);
    coeffs[0] = fmpz_mod_ui(x, fmpz_mat_entry(w, 0, 0), 2);
    coeffs[1] = fmpz_mod_ui(x, fmpz_mat_entry(w, 0, 1), 2);
    coeffs[2] = fmpz_mod_ui(x, fmpz_mat_entry(w, 1, 0), 2);
    coeffs[3] = fmpz_mod_ui(x, fmpz_mat_entry(w, 1, 1), 2);
    fmpz_clear(x);
}

static slong
g2_block_det_mod_2(slong * coeffs)
{
    return (coeffs[0] * coeffs[3] + coeffs[1] * coeffs[2]) % 2;
}

static slong
g2_character_formula(slong * a, slong * b, slong * c, slong * d)
{
    return (a[0] * c[0] + a[1] * c[0] + a[1] * c[1] + a[2] * c[2] + a[3] * c[2]
        + a[3] * c[3] + c[0] * c[1] + c[1] * c[2] + c[2] * c[3] + c[0] * d[3]
        + c[1] * d[2] + c[1] * d[3] + c[2] * d[1] + c[3] * d[0] + c[3] * d[1]) % 2;
}

static slong
g2_character_switch(slong * a, slong * b, slong * c, slong * d, int twice)
{
    slong row[4];

    if (g2_block_det_mod_2(c) == 1)
    {
        return g2_character_formula(a, b, c, d);
    }
    if (g2_block_det_mod_2(a) == 1)
    {
        return g2_character_formula(c, d, a, b);
    }
    if (g2_block_det_mod_2(d) == 1)
    {
        return g2_character_formula(b, a, d, c);
    }
    if (g2_block_det_mod_2(b) == 1)
    {
        return g2_character_formula(d, c, b, a);
    }

    if (twice)
    {
        flint_throw(FLINT_ERROR, "error: went through g2_character_switch twice\n");
    }
    row[0] = a[0];
    row[1] = a[1];
    row[2] = b[0];
    row[3] = b[1];
    a[0] = c[0];
    a[1] = c[1];
    b[0] = d[0];
    b[1] = d[1];
    c[0] = row[0];
    c[1] = row[1];
    d[0] = row[2];
    d[1] = row[3];
    return 1 - g2_character_switch(a, b, c, d, 1);
}

slong
acb_theta_g2_character(const fmpz_mat_t mat)
{
    fmpz_mat_t w;
    slong coeffs[16];
    slong j, k;

    for (j = 0; j < 2; j++)
    {
        for (k = 0; k < 2; k++)
        {
            fmpz_mat_window_init(w, mat, 2 * j, 2 * k, 2 * j + 2, 2 * k + 2);
            g2_block_coeffs_mod_2(coeffs + 4 * (2 * j + k), w);
            fmpz_mat_window_clear(w);
        }
    }
    return g2_character_switch(coeffs, coeffs + 4, coeffs + 8, coeffs + 12, 0);
}
