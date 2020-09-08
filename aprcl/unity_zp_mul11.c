/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

/*
    Computes f = g * h for p = 11. 
    g and h must be reduced by F_11 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 60.
    Resulting f reduced by F_11 cyclotomic polynomial.
*/
void
unity_zp_mul11(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    int i;

    /*
        g = (x0, ... , x9);
        h = (y0, ... , y9);
        f = (z0, ... , z9);

        x0 = t[40]; ... ; x9 = t[49];
        y0 = t[50]; ... ; y9 = t[59];
        
        for auxiliary routine 4:
        a0 = t[0]; ... ; a4 = t[4];
        c0 = t[5]; ... ; c8 = t[13];

        for auxiliary routine 3:
        a0 = t[0]; ... ; a4 = t[4];
        b0 = t[5]; ... ; a4 = b[9];
        c0 = t[10]; ... ; t[18];

        d_{1, i} = t[50 + i] for i in [0, 9];
        d_{2, i} = t[10 + i] for i in [0, 9];
        d_{3, i} = t[40 + i] for i in [0, 9].
    */

    /* set xi */
    for (i = 0; i < 10; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[40 + i], g->poly, i, g->ctx);

    /* set yi */
    for (i = 0; i < 10; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[50 + i], h->poly, i, h->ctx);

    /* set ai = xi + x{i + 5} and bi = yi + y{i + 5} */
    for (i = 0; i < 5; i++)
    {
        fmpz_add(t[i], t[40 + i], t[45 + i]);
        fmpz_add(t[5 + i], t[50 + i], t[55 + i]);
    }

    /* 
        apply auxiliary routine 3 with (a0, ... , a4) and (b0, ... , b4)
        store result in (c0, ... , c8)
    */
    unity_zp_ar3(t);

    /* set d_{3, i} = c_i */
    for (i = 0; i < 9; i++)
        fmpz_set(t[40 + i], t[10 + i]);

    /* set ai = xi and bi = yi */
    for (i = 0; i < 5; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(t[i], g->poly, i, g->ctx);
        fmpz_mod_poly_get_coeff_fmpz(t[5 + i], h->poly, i, h->ctx);
    }

    /* 
        apply auxiliary routine 3 with (a0, ... , a4) and (b0, ... , b4)
        store result in (c0, ... , c8)
    */
    unity_zp_ar3(t);

    /* set d_{1, i} = c_i */
    for (i = 0; i < 9; i++)
        fmpz_set(t[50 + i], t[10 + i]);

    /* set ai = x{i + 5} and bi = y{i + 5} */ 
    for (i = 0; i < 5; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(t[i], g->poly, 5 + i, g->ctx);
        fmpz_mod_poly_get_coeff_fmpz(t[5 + i], h->poly, 5 + i, h->ctx);
    }

    /* 
        apply auxiliary routine 3 with (a0, ... , a4) and (b0, ... , b4)
        store result in (c0, ... , c8)
    */
    unity_zp_ar3(t);

    /* now we call c_i as d_{2, i} */

    /* set d_{3, i} = d_{3, i} - d_{2, i} - d_{1, i} */
    for (i = 0; i < 9; i++)
    {
        fmpz_sub(t[40 + i], t[40 + i], t[10 + i]);
        fmpz_sub(t[40 + i], t[40 + i], t[50 + i]);
    }

    /* a1 = d_{2, 0} + d_{3, 5} */
    fmpz_add(t[1], t[10], t[45]);

    /* d_{1, i} += d_{2, i} */
    for (i = 0; i < 8; i++)
        fmpz_add(t[50 + i], t[50 + i], t[11 + i]);

    /* d_{1, i} += d_{3, i + 6}, i in 0, 1, 2 */
    for (i = 0; i < 3; i++)
        fmpz_add(t[50 + i], t[50 + i], t[46 + i]);

    /* d_{1, i} += d_{3, i - 5}, i in 5, 6, 7, 8 */
    for (i = 5; i < 9; i++)
        fmpz_add(t[50 + i], t[50 + i], t[35 + i]);

    /* yi = d_{1, i} - a1 */
    for (i = 0; i < 9; i++)
    {
        fmpz_sub(t[0], t[50 + i], t[1]);
        unity_zp_coeff_set_fmpz(f, i, t[0]);
    }
    /* a0 = d_{3, 4} - a1 */
    fmpz_sub(t[0], t[44], t[1]);
    unity_zp_coeff_set_fmpz(f, 9, t[0]);    /*  y9 = a0 mod n   */
}

