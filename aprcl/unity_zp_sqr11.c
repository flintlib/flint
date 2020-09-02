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
    Computes f = g * g for p = 11. 
    g must be reduced by F_11 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 69.
    Resulting f reduced by F_11 cyclotomic polynomial.
*/
void
unity_zp_sqr11(unity_zp f, const unity_zp g, fmpz_t * t)
{
    int i;

    /*
        g = (x0, ... , x9);
        f = (y0, ... , y9);

        x0 = t[30]; ... ; x9 = t[39];
        
        for auxiliary routine 4:
        a0 = t[0]; ... ; a4 = t[4];
        c0 = t[5]; ... ; c8 = t[13];

        for auxiliary routine 3:
        a0 = t[0]; ... ; a4 = t[4];
        b0 = t[5]; ... ; a4 = b[9];
        c0 = t[10]; ... ; t[18];

        d_{1, i} = t[50 + i] for i in [0, 9];
        d_{2, i} = t[60 + i] for i in [0, 9];
        d_{3, i} = t[10 + i] for i in [0, 9].
    */

    /* set xi */
    for (i = 0; i < 10; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[30 + i], g->poly, i, g->ctx);

    fmpz_set(t[0], t[30]);                  /*  set a0 = x0     */
    fmpz_set(t[1], t[31]);                  /*  set a1 = x1     */
    fmpz_set(t[2], t[32]);                  /*  set a2 = x2     */
    fmpz_set(t[3], t[33]);                  /*  set a3 = x3     */
    fmpz_set(t[4], t[34]);                  /*  set a4 = x4     */

    /* 
        apply auxiliary routine 4 with (a0, .. , a4)
        store result in (c0, .. , c8)
    */
    unity_zp_ar4(t);

    /* set d_{1, i} = c_i */
    for (i = 0; i < 9; i++)
        fmpz_set(t[50 + i], t[5 + i]);

    fmpz_set(t[0], t[35]);                  /*  set a0 = x5     */
    fmpz_set(t[1], t[36]);                  /*  set a1 = x6     */
    fmpz_set(t[2], t[37]);                  /*  set a2 = x7     */
    fmpz_set(t[3], t[38]);                  /*  set a3 = x8     */
    fmpz_set(t[4], t[39]);                  /*  set a4 = x9     */

    /* 
        apply auxiliary routine 4 with (a0, ... , a4)
        store result in (c0, ... , c8)
    */
    unity_zp_ar4(t);

    /* set d_{2, i} = c_i */
    for (i = 0; i < 9; i++)
        fmpz_set(t[60 + i], t[5 + i]);

    fmpz_set(t[0], t[35]);                  /*  set a0 = x5     */
    fmpz_set(t[1], t[36]);                  /*  set a1 = x6     */
    fmpz_set(t[2], t[37]);                  /*  set a2 = x7     */
    fmpz_set(t[3], t[38]);                  /*  set a3 = x8     */
    fmpz_set(t[4], t[39]);                  /*  set a4 = x9     */

    fmpz_mul_2exp(t[5], t[30], 1);          /*  b0 = 2 * x0     */
    fmpz_mul_2exp(t[6], t[31], 1);          /*  b1 = 2 * x1     */
    fmpz_mul_2exp(t[7], t[32], 1);          /*  b2 = 2 * x2     */
    fmpz_mul_2exp(t[8], t[33], 1);          /*  b3 = 2 * x3     */
    fmpz_mul_2exp(t[9], t[34], 1);          /*  b4 = 2 * x4     */

    /* 
        apply auxiliary routine 3 with (a0, ... , a4) and (b0, ... , b4)
        store result in (c0, ... , c8)
    */
    unity_zp_ar3(t);

    /* now we call c_i as d_{3, i} */

    /* a1 = d_{2, 0} + d_{3, 5} */
    fmpz_add(t[1], t[60], t[15]);

    /* d_{1, i} += d_{2, i} */
    for (i = 0; i < 8; i++)
        fmpz_add(t[50 + i], t[50 + i], t[61 + i]);

    /* d_{1, i} += d_{3, i + 6}, i in 0, 1, 2 */
    for (i = 0; i < 3; i++)
        fmpz_add(t[50 + i], t[50 + i], t[16 + i]);

    /* d_{1, i} += d_{3, i - 5}, i in 5, 6, 7, 8 */
    for (i = 5; i < 9; i++)
        fmpz_add(t[50 + i], t[50 + i], t[5 + i]);

    /* yi = d_{1, i} - a1 */
    for (i = 0; i < 9; i++)
    {
        fmpz_sub(t[0], t[50 + i], t[1]);
        unity_zp_coeff_set_fmpz(f, i, t[0]);
    }

    /* a0 = d_{3, 4} - a1 */
    fmpz_sub(t[0], t[14], t[1]);
    unity_zp_coeff_set_fmpz(f, 9, t[0]);    /*  y9 = a0 mod n   */
}

