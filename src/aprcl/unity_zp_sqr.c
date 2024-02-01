/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "aprcl.h"

void
unity_zp_sqr(unity_zp f, const unity_zp g)
{
    if (g->poly->length == 0)
    {
        fmpz_mod_poly_zero(f->poly, f->ctx);
        return;
    }

    fmpz_mod_poly_fit_length(f->poly, g->poly->length * 2 - 1, f->ctx);

    _fmpz_poly_sqr(f->poly->coeffs, g->poly->coeffs, g->poly->length);
    _fmpz_mod_poly_set_length(f->poly, 2 * g->poly->length - 1);

    _unity_zp_reduce_cyclotomic_divmod(f);
}

void
unity_zp_sqr_inplace(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /* squaring for p^k = 4 */
    if (f->p == 2 && f->exp == 2)
    {
        unity_zp_sqr4(f, g, t);
        return;
    }

    /* squaring for p^k = 8 */
    if (f->p == 2 && f->exp == 3)
    {
        unity_zp_sqr8(f, g, t);
        return;
    }

    /* squaring for p^k = 16 */
    if (f->p == 2 && f->exp == 4)
    {
        unity_zp_sqr16(f, g, t);
        return;
    }

    /* squaring for p^k = 3 */
    if (f->p == 3 && f->exp == 1)
    {
        unity_zp_sqr3(f, g, t);
        return;
    }

    /* squaring for p^k = 9 */
    if (f->p == 3 && f->exp == 2)
    {
        unity_zp_sqr9(f, g, t);
        return;
    }

    /* squaring for p^k = 5 */
    if (f->p == 5 && f->exp == 1)
    {
        unity_zp_sqr5(f, g, t);
        return;
    }

    /* squaring for p^k = 7 */
    if (f->p == 7 && f->exp == 1)
    {
        unity_zp_sqr7(f, g, t);
        return;
    }

    /* squaring for p^k = 11 */
    if (f->p == 11 && f->exp == 1)
    {
        unity_zp_sqr11(f, g, t);
        return;
    }

    /* traditional squaring */
    unity_zp_sqr(f, g);
}

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

/*
    Computes f = g * g for p = 2^2.
    g must be reduced by F_4 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 3.
    Resulting f reduced by F_4 cyclotomic polynomial.
*/
void
unity_zp_sqr4(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1);
        f = (y0, y1);

        x0 = t[0]; x1 = t[1];
        m1 = t[2]; m2 = t[3];
        d1 = t[4].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);

    fmpz_sub(t[2], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_add(t[3], t[0], t[1]);             /*  m2 = x0 + x1    */
    fmpz_mul(t[4], t[2], t[3]);             /*  d1 = m1 * m2    */
    fmpz_add(t[2], t[0], t[0]);             /*  m1 = x0 + x0    */
    unity_zp_coeff_set_fmpz(f, 0, t[4]);    /*  y0 = d1 mod n   */
    fmpz_mul(t[4], t[2], t[1]);             /*  d1 = m1 * x1    */
    unity_zp_coeff_set_fmpz(f, 1, t[4]);    /*  y1 = d1 mod n   */
}

/*
    Computes f = g * g for p = 2^3.
    g must be reduced by F_8 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 16.
    Resulting f reduced by F_8 cyclotomic polynomial.
*/
void
unity_zp_sqr8(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3);
        f = (y0, y1, y2, y3);

        x0 = t[0]; x1 = t[1]; x2 = t[2]; x3 = t[3];
        m1 = t[4]; m2 = t[5]; m3 = t[6]; m4 = t[7];
        m5 = t[8]; m6 = t[9]; m7 = t[10]; m8 = t[11];
        d1 = t[12]; d2 = t[13]; d3 = t[14]; d4 = t[15];
        d5 = t[16]; d6 = t[17].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);

    fmpz_sub(t[4], t[0], t[2]);             /*  m1 = x0 - x2    */
    fmpz_add(t[5], t[0], t[2]);             /*  m2 = x0 + x2    */
    fmpz_sub(t[6], t[1], t[3]);             /*  m3 = x1 - x3    */
    fmpz_add(t[7], t[1], t[3]);             /*  m4 = x1 + x3    */
    fmpz_add(t[8], t[0], t[0]);             /*  m5 = x0 + x0    */
    fmpz_add(t[9], t[1], t[1]);             /*  m6 = x1 + x1    */

    fmpz_add(t[10], t[4], t[6]);            /*  m7 = m1 + m3    */
    fmpz_add(t[11], t[5], t[7]);            /*  m8 = m2 + m4    */
    fmpz_mul(t[12], t[4], t[5]);            /*  d1 = m1 * m2    */
    fmpz_mul(t[13], t[6], t[7]);            /*  d2 = m3 * m4    */
    fmpz_mul(t[14], t[9], t[3]);            /*  d3 = m6 * x3    */
    fmpz_mul(t[15], t[8], t[2]);            /*  d4 = m5 * x2    */

    fmpz_add(t[5], t[2], t[3]);             /*  m2 = x2 + x3    */
    fmpz_sub(t[16], t[12], t[14]);          /*  d5 = d1 - d3    */
    unity_zp_coeff_set_fmpz(f, 0, t[16]);   /*  y0 = d5 mod n   */
    fmpz_add(t[17], t[13], t[15]);          /*  d6 = d2 + d4    */
    unity_zp_coeff_set_fmpz(f, 2, t[17]);   /*  y2 = d6 mod n   */
    fmpz_mul(t[16], t[10], t[11]);          /*  d5 = m7 * m8    */

    fmpz_add(t[17], t[12], t[13]);          /*  d6 = d1 + d2    */
    fmpz_sub(t[13], t[16], t[17]);          /*  d2 = d5 - d6    */
    unity_zp_coeff_set_fmpz(f, 1, t[13]);   /*  y1 = d2 mod n   */
    fmpz_add(t[4], t[8], t[9]);             /*  m1 = m5 + m6    */
    fmpz_mul(t[12], t[4], t[5]);            /*  d1 = m1 * m2    */
    fmpz_add(t[17], t[14], t[15]);          /*  d6 = d3 + d4    */

    fmpz_sub(t[13], t[12], t[17]);          /*  d2 = d1 - d6    */
    unity_zp_coeff_set_fmpz(f, 3, t[13]);   /*  y3 = d2 mod n   */
}

/*
    Computes f = g * g for p = 2^4.
    g must be reduced by F_16 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 44.
    Resulting f reduced by F_16 cyclotomic polynomial.
*/
void
unity_zp_sqr16(unity_zp f, const unity_zp g, fmpz_t * t)
{
    ulong i;

    /*
        g = (x0, x1, x2, x3, x4, x5, x6, x7);
        f = (y0, y1, y2, y3, y4, y5, y6, y7);

        x0 = t[30]; x1 = t[31]; x2 = t[32]; x3 = t[33];
        x4 = t[34]; x5 = t[35]; x6 = t[37];

        a0 = t[0]; a1 = t[1]; a2 = t[2]; a3 = t[3];
        b0 = t[4]; b1 = t[5]; b2 = t[6]; b3 = t[7];
        c0 = t[8]; c1 = t[9]; c2 = t[10]; c3 = t[11];
        c4 = t[12]; c5 = t[13]; c6 = t[14];

        d0 = t[38]; d1 = t[39]; d2 = t[40]; d3 = t[41];
        d4 = t[42]; d5 = t[43]; d6 = t[44]; d7 = t[45].
    */

    /* set xi */
    for (i = 0; i < 8; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[30 + i], g->poly, i, g->ctx);

    fmpz_add(t[0], t[30], t[34]);           /*  a0 = x0 + x4    */
    fmpz_add(t[1], t[31], t[35]);           /*  a1 = x1 + x5    */
    fmpz_add(t[2], t[32], t[36]);           /*  a2 = x2 + x6    */
    fmpz_add(t[3], t[33], t[37]);           /*  a3 = x3 + x7    */
    fmpz_sub(t[4], t[30], t[34]);           /*  b0 = x0 - x4    */
    fmpz_sub(t[5], t[31], t[35]);           /*  b1 = x1 - x5    */
    fmpz_sub(t[6], t[32], t[36]);           /*  b2 = x2 - x6    */
    fmpz_sub(t[7], t[33], t[37]);           /*  b3 = x3 - x7    */

    /*
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    /* set d_i = c_i */
    for (i = 8; i < 15; i++)
        fmpz_set(t[30 + i], t[i]);

    fmpz_add(t[0], t[30], t[30]);           /*  a0 = x0 + x0    */
    fmpz_add(t[1], t[31], t[31]);           /*  a1 = x1 + x1    */
    fmpz_add(t[2], t[32], t[32]);           /*  a2 = x2 + x2    */
    fmpz_add(t[3], t[33], t[33]);           /*  a3 = x3 + x3    */
    fmpz_set(t[4], t[34]);                  /*  set b0 = x4     */
    fmpz_set(t[5], t[35]);                  /*  set b1 = x5     */
    fmpz_set(t[6], t[36]);                  /*  set b2 = x6     */
    fmpz_set(t[7], t[37]);                  /*  set b3 = x7     */

    /*
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    fmpz_sub(t[16], t[38], t[12]);          /*  d7 = d0 - c4    */
    unity_zp_coeff_set_fmpz(f, 0, t[16]);   /*  y0 = d7 mod n   */
    fmpz_sub(t[16], t[39], t[13]);          /*  d7 = d1 - c5    */
    unity_zp_coeff_set_fmpz(f, 1, t[16]);   /*  y1 = d7 mod n   */
    fmpz_sub(t[16], t[40], t[14]);          /*  d7 = d2 - c6    */
    unity_zp_coeff_set_fmpz(f, 2, t[16]);   /*  y2 = d7 mod n   */
    unity_zp_coeff_set_fmpz(f, 3, t[41]);   /*  y3 = d3 mod n   */
    fmpz_add(t[16], t[42], t[8]);           /*  d7 = d4 + c0    */
    unity_zp_coeff_set_fmpz(f, 4, t[16]);   /*  y4 = d7 mod n   */
    fmpz_add(t[16], t[43], t[9]);           /*  d7 = d5 + c1    */
    unity_zp_coeff_set_fmpz(f, 5, t[16]);   /*  y5 = d7 mod n   */
    fmpz_add(t[16], t[44], t[10]);          /*  d7 = d6 + c2    */
    unity_zp_coeff_set_fmpz(f, 6, t[16]);   /*  y6 = d7 mod n   */
    unity_zp_coeff_set_fmpz(f, 7, t[11]);   /*  y7 = c3 mod n   */
}

/*
    Computes f = g * g for p = 3.
    g must be reduced by F_3 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 3.
    Resulting f reduced by F_3 cyclotomic polynomial.
*/
void
unity_zp_sqr3(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1);
        f = (y0, y1);

        x0 = t[0]; x1 = t[1];
        m1 = t[2]; m2 = t[3];
        d1 = t[4].
    */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);

    fmpz_sub(t[2], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_add(t[3], t[0], t[1]);             /*  m2 = x0 + x1    */
    fmpz_mul(t[4], t[2], t[3]);             /*  d1 = m1 * m2    */
    fmpz_add(t[3], t[2], t[0]);             /*  m2 = m1 + m0    */
    unity_zp_coeff_set_fmpz(f, 0, t[4]);    /*  y0 = d1 mod n   */
    fmpz_mul(t[4], t[1], t[3]);             /*  d1 = x1 * m2    */
    unity_zp_coeff_set_fmpz(f, 1, t[4]);    /*  y1 = d1 mod n   */
}

/*
    Computes f = g * g for p = 3.
    g must be reduced by F_3 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 30.
    Resulting f reduced by F_3 cyclotomic polynomial.
*/
void
unity_zp_sqr9(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3, x4, x5);
        f = (y0, y1, y2, y3, y4, y5);

        x0 = t[20]; x1 = t[21]; x2 = t[22];
        x3 = t[23]; x4 = t[24]; x5 = t[25];

        a0 = t[0]; a1 = t[1]; a2 = t[2];
        b0 = t[3]; b1 = t[4]; b2 = t[5];
        c0 = t[6]; c1 = t[7]; c2 = t[8];
        c3 = t[9]; c4 = t[10];

        d0 = t[26]; d1 = t[27]; d2 = t[28];
        d3 = t[29]; d4 = t[30].
    */
    fmpz_mod_poly_get_coeff_fmpz(t[20], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[21], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[22], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[23], g->poly, 3, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[24], g->poly, 4, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[25], g->poly, 5, g->ctx);

    fmpz_sub(t[0], t[20], t[23]);           /*  a0 = x0 - x3    */
    fmpz_sub(t[1], t[21], t[24]);           /*  a1 = x1 - x4    */
    fmpz_sub(t[2], t[22], t[25]);           /*  a2 = x2 - x5    */
    fmpz_add(t[3], t[20], t[23]);           /*  b0 = x0 + x3    */
    fmpz_add(t[4], t[21], t[24]);           /*  b1 = x1 + x4    */
    fmpz_add(t[5], t[22], t[25]);           /*  b2 = x2 + x5    */

    /*
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[26], t[6]);                  /*  set d0 = c0     */
    fmpz_set(t[27], t[7]);                  /*  set d1 = c1     */
    fmpz_set(t[28], t[8]);                  /*  set d2 = c2     */
    fmpz_set(t[29], t[9]);                  /*  set d3 = c3     */
    fmpz_set(t[30], t[10]);                 /*  set d4 = c4     */

    fmpz_add(t[3], t[20], t[0]);            /*  b0 = x0 + a0    */
    fmpz_add(t[4], t[21], t[1]);            /*  b1 = x1 + a1    */
    fmpz_add(t[5], t[22], t[2]);            /*  b2 = x2 + a2    */
    fmpz_set(t[0], t[23]);                  /*  set a0 = x3     */
    fmpz_set(t[1], t[24]);                  /*  set a1 = x4     */
    fmpz_set(t[2], t[25]);                  /*  set a2 = x5     */

    /*
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_sub(t[0], t[26], t[9]);            /*  a0 = d0 - c3    */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  y0 = a0 mod n   */
    fmpz_sub(t[0], t[27], t[10]);           /*  a0 = d1 - c4    */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  y1 = a0 mod n   */
    unity_zp_coeff_set_fmpz(f, 2, t[28]);   /*  y2 = d2 mod n   */
    fmpz_add(t[0], t[29], t[6]);            /*  a0 = d3 + c0    */
    fmpz_sub(t[1], t[0], t[9]);             /*  a1 = a0 - c3    */
    unity_zp_coeff_set_fmpz(f, 3, t[1]);    /*  y3 = a1 mod n   */
    fmpz_add(t[0], t[30], t[7]);            /*  a0 = d4 + c1    */
    fmpz_sub(t[1], t[0], t[10]);            /*  a1 = a0 - c4    */
    unity_zp_coeff_set_fmpz(f, 4, t[1]);    /*  y4 = a1 mod n   */
    unity_zp_coeff_set_fmpz(f, 5, t[8]);    /*  y5 = c2 mod n   */
}

/*
    Computes f = g * g for p = 5.
    g must be reduced by F_5 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 14.
    Resulting f reduced by F_5 cyclotomic polynomial.
*/
void
unity_zp_sqr5(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3);
        f = (y0, y1, y2, y3);

        x0 = t[0]; x1 = t[1]; x2 = t[2]; x3 = t[3];
        m1 = t[4]; m2 = t[5]; m3 = t[6]; m4 = t[7];
        m5 = t[8]; m6 = t[9]; m7 = t[10]; m8 = t[11];
        d1 = t[12]; d2 = t[13]; d3 = t[14]; d4 = t[15].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);

    fmpz_sub(t[4], t[0], t[2]);             /*  m1 = x0 - x2    */
    fmpz_add(t[5], t[0], t[2]);             /*  m2 = x0 + x2    */
    fmpz_sub(t[6], t[2], t[1]);             /*  m3 = x2 - x1    */
    fmpz_sub(t[7], t[0], t[3]);             /*  m4 = x0 - x3    */
    fmpz_sub(t[8], t[1], t[0]);             /*  m5 = x1 - x0    */
    fmpz_sub(t[9], t[2], t[3]);             /*  m6 = x2 - x3    */

    fmpz_sub(t[10], t[1], t[3]);            /*  m7 = x1 - x3    */
    fmpz_add(t[11], t[3], t[3]);            /*  m8 = x3 + x3    */
    fmpz_mul(t[12], t[4], t[5]);            /*  d1 = m1 * m2    */
    fmpz_mul(t[13], t[6], t[11]);           /*  d2 = m3 * m8    */
    fmpz_add(t[14], t[12], t[13]);          /*  d3 = d1 + d2    */
    unity_zp_coeff_set_fmpz(f, 0, t[14]);   /*  y0 = d3 mod n   */
    fmpz_add(t[11], t[8], t[10]);           /*  m8 = m5 + m7    */

    fmpz_mul(t[13], t[7], t[11]);           /*  d2 = m4 * m8    */
    fmpz_add(t[15], t[12], t[13]);          /*  d4 = d1 + d2    */
    unity_zp_coeff_set_fmpz(f, 1, t[15]);   /*  y1 = d4 mod n   */
    fmpz_add(t[6], t[4], t[0]);             /*  m3 = m1 + x0    */
    fmpz_mul(t[12], t[2], t[6]);            /*  d1 = x2 * m3    */
    fmpz_sub(t[5], t[10], t[3]);            /*  m2 = m7 - x3    */
    fmpz_mul(t[13], t[5], t[1]);            /*  d2 = m2 * x1    */

    fmpz_add(t[14], t[12], t[13]);          /*  d3 = d1 + d2    */
    unity_zp_coeff_set_fmpz(f, 2, t[14]);   /*  y2 = d3 mod n   */
    fmpz_add(t[10], t[9], t[9]);            /*  m7 = m6 + m6    */
    fmpz_mul(t[13], t[10], t[8]);           /*  d2 = m7 * m5    */
    fmpz_add(t[14], t[12], t[13]);          /*  d3 = d1 + d2    */
    unity_zp_coeff_set_fmpz(f, 3, t[14]);   /*  y3 = d3 mod n   */
}

/*
    Computes f = g * g for p = 7.
    g must be reduced by F_7 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 31.
    Resulting f reduced by F_7 cyclotomic polynomial.
*/
void
unity_zp_sqr7(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, ... , x5);
        f = (y0, ... , y5);

        x0 = t[0]; x1 = t[1]; x2 = t[2];
        x3 = t[3]; x4 = t[4]; x5 = t[5];

        m1 = t[6]; m2 = t[7]; ... ; m18 = t[23];
        d1 = t[24]; ... ; d8 = t[31].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[4], g->poly, 4, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[5], g->poly, 5, g->ctx);

    fmpz_sub(t[6], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_sub(t[7], t[1], t[2]);             /*  m2 = x1 - x2    */
    fmpz_sub(t[8], t[2], t[3]);             /*  m3 = x2 - x3    */
    fmpz_sub(t[9], t[3], t[4]);             /*  m4 = x3 - x4    */
    fmpz_sub(t[10], t[5], t[4]);            /*  m5 = x5 - x4    */
    fmpz_add(t[11], t[6], t[7]);            /*  m6 = x6 + x7    */

    fmpz_add(t[12], t[7], t[8]);            /*  m7 = m2 + m3    */
    fmpz_add(t[13], t[8], t[9]);            /*  m8 = m3 + m4    */
    fmpz_sub(t[14], t[3], t[5]);            /*  m9 = x3 - x5    */
    fmpz_add(t[15], t[8], t[11]);           /*  m10 = m3 + m6   */
    fmpz_add(t[16], t[9], t[12]);           /*  m11 = m4 + m7   */

    fmpz_add(t[18], t[11], t[13]);          /*  m13 = m6 + m8   */
    fmpz_add(t[19], t[12], t[14]);          /*  m14 = m7 + m9   */
    fmpz_add(t[21], t[0], t[1]);            /*  m16 = x0 + x1   */
    fmpz_add(t[22], t[0], t[15]);           /*  m17 = m8 + m2   */
    fmpz_mul(t[24], t[3], t[22]);           /*  d1 = x3 * m17   */

    fmpz_sub(t[22], t[19], t[4]);           /*  m17 = m14 - m9  */
    fmpz_add(t[23], t[19], t[4]);           /*  m18 = m14 + m9  */
    fmpz_mul(t[25], t[22], t[23]);          /*  d2 = m17 * m18  */
    fmpz_sub(t[22], t[13], t[7]);           /*  m17 = m8 - m2   */
    fmpz_mul(t[26], t[16], t[22]);          /*  d3 = m11 * m17  */

    fmpz_add(t[22], t[19], t[14]);          /*  m17 = m14 + m9  */
    fmpz_mul(t[27], t[22], t[12]);          /*  d4 = m17 * m7   */
    fmpz_add(t[22], t[1], t[1]);            /*  m17 = x1 + x1   */
    fmpz_mul(t[28], t[22], t[11]);          /*  d5 = m17 * m6   */
    fmpz_mul(t[29], t[6], t[21]);           /*  d6 = m1 * m16   */
    fmpz_add(t[22], t[8], t[8]);            /*  m17 = m3 + m3   */

    fmpz_add(t[7], t[0], t[18]);            /*  m2 = x0 + m13   */
    fmpz_mul(t[30], t[22], t[10]);          /*  d7 = m17 * m5   */
    fmpz_add(t[31], t[24], t[25]);          /*  d8 = d1 + d2    */
    fmpz_add(t[24], t[31], t[26]);          /*  d1 = d8 + d5    */
    unity_zp_coeff_set_fmpz(f, 3, t[24]);   /*  y3 = d1 mod n   */
    fmpz_add(t[31], t[26], t[27]);          /*  d8 = d3 + d4    */

    fmpz_add(t[24], t[31], t[28]);          /*  d1 = d8 + d5    */
    unity_zp_coeff_set_fmpz(f, 1, t[24]);   /*  y1 = d1 mod n   */
    fmpz_add(t[31], t[27], t[29]);          /*  d8 = d4 + d6    */
    fmpz_add(t[24], t[31], t[30]);          /*  d1 = d8 + d7    */
    unity_zp_coeff_set_fmpz(f, 0, t[24]);   /*  y0 = d1 mod n   */
    fmpz_add(t[22], t[12], t[19]);          /*  m17 = m7 + m14  */

    fmpz_mul(t[24], t[14], t[22]);          /*  d1 = m9 * m17   */
    fmpz_sub(t[22], t[13], t[5]);           /*  m17 = m8 - x5   */
    fmpz_add(t[23], t[2], t[10]);           /*  m18 = x2 + m5   */
    fmpz_mul(t[25], t[22], t[23]);          /*  d2 = m17 * m18  */
    fmpz_mul(t[26], t[7], t[4]);            /*  d3 = m2 * x4    */
    fmpz_add(t[22], t[8], t[13]);           /*  m17 = m3 + m8   */

    fmpz_mul(t[27], t[22], t[9]);           /*  d4 = m17 * m4   */
    fmpz_add(t[22], t[6], t[6]);            /*  m17 = m1 + m1   */
    fmpz_mul(t[28], t[22], t[10]);          /*  d5 = m17 * m5   */
    fmpz_sub(t[22], t[19], t[10]);          /*  m17 = m14 - m5  */
    fmpz_mul(t[29], t[22], t[16]);          /*  d6 = m17 * m11  */

    fmpz_add(t[22], t[2], t[2]);            /*  m17 = x2 + x2   */
    fmpz_mul(t[30], t[22], t[15]);          /*  d7 = m17 * m10  */
    fmpz_add(t[31], t[24], t[25]);          /*  d8 = d1 + d2    */
    fmpz_add(t[24], t[31], t[26]);          /*  d1 = d8 + d3    */
    unity_zp_coeff_set_fmpz(f, 4, t[24]);   /*  y4 = d1 mod n   */
    fmpz_add(t[31], t[26], t[27]);          /*  d8 = d3 + d4    */

    fmpz_add(t[24], t[31], t[28]);          /*  d1 = d8 + d5    */
    unity_zp_coeff_set_fmpz(f, 5, t[24]);   /*  y5 = d1 mod n   */
    fmpz_add(t[31], t[27], t[29]);          /*  d8 = d4 + d6    */
    fmpz_add(t[24], t[31], t[30]);          /*  d1 = d8 + d7    */
    unity_zp_coeff_set_fmpz(f, 2, t[24]);   /*  y2 = d1 mod n   */
}
