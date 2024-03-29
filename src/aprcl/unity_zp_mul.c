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

#if FLINT_WANT_ASSERT
# include "fmpz_mod.h"
#endif

void
unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h)
{
    slong glen, hlen;

    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(h->ctx)));

    glen = g->poly->length;
    hlen = h->poly->length;

    if (glen == 0 || hlen == 0)
    {
        fmpz_mod_poly_zero(f->poly, f->ctx);
        return;
    }

    fmpz_mod_poly_fit_length(f->poly, glen + hlen - 1, f->ctx);

    if (glen >= hlen)
        _fmpz_poly_mul(f->poly->coeffs, g->poly->coeffs, glen, h->poly->coeffs, hlen);
    else
        _fmpz_poly_mul(f->poly->coeffs, h->poly->coeffs, hlen, g->poly->coeffs, glen);

    _fmpz_mod_poly_set_length(f->poly, glen + hlen - 1);

    _unity_zp_reduce_cyclotomic_divmod(f);
}

void
unity_zp_mul_inplace(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /* multiplication for p^k = 4 */
    if (f->p == 2 && f->exp == 2)
    {
        unity_zp_mul4(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 8 */
    if (f->p == 2 && f->exp == 3)
    {
        unity_zp_mul8(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 16 */
    if (f->p == 2 && f->exp == 4)
    {
        unity_zp_mul16(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 3 */
    if (f->p == 3 && f->exp == 1)
    {
        unity_zp_mul3(f, g, h, t);
        return;
    }

    /* multiplicatiom for p^k = 9 */
    if (f->p == 3 && f->exp == 2)
    {
        unity_zp_mul9(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 5 */
    if (f->p == 5 && f->exp == 1)
    {
        unity_zp_mul5(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 7 */
    if (f->p == 7 && f->exp == 1)
    {
        unity_zp_mul7(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 11 */
    if (f->p == 11 && f->exp == 1)
    {
        unity_zp_mul11(f, g, h, t);
        return;
    }

    /* traditional multiplication */
    unity_zp_mul(f, g, h);
}

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

/*
    Computes f = g * h for p = 2^2.
    g and h must be reduced by F_4 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 9.
    Resulting f reduced by F_4 cyclotomic polynomial.
*/
void
unity_zp_mul4(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, x1);
        h = (y0, y1);
        f = (z0, z1);

        x0 = t[0]; x1 = t[1];
        y0 = t[2]; y1 = t[3];

        m1 = t[4]; m2 = t[5]; m3 = t[6];
        d1 = t[7]; d2 = t[8]; d3 = t[9].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[2], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], h->poly, 1, h->ctx);

    fmpz_add(t[4], t[0], t[1]);             /*  m1 = x0 + x1    */
    fmpz_add(t[5], t[2], t[3]);             /*  m2 = y0 + y1    */
    fmpz_sub(t[6], t[3], t[2]);             /*  m3 = y1 - y0    */
    fmpz_mul(t[7], t[4], t[2]);             /*  d1 = m1 * y0    */
    fmpz_mul(t[8], t[5], t[1]);             /*  d2 = m2 * x1    */
    fmpz_mul(t[9], t[6], t[0]);             /*  d3 = m3 * x0    */
    fmpz_sub(t[0], t[7], t[8]);             /*  t[0] = d1 - d2  */

    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = t[0] mod n */
    fmpz_add(t[0], t[7], t[9]);             /*  t[0] = d1 + d3  */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = t[0] mod n */
}

/*
    Computes f = g * h for p = 2^3.
    g and h must be reduced by F_8 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 27.
    Resulting f reduced by F_8 cyclotomic polynomial.
*/
void
unity_zp_mul8(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3);
        h = (y0, y1, y2, y3);
        f = (z0, z1, z2, z3);

        x0 = t[0]; x1 = t[1]; x2 = t[2]; x3 = t[3];
        y0 = t[4]; y1 = t[5]; y2 = t[6]; y3 = t[7];

        m1 = t[8]; ... ; m8 = t[15];
        d1 = t[16]; ... ; d12 = t[28].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[4], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[5], h->poly, 1, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[6], h->poly, 2, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[7], h->poly, 3, h->ctx);

    fmpz_add(t[8], t[1], t[3]);             /*  m1 = x1 + x3    */
    fmpz_add(t[9], t[5], t[7]);             /*  m2 = y1 + y3    */
    fmpz_add(t[10], t[2], t[3]);            /*  m3 = x2 + x3    */
    fmpz_add(t[11], t[6], t[7]);            /*  m4 = y2 + y3    */
    fmpz_add(t[12], t[0], t[1]);            /*  m5 = x0 + x1    */
    fmpz_add(t[13], t[4], t[5]);            /*  m6 = y0 + y1    */

    fmpz_add(t[14], t[0], t[2]);            /*  m7 = x0 + x2    */
    fmpz_add(t[15], t[4], t[6]);            /*  m8 = y0 + y2    */
    fmpz_mul(t[16], t[0], t[4]);            /*  d0 = x0 * y0    */
    fmpz_mul(t[17], t[1], t[5]);            /*  d1 = x1 * y1    */
    fmpz_mul(t[18], t[2], t[6]);            /*  d2 = x2 * y2    */
    fmpz_mul(t[19], t[3], t[7]);            /*  d3 = x3 * y3    */
    fmpz_mul(t[22], t[12], t[13]);          /*  d6 = m5 * m6    */

    fmpz_mul(t[23], t[14], t[15]);          /*  d7 = m7 * m8    */
    fmpz_mul(t[24], t[8], t[9]);            /*  d8 = m1 * m2    */
    fmpz_mul(t[25], t[10], t[11]);          /*  d9 = m3 * m4    */
    fmpz_add(t[10], t[8], t[14]);           /*  m3 = m1 + m7    */
    fmpz_add(t[11], t[9], t[15]);           /*  m4 = m2 + m8    */
    fmpz_mul(t[20], t[10], t[11]);          /*  d4 = m3 * m4    */

    fmpz_add(t[26], t[16], t[17]);          /*  d10 = d0 + d1   */
    fmpz_add(t[27], t[18], t[19]);          /*  d11 = d2 + d3   */
    fmpz_add(t[28], t[26], t[19]);          /*  d12 = d10 + d3  */
    fmpz_add(t[21], t[24], t[18]);          /*  d5 = d8 + d2    */
    fmpz_sub(t[0], t[28], t[21]);           /*  t[0] = d12 - d5 */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = t[0] mod n */

    fmpz_add(t[28], t[22], t[27]);          /*  d12 = d6 + d11  */
    fmpz_add(t[21], t[26], t[25]);          /*  d5 = d10 + d9   */
    fmpz_sub(t[0], t[28], t[21]);           /*  t[0] = d12 - d5 */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = t[0] mod n */
    fmpz_add(t[28], t[17], t[23]);          /*  d12 = d1 + d7   */
    fmpz_add(t[21], t[16], t[27]);          /*  d5 = d0 + d11   */

    fmpz_sub(t[0], t[28], t[21]);           /*  t[0] = d12 - d5 */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /*  z2 = t[0] mod n */
    fmpz_add(t[28], t[23], t[22]);          /*  d12 = d7 + d6   */
    fmpz_add(t[21], t[28], t[24]);          /*  d5 = d12 + d8   */
    fmpz_add(t[28], t[21], t[25]);          /*  d12 = d5 + d8   */
    fmpz_add(t[19], t[26], t[27]);          /*  d3 = d10 + d11  */

    fmpz_add(t[21], t[19], t[20]);          /*  d5 = d3 + d4    */
    fmpz_sub(t[0], t[21], t[28]);           /*  t[0] = d5 - d12 */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /*  z3 = t[0] mod n */
}

/*
    Computes f = g * h for p = 2^4.
    g and h must be reduced by F_16 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 63.
    Resulting f reduced by F_16 cyclotomic polynomial.
*/
void
unity_zp_mul16(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    int i;

    /*
        g = (x0, ... , x7);
        h = (y0, ... , y7);
        f = (z0, ... , z7);

        x0 = t[30]; x1 = t[31]; x2 = t[32]; x3 = t[33];
        x4 = t[34]; x5 = t[35]; x6 = t[36]; x7 = t[37];

        y0 = t[40]; y1 = t[41]; y2 = t[42]; y3 = t[43];
        y4 = t[44]; y5 = t[45]; y6 = t[46]; y7 = t[47];

        a0 = t[0]; ... ; a3 = t[3];
        b0 = t[4]; ... ; b3 = t[7];
        c0 = t[8]; ... ; t[14];

        d0 = t[50]; ... ; d12 = t[62];
    */

    /* set xi */
    for (i = 0; i < 8; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[30 + i], g->poly, i, g->ctx);

    /* set yi */
    for (i = 0; i < 8; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[40 + i], h->poly, i, h->ctx);

    fmpz_add(t[0], t[30], t[34]);           /*  a0 = x0 + x4    */
    fmpz_add(t[1], t[31], t[35]);           /*  a1 = x1 + x5    */
    fmpz_add(t[2], t[32], t[36]);           /*  a2 = x2 + x6    */
    fmpz_add(t[3], t[33], t[37]);           /*  a3 = x3 + x7    */
    fmpz_set(t[4], t[40]);                  /*  set b0 = y0     */
    fmpz_set(t[5], t[41]);                  /*  set b1 = y1     */
    fmpz_set(t[6], t[42]);                  /*  set b2 = y2     */
    fmpz_set(t[7], t[43]);                  /*  set b3 = y3     */

    /*
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    /* set d_i = c_i */
    for (i = 0; i < 7; i++)
        fmpz_set(t[50 + i], t[8 + i]);

    fmpz_add(t[0], t[40], t[44]);           /*  a0 = y0 + y4    */
    fmpz_add(t[1], t[41], t[45]);           /*  a1 = y1 + y5    */
    fmpz_add(t[2], t[42], t[46]);           /*  a2 = y2 + y6    */
    fmpz_add(t[3], t[43], t[47]);           /*  a3 = y3 + y7    */
    fmpz_set(t[4], t[34]);                  /*  set b0 = x4     */
    fmpz_set(t[5], t[35]);                  /*  set b1 = x5     */
    fmpz_set(t[6], t[36]);                  /*  set b2 = x6     */
    fmpz_set(t[7], t[37]);                  /*  set b3 = x7     */

    /*
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    /* set d_{7 + i} = c_i */
    for (i = 0; i < 7; i++)
        fmpz_set(t[57 + i], t[8 + i]);

    fmpz_sub(t[0], t[44], t[40]);           /*  a0 = y4 - y0    */
    fmpz_sub(t[1], t[45], t[41]);           /*  a1 = y5 - y1    */
    fmpz_sub(t[2], t[46], t[42]);           /*  a2 = y6 - y2    */
    fmpz_sub(t[3], t[47], t[43]);           /*  a3 = y7 - y3    */
    fmpz_set(t[4], t[30]);                  /*  set b0 = x0     */
    fmpz_set(t[5], t[31]);                  /*  set b1 = x1     */
    fmpz_set(t[6], t[32]);                  /*  set b2 = x2     */
    fmpz_set(t[7], t[33]);                  /*  set b3 = x3     */

    /*
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    fmpz_add(t[1], t[54], t[57]);           /*  a1 = d4 + d7    */
    fmpz_add(t[2], t[1], t[12]);            /*  a2 = a1 + c3    */
    fmpz_sub(t[0], t[50], t[2]);            /*  a0 = d0 - a2    */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = a0 mod n   */

    fmpz_add(t[1], t[55], t[58]);           /*  a1 = d5 + d8    */
    fmpz_add(t[2], t[1], t[13]);            /*  a2 = a1 + c4    */
    fmpz_sub(t[0], t[51], t[2]);            /*  a0 = d1 - a2    */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = a0 mod n   */

    fmpz_add(t[1], t[56], t[59]);           /*  a1 = d6 + d9    */
    fmpz_add(t[2], t[1], t[14]);            /*  a2 = a1 + c5    */
    fmpz_sub(t[0], t[52], t[2]);            /*  a0 = d2 - a2    */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /*  z2 = a0 mod n   */

    fmpz_sub(t[0], t[53], t[60]);           /*  a0 = d3 - d10   */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /*  z3 = a0 mod n   */

    fmpz_add(t[1], t[54], t[50]);           /*  a1 = d4 + d0    */
    fmpz_add(t[2], t[1], t[8]);             /*  a2 = a1 + c0    */
    fmpz_sub(t[0], t[2], t[61]);            /*  a0 = a2 - d11   */
    unity_zp_coeff_set_fmpz(f, 4, t[0]);    /*  z4 = a0 mod n   */

    fmpz_add(t[1], t[55], t[51]);           /*  a1 = d5 + d1    */
    fmpz_add(t[2], t[1], t[9]);             /*  a2 = a1 + c0    */
    fmpz_sub(t[0], t[2], t[62]);            /*  a0 = a2 - d12   */
    unity_zp_coeff_set_fmpz(f, 5, t[0]);    /*  z5 = a0 mod n   */

    fmpz_add(t[1], t[56], t[52]);           /*  a1 = d6 + d2    */
    fmpz_add(t[2], t[1], t[10]);            /*  a2 = a1 + c1    */
    fmpz_sub(t[0], t[2], t[63]);            /*  a0 = a2 - a13   */
    unity_zp_coeff_set_fmpz(f, 6, t[0]);    /*  z6 = a0 mod n   */

    fmpz_add(t[1], t[53], t[11]);           /*  a1 = d3 + c2    */
    unity_zp_coeff_set_fmpz(f, 7, t[1]);    /*  z7 = a1 mod n   */
}

/*
    Computes f = g * h for p = 3.
    g and h must be reduced by F_3 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 3.
    Resulting f reduced by F_3 cyclotomic polynomial.
*/
void
unity_zp_mul3(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, x1);
        h = (y0, y1);
        f = (z0, z1);

        x0 = t[0]; x1 = t[1];
        y0 = t[2]; y1 = t[3];
        m1 = t[4]; m2 = t[5];
        d1 = t[6]; d2 = t[7]; d3 = t[8].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[2], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], h->poly, 1, h->ctx);

    fmpz_mul(t[6], t[0], t[2]);             /*  d1 = x0 * y0    */
    fmpz_mul(t[7], t[1], t[3]);             /*  d2 = x1 * y1    */
    fmpz_sub(t[4], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_sub(t[5], t[3], t[2]);             /*  m2 = y1 - y0    */
    fmpz_mul(t[8], t[4], t[5]);             /*  d3 = m1 * m2    */
    fmpz_add(t[8], t[8], t[6]);             /*  d3 = d3 + d1    */

    unity_zp_coeff_set_fmpz(f, 1, t[8]);    /*  z1 = d3 mod n   */
    fmpz_sub(t[0], t[6], t[7]);             /*  x0 = d1 - d2    */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = x0 mod n   */
}

/*
    Computes f = g * h for p = 3^2.
    g and h must be reduced by F_9 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 50.
    Resulting f reduced by F_9 cyclotomic polynomial.
*/
void
unity_zp_mul9(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, ... , x5);
        h = (y0, ... , y5);
        f = (z0, ... , z5);

        x0 = t[20]; ... ; x5 = t[25];
        y0 = t[26]; ... ; y5 = t[31];

        a0 = t[0]; a1 = t[1]; a2 = t[2];
        b0 = t[3]; b1 = t[4]; b2 = t[5];
        c0 = t[6]; ... ; c4 = t[10];

        d0 = t[32]; ... ; d19 = t[51].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[20], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[21], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[22], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[23], g->poly, 3, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[24], g->poly, 4, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[25], g->poly, 5, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[26], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[27], h->poly, 1, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[28], h->poly, 2, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[29], h->poly, 3, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[30], h->poly, 4, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[31], h->poly, 5, h->ctx);

    fmpz_set(t[0], t[20]);                  /*  set a0 = x0     */
    fmpz_set(t[1], t[21]);                  /*  set a1 = x1     */
    fmpz_set(t[2], t[22]);                  /*  set a2 = x2     */
    fmpz_set(t[3], t[26]);                  /*  set b0 = y0     */
    fmpz_set(t[4], t[27]);                  /*  set b1 = y1     */
    fmpz_set(t[5], t[28]);                  /*  set b2 = y2     */

    /*
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[32], t[6]);                  /*  set d0 = c0     */
    fmpz_set(t[33], t[7]);                  /*  set d1 = c1     */
    fmpz_set(t[34], t[8]);                  /*  set d2 = c2     */
    fmpz_set(t[35], t[9]);                  /*  set d3 = c3     */
    fmpz_set(t[36], t[10]);                 /*  set d4 = c4     */

    fmpz_set(t[0], t[23]);                  /*  set a0 = x3     */
    fmpz_set(t[1], t[24]);                  /*  set a1 = x4     */
    fmpz_set(t[2], t[25]);                  /*  set a2 = x5     */
    fmpz_set(t[3], t[29]);                  /*  set b0 = y3     */
    fmpz_set(t[4], t[30]);                  /*  set b1 = y4     */
    fmpz_set(t[5], t[31]);                  /*  set b2 = y5     */

    /*
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[38], t[6]);                  /*  set d6 = c0     */
    fmpz_set(t[39], t[7]);                  /*  set d7 = c1     */
    fmpz_set(t[40], t[8]);                  /*  set d8 = c2     */
    fmpz_set(t[41], t[9]);                  /*  set d9 = c3     */
    fmpz_set(t[42], t[10]);                 /*  set d10 = c4    */

    fmpz_sub(t[0], t[20], t[23]);           /*  a0 = x0 - x3    */
    fmpz_sub(t[1], t[21], t[24]);           /*  a1 = x1 - x4    */
    fmpz_sub(t[2], t[22], t[25]);           /*  a2 = x2 - x5    */
    fmpz_sub(t[3], t[29], t[26]);           /*  b0 = y3 - y0    */
    fmpz_sub(t[4], t[30], t[27]);           /*  b1 = y4 - y1    */
    fmpz_sub(t[5], t[31], t[28]);           /*  b2 = y5 - y2    */

    /*
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[43], t[6]);                  /*  set d11 = c0    */
    fmpz_set(t[44], t[7]);                  /*  set d12 = c1    */
    fmpz_set(t[45], t[8]);                  /*  set d13 = c2    */
    fmpz_set(t[46], t[9]);                  /*  set d14 = c3    */
    fmpz_set(t[47], t[10]);                 /*  set d15 = c4    */

    fmpz_add(t[50], t[38], t[46]);          /*  d18 = d6 + d14  */
    fmpz_add(t[48], t[50], t[35]);          /*  d16 = d18 + d3  */
    fmpz_add(t[50], t[39], t[47]);          /*  d18 = d7 + d15  */
    fmpz_add(t[49], t[50], t[36]);          /*  d17 = d18 + d4  */
    fmpz_add(t[50], t[35], t[43]);          /*  d18 = d3 + d11  */
    fmpz_add(t[35], t[50], t[32]);          /*  d3 = d18 + d0   */
    fmpz_add(t[50], t[36], t[44]);          /*  d18 = d4 + d12  */

    fmpz_add(t[36], t[50], t[33]);          /*  d4 = d18 + d1   */
    fmpz_add(t[37], t[34], t[45]);          /*  d5 = d2 + d13   */
    fmpz_sub(t[0], t[32], t[48]);           /*  a0 = d0 - d16   */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = a0 mod n   */

    fmpz_sub(t[0], t[33], t[49]);           /*  a0 = d1 - d17   */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = a0 mod n   */

    fmpz_sub(t[0], t[34], t[40]);           /*  a0 = d2 - d8    */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /*  z2 = a0 mod n   */

    unity_zp_coeff_set_fmpz(f, 5, t[37]);   /*  z5 = d5 mod n   */

    fmpz_add(t[50], t[35], t[38]);          /*  d18 = d3 + d6   */
    fmpz_add(t[51], t[48], t[41]);          /*  d19 = d16 + d9  */
    fmpz_sub(t[0], t[50], t[51]);           /*  a0 = d18 - d19  */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /*  z3 = a0 mod n   */

    fmpz_add(t[50], t[36], t[39]);          /*  d18 = d4 + d7   */
    fmpz_add(t[51], t[42], t[49]);          /*  d19 = d10 + d17 */
    fmpz_sub(t[0], t[50], t[51]);           /*  a0 = d18 - d19  */
    unity_zp_coeff_set_fmpz(f, 4, t[0]);    /*  z4 = a0 mod n   */
}

/*
    Computes f = g * h for p = 5.
    g and h must be reduced by F_5 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 24.
    Resulting f reduced by F_5 cyclotomic polynomial.
*/
void
unity_zp_mul5(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3);
        h = (y0, y1, y2, y3);
        f = (z0, z1, z2, z3);

        x0 = t[0]; x1 = t[1]; x2 = t[2]; x3 = t[3];
        y0 = t[4]; y1 = t[5]; y2 = t[6]; y3 = t[7];

        m1 = t[8]; m2 = t[9]; m3 = t[10]; m4 = t[11];
        m5 = t[12]; m6 = t[13]; m7 = t[14]; m8 = t[15];

        d0 = t[16]; d1 = t[17]; d2 = t[18]; d3 = t[19];
        d4 = t[20]; d5 = t[21]; d6 = t[22]; d7 = t[23];
        d8 = t[24].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[4], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[5], h->poly, 1, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[6], h->poly, 2, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[7], h->poly, 3, h->ctx);

    fmpz_sub(t[8], t[1], t[3]);             /* m1 = x1 - x3  */
    fmpz_sub(t[9], t[5], t[7]);             /* m2 = y1 - y3  */
    fmpz_sub(t[10], t[2], t[3]);            /* m3 = x2 - x3  */
    fmpz_sub(t[11], t[7], t[6]);            /* m4 = y3 - y2  */
    fmpz_sub(t[12], t[0], t[1]);            /* m5 = x0 - x1  */
    fmpz_sub(t[13], t[5], t[4]);            /* m6 = y1 - y0  */

    fmpz_sub(t[14], t[0], t[2]);            /* m7 = x0 - x2  */
    fmpz_sub(t[15], t[6], t[4]);            /* m8 = y2 - y0  */
    fmpz_mul(t[16], t[0], t[4]);            /* d0 = x0 * y0  */
    fmpz_mul(t[18], t[8], t[9]);            /* d2 = m1 * m2  */
    fmpz_add(t[17], t[16], t[18]);          /* d1 = d0 + d2  */
    fmpz_mul(t[18], t[10], t[11]);          /* d2 = m3 * m4  */

    fmpz_mul(t[19], t[12], t[13]);          /* d3 = m5 * m6  */
    fmpz_mul(t[20], t[14], t[15]);          /* d4 = m7 * m8  */
    fmpz_mul(t[21], t[1], t[5]);            /* d5 = x1 * y1  */
    fmpz_mul(t[22], t[2], t[6]);            /* d6 = x2 * y2  */
    fmpz_mul(t[23], t[3], t[7]);            /* d7 = x3 * y3  */
    fmpz_add(t[24], t[17], t[18]);          /* d8 = d1 + d2  */

    fmpz_sub(t[0], t[24], t[21]);           /* x0 = d8 - d5  */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /* z0 = x0 mod n */

    fmpz_add(t[24], t[17], t[19]);          /* d8 = d1 + d3  */
    fmpz_sub(t[0], t[24], t[22]);           /* x0 = d8 - d6  */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /* z1 = x0 mod n */

    fmpz_add(t[24], t[17], t[20]);          /* d8 = d1 + d4  */
    fmpz_sub(t[0], t[24], t[23]);           /* x0 = d8 - d7  */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /* z2 = x0 mod n */

    fmpz_sub(t[10], t[8], t[14]);           /* m3 = m1 - m7  */
    fmpz_add(t[11], t[9], t[15]);           /* m4 = m2 + m8  */
    fmpz_mul(t[17], t[10], t[11]);          /* d1 = m3 * m4  */
    fmpz_add(t[24], t[16], t[17]);          /* d8 = d0 + d1  */
    fmpz_add(t[23], t[24], t[18]);          /* d7 = d8 + d2  */
    fmpz_add(t[24], t[23], t[19]);          /* d8 = d7 + d3  */
    fmpz_add(t[0], t[24], t[20]);           /* d7 = d8 + d4  */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /* z3 = d7 mod n */
}

/*
    Computes f = g * h for p = 7.
    g and h must be reduced by F_7 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 68.
    Resulting f reduced by F_7 cyclotomic polynomial.
*/
void
unity_zp_mul7(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, ... , x5);
        h = (y0, ... , y5);
        f = (z0, ... , z5);

        x0 = t[30]; x1 = t[31]; x2 = t[32]; x3 = t[33];
        x4 = t[34]; x5 = t[35];

        y0 = t[40]; y1 = t[41]; y2 = t[42]; y3 = t[43];
        y4 = t[44]; y5 = t[45];

        a0 = t[0]; a1 = t[1] ; a2 = t[2];
        b0 = t[3]; b1 = t[4] ; b2 = t[5];
        c0 = t[6]; ... ; c4 = t[10];

        d0 = t[50]; ... ; d18 = t[68];
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[30], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[31], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[32], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[33], g->poly, 3, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[34], g->poly, 4, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[35], g->poly, 5, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[40], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[41], h->poly, 1, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[42], h->poly, 2, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[43], h->poly, 3, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[44], h->poly, 4, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[45], h->poly, 5, h->ctx);

    fmpz_set(t[0], t[30]);                  /*  set a0 = x0     */
    fmpz_set(t[1], t[31]);                  /*  set a1 = x1     */
    fmpz_set(t[2], t[32]);                  /*  set a2 = x2     */
    fmpz_set(t[3], t[40]);                  /*  set b0 = y0     */
    fmpz_set(t[4], t[41]);                  /*  set b1 = y1     */
    fmpz_set(t[5], t[42]);                  /*  set b2 = y2     */

    /*
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[50], t[6]);                  /*  set d0 = c0     */
    fmpz_set(t[51], t[7]);                  /*  set d1 = c1     */
    fmpz_set(t[52], t[8]);                  /*  set d2 = c2     */
    fmpz_set(t[53], t[9]);                  /*  set d3 = c3     */
    fmpz_set(t[54], t[10]);                 /*  set d4 = c4     */

    fmpz_set(t[0], t[33]);                  /*  set a0 = x3     */
    fmpz_set(t[1], t[34]);                  /*  set a1 = x4     */
    fmpz_set(t[2], t[35]);                  /*  set a2 = x5     */
    fmpz_set(t[3], t[43]);                  /*  set b0 = y3     */
    fmpz_set(t[4], t[44]);                  /*  set b1 = y4     */
    fmpz_set(t[5], t[45]);                  /*  set b2 = y5     */

    /*
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[56], t[6]);                  /*  set d6 = c0     */
    fmpz_set(t[57], t[7]);                  /*  set d7 = c1     */
    fmpz_set(t[58], t[8]);                  /*  set d8 = c2     */
    fmpz_set(t[59], t[9]);                  /*  set d9 = c3     */
    fmpz_set(t[60], t[10]);                 /*  set d10 = c4    */

    fmpz_sub(t[0], t[30], t[33]);           /*  a0 = x0 - x3    */
    fmpz_sub(t[1], t[31], t[34]);           /*  a1 = x1 - x4    */
    fmpz_sub(t[2], t[32], t[35]);           /*  a2 = x2 - x5    */
    fmpz_sub(t[3], t[43], t[40]);           /*  b0 = y3 - y0    */
    fmpz_sub(t[4], t[44], t[41]);           /*  b1 = y4 - y1    */
    fmpz_sub(t[5], t[45], t[42]);           /*  b2 = y5 - y2    */

    /*
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[61], t[6]);                  /*  set d11 = c0    */
    fmpz_set(t[62], t[7]);                  /*  set d12 = c1    */
    fmpz_set(t[63], t[8]);                  /*  set d13 = c2    */
    fmpz_set(t[64], t[9]);                  /*  set d14 = c3    */
    fmpz_set(t[65], t[10]);                 /*  set d15 = c4    */

    fmpz_add(t[68], t[56], t[64]);          /*  d18 = d6 + d14  */
    fmpz_add(t[66], t[68], t[53]);          /*  d16 = d18 + d3  */
    fmpz_add(t[68], t[57], t[65]);          /*  d18 = d7 + d15  */
    fmpz_add(t[67], t[68], t[54]);          /*  d17 = d18 + d4  */
    fmpz_add(t[68], t[53], t[61]);          /*  d18 = d3 + d11  */
    fmpz_add(t[53], t[68], t[50]);          /*  d3 = d18 + d0   */
    fmpz_add(t[68], t[54], t[62]);          /*  d18 = d4 + d12  */
    fmpz_add(t[54], t[68], t[51]);          /*  d4 = d18 + d1   */

    fmpz_add(t[55], t[52], t[63]);          /*  d5 = d2 + d13   */
    fmpz_add(t[63], t[53], t[56]);          /*  d13 = d3 + d6   */
    fmpz_add(t[64], t[54], t[57]);          /*  d14 = d4 + d7   */
    fmpz_add(t[65], t[55], t[58]);          /*  d15 = d5 + d8   */
    fmpz_add(t[56], t[66], t[59]);          /*  d6 = d16 + d9   */
    fmpz_add(t[57], t[67], t[60]);          /*  d7 = d17 + d10  */

    fmpz_add(t[68], t[50], t[57]);          /*  d18 = d10 + d7  */
    fmpz_sub(t[0], t[68], t[56]);           /*  a0 = d18 - d6   */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = a0 mod n   */

    fmpz_add(t[68], t[51], t[58]);          /*  d18 = d1 + d8   */
    fmpz_sub(t[0], t[68], t[56]);           /*  a0 = d18 - d6   */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = a0 mod n   */

    fmpz_add(t[68], t[52], t[59]);          /*  d18 = d2 + d9   */
    fmpz_sub(t[0], t[68], t[56]);           /*  a0 = d18 - d6   */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /*  z2 = a0 mod n   */

    fmpz_add(t[68], t[63], t[60]);          /*  d18 = d13 + d10 */
    fmpz_sub(t[0], t[68], t[56]);           /*  a0 = d18 - d6   */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /*  z3 = a0 mod n   */

    fmpz_sub(t[0], t[64], t[56]);           /*  a0 = d14 - d6   */
    unity_zp_coeff_set_fmpz(f, 4, t[0]);    /*  z4 = a0 mod n   */

    fmpz_sub(t[0], t[65], t[56]);           /*  a0 = d15 - d6   */
    unity_zp_coeff_set_fmpz(f, 5, t[0]);    /*  z5 = a0 mod n   */
}

void
unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s)
{
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));

    fmpz_mod_poly_scalar_mul_ui(f->poly, g->poly, s, f->ctx);
}
