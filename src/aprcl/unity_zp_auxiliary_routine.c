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
    Input a0, a1, a2 store in t[0 .. 2] and b0, b1, b2 store in t[3 .. 5].
    Output c0, ..., c4 store in t[6 .. 10].

    Used t[0 .. 17].

    c0 = a0 * b0;
    c1 = a0 * b1 + b0 * a1;
    c2 = a0 * b2 + a1 * b1 + a2 * b0;
    c3 = a1 * b2 + a2 * b1;
    c4 = a2 * b2.
*/
void
unity_zp_ar1(fmpz_t * t)
{
    /*
        a0 = t[0]; a1 = t[1]; a2 = t[2];
        b0 = t[3]; b1 = t[4]; b2 = t[5];
        
        c0 = t[6]; c1 = t[7]; c2 = t[8];
        c3 = t[9]; c4 = t[10];

        m1 = t[11]; m2 = t[12];

        d1 = t[13]; d2 = t[14]; d3 = t[15]; d4 = t[16]; d5 = t[17].
    */

    fmpz_mul(t[6], t[0], t[3]);         /* c0 = a0 * b0 */
    fmpz_mul(t[13], t[1], t[4]);        /* d1 = a1 * b1 */
    fmpz_mul(t[10], t[2], t[5]);        /* c4 = a2 * b2 */
    fmpz_add(t[11], t[0], t[1]);        /* m1 = a0 + a1 */
    fmpz_add(t[12], t[3], t[4]);        /* m2 = b0 + b1 */
    fmpz_mul(t[15], t[11], t[12]);      /* d3 = m1 * m2 */
    
    fmpz_add(t[11], t[0], t[2]);        /* m1 = a0 + a2 */
    fmpz_add(t[12], t[3], t[5]);        /* m2 = b0 + b2 */
    fmpz_mul(t[16], t[11], t[12]);      /* d4 = m1 * m2 */
    fmpz_add(t[11], t[1], t[2]);        /* m1 = a1 + a2 */
    fmpz_add(t[12], t[4], t[5]);        /* m2 = b1 + b2 */
    fmpz_mul(t[17], t[11], t[12]);      /* d5 = m1 * m2 */

    fmpz_add(t[14], t[6], t[13]);       /* d2 = c0 + d1 */
    fmpz_sub(t[7], t[15], t[14]);       /* c1 = d3 - d2 */
    fmpz_add(t[14], t[16], t[13]);      /* d2 = d4 + d1 */
    fmpz_add(t[16], t[6], t[10]);       /* d4 = c0 + c4 */
    fmpz_sub(t[8], t[14], t[16]);       /* c2 = d2 - d4 */
    fmpz_add(t[14], t[13], t[10]);      /* d2 = d1 + c4 */

    fmpz_sub(t[9], t[17], t[14]);       /* c3 = d5 - d2 */
}

/*
    Input a0, ..., a3 store in t[0 .. 3] and b0, ..., b3 store in t[4 .. 7].
    Output c0, ..., c6 store in t[8 .. 14].

    Used t[0 .. 27].

    c0 = a0 * b0;
    c1 = a0 * b1 + b0 * a1;
    c2 = a0 * b2 + a1 * b1 + a2 * b0;
    c3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0;
    c4 = a1 * b3 + a2 * b2 + a3 * b1;
    c5 = a2 * b3 + a3 * b2;
    c6 = a3 * b3.
*/
void
unity_zp_ar2(fmpz_t * t)
{
    /*
        a0 = t[0]; a1 = t[1]; a2 = t[2]; a3 = t[3];
        b0 = t[4]; b1 = t[5]; b2 = t[6]; b3 = t[7];

        c0 = t[8]; c1 = t[9]; c2 = t[10]; c3 = t[11];
        c4 = t[12]; c5 = t[13]; c6 = t[14];

        m1 = t[15]; m2 = t[16]; m3 = t[17]; m4 = t[18];
        m5 = t[19];

        d1 = t[20]; d2 = t[21]; d3 = t[22]; d4 = t[23];
        d5 = t[24]; d6 = t[25]; d7 = t[26]; d8 = t[27];
        d9 = t[28].
    */

    fmpz_mul(t[8], t[0], t[4]);         /* c0 = a0 * b0 */
    fmpz_mul(t[20], t[1], t[5]);        /* d1 = a1 * b1 */
    fmpz_mul(t[21], t[2], t[6]);        /* d2 = a2 * b2 */
    fmpz_mul(t[14], t[3], t[7]);        /* c6 = a3 * b3 */
    fmpz_add(t[15], t[0], t[1]);        /* m1 = a0 + a1 */
    fmpz_add(t[16], t[4], t[5]);        /* m2 = b0 + b1 */
    fmpz_mul(t[22], t[15], t[16]);      /* d3 = m1 * m2 */

    fmpz_add(t[15], t[0], t[2]);        /* m1 = a0 + a2 */
    fmpz_add(t[16], t[4], t[6]);        /* m2 = b0 + b2 */
    fmpz_mul(t[23], t[15], t[16]);      /* d4 = m1 * m2 */
    fmpz_add(t[17], t[2], t[3]);        /* m3 = a2 + a3 */
    fmpz_add(t[18], t[6], t[7]);        /* m4 = b2 + b3 */
    fmpz_mul(t[24], t[17], t[18]);      /* d5 = m3 * m4 */

    fmpz_add(t[17], t[1], t[3]);        /* m3 = a1 + a3 */
    fmpz_add(t[18], t[5], t[7]);        /* m4 = b1 + b3 */
    fmpz_mul(t[25], t[17], t[18]);      /* d6 = m3 * m4 */
    fmpz_add(t[26], t[8], t[20]);       /* d7 = c0 + d1 */
    fmpz_sub(t[9], t[22], t[26]);       /* c1 = d3 - d7 */
    fmpz_add(t[26], t[8], t[21]);       /* d7 = c0 + d2 */

    fmpz_add(t[27], t[20], t[23]);      /* d8 = d1 + d4 */
    fmpz_sub(t[10], t[27], t[26]);      /* c2 = d8 - d7 */
    fmpz_add(t[19], t[15], t[17]);      /* m5 = m1 + m3 */
    fmpz_add(t[17], t[16], t[18]);      /* m3 = m2 + m4 */
    fmpz_add(t[26], t[21], t[14]);      /* d7 = d2 + c6 */
    fmpz_sub(t[13], t[24], t[26]);      /* c5 = d5 - d7 */

    fmpz_mul(t[26], t[17], t[19]);      /* d7 = m3 * m5 */
    fmpz_add(t[27], t[9], t[13]);       /* d8 = c1 + c5 */
    fmpz_add(t[28], t[27], t[25]);      /* d9 = d8 + d6 */
    fmpz_add(t[27], t[28], t[23]);      /* d8 = d9 + d4 */
    fmpz_sub(t[11], t[26], t[27]);      /* c3 = d7 - d8 */
    fmpz_add(t[26], t[25], t[21]);      /* d7 = d6 + d2 */

    fmpz_add(t[27], t[20], t[14]);      /* d8 = d1 + c6 */
    fmpz_sub(t[12], t[26], t[27]);      /* c4 = d7 - d8 */
}

/*
    Input a0, ... , a4 store in t[0 .. 4] and b0, ..., b4 in t[5 .. 9].
    
    Output c0, ..., c8 store in t[10 .. 18].

    Used t[0 .. 38].

    c0 = a0 * b0;
    c1 = a0 * b1 + a1 * b0;
    c2 = a0 * b2 + a1 * b1 + a2 * b0;
    c3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0;
    c4 = a0 * b4 + a1 * b3 + a2 * b2 + a3 * b1 + a4 * b0;
    c5 = a1 * b4 + a2 * b3 + a3 * b2 + a4 * b1;
    c6 = a2 * b4 + a3 * b3 + a4 * b2;
    c7 = a3 * b4 + a4 * b3;
    c8 = a4 * b4.
*/
void
unity_zp_ar3(fmpz_t * t)
{
    /*
        a0 = t[0] = t[20]; a1 = t[1] = t[21]; a2 = t[2] = t[22];
        a3 = t[3] = t[23]; a4 = t[4] = t[24];

        b0 = t[5] = t[25]; b1 = t[6] = t[26]; b2 = t[7] = t[27];
        b3 = t[8] = t[28]; b4 = t[9] = t[29];

        c0 = t[30] = t[10]; c1 = t[31] = t[11]; c2 = t[31] = t[12];
        c3 = t[33] = t[13]; c4 = t[34] = t[14]; c5 = t[35] = t[15];
        c6 = t[36] = t[16]; c7 = t[37] = t[17]; c8 = t[38] = t[18];

        m0 = t[0]; m1 = t[1]; m2 = t[2]; m3 = t[3]; m4 = t[15]; m5 = t[16];

        d0 = t[6]; d1 = t[7]; d2 = t[8]; d3 = t[9]; d4 = t[10];
        d5 = t[11]; d6 = t[12]; d7 = t[13].
    */

    fmpz_set(t[20], t[0]);              /* store a0 */
    fmpz_set(t[21], t[1]);              /* store a1 */
    fmpz_set(t[22], t[2]);              /* store a2 */
    fmpz_set(t[23], t[3]);              /* store a3 */
    fmpz_set(t[24], t[4]);              /* store a4 */
    fmpz_set(t[25], t[5]);              /* store b0 */
    fmpz_set(t[26], t[6]);              /* store b1 */
    fmpz_set(t[27], t[7]);              /* store b2 */
    fmpz_set(t[28], t[8]);              /* store b3 */
    fmpz_set(t[29], t[9]);              /* store b4 */

    fmpz_set(t[3], t[25]);              /* set t[3] = b0 */
    fmpz_set(t[4], t[26]);              /* set t[4] = b1 */
    fmpz_set(t[5], t[27]);              /* set t[5] = b2 */

    /* apply ar1 to a0, a1, a2 and b0, b1, b2 */
    unity_zp_ar1(t);

    fmpz_set(t[30], t[6]);              /* store c0 */
    fmpz_set(t[31], t[7]);              /* store c1 */
    fmpz_set(t[32], t[8]);              /* store c2 */
    fmpz_set(t[33], t[9]);              /* store c3 */
    fmpz_set(t[34], t[10]);             /* store c4 */

    fmpz_add(t[0], t[20], t[23]);       /* m0 = a0 + a3 */
    fmpz_add(t[1], t[21], t[24]);       /* m1 = a1 + a4 */
    fmpz_add(t[3], t[25], t[28]);       /* m2 = b0 + b3 */
    fmpz_add(t[4], t[26], t[29]);       /* m3 = b1 + b4 */

    /* apply ar1 to m0, m1, a2 and m2, m3, b2 */
    unity_zp_ar1(t);

    fmpz_add(t[15], t[23], t[24]);      /* m4 = a3 + a4 */
    fmpz_add(t[16], t[28], t[29]);      /* m5 = b3 + b4 */
    
    fmpz_mul(t[11], t[23], t[28]);      /* d5 = a3 * b3 */
    fmpz_mul(t[38], t[24], t[29]);      /* c8 = a4 * b4 */
    fmpz_mul(t[12], t[15], t[16]);      /* d6 = m4 * m5 */
    fmpz_add(t[13], t[11], t[38]);      /* d7 = d5 + c8 */
    fmpz_sub(t[37], t[12], t[13]);      /* c7 = d6 - d7 */
    fmpz_add(t[12], t[9], t[11]);       /* d6 = d3 + d5 */

    fmpz_sub(t[36], t[12], t[33]);      /* c6 = d6 - c3 */
    fmpz_add(t[12], t[30], t[11]);      /* d6 = c0 + d5 */
    fmpz_add(t[13], t[33], t[6]);       /* d7 = c3 + d0 */
    fmpz_sub(t[33], t[13], t[12]);      /* c3 = d7 - d6 */
    fmpz_add(t[12], t[31], t[37]);      /* d6 = c1 + c7 */
    fmpz_add(t[13], t[34], t[7]);       /* d7 = c4 + d1 */

    fmpz_sub(t[34], t[13], t[12]);      /* c4 = d7 - d6 */
    fmpz_add(t[12], t[32], t[38]);      /* d6 = c2 + c8 */
    fmpz_sub(t[35], t[8], t[12]);       /* c5 = d2 - d6 */

    fmpz_set(t[10], t[30]);             /* store c0 */
    fmpz_set(t[11], t[31]);             /* store c1 */
    fmpz_set(t[12], t[32]);             /* store c2 */
    fmpz_set(t[13], t[33]);             /* store c3 */
    fmpz_set(t[14], t[34]);             /* store c4 */
    fmpz_set(t[15], t[35]);             /* store c5 */
    fmpz_set(t[16], t[36]);             /* store c6 */
    fmpz_set(t[17], t[37]);             /* store c7 */
    fmpz_set(t[18], t[38]);             /* store c8 */
}

/*
    Input a0, ... , a4 store in t[0 .. 4].
    Output c0, ..., c8 store in t[5 .. 13].

    Used t[0 .. 29].

    c0 = a0 * a0;
    c1 = 2 * a0 * a1;
    c2 = 2 * a0 * a2 + a1 * a1;
    c3 = 2 * a0 * a3 + a * a1 * a2;
    c4 = 2 * a0 * a4 + 2 * a1 * a3 + a2 * a2;
    c5 = 2 * a1 * a4 + 2 * a2 * a3;
    c6 = 2 * a2 * a4 + a3 * a3;
    c7 = 2 * a3 * a4;
    c8 = a4 * a4.
*/
void
unity_zp_ar4(fmpz_t * t)
{
    /*
        a0 = t[0]; a1 = t[1]; a2 = t[2]; a3 = t[3]; a4 = t[4];

        c0 = t[5]; c1 = t[6]; c2 = t[7]; c3 = t[8]; c4 = t[9];
        c5 = t[10]; c6 = t[11]; c7 = t[12]; c8 = t[13];

        m1 = t[14]; m2 = t[15]; m3 = t[16]; m4 = t[17];
        m5 = t[18]; m6 = t[19]; m7 = t[20]; m8 = t[21];
        m9 = t[22]; m10 = t[23];

        d1 = t[24]; d2 = t[25]; d3 = t[26]; d4 = t[27];
        d5 = t[28]; d6 = t[29].
    */

    fmpz_add(t[14], t[2], t[2]);        /* m1 = a2 + a2  */
    fmpz_add(t[15], t[0], t[1]);        /* m2 = a0 + a1  */
    fmpz_add(t[16], t[1], t[14]);       /* m3 = a1 + m1  */
    fmpz_add(t[17], t[3], t[4]);        /* m4 = a3 + a4  */
    fmpz_add(t[18], t[3], t[14]);       /* m5 = a3 + m1  */
    fmpz_add(t[19], t[0], t[0]);        /* m6 = a0 + a0  */

    fmpz_add(t[19], t[19], t[14]);      /* m6 = m6 + m1  */
    fmpz_add(t[20], t[1], t[3]);        /* m7 = a1 + a3  */
    fmpz_add(t[21], t[4], t[4]);        /* m8 = a4 + a4  */
    fmpz_add(t[21], t[21], t[14]);      /* m8 = m8 + m1  */
    fmpz_add(t[22], t[0], t[3]);        /* m9 = a0 + a3  */
    fmpz_add(t[23], t[1], t[4]);        /* m10 = a1 + a4 */

    fmpz_mul(t[5], t[0], t[0]);         /* c0 = a0 * a0  */
    fmpz_mul(t[24], t[0], t[1]);        /* d1 = a0 * a1  */
    fmpz_mul(t[13], t[4], t[4]);        /* c8 = a4 * a4  */
    fmpz_mul(t[25], t[3], t[4]);        /* d2 = a3 * a4  */
    fmpz_mul(t[26], t[1], t[14]);       /* d3 = a1 * m1  */
    fmpz_mul(t[27], t[3], t[14]);       /* d4 = a3 * m1  */
    fmpz_add(t[6], t[24], t[24]);       /* c1 = d1 + d1  */

    fmpz_add(t[12], t[25], t[25]);      /* c7 = d2 + d2  */
    fmpz_mul(t[28], t[15], t[16]);      /* d5 = m2 * m3  */
    fmpz_add(t[29], t[24], t[26]);      /* d6 = d1 + d3  */
    fmpz_sub(t[7], t[28], t[29]);       /* c2 = d5 - d6  */
    fmpz_mul(t[28], t[17], t[18]);      /* d5 = m4 * m5  */
    fmpz_add(t[29], t[25], t[27]);      /* d6 = d2 + d4  */

    fmpz_sub(t[11], t[28], t[29]);      /* c6 = d5 - d6  */
    fmpz_mul(t[28], t[19], t[20]);      /* d5 = m6 * m7  */
    fmpz_add(t[29], t[6], t[27]);       /* d6 = c1 + d4  */
    fmpz_sub(t[8], t[28], t[29]);       /* c3 = d5 - d6  */
    fmpz_mul(t[28], t[20], t[21]);      /* d5 = m7 * m8  */
    fmpz_add(t[29], t[12], t[26]);      /* d6 = c7 + d3  */

    fmpz_sub(t[10], t[28], t[29]);      /* c5 = d5 - d6  */
    fmpz_mul(t[28], t[22], t[23]);      /* d5 = m9 * m10 */
    fmpz_add(t[29], t[24], t[25]);      /* d6 = d1 + d2  */
    fmpz_sub(t[28], t[28], t[29]);      /* d6 = d5 - d6  */
    fmpz_add(t[29], t[28], t[28]);      /* d6 = d5 + d5  */
    fmpz_mul(t[28], t[2], t[2]);        /* d5 = a2 * a2  */
    fmpz_add(t[9], t[28], t[29]);       /* c4 = d5 + d6  */
}

