/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Vladimir Glazachev
   
******************************************************************************/

#include "aprcl.h"

void
unity_zp_ar1(fmpz_t * t)
{
    /*
        a0 = t[0]; a1 = t[1]; a2 = t[2];
        b0 = t[3]; b1 = t[4]; b2 = t[5];
        
        c0 = t[6]; c1 = t[7]; c2 = t[8];
        c3 = t[9]; c4 = t[10];
    */

    fmpz_mul(t[6], t[0], t[3]);
    fmpz_mul(t[13], t[1], t[4]);
    fmpz_mul(t[10], t[2], t[5]);
    fmpz_add(t[11], t[0], t[1]);
    fmpz_add(t[12], t[3], t[4]);
    fmpz_mul(t[15], t[11], t[12]);
    
    fmpz_add(t[11], t[0], t[2]);
    fmpz_add(t[12], t[3], t[5]);
    fmpz_mul(t[16], t[11], t[12]);
    fmpz_add(t[11], t[1], t[2]);
    fmpz_add(t[12], t[4], t[5]);
    fmpz_mul(t[17], t[11], t[12]);

    fmpz_add(t[14], t[6], t[13]);
    fmpz_sub(t[7], t[15], t[14]);
    fmpz_add(t[14], t[16], t[13]);
    fmpz_add(t[16], t[6], t[10]);
    fmpz_sub(t[8], t[14], t[16]);
    fmpz_add(t[14], t[13], t[10]);

    fmpz_sub(t[9], t[17], t[14]);
}

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

    fmpz_mul(t[8], t[0], t[4]);
    fmpz_mul(t[20], t[1], t[5]);
    fmpz_mul(t[21], t[2], t[6]);
    fmpz_mul(t[14], t[3], t[7]);
    fmpz_add(t[15], t[0], t[1]);
    fmpz_add(t[16], t[4], t[5]);
    fmpz_mul(t[22], t[15], t[16]);

    fmpz_add(t[15], t[0], t[2]);
    fmpz_add(t[16], t[4], t[6]);
    fmpz_mul(t[23], t[15], t[16]);
    fmpz_add(t[17], t[2], t[3]);
    fmpz_add(t[18], t[6], t[7]);
    fmpz_mul(t[24], t[17], t[18]);

    fmpz_add(t[17], t[1], t[3]);
    fmpz_add(t[18], t[5], t[7]);
    fmpz_mul(t[25], t[17], t[18]);
    fmpz_add(t[26], t[8], t[20]);
    fmpz_sub(t[9], t[22], t[26]);
    fmpz_add(t[26], t[8], t[21]);

    fmpz_add(t[27], t[20], t[23]);
    fmpz_sub(t[10], t[27], t[26]);
    fmpz_add(t[19], t[15], t[17]);
    fmpz_add(t[17], t[16], t[18]);
    fmpz_add(t[26], t[21], t[14]);
    fmpz_sub(t[13], t[24], t[26]);

    fmpz_mul(t[26], t[17], t[19]);
    fmpz_add(t[27], t[9], t[13]);
    fmpz_add(t[28], t[27], t[25]);
    fmpz_add(t[27], t[28], t[23]);
    fmpz_sub(t[11], t[26], t[27]);
    fmpz_add(t[26], t[25], t[21]);

    fmpz_add(t[27], t[20], t[14]);
    fmpz_sub(t[12], t[26], t[27]);
}

void
unity_zp_ar3(fmpz_t * t)
{
    fmpz_set(t[30], t[0]);              /* store a0 */
    fmpz_set(t[31], t[1]);              /* store a1 */
    fmpz_set(t[32], t[2]);              /* store a2 */
    fmpz_set(t[33], t[3]);              /* store a3 */
    fmpz_set(t[34], t[4]);              /* store a4 */
    fmpz_set(t[35], t[5]);              /* store b0 */
    fmpz_set(t[36], t[6]);              /* store b1 */
    fmpz_set(t[37], t[7]);              /* store b2 */
    fmpz_set(t[38], t[8]);              /* store b3 */
    fmpz_set(t[39], t[9]);              /* store b4 */

    fmpz_set(t[3], t[35]);              /* set t[3] = b0 */
    fmpz_set(t[4], t[36]);              /* set t[4] = b1 */
    fmpz_set(t[5], t[37]);              /* set t[5] = b2 */

    /* apply ar1 to t[30], t[31], a2 and t[35], t[36], b2 */
    unity_zp_ar1(t);

    fmpz_set(t[40], t[6]);              /* store c0 */
    fmpz_set(t[41], t[7]);              /* store c1 */
    fmpz_set(t[42], t[8]);              /* store c2 */
    fmpz_set(t[43], t[9]);              /* store c3 */
    fmpz_set(t[44], t[10]);             /* store c4 */

    fmpz_add(t[0], t[30], t[33]);       /* m0 = a0 + a3 */
    fmpz_add(t[1], t[31], t[34]);       /* m1 = a1 + a4 */
    fmpz_add(t[3], t[35], t[38]);       /* m2 = b0 + b3 */
    fmpz_add(t[4], t[36], t[39]);       /* m3 = b1 + b4 */

    /* apply ar1 to m0, m1, a2 and m2, m3, b2 */
    unity_zp_ar1(t);

    fmpz_add(t[15], t[33], t[34]);      /* m4 = a3 + a4 */
    fmpz_add(t[16], t[38], t[39]);      /* m5 = b3 + b4 */
    
    fmpz_mul(t[11], t[33], t[38]);      /* d5 = a3 * b3 */
    fmpz_mul(t[48], t[34], t[39]);      /* c8 = a4 * b4 */
    fmpz_mul(t[12], t[15], t[16]);      /* d6 = m4 * m5 */
    fmpz_add(t[13], t[11], t[48]);      /* d7 = d5 + c8 */
    fmpz_sub(t[47], t[12], t[13]);      /* c7 = d6 - d7 */
    fmpz_add(t[12], t[9], t[11]);       /* d6 = d3 + d5 */

    fmpz_sub(t[46], t[12], t[43]);      /* c6 = d6 - c3 */
    fmpz_add(t[12], t[40], t[11]);      /* d6 = c0 + d5 */
    fmpz_add(t[13], t[43], t[6]);       /* d7 = c3 + d0 */
    fmpz_sub(t[43], t[13], t[12]);      /* c3 = d7 - d6 */
    fmpz_add(t[12], t[41], t[47]);      /* d6 = c1 + c7 */
    fmpz_add(t[13], t[44], t[7]);       /* d7 = c4 + d1 */

    fmpz_sub(t[44], t[13], t[12]);      /* c4 = d7 - d6 */
    fmpz_add(t[12], t[42], t[48]);      /* d6 = c2 + c8 */
    fmpz_sub(t[45], t[8], t[12]);       /* c5 = d2 - d6 */

    fmpz_set(t[10], t[40]);             /* store c0 */
    fmpz_set(t[11], t[41]);             /* store c1 */
    fmpz_set(t[12], t[42]);             /* store c2 */
    fmpz_set(t[13], t[43]);             /* store c3 */
    fmpz_set(t[14], t[44]);             /* store c4 */
    fmpz_set(t[15], t[45]);             /* store c5 */
    fmpz_set(t[16], t[46]);             /* store c6 */
    fmpz_set(t[17], t[47]);             /* store c7 */
    fmpz_set(t[18], t[48]);             /* store c8 */
}

void
unity_zp_ar4(fmpz_t * t)
{
    fmpz_add(t[14], t[2], t[2]);
    fmpz_add(t[15], t[0], t[1]);
    fmpz_add(t[16], t[1], t[14]);
    fmpz_add(t[17], t[3], t[4]);
    fmpz_add(t[18], t[3], t[14]);
    fmpz_add(t[19], t[0], t[0]);

    fmpz_add(t[19], t[19], t[14]);
    fmpz_add(t[20], t[1], t[3]);
    fmpz_add(t[21], t[4], t[4]);
    fmpz_add(t[21], t[21], t[14]);
    fmpz_add(t[22], t[0], t[3]);
    fmpz_add(t[23], t[1], t[4]);

    fmpz_mul(t[5], t[0], t[0]);
    fmpz_mul(t[24], t[0], t[1]);
    fmpz_mul(t[13], t[4], t[4]);
    fmpz_mul(t[25], t[3], t[4]);
    fmpz_mul(t[26], t[1], t[14]);
    fmpz_mul(t[27], t[3], t[14]);
    fmpz_add(t[6], t[24], t[24]);

    fmpz_add(t[12], t[25], t[25]);
    fmpz_mul(t[28], t[15], t[16]);
    fmpz_add(t[29], t[24], t[26]);
    fmpz_sub(t[7], t[28], t[29]);
    fmpz_mul(t[28], t[17], t[18]);
    fmpz_add(t[29], t[25], t[27]);

    fmpz_sub(t[11], t[28], t[29]);
    fmpz_mul(t[28], t[19], t[20]);
    fmpz_add(t[29], t[6], t[27]);
    fmpz_sub(t[8], t[28], t[29]);
    fmpz_mul(t[28], t[20], t[21]);
    fmpz_add(t[29], t[12], t[26]);

    fmpz_sub(t[10], t[28], t[29]);
    fmpz_mul(t[28], t[22], t[23]);
    fmpz_add(t[29], t[24], t[25]);
    fmpz_sub(t[28], t[28], t[29]);
    fmpz_add(t[29], t[28], t[28]);
    fmpz_mul(t[28], t[2], t[2]);
    fmpz_add(t[9], t[28], t[29]);
}

