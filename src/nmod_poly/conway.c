/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "ulong_extras.h"
#include "nmod_poly.h"

extern uint8_t __nmod_poly_cp_primes0[];
extern uint16_t __nmod_poly_cp_degrees0[];
extern uint8_t __nmod_poly_numntcoeffs0[];
extern uint8_t __nmod_poly_ntcoeffs0[];

#define primes __nmod_poly_cp_primes0
#define degrees __nmod_poly_cp_degrees0
#define num_nontrivialcoeffs __nmod_poly_numntcoeffs0
#define nontrivialcoeffs __nmod_poly_ntcoeffs0
static int
conway_polynomial_lt_260(nn_ptr op, ulong prime, ulong deg)
{
    ulong ix, jx, kx;
    ulong numnt, sum;
    const uint8_t * ntcoeffs;

    /* primes in __nmod_poly_cp_primes0 are offset by 2 */
    prime -= 2;

    /* Find corresponding index. As we offset the primes by two in order to fit
     * 257 into a byte, we cannot check if the first "prime" is zero. */
    ix = 0;
    if (primes[ix] == prime)
        /* Do nothing */;
    else
    {
        do
        {
            ix++;
            if (primes[ix] == prime)
                break;
            else if (primes[ix] == 0)
                return 0;
        } while (1);
    }

    /* Degrees are stored as

       [{degrees for prime = 2}, {degrees for prime = 3}, ...]

       in the one-dimensional array __nmod_poly_cp_degrees0. Hence, we first has to
       get to the correct place in this array where the degrees for `prime' are
       stored. */
    kx = 0;
    jx = 0;
    while (jx != ix)
    {
        while (degrees[kx] < degrees[kx + 1])
            kx++;

        kx++;
        jx++;
    }

    if (deg == 1)
        /* Do nothing */;
    else
    {
        do
            kx++;
        while (degrees[kx] < deg && degrees[kx] != 1 && degrees[kx] != 0);
        /* degrees[kx] == 1 means that we've gone to the next prime, and
           degrees[kx] == 0 means that we've reached the end of the list. */
    }

    if (degrees[kx] != deg)
        return 0;

    /* Sum up all numbers of non-trivial coefficients before this one to get the
       index where our non-trivial coefficients appear. */
    sum = 0;
    for (jx = 0; jx < kx; jx++)
        sum += num_nontrivialcoeffs[jx];

    numnt = num_nontrivialcoeffs[kx];

    /* Set trivial coefficients */
    for (jx = 1; jx < deg; jx++)
        op[jx] = 0;
    op[deg] = 1;

    ntcoeffs = nontrivialcoeffs + sum;
    for (jx = 0; jx < numnt; jx++)
        op[jx] = ntcoeffs[jx];

    return 1;
}
#undef primes
#undef degrees
#undef num_nontrivialcoeffs
#undef nontrivialcoeffs

extern uint16_t __nmod_poly_cp_primes1[];
extern uint8_t __nmod_poly_cp_sm_coeffs1[];
extern uint16_t __nmod_poly_cp_md_coeffs1[];

#define primes __nmod_poly_cp_primes1
#define small_coeffs __nmod_poly_cp_sm_coeffs1
#define big_coeffs __nmod_poly_cp_md_coeffs1
static int
conway_polynomial_lt_300(nn_ptr op, ulong prime, ulong deg)
{
    ulong ix = 0;
    const uint8_t * ap;
    const uint16_t * bp;

    if (deg > 12)
        return 0;

    /* Find corresponding index */
    do
    {
        if (primes[ix] == prime)
            break;
        else if (primes[ix] == 0)
            return 0;
        ix++;
    } while (1);

    ap = small_coeffs + 23 * ix;
    bp = big_coeffs + 11 * ix;

    for (ix = 1; ix < deg; ix++)
        op[ix] = 0;
    op[deg] = 1;

    switch (deg)
    {
        case 1:
            op[0] = bp[0];
            break;
        case 2:
            op[0] = ap[0];
            op[1] = bp[1];
            break;
        case 3:
            op[0] = bp[0];
            op[1] = ap[1];
            break;
        case 4:
            op[0] = ap[0];
            op[1] = bp[2];
            op[2] = ap[2];
            break;
        case 5:
            op[0] = bp[0];
            op[1] = ap[3];
            break;
        case 6:
            op[0] = ap[0];
            op[1] = bp[3];
            op[2] = ap[4];
            op[3] = ap[5];
            op[4] = ap[6];
            break;
        case 7:
            op[0] = bp[0];
            op[1] = ap[7];
            break;
        case 8:
            op[0] = ap[0];
            op[1] = ap[8];
            op[2] = bp[4];
            op[3] = ap[9];
            op[4] = ap[10];
            break;
        case 9:
            op[0] = bp[0];
            op[1] = bp[5];
            op[2] = bp[6];
            op[3] = ap[11];
            break;
        case 10:
            op[0] = ap[0];
            op[1] = bp[7];
            op[2] = ap[12];
            op[3] = bp[8];
            op[4] = ap[13];
            op[5] = bp[9];
            op[6] = ap[14];
            break;
        case 11:
            op[0] = bp[0];
            op[1] = ap[15];
            break;
        case 12:
            op[0] = ap[0];
            op[1] = ap[16];
            op[2] = bp[10];
            op[3] = ap[17];
            op[4] = ap[18];
            op[5] = ap[19];
            op[6] = ap[20];
            op[7] = ap[21];
            op[8] = ap[22];
            break;
        default:
            FLINT_UNREACHABLE;
    }

    return 1;
}
#undef primes
#undef small_coeffs
#undef big_coeffs

extern uint16_t __nmod_poly_cp_primes2[];
extern uint8_t __nmod_poly_cp_sm_coeffs2[];
extern uint16_t __nmod_poly_cp_md_coeffs2[];

#define primes __nmod_poly_cp_primes2
#define small_coeffs __nmod_poly_cp_sm_coeffs2
#define big_coeffs __nmod_poly_cp_md_coeffs2
static int
conway_polynomial_lt_1000(nn_ptr op, ulong prime, ulong deg)
{
    ulong ix = 0;
    const uint8_t * ap;
    const uint16_t * bp;

    if (deg > 9)
        return 0;

    /* Find corresponding index */
    do
    {
        if (primes[ix] == prime)
            break;
        else if (primes[ix] == 0)
            return 0;
        ix++;
    } while (1);

    ap = small_coeffs + 8 * ix;
    bp = big_coeffs + 11 * ix;

    for (ix = 1; ix < deg; ix++)
        op[ix] = 0;
    op[deg] = 1;

    switch (deg)
    {
        case 1:
            op[0] = bp[0];
            break;
        case 2:
            op[0] = ap[0];
            op[1] = bp[1];
            break;
        case 3:
            op[0] = bp[0];
            op[1] = ap[1];
            break;
        case 4:
            op[0] = ap[0];
            op[1] = bp[2];
            op[2] = ap[2];
            break;
        case 5:
            op[0] = bp[0];
            op[1] = ap[3];
            break;
        case 6:
            op[0] = ap[0];
            op[1] = bp[3];
            op[2] = bp[4];
            op[3] = bp[5];
            op[4] = ap[4];
            break;
        case 7:
            op[0] = bp[0];
            op[1] = ap[5];
            break;
        case 8:
            op[0] = ap[0];
            op[1] = bp[6];
            op[2] = bp[7];
            op[3] = bp[8];
            op[4] = ap[6];
            break;
        case 9:
            op[0] = bp[0];
            op[1] = bp[9];
            op[2] = bp[10];
            op[3] = ap[7];
            break;
        default:
            FLINT_UNREACHABLE;
    }

    return 1;
}
#undef primes
#undef small_coeffs
#undef big_coeffs

extern uint16_t __nmod_poly_cp_primes3[];
extern uint8_t __nmod_poly_cp_sm_coeffs3[];
extern uint16_t __nmod_poly_cp_md_coeffs3[];

#define primes __nmod_poly_cp_primes3
#define small_coeffs __nmod_poly_cp_sm_coeffs3
#define big_coeffs __nmod_poly_cp_md_coeffs3
static int
conway_polynomial_lt_3371(nn_ptr op, ulong prime, ulong deg)
{
    ulong ix = 0;
    const uint8_t * ap;
    const uint16_t * bp;

    if (deg > 9 || deg == 8)
        return 0;
    else if (deg > 6 && (prime == 2689 || prime == 2797 || prime == 2833
                || prime == 3019 || prime == 3163 || prime == 3209
                || prime == 3331))
        return 0;

    /* Find corresponding index */
    do
    {
        if (primes[ix] == prime)
            break;
        else if (primes[ix] == 0)
            return 0;
        ix++;
    } while (1);

    ap = small_coeffs + 7 * ix;
    bp = big_coeffs + 8 * ix;

    for (ix = 1; ix < deg; ix++)
        op[ix] = 0;
    op[deg] = 1;

    switch (deg)
    {
        case 1:
            op[0] = bp[0];
            break;
        case 2:
            op[0] = ap[0];
            op[1] = bp[1];
            break;
        case 3:
            op[0] = bp[0];
            op[1] = ap[1];
            break;
        case 4:
            op[0] = ap[0];
            op[1] = bp[2];
            op[2] = ap[2];
            break;
        case 5:
            op[0] = bp[0];
            op[1] = ap[3];
            break;
        case 6:
            op[0] = ap[0];
            op[1] = bp[3];
            op[2] = bp[4];
            op[3] = bp[5];
            op[4] = ap[4];
            break;
        case 7:
            op[0] = bp[0];
            op[1] = ap[5];
            break;
        case 9:
            op[0] = bp[0];
            op[1] = bp[6];
            op[2] = bp[7];
            op[3] = ap[6];
            break;
        default:
            FLINT_UNREACHABLE;
    }

    return 1;
}
#undef primes
#undef small_coeffs
#undef big_coeffs

extern uint16_t __nmod_poly_cp_primes4[];
extern uint8_t __nmod_poly_cp_sm_coeffs4[];
extern uint16_t __nmod_poly_cp_md_coeffs4[];

#define primes __nmod_poly_cp_primes4
#define small_coeffs __nmod_poly_cp_sm_coeffs4
#define big_coeffs __nmod_poly_cp_md_coeffs4
static int
conway_polynomial_lt_11000(nn_ptr op, ulong prime, ulong deg)
{
    ulong ix = 0;
    const uint8_t * ap;
    const uint16_t * bp;

    if (deg > 6)
        return 0;

    /* Find corresponding index */
    do
    {
        if (primes[ix] == prime)
            break;
        else if (primes[ix] == 0)
            return 0;
        ix++;
    } while (1);

    ap = small_coeffs + 5 * ix;
    bp = big_coeffs + 6 * ix;

    for (ix = 1; ix < deg; ix++)
        op[ix] = 0;
    op[deg] = 1;

    switch (deg)
    {
        case 1:
            op[0] = bp[0];
            break;
        case 2:
            op[0] = ap[0];
            op[1] = bp[1];
            break;
        case 3:
            op[0] = bp[0];
            op[1] = ap[1];
            break;
        case 4:
            op[0] = ap[0];
            op[1] = bp[2];
            op[2] = ap[2];
            break;
        case 5:
            op[0] = bp[0];
            op[1] = ap[3];
            break;
        case 6:
            op[0] = ap[0];
            op[1] = bp[3];
            op[2] = bp[4];
            op[3] = bp[5];
            op[4] = ap[4];
            break;
        default:
            FLINT_UNREACHABLE;
    }

    return 1;
}
#undef primes
#undef small_coeffs
#undef big_coeffs

extern uint16_t __nmod_poly_cp_primes5[];
extern uint8_t __nmod_poly_cp_sm_coeffs5[];
extern uint16_t __nmod_poly_cp_md_coeffs5[];

#define primes __nmod_poly_cp_primes5
#define small_coeffs __nmod_poly_cp_sm_coeffs5
#define big_coeffs __nmod_poly_cp_md_coeffs5
static int
conway_polynomial_lt_65536(nn_ptr op, ulong prime, ulong deg)
{
    ulong ix = 0;
    const uint8_t * ap;
    const uint16_t * bp;

    if (deg > 4)
        return 0;

    /* Find corresponding index */
    do
    {
        if (primes[ix] == prime)
            break;
        else if (primes[ix] == 0)
            return 0;
        ix++;
    } while (1);

    ap = small_coeffs + 3 * ix;
    bp = big_coeffs + 3 * ix;

    for (ix = 1; ix < deg; ix++)
        op[ix] = 0;
    op[deg] = 1;

    switch (deg)
    {
        case 1:
            op[0] = bp[0];
            break;
        case 2:
            op[0] = ap[0];
            op[1] = bp[1];
            break;
        case 3:
            op[0] = bp[0];
            op[1] = ap[1];
            break;
        case 4:
            op[0] = ap[0];
            op[1] = bp[2];
            op[2] = ap[2];
            break;
        default:
            FLINT_UNREACHABLE;
    }

    return 1;
}
#undef primes
#undef small_coeffs
#undef big_coeffs

extern uint16_t __nmod_poly_cp_primes6[];
extern uint8_t __nmod_poly_cp_sm_coeffs6[];
extern uint32_t __nmod_poly_cp_lg_coeffs6[];

#define primes __nmod_poly_cp_primes6
#define small_coeffs __nmod_poly_cp_sm_coeffs6
#define big_coeffs __nmod_poly_cp_lg_coeffs6
static int
conway_polynomial_lt_109988(nn_ptr op, ulong prime, ulong deg)
{
    ulong ix = 0;
    const uint8_t * ap;
    const uint32_t * bp;

    if (deg != 4)
        return 0;

    /* primes in __nmod_poly_cp_primes6 are offset by 2^16 */
    prime -= 1 << 16;

    /* Find corresponding index */
    do
    {
        if (primes[ix] == prime)
            break;
        else if (primes[ix] == 0)
            return 0;
        ix++;
    } while (1);

    ap = small_coeffs + 2 * ix;
    bp = big_coeffs + 1 * ix;

    op[0] = ap[0];
    op[1] = bp[0];
    op[2] = ap[1];
    op[3] = 0;
    op[4] = 1;

    return 1;
}
#undef primes
#undef small_coeffs
#undef big_coeffs

int
_nmod_poly_conway(nn_ptr op, ulong prime, slong deg)
{
    if (deg <= 0)
        return 0;
    else if (prime < 2)
        return 0;

    if (prime < 260)
        return conway_polynomial_lt_260(op, prime, deg);
    else if (prime % 2 == 0)
        return 0;
    else if (prime < 300)
        return conway_polynomial_lt_300(op, prime, deg);
    else if (prime < 1000)
        return conway_polynomial_lt_1000(op, prime, deg);
    else if (prime < 3371)
        return conway_polynomial_lt_3371(op, prime, deg);
    else if (prime < 11000)
        return conway_polynomial_lt_11000(op, prime, deg);
    else if (prime < 65536)
        return conway_polynomial_lt_65536(op, prime, deg);
    else if (prime < 109988)
        return conway_polynomial_lt_109988(op, prime, deg);
    else
        return 0;
}

ulong
_nmod_poly_conway_rand(slong * degree, flint_rand_t state, int type)
{
    ulong prime;

    switch (type)
    {
        case 0: /* prime whatever, degree whatever */
        case 1: /* degree < 15 */
            do
                prime = n_randprime(state, 2 + n_randint(state, 16), 1);
            while (prime > 109987);
            break;
        case 2: /* prime < 2^10 */
        case 3: /* prime < 2^10 and degree < 15 */
            prime = n_randprime(state, 2 + n_randint(state, 9), 1);
            break;

        default: flint_throw(FLINT_ERROR, "wrong type in %s", __func__);
    }

    if (prime < 260)
    {
        /* Find the position it corresponds to in __nmod_poly_cp_primes0, and
         * then generate degree based on __nmod_poly_cp_degrees0. */
#define primes __nmod_poly_cp_primes0
#define degrees __nmod_poly_cp_degrees0
        slong ix = 0, jx, kx;

        /* Primes are offset by 2 in table */
        prime -= 2;

        /* Search for prime's index */
        while (prime != primes[ix])
            ix++;

        /* Reset prime */
        prime += 2;

        /* Search for starting index of degrees for prime */
        kx = 0;
        jx = 0;
        while (jx != ix)
        {
            while (degrees[kx] < degrees[kx + 1])
                kx++;

            kx++;
            jx++;
        }

        /* Now kx points to the starting degree for prime */
        jx = 0;
        if (type % 2 == 1)
        {
            do
                jx++;
            while (degrees[kx + jx] < 15);
        }
        else
        {
            do
                jx++;
            while (degrees[kx + jx] != 1 && degrees[kx + jx] != 0);
        }

        /* The degrees we will be using are degrees[kx + 0], ...,
         * degrees[kx + jx - 1]. */
        *degree = degrees[kx + n_randint(state, jx)];
#undef primes
#undef degrees
    }
    else if (prime < 300)
    {
        *degree = 1 + n_randint(state, 12);
    }
    else if (prime < 1000)
    {
        *degree = 1 + n_randint(state, 9);
    }
    else if (prime < 3371)
    {
        switch (prime)
        {
            case 2689:
            case 2797:
            case 2833:
            case 3019:
            case 3163:
            case 3209:
            case 3331: *degree = 1 + n_randint(state, 6);
                       break;

            default: *degree = 1 + n_randint(state, 8);
                     if (*degree == 8)
                         *degree += 1;
        }
    }
    else if (prime < 11000)
        *degree = 1 + n_randint(state, 6);
    else if (prime < 65536)
        *degree = 1 + n_randint(state, 4);
    else
        *degree = 4;

    return prime;
}
