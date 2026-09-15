/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Opus 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#define TOOM_SCALAR_N_MAX FMPZ_POLY_TOOM_SCALAR_N_MAX

static const int32_t toom_m[2][9][9] = {
    {
        { 0 },
        { -1, 1 },
        { 15, -24, 1 },
        { -7, 14, -1, 2 },
        { 210, -480, 45, -160, 1 },
        { -198, 495, -55, 264, -3, 1 },
        { 3003, -8008, 1001, -5824, 91, -56, 1 },
        { -715, 2002, -273, 1820, -35, 30, -1, 4 },
        { 43758, -127296, 18564, -137088, 3060, -3264, 153, -1152, 1 },
    },
    {
        { 0 },
        { -1, 1 },
        { 5, -8, 1 },
        { -7, 14, -3, 2 },
        { 42, -96, 27, -32, 1 },
        { -66, 165, -55, 88, -5, 1 },
        { 429, -1144, 429, -832, 65, -24, 1 },
        { -715, 2002, -819, 1820, -175, 90, -7, 4 },
        { 4862, -14144, 6188, -15232, 1700, -1088, 119, -128, 1 },
    },
};

static const uint32_t toom_h[8][7] = {
    { 0 },
    { 5, 21, 85, 341, 1365, 5461, 21845 },
    { 14, 147, 1408, 13013, 118482, 1071799, 0 },
    { 30, 627, 11440, 196053, 3255330, 0, 0 },
    { 55, 2002, 61490, 1733303, 0, 0, 0 },
    { 91, 5278, 251498, 0, 0, 0, 0 },
    { 140, 12138, 0, 0, 0, 0, 0 },
    { 204, 0, 0, 0, 0, 0, 0 },
};

static const uint16_t tmid_r0[21] = {
    0, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
};

static const uint16_t tmid_rinf[21] = {
    0, 0, 0, 1, 2, 2, 8, 8, 24, 8, 32, 32, 160, 32, 96, 32, 224, 32, 128,
    128, 1152
};

static const uint16_t tmid_rpair[21][10] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 4, 1, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 4, 1, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 4, 1, 4, 0, 0, 0, 0, 0, 0 },
    { 0, 4, 1, 4, 0, 0, 0, 0, 0, 0 },
    { 0, 16, 4, 16, 1, 0, 0, 0, 0, 0 },
    { 0, 16, 4, 16, 1, 0, 0, 0, 0, 0 },
    { 0, 16, 4, 16, 1, 16, 0, 0, 0, 0 },
    { 0, 16, 4, 16, 1, 16, 0, 0, 0, 0 },
    { 0, 16, 4, 16, 1, 16, 4, 0, 0, 0 },
    { 0, 16, 4, 16, 1, 16, 4, 0, 0, 0 },
    { 0, 16, 4, 16, 1, 16, 4, 16, 0, 0 },
    { 0, 16, 4, 16, 1, 16, 4, 16, 0, 0 },
    { 0, 64, 16, 64, 4, 64, 16, 64, 1, 0 },
    { 0, 64, 16, 64, 4, 64, 16, 64, 1, 0 },
    { 0, 64, 16, 64, 4, 64, 16, 64, 1, 64 },
};

#if FLINT_BITS == 64

static const ulong toom_D[2][9] = {
    {
        0, 3, 360, 2520, 1814400, 59875200, UWORD(43589145600),
        UWORD(653837184000), UWORD(3201186852864000)
    },
    {
        0, 3, 120, 2520, 362880, 19958400, UWORD(6227020800),
        UWORD(653837184000), UWORD(355687428096000)
    },
};

static const ulong tmid_Dtop[2][10] = {
    {
        1, 1, 3, 360, 2520, 1814400, 59875200, UWORD(43589145600),
        UWORD(653837184000), UWORD(3201186852864000)
    },
    {
        1, 1, 3, 120, 2520, 362880, 19958400, UWORD(6227020800),
        UWORD(653837184000), UWORD(355687428096000)
    },
};

static const ulong tmid_A[2][10][9] = {
    {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 120, 1, 0, 0, 0, 0, 0, 0 },
        { 0, 840, 7, 1, 0, 0, 0, 0, 0 },
        { 0, 604800, 5040, 720, 1, 0, 0, 0, 0 },
        { 0, 19958400, 166320, 23760, 33, 1, 0, 0, 0 },
        { 0, UWORD(14529715200), 121080960, 17297280, 24024, 728, 1, 0, 0 },
        { 0, UWORD(217945728000), 1816214400, 259459200, 360360, 10920, 15, 1, 0 },
        { 0, UWORD(1067062284288000), UWORD(8892185702400), UWORD(1270312243200), 1764322560, 53464320, 73440, 4896, 1 },
    },
    {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 40, 1, 0, 0, 0, 0, 0, 0 },
        { 0, 840, 21, 1, 0, 0, 0, 0, 0 },
        { 0, 120960, 3024, 144, 1, 0, 0, 0, 0 },
        { 0, 6652800, 166320, 7920, 55, 1, 0, 0, 0 },
        { 0, 2075673600, 51891840, 2471040, 17160, 312, 1, 0, 0 },
        { 0, UWORD(217945728000), UWORD(5448643200), 259459200, 1801800, 32760, 105, 1, 0 },
        { 0, UWORD(118562476032000), UWORD(2964061900800), UWORD(141145804800), 980179200, 17821440, 57120, 544, 1 },
    },
};

static const ulong tmid_G1[21] = {
    0, 0, 0, 1, 2, 12, 24, 2880, 2880, 29030400, 80640,
    UWORD(29262643200), 58060800, UWORD(2317601341440000), 1916006400,
    UWORD(11931011705733120000), UWORD(1394852659200),
    UWORD(3487131648000), UWORD(83691159552000), UWORD(127682666496000),
    UWORD(409751917166592000)
};

static const ulong tmid_G2[21] = {
    0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    UWORD(2092278988800), 1, UWORD(233139658752000), 1
};

static const ulong tmid_rT[21] = {
    0, 0, 0, 1, 0, 2, 0, 24, 0, 2880, 0, 80640, 0, 58060800, 0,
    1916006400, 0, UWORD(1394852659200), 0, UWORD(83691159552000), 0
};

static const ulong tmid_gT[21] = {
    0, 0, 0, 1, 0, 6, 0, 120, 0, 10080, 0, 362880, 0, 39916800, 0,
    UWORD(6227020800), 0, UWORD(5230697472000), 0,
    UWORD(355687428096000), 0
};

static const ulong tmid_sc0[21] = {
    0, 0, 0, 1, 1, 1, 12, 12, 1440, 1440, 40320, 40320, 29030400,
    29030400, 958003200, 958003200, UWORD(697426329600),
    UWORD(697426329600), UWORD(41845579776000), UWORD(41845579776000),
    UWORD(204875958583296000)
};

static const ulong tmid_scinf[21] = {
    0, 0, 0, 1, 0, 1, 0, 3, 0, 360, 0, 2520, 0, 1814400, 0, 59875200, 0,
    UWORD(43589145600), 0, UWORD(653837184000), 0
};

static const ulong tmid_qN[21] = {
    0, 0, 0, 1, 0, 16, 0, 729, 0, 65536, 0, 9765625, 0,
    UWORD(2176782336), 0, UWORD(678223072849), 0, UWORD(281474976710656),
    0, UWORD(150094635296999121), 0
};

static const ulong tmid_fA[21][10] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 2, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 2, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 3, 6, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 3, 6, 3, 0, 0, 0, 0, 0, 0 },
    { 0, 1440, 2880, 1440, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 2, 1, 4, 0, 0, 0, 0, 0 },
    { 0, 2520, 5040, 2520, 10080, 0, 0, 0, 0, 0 },
    { 0, 5, 10, 5, 20, 5, 0, 0, 0, 0 },
    { 0, 3628800, 7257600, 3628800, 14515200, 3628800, 0, 0, 0, 0 },
    { 0, 3, 6, 3, 12, 3, 6, 0, 0, 0 },
    {
        0, 59875200, 119750400, 59875200, 239500800, 59875200, 119750400,
        0, 0, 0
    },
    { 0, 7, 14, 7, 28, 7, 14, 7, 0, 0 },
    {
        0, UWORD(348713164800), UWORD(697426329600), UWORD(348713164800),
        UWORD(1394852659200), UWORD(348713164800), UWORD(697426329600),
        UWORD(348713164800), 0, 0
    },
    { 0, 1, 2, 1, 4, 1, 2, 1, 8, 0 },
    {
        0, UWORD(653837184000), UWORD(1307674368000),
        UWORD(653837184000), UWORD(2615348736000), UWORD(653837184000),
        UWORD(1307674368000), UWORD(653837184000), UWORD(5230697472000), 0
    },
    { 0, 9, 18, 9, 36, 9, 18, 9, 72, 9 },
};

static const ulong tmid_pwodd[21][10] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 16, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 64, 2187, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 256, 19683, 65536, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 1024, 177147, 1048576, 48828125, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    {
        0, 1, 4096, 1594323, 16777216, 1220703125, UWORD(6530347008), 0,
        0, 0
    },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    {
        0, 1, 16384, 14348907, 268435456, UWORD(30517578125),
        UWORD(235092492288), UWORD(4747561509943), 0, 0
    },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    {
        0, 1, 65536, 129140163, UWORD(4294967296), UWORD(762939453125),
        UWORD(8463329722368), UWORD(232630513987207),
        UWORD(281474976710656), 0
    },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    {
        0, 1, 262144, 1162261467, UWORD(68719476736),
        UWORD(19073486328125), UWORD(304679870005248),
        UWORD(11398895185373143), UWORD(18014398509481984),
        UWORD(1350851717672992089)
    },
};

static const ulong tmid_pwev[21][10] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 16, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 64, 6561, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 256, 59049, 65536, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 1024, 531441, 1048576, 244140625, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    {
        0, 1, 4096, 4782969, 16777216, UWORD(6103515625),
        UWORD(19591041024), 0, 0, 0
    },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    {
        0, 1, 16384, 43046721, 268435456, UWORD(152587890625),
        UWORD(705277476864), UWORD(33232930569601), 0, 0
    },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    {
        0, 1, 65536, 387420489, UWORD(4294967296), UWORD(3814697265625),
        UWORD(25389989167104), UWORD(1628413597910449),
        UWORD(281474976710656), 0
    },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

#else

static const ulong toom_D[2][9] = {
    { 0, 3, 360, 2520, 1814400, 59875200, 0, 0, 0 },
    { 0, 3, 120, 2520, 362880, 19958400, 0, 0, 0 },
};

static const ulong tmid_Dtop[2][10] = {
    { 1, 1, 3, 360, 2520, 1814400, 59875200, 0, 0, 0 },
    { 1, 1, 3, 120, 2520, 362880, 19958400, 0, 0, 0 },
};

static const ulong tmid_A[2][10][9] = {
    {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 120, 1, 0, 0, 0, 0, 0, 0 },
        { 0, 840, 7, 1, 0, 0, 0, 0, 0 },
        { 0, 604800, 5040, 720, 1, 0, 0, 0, 0 },
        { 0, 19958400, 166320, 23760, 33, 1, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    },
    {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 40, 1, 0, 0, 0, 0, 0, 0 },
        { 0, 840, 21, 1, 0, 0, 0, 0, 0 },
        { 0, 120960, 3024, 144, 1, 0, 0, 0, 0 },
        { 0, 6652800, 166320, 7920, 55, 1, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    },
};

static const ulong tmid_G1[21] = {
    0, 0, 0, 1, 2, 12, 24, 2880, 2880, 29030400, 80640, 241920, 58060800,
    34214400, 0, 0, 0, 0, 0, 0, 0
};

static const ulong tmid_G2[21] = {
    0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 120960, 1, 67737600, 0, 0, 0, 0, 0,
    0, 0
};

static const ulong tmid_rT[21] = {
    0, 0, 0, 1, 0, 2, 0, 24, 0, 2880, 0, 80640, 0, 58060800, 0, 0, 0, 0,
    0, 0, 0
};

static const ulong tmid_gT[21] = {
    0, 0, 0, 1, 0, 6, 0, 120, 0, 10080, 0, 362880, 0, 39916800, 0, 0, 0,
    0, 0, 0, 0
};

static const ulong tmid_sc0[21] = {
    0, 0, 0, 1, 1, 1, 12, 12, 1440, 1440, 40320, 40320, 29030400,
    29030400, 0, 0, 0, 0, 0, 0, 0
};

static const ulong tmid_scinf[21] = {
    0, 0, 0, 1, 0, 1, 0, 3, 0, 360, 0, 2520, 0, 1814400, 0, 0, 0, 0, 0,
    0, 0
};

static const ulong tmid_qN[21] = {
    0, 0, 0, 1, 0, 16, 0, 729, 0, 65536, 0, 9765625, 0,
    UWORD(2176782336), 0, 0, 0, 0, 0, 0, 0
};

static const ulong tmid_fA[21][10] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 2, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 2, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 3, 6, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 3, 6, 3, 0, 0, 0, 0, 0, 0 },
    { 0, 1440, 2880, 1440, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 2, 1, 4, 0, 0, 0, 0, 0 },
    { 0, 2520, 5040, 2520, 10080, 0, 0, 0, 0, 0 },
    { 0, 5, 10, 5, 20, 5, 0, 0, 0, 0 },
    { 0, 3628800, 7257600, 3628800, 14515200, 3628800, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static const ulong tmid_pwodd[21][10] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 16, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 64, 2187, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 256, 19683, 65536, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 1024, 177147, 1048576, 48828125, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static const ulong tmid_pwev[21][10] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 16, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 64, 6561, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 256, 59049, 65536, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 1024, 531441, 1048576, 244140625, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

#endif

static ulong
toom_pow_ui(ulong k, ulong e)
{
    ulong r = 1;
    while (e--) r *= k;
    return r;
}

/* ------------------------------------------------------------------ */
/* evaluation                                                         */
/* ------------------------------------------------------------------ */

/*
    E = even part of a at k: a_0 + a_2 k^2 + ...  Requires len >= 3,
    k >= 1, E distinct from a.
*/
static void
toom_eval_even(fmpz_t E, const fmpz * a, slong len, ulong k, slong st)
{
    slong i;

    if (k == 1)
    {
        fmpz_add(E, a + 0, a + 2 * st);

        for (i = 4; i < len; i += 2)
            fmpz_add(E, E, a + i * st);
    }
    else
    {
        ulong kk = k * k, pw;

        fmpz_mul_ui(E, a + 2 * st, kk);

        for (pw = kk, i = 4; i < len; i += 2)
        {
            pw *= kk;
            fmpz_addmul_ui(E, a + i * st, pw);
        }

        fmpz_add(E, E, a + 0);
    }
}

/*
    O = odd part of a at k: a_1 k + a_3 k^3 + ...  Requires len >= 2,
    k >= 1, and NOT (k == 1 && len <= 3) -- that case is a bare copy of
    a_1, which callers fold directly instead.
*/
static void
toom_eval_odd(fmpz_t O, const fmpz * a, slong len, ulong k, slong st)
{
    slong i;

    if (k == 1)
    {
        fmpz_add(O, a + st, a + 3 * st);

        for (i = 5; i < len; i += 2)
            fmpz_add(O, O, a + i * st);
    }
    else
    {
        ulong kk = k * k, pw;

        fmpz_mul_ui(O, a + st, k);

        for (pw = k, i = 3; i < len; i += 2)
        {
            pw *= kk;
            fmpz_addmul_ui(O, a + i * st, pw);
        }
    }
}

/* r = a(k), len >= 2, k >= 1 */
static void
toom_eval(fmpz_t r, const fmpz * a, slong len, ulong k, slong st)
{
    slong i;

    if (k == 1)
    {
        fmpz_add(r, a + 0, a + st);

        for (i = 2; i < len; i++)
            fmpz_add(r, r, a + i * st);
    }
    else
    {
        ulong pw;

        fmpz_mul_ui(r, a + st, k);

        for (pw = k, i = 2; i < len; i++)
        {
            pw *= k;
            fmpz_addmul_ui(r, a + i * st, pw);
        }

        fmpz_add(r, r, a + 0);
    }
}

/* ------------------------------------------------------------------ */
/* interpolation                                                      */
/* ------------------------------------------------------------------ */

/*
    v[0], v[2], ..., v[2(n-1)] hold the scaled samples of
    x_0 + x_1 u + ... + x_(n-1) u^(n-1) at u = 1, 4, ..., n^2;
    on exit v[2i] = x_i.  In place, no scratch: row r accumulates
    into v[2r], which is dead as a sample once rows above are done.
    Handles n <= 1 as a no-op.
*/
static void
toom_interp(fmpz * v, slong n, int par, slong r0)
{
    slong r, i, m;

    for (r = n - 1; r >= FLINT_MAX(r0, 1); r--)
    {
        m = toom_m[par][r][r];      /* always 1, 2 or 4 */
        if (m != 1)
            fmpz_mul_2exp(v + 2 * r, v + 2 * r, m == 2 ? 1 : 2);

        for (i = 0; i < r; i++)
        {
            m = toom_m[par][r][i];

            if      (m ==  1)
                fmpz_add(v + 2 * r, v + 2 * r, v + 2 * i);
            else if (m == -1)
                fmpz_sub(v + 2 * r, v + 2 * r, v + 2 * i);
            else if (m > 0)
                fmpz_addmul_ui(v + 2 * r, v + 2 * i, (ulong)  m);
            else
                fmpz_submul_ui(v + 2 * r, v + 2 * i, (ulong) -m);
        }

        fmpz_divexact_ui(v + 2 * r, v + 2 * r, toom_D[par][r]);

        for (i = r + 1; i < n; i++)
            fmpz_submul_ui(v + 2 * r, v + 2*i, toom_h[r][i - r - 1]);
    }

    if (r0 == 0)
        for (i = 1; i < n; i++)
            fmpz_sub(v, v, v + 2 * i);
}

/* ------------------------------------------------------------------ */
/* hardcoded smallest shapes (measurably faster than the generic
   driver at the bit sizes where these shapes are used)             */
/* ------------------------------------------------------------------ */

/* c = a * b, len 2 x 2; points 0, 1, inf */
static void
toom_mul_2_2(fmpz * c, const fmpz * a, const fmpz * b, int sq)
{
    fmpz_add(c + 1, a + 0, a + 1);              /* a(1)          */
    if (sq)
        fmpz_mul(c + 1, c + 1, c + 1);
    else
    {
        fmpz_add(c + 2, b + 0, b + 1);          /* b(1)          */
        fmpz_mul(c + 1, c + 1, c + 2);          /* v1            */
    }
    if (sq)
    {
        fmpz_mul(c + 0, a + 0, a + 0);
        fmpz_mul(c + 2, a + 1, a + 1);
    }
    else
    {
        fmpz_mul(c + 0, a + 0, b + 0);          /* c0, final     */
        fmpz_mul(c + 2, a + 1, b + 1);          /* c2, final     */
    }
    fmpz_sub(c + 1, c + 1, c + 0);
    fmpz_sub(c + 1, c + 1, c + 2);              /* c1, final     */
}

/* c = a * b, len 3 x 2; points 0, +-1, inf */
static void
toom_mul_3_2(fmpz * c, const fmpz * a, const fmpz * b)
{
    fmpz_add(c + 1, a + 0, a + 2);
    fmpz_add(c + 2, c + 1, a + 1);              /* a(1)          */
    fmpz_sub(c + 1, c + 1, a + 1);              /* a(-1)         */
    fmpz_add(c + 3, b + 0, b + 1);              /* b(1)          */
    fmpz_mul(c + 2, c + 2, c + 3);              /* v1            */
    fmpz_sub(c + 3, b + 0, b + 1);              /* b(-1)         */
    fmpz_mul(c + 1, c + 1, c + 3);              /* vm1           */
    fmpz_sub(c + 3, c + 2, c + 1);
    fmpz_add(c + 2, c + 2, c + 1);
    fmpz_tdiv_q_2exp(c + 1, c + 3, 1);          /* J1 = c1 + c3  */
    fmpz_tdiv_q_2exp(c + 2, c + 2, 1);          /* E1 = c0 + c2  */
    fmpz_mul(c + 0, a + 0, b + 0);              /* c0, final     */
    fmpz_mul(c + 3, a + 2, b + 1);              /* c3, final     */
    fmpz_sub(c + 2, c + 2, c + 0);              /* c2, final     */
    fmpz_sub(c + 1, c + 1, c + 3);              /* c1, final     */
}

/* c = a * b, len 3 x 3; points 0, +-1, 2, inf */
static void
toom_mul_3_3(fmpz * c, const fmpz * a, const fmpz * b, int sq)
{
    fmpz_add(c + 1, a + 0, a + 2);
    fmpz_add(c + 2, c + 1, a + 1);              /* a(1)          */
    fmpz_sub(c + 1, c + 1, a + 1);              /* a(-1)         */
    if (sq)
    {
        fmpz_mul(c + 2, c + 2, c + 2);          /* v1            */
        fmpz_mul(c + 1, c + 1, c + 1);          /* vm1           */
    }
    else
    {
        fmpz_add(c + 3, b + 0, b + 2);
        fmpz_add(c + 4, c + 3, b + 1);          /* b(1)          */
        fmpz_sub(c + 3, c + 3, b + 1);          /* b(-1)         */
        fmpz_mul(c + 2, c + 2, c + 4);          /* v1            */
        fmpz_mul(c + 1, c + 1, c + 3);          /* vm1           */
    }
    fmpz_sub(c + 3, c + 2, c + 1);
    fmpz_add(c + 2, c + 2, c + 1);
    fmpz_tdiv_q_2exp(c + 1, c + 3, 1);          /* J1 = c1 + c3  */
    fmpz_tdiv_q_2exp(c + 2, c + 2, 1);          /* E1            */

    fmpz_mul_2exp(c + 3, a + 2, 2);
    fmpz_addmul_ui(c + 3, a + 1, 2);
    fmpz_add(c + 3, c + 3, a + 0);              /* a(2)          */
    if (sq)
        fmpz_mul(c + 3, c + 3, c + 3);          /* v2            */
    else
    {
        fmpz_mul_2exp(c + 4, b + 2, 2);
        fmpz_addmul_ui(c + 4, b + 1, 2);
        fmpz_add(c + 4, c + 4, b + 0);          /* b(2)          */
        fmpz_mul(c + 3, c + 3, c + 4);          /* v2            */
    }

    if (sq)
    {
        fmpz_mul(c + 0, a + 0, a + 0);
        fmpz_mul(c + 4, a + 2, a + 2);
    }
    else
    {
        fmpz_mul(c + 0, a + 0, b + 0);          /* c0, final     */
        fmpz_mul(c + 4, a + 2, b + 2);          /* c4, final     */
    }

    fmpz_sub(c + 2, c + 2, c + 0);
    fmpz_sub(c + 2, c + 2, c + 4);              /* c2, final     */
    fmpz_sub(c + 3, c + 3, c + 0);
    fmpz_submul_ui(c + 3, c + 2, 4);
    fmpz_submul_ui(c + 3, c + 4, 16);           /* 2(c1 + 4c3)   */
    fmpz_submul_ui(c + 3, c + 1, 2);            /* 6 c3          */
    fmpz_divexact_ui(c + 3, c + 3, 6);          /* c3, final     */
    fmpz_sub(c + 1, c + 1, c + 3);              /* c1, final     */
}

/*
    c = a * b (or a^2 if sq), la = lb = 4, N = 7.  Nodes 0, +-1, +-2, 3,
    infinity.  Slots: E_k -> c[2k], J_k -> c[2k-1], v_3 -> c[5],
    endpoints c[0], c[6]; the endpoint slots double as scratch until the
    pair products are folded.
*/
static void
toom_mul_4_4(fmpz * c, const fmpz * a, const fmpz * b, int sq)
{
    /* ---- k = 1 ---------------------------------------------------- */

    fmpz_add(c + 6, a + 0, a + 2);
    fmpz_add(c + 0, a + 1, a + 3);
    fmpz_add(c + 2, c + 6, c + 0);              /* a(1)          */
    fmpz_sub(c + 1, c + 6, c + 0);              /* a(-1)         */
    if (sq)
    {
        fmpz_mul(c + 2, c + 2, c + 2);          /* v1            */
        fmpz_mul(c + 1, c + 1, c + 1);          /* vm1           */
    }
    else
    {
        fmpz_add(c + 3, b + 0, b + 2);
        fmpz_add(c + 4, b + 1, b + 3);
        fmpz_add(c + 5, c + 3, c + 4);          /* b(1)          */
        fmpz_sub(c + 3, c + 3, c + 4);          /* b(-1)         */
        fmpz_mul(c + 2, c + 2, c + 5);          /* v1            */
        fmpz_mul(c + 1, c + 1, c + 3);          /* vm1           */
    }
    fmpz_sub(c + 0, c + 2, c + 1);
    fmpz_add(c + 2, c + 2, c + 1);
    fmpz_tdiv_q_2exp(c + 2, c + 2, 1);          /* E1            */
    fmpz_tdiv_q_2exp(c + 1, c + 0, 1);          /* J1            */

    /* ---- k = 2 ---------------------------------------------------- */

    fmpz_mul_2exp(c + 0, a + 2, 2);
    fmpz_add(c + 0, c + 0, a + 0);              /* a0 + 4 a2     */
    fmpz_mul_2exp(c + 6, a + 3, 3);
    fmpz_addmul_ui(c + 6, a + 1, 2);            /* 2 a1 + 8 a3   */
    fmpz_add(c + 4, c + 0, c + 6);              /* a(2)          */
    fmpz_sub(c + 3, c + 0, c + 6);              /* a(-2)         */
    if (sq)
    {
        fmpz_mul(c + 4, c + 4, c + 4);          /* v2            */
        fmpz_mul(c + 3, c + 3, c + 3);          /* vm2           */
    }
    else
    {
        fmpz_mul_2exp(c + 0, b + 2, 2);
        fmpz_add(c + 0, c + 0, b + 0);
        fmpz_mul_2exp(c + 6, b + 3, 3);
        fmpz_addmul_ui(c + 6, b + 1, 2);
        fmpz_add(c + 5, c + 0, c + 6);          /* b(2)          */
        fmpz_sub(c + 0, c + 0, c + 6);          /* b(-2)         */
        fmpz_mul(c + 4, c + 4, c + 5);          /* v2            */
        fmpz_mul(c + 3, c + 3, c + 0);          /* vm2           */
    }
    fmpz_sub(c + 6, c + 4, c + 3);
    fmpz_add(c + 4, c + 4, c + 3);
    fmpz_tdiv_q_2exp(c + 4, c + 4, 1);          /* E2            */
    fmpz_tdiv_q_2exp(c + 3, c + 6, 2);          /* J2            */

    /* ---- extra point 3, endpoints --------------------------------- */

    fmpz_mul_ui(c + 0, a + 3, 27);
    fmpz_addmul_ui(c + 0, a + 2, 9);
    fmpz_addmul_ui(c + 0, a + 1, 3);
    fmpz_add(c + 0, c + 0, a + 0);              /* a(3)          */
    if (sq)
        fmpz_mul(c + 5, c + 0, c + 0);          /* v3            */
    else
    {
        fmpz_mul_ui(c + 6, b + 3, 27);
        fmpz_addmul_ui(c + 6, b + 2, 9);
        fmpz_addmul_ui(c + 6, b + 1, 3);
        fmpz_add(c + 6, c + 6, b + 0);          /* b(3)          */
        fmpz_mul(c + 5, c + 0, c + 6);          /* v3            */
    }

    if (sq)
    {
        fmpz_mul(c + 0, a + 0, a + 0);
        fmpz_mul(c + 6, a + 3, a + 3);
    }
    else
    {
        fmpz_mul(c + 0, a + 0, b + 0);          /* c0, final     */
        fmpz_mul(c + 6, a + 3, b + 3);          /* c6, final     */
    }

    /* ---- even half: samples of c2 + c4 u at u = 1, 4 --------------- */

    fmpz_sub(c + 2, c + 2, c + 0);
    fmpz_sub(c + 2, c + 2, c + 6);              /* c2 + c4       */
    fmpz_sub(c + 4, c + 4, c + 0);
    fmpz_submul_ui(c + 4, c + 6, 64);
    fmpz_tdiv_q_2exp(c + 4, c + 4, 2);          /* c2 + 4 c4     */

    fmpz_sub(c + 4, c + 4, c + 2);
    fmpz_divexact_ui(c + 4, c + 4, 3);          /* c4, final     */
    fmpz_sub(c + 2, c + 2, c + 4);              /* c2, final     */

    /* ---- odd half -------------------------------------------------- */

    fmpz_sub(c + 5, c + 5, c + 0);
    fmpz_submul_ui(c + 5, c + 2, 9);
    fmpz_submul_ui(c + 5, c + 4, 81);
    fmpz_submul_ui(c + 5, c + 6, 729);          /* 3(c1 + 9c3 + 81c5) */

    fmpz_addmul_ui(c + 5, c + 1, 5);
    fmpz_submul_ui(c + 5, c + 3, 8);
    fmpz_divexact_ui(c + 5, c + 5, 120);        /* c5, final     */

    fmpz_sub(c + 3, c + 3, c + 1);
    fmpz_divexact_ui(c + 3, c + 3, 3);
    fmpz_submul_ui(c + 3, c + 5, 5);            /* c3, final     */

    fmpz_sub(c + 1, c + 1, c + 3);
    fmpz_sub(c + 1, c + 1, c + 5);              /* c1, final     */
}

/*
    c = a * b (or a^2 if sq), la >= lb >= 2, N = la + lb - 1 >= 3,
    no aliasing.  Uses the N output slots only:

        E_k -> c[2k], J_k -> c[2k-1], k = 1..p;
        extra product -> c[N-2] (N odd);
        c[0], c[N-1]: scratch during the pairs (they are the last two
        slots still free), endpoint products afterwards.
*/
static void
toom_mul_generic(fmpz * c, const fmpz * a, slong la,
                 const fmpz * b, slong lb, int sq,
                 slong wlo, int rev)
{
    slong N = la + lb - 1;
    slong p = (N - 2) >> 1;
    int extra = N & 1;
    ulong q = (ulong) (p + 1);
    slong st = rev ? -1 : 1;
    const fmpz * ap = rev ? a + la - 1 : a;
    const fmpz * bp = rev ? b + lb - 1 : b;
    fmpz * X = c + 0, * Y = c + N - 1;
    slong i, k;
    slong r0_e, r0_o;

    r0_e = (extra || wlo <= 2) ? 0 : (wlo - 1) / 2;
    r0_o = wlo / 2;

    /* ---- pointwise products -------------------------------------- */

    for (k = 1; k <= p; k++)
    {
        fmpz * A = c + 2*k - 1, * B = c + 2*k;

        /* a(k) -> B, a(-k) -> A  (la >= 3 whenever p >= 1) */
        toom_eval_even(B, ap, la, (ulong) k, st);
        if (k == 1 && la == 3)
        {
            fmpz_sub(A, B, ap + st);
            fmpz_add(B, B, ap + st);
        }
        else
        {
            toom_eval_odd(X, ap, la, (ulong) k, st);
            fmpz_sub(A, B, X);
            fmpz_add(B, B, X);
        }

        if (sq)
        {
            fmpz_mul(B, B, B);                      /* v_k    */
            fmpz_mul(A, A, A);                      /* v_{-k} */
        }
        else if (lb == 2)
        {
            if (k == 1)
            {
                fmpz_sub(Y, bp + 0, bp + st);       /* b(-1)  */
                fmpz_mul(A, A, Y);
                fmpz_add(Y, Y, bp + st);
                fmpz_add(Y, Y, bp + st);            /* b(1)   */
                fmpz_mul(B, B, Y);
            }
            else
            {
                fmpz_mul_ui(X, bp + st, (ulong) k);
                fmpz_sub(Y, bp + 0, X);             /* b(-k)  */
                fmpz_mul(A, A, Y);
                fmpz_add(Y, Y, X);
                fmpz_add(Y, Y, X);                  /* b(k)   */
                fmpz_mul(B, B, Y);
            }
        }
        else if (k == 1 && lb == 3)
        {
            toom_eval_even(Y, bp, lb, 1, st);
            fmpz_sub(Y, Y, bp + st);                /* b(-1)  */
            fmpz_mul(A, A, Y);
            fmpz_add(Y, Y, bp + st);
            fmpz_add(Y, Y, bp + st);                /* b(1)   */
            fmpz_mul(B, B, Y);
        }
        else
        {
            toom_eval_even(Y, bp, lb, (ulong) k, st);
            toom_eval_odd(X, bp, lb, (ulong) k, st);
            fmpz_sub(Y, Y, X);                      /* b(-k)  */
            fmpz_mul(A, A, Y);                      /* v_{-k} */
            fmpz_add(Y, Y, X);
            fmpz_add(Y, Y, X);                      /* b(k)   */
            fmpz_mul(B, B, Y);                      /* v_k    */
        }

        /* fold: E_k -> B, J_k -> A */
        fmpz_sub(X, B, A);
        fmpz_add(B, B, A);
        fmpz_tdiv_q_2exp(B, B, 1);
        fmpz_tdiv_q_2exp(A, X, 1 + flint_ctz((ulong) k));
    }

    if (extra)
    {
        toom_eval(X, ap, la, q, st);
        if (sq)
            fmpz_mul(c + N - 2, X, X);
        else
        {
            toom_eval(Y, bp, lb, q, st);
            fmpz_mul(c + N - 2, X, Y);              /* v_q */
        }
    }

    if (sq)
    {
        fmpz_mul(c + 0,     ap + 0,            ap + 0);
        fmpz_mul(c + N - 1, ap + (la - 1)*st,  ap + (la - 1)*st);
    }
    else
    {
        fmpz_mul(c + 0,     ap + 0,            bp + 0);
        fmpz_mul(c + N - 1, ap + (la - 1)*st,  bp + (lb - 1)*st);
    }

    /* ---- even half ------------------------------------------------ */

    for (k = 1; k <= p; k++)
    {
        flint_bitcnt_t r = 2 * flint_ctz((ulong) k);

        fmpz_sub(c + 2*k, c + 2*k, c + 0);
        if (extra)          /* N odd: top even index N-1 is known */
        {
            ulong pw = toom_pow_ui((ulong) k, (ulong) (N - 1));
            if (pw == 1)
                fmpz_sub(c + 2*k, c + 2*k, c + N - 1);
            else
                fmpz_submul_ui(c + 2*k, c + N - 1, pw);
        }
        if (r)
            fmpz_tdiv_q_2exp(c + 2*k, c + 2*k, r);
    }

    toom_interp(c + 2, p, 0, r0_e);

    /* ---- odd half ------------------------------------------------- */

    if (!extra)             /* N even: top odd index N-1 is known */
    {
        for (k = 1; k <= p; k++)
        {
            ulong pw = ((ulong) k >> flint_ctz((ulong) k))
                       * toom_pow_ui((ulong) k, (ulong) (N - 2));
            if (pw == 1)
                fmpz_sub(c + 2*k - 1, c + 2*k - 1, c + N - 1);
            else
                fmpz_submul_ui(c + 2*k - 1, c + N - 1, pw);
        }
        toom_interp(c + 1, p, 1, r0_o);
    }
    else                    /* fold the extra point into an odd sample */
    {
        flint_bitcnt_t r = flint_ctz(q);

        fmpz_sub(c + N - 2, c + N - 2, c + 0);
        for (i = 1; i <= (N - 1) / 2; i++)
            fmpz_submul_ui(c + N - 2, c + 2*i, toom_pow_ui(q, 2 * (ulong) i));
        if (r)
            fmpz_tdiv_q_2exp(c + N - 2, c + N - 2, r);   /* T_q */

        toom_interp(c + 1, p + 1, 1, r0_o);
    }
}

/* nmax = -1 selects the bit-size dependent limit */
int
_fmpz_poly_mul_toom_scalar_nmax(fmpz * res, const fmpz * poly1, slong len1,
                                const fmpz * poly2, slong len2, slong nmax)
{
    nmax = FLINT_MIN(nmax, TOOM_SCALAR_N_MAX);

    if (len1 + len2 - 1 > nmax)
        return 0;

    if (len2 == 1)
    {
        slong i;
        for (i = 0; i < len1; i++)
            fmpz_mul(res + i, poly1 + i, poly2 + 0);
        return 1;
    }

    if (len1 == 2)              /* len2 == 2 */
    {
        toom_mul_2_2(res, poly1, poly2, poly1 == poly2);
        return 1;
    }
    if (len1 == 3)
    {
        if (len2 == 2)
            toom_mul_3_2(res, poly1, poly2);
        else
            toom_mul_3_3(res, poly1, poly2, poly1 == poly2);
        return 1;
    }
    if (len1 == 4 && len2 == 4)
    {
        toom_mul_4_4(res, poly1, poly2, poly1 == poly2);
        return 1;
    }

    toom_mul_generic(res, poly1, len1, poly2, len2,
                     poly1 == poly2 && len1 == len2, 0, 0);
    return 1;
}

int
_fmpz_poly_mul_toom_scalar(fmpz * res, const fmpz * poly1, slong len1,
                           const fmpz * poly2, slong len2)
{
    return _fmpz_poly_mul_toom_scalar_nmax(res, poly1, len1, poly2, len2, TOOM_SCALAR_N_MAX);
}

int
fmpz_poly_mul_toom_scalar(fmpz_poly_t res, const fmpz_poly_t poly1,
                          const fmpz_poly_t poly2)
{
    slong len1 = poly1->length, len2 = poly2->length;
    slong N = len1 + len2 - 1;
    int ok;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return 1;
    }

    if (N > TOOM_SCALAR_N_MAX)
        return 0;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, N);
        ok = (len1 >= len2)
           ? _fmpz_poly_mul_toom_scalar(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2)
           : _fmpz_poly_mul_toom_scalar(t->coeffs, poly2->coeffs, len2,
                                        poly1->coeffs, len1);
        if (ok)
        {
            _fmpz_poly_set_length(t, N);
            fmpz_poly_swap(res, t);
        }
        fmpz_poly_clear(t);
        return ok;
    }

    fmpz_poly_fit_length(res, N);
    ok = (len1 >= len2)
       ? _fmpz_poly_mul_toom_scalar(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2)
       : _fmpz_poly_mul_toom_scalar(res->coeffs, poly2->coeffs, len2,
                                    poly1->coeffs, len1);
    if (ok)
        _fmpz_poly_set_length(res, N);
    return ok;
}

/*
    General middle products.

    Computes coefficients [nlo, nhi) of poly1 * poly2 with the minimal
    number of pointwise multiplications reachable by either of the two
    evaluation-style routes:

      full route:  trim unused coefficients of both operands, run the
                   full product engine on the trimmed shape
                   (m' + n' - 1 products), keep the window;
                                                                       
      mid route:   realise the window as the middle range of a padded
                   <w, min(m',n')> product and run that product
                   *transposed* (w + min(m',n') - 1 products).
                                                                       
    After trimming, (m', n', w) always satisfies the triangle
    conditions max <= sum of others - 1, and the combined cost is

         (sum of the two smallest of m', n', w) - 1,

    symmetric in all three parameters as transposition invariance
    requires.  This is provably optimal whenever the window touches an
    end of the product (Alder-Strassen for mullow/mulhigh windows,
    transposition of the full product for exact middle windows).  For
    doubly punctured windows -- lo > 0 and hi < m' + n' - 2
    simultaneously -- strictly smaller ranks can exist: the window
    [1,3] of a 3x3 product has rank 4 (Waring decomposition of
    x^3 + 6xyz) versus 5 here.  Such decompositions have inconvenient
    rational constants and are not used.

    The default mid-route engine transposes the forward pipeline of
    this file literally (Tellegen): reverse the op sequence of
    fold -> setup -> Newton and transpose each elementary op.  Exact
    divisions do not transpose (the divisibility argument lives on the
    forward side), so the Newton recurrences are rewritten division-
    free with per-row denominators delta_r = lcm(D_r .. D_{n-1}) and
    every power-of-two shift is absorbed into a static per-slot scale.
    All resulting constants depend only on N = w + min(m',n') - 1 and
    are precomputed into the tmid_* tables. The merged divisors
    form a divisibility chain.
*/

#define TMID_ADDMUL_SI(r, x, m)                       \
    do {                                              \
        slong m__ = (m);                              \
        if (m__ >= 0)                                 \
            fmpz_addmul_ui((r), (x), (ulong) m__);    \
        else                                          \
            fmpz_submul_ui((r), (x), (ulong) (-m__)); \
    } while (0)

/* y_t = Along[ybase + t], zero outside [0, L) */
#define TMID_Y(t) (((ybase + (t)) >= 0 && (ybase + (t)) < L) \
                   ? Along + (ybase + (t)) : NULL)

/* Brev(u) for Brev_j = B[S-1-j]; pair-shared evaluation. */
static void
tmid_eval_brev_pair(fmpz_t vp, fmpz_t vm, fmpz_t t,
                    const fmpz * B, slong S, ulong k)
{
    ulong k2 = k * k;
    slong j;

    {
        slong top = (S - 1) & ~(slong) 1;         /* largest even j < S */
        fmpz_set(vp, B + (S - 1 - top));
        for (j = top - 2; j >= 0; j -= 2)
        {
            fmpz_mul_ui(vp, vp, k2);
            fmpz_add(vp, vp, B + (S - 1 - j));
        }
    }
    if (S >= 2)
    {
        slong topo = ((S - 2) & ~(slong) 1) + 1;  /* largest odd j < S */
        fmpz_set(t, B + (S - 1 - topo));
        for (j = topo - 2; j >= 1; j -= 2)
        {
            fmpz_mul_ui(t, t, k2);
            fmpz_add(t, t, B + (S - 1 - j));
        }
        fmpz_mul_ui(t, t, k);                     /* k * O(k^2) */
    }
    else
        fmpz_zero(t);

    fmpz_sub(vm, vp, t);
    fmpz_add(vp, vp, t);
}

static void
tmid_eval_brev(fmpz_t v, const fmpz * B, slong S, ulong u)
{
    slong j;
    fmpz_set(v, B + 0);                           /* Brev_{S-1} */
    for (j = S - 2; j >= 0; j--)
    {
        fmpz_mul_ui(v, v, u);
        fmpz_add(v, v, B + (S - 1 - j));
    }
}

/* ------------------------------------------------------------------ */
/* transposed-Newton engine (table driven)                            */
/* ------------------------------------------------------------------ */

static void
toom_tinterp(fmpz * ys, slong base, slong n, int par)
{
    const ulong * Aq = tmid_A[par][n];
    slong r, i, j;

    if (n <= 1)
        return;

    for (j = 1; j < n; j++)
        fmpz_sub(ys + base + 2*j, ys + base + 2*j, ys + base);
    fmpz_mul_ui(ys + base, ys + base, tmid_Dtop[par][n]);

    for (r = 1; r < n; r++)
    {
        for (j = r + 1; j < n; j++)
            fmpz_submul_ui(ys + base + 2*j, ys + base + 2*r,
                           toom_h[r][j - r - 1]);
        for (i = 0; i < r; i++)
            TMID_ADDMUL_SI(ys + base + 2*i, ys + base + 2*r,
                           (slong) Aq[r] * toom_m[par][r][i]);
        fmpz_mul_ui(ys + base + 2*r, ys + base + 2*r,
                    Aq[r] * (ulong) toom_m[par][r][r]);
    }
}

static void
toom_mulmid_engine(fmpz * res, const fmpz * Along, slong L,
                   slong lo, const fmpz * B, slong S, slong w)
{
    slong N = w + S - 1;
    slong p = (N - 2) >> 1;
    int extra = N & 1;
    ulong q = (ulong) (p + 1);
    slong n_e = p, n_o = p + extra;
    ulong Go, gT, sc, wt, pw, qq;
    ulong pwk[10];
    slong ybase = lo - S + 1;
    slong t, k, i, M2;
    fmpz * ys, * T1, * T2, * T3;

    if (N == 2)                       /* S = 2, w = 1 */
    {
        fmpz_zero(res);
        if (TMID_Y(0) != NULL) fmpz_addmul(res, TMID_Y(0), B + 1);
        if (TMID_Y(1) != NULL) fmpz_addmul(res, TMID_Y(1), B + 0);
        return;
    }

    ys = _fmpz_vec_init(N + 3);
    T1 = ys + N; T2 = ys + N + 1; T3 = ys + N + 2;

    for (t = 0; t < N; t++)
        if (TMID_Y(t) != NULL)
            fmpz_set(ys + t, TMID_Y(t));

    Go = tmid_Dtop[1][n_o];
    (void) n_e;
    M2 = 0;
    for (k = 1; k <= p; k++)
        M2 = FLINT_MAX(M2, 2 * (slong) flint_ctz((ulong) k));

    /* ---- transposed odd Newton ------------------------------------ */

    toom_tinterp(ys, 1, n_o, 1);

    /* ---- transposed odd setup / T stage --------------------------- */

    if (extra)
    {
        gT = tmid_gT[N];
        if (gT != 1) fmpz_mul_ui(ys + 0, ys + 0, gT);
        fmpz_sub(ys + 0, ys + 0, ys + N - 2);
        qq = q * q;
        pw = 1;
        for (i = 1; i <= p; i++)
        {
            pw *= qq;                              /* q^{2i} */
            if (gT != 1) fmpz_mul_ui(ys + 2*i, ys + 2*i, gT);
            fmpz_submul_ui(ys + 2*i, ys + N - 2, pw);
        }
        if (gT != 1) fmpz_mul_ui(ys + N - 1, ys + N - 1, gT);
        fmpz_submul_ui(ys + N - 1, ys + N - 2, tmid_qN[N]);
    }
    else
    {
        if (Go != 1) fmpz_mul_ui(ys + N - 1, ys + N - 1, Go);
        for (k = 1; k <= p; k++)
            fmpz_submul_ui(ys + N - 1, ys + 2*k - 1, tmid_pwodd[N][k]);
    }

    /* ---- transposed even Newton ----------------------------------- */

    toom_tinterp(ys, 2, n_e, 0);

    /* ---- transposed even setup ------------------------------------ */

    sc = tmid_sc0[N];                 /* Ge << M2 */
    if (sc != 1) fmpz_mul_ui(ys + 0, ys + 0, sc);
    for (k = 1; k <= p; k++)
    {
        wt = UWORD(1) << (M2 - 2 * flint_ctz((ulong) k));
        if (wt == 1)
            fmpz_sub(ys + 0, ys + 0, ys + 2*k);
        else
            fmpz_submul_ui(ys + 0, ys + 2*k, wt);
    }
    if (extra)
    {
        sc = tmid_scinf[N];           /* Ge */
        if (sc != 1) fmpz_mul_ui(ys + N - 1, ys + N - 1, sc);
        for (k = 1; k <= p; k++)
            fmpz_submul_ui(ys + N - 1, ys + 2*k, tmid_pwev[N][k]);
    }

    /* ---- transposed folds (fold correction factor is always 1) ---- */

    for (k = 1; k <= p; k++)
    {
        fmpz_mul_ui(T1, ys + 2*k - 1, tmid_fA[N][k]);
        fmpz_sub(ys + 2*k - 1, ys + 2*k, T1);
        fmpz_add(ys + 2*k,     ys + 2*k, T1);
    }

    /* ---- scale each slot by Gamma_tot / gamma_slot ----------------- */

    if ((wt = tmid_r0[N]) != 1)
        fmpz_mul_ui(ys + 0, ys + 0, wt);
    if ((wt = tmid_rinf[N]) != 1)
        fmpz_mul_ui(ys + N - 1, ys + N - 1, wt);
    if (extra && (wt = tmid_rT[N]) != 1)
        fmpz_mul_ui(ys + N - 2, ys + N - 2, wt);
    for (k = 1; k <= p; k++)
        if ((wt = tmid_rpair[N][k]) != 1)
        {
            fmpz_mul_ui(ys + 2*k - 1, ys + 2*k - 1, wt);
            fmpz_mul_ui(ys + 2*k,     ys + 2*k,     wt);
        }

    /* ---- pointwise: N products ------------------------------------ */

    for (k = 1; k <= p; k++)
    {
        tmid_eval_brev_pair(T1, T2, T3, B, S, (ulong) k);
        fmpz_mul(ys + 2*k,     ys + 2*k,     T1);    /* * Brev(+k) */
        fmpz_mul(ys + 2*k - 1, ys + 2*k - 1, T2);    /* * Brev(-k) */
    }
    fmpz_mul(ys + 0, ys + 0, B + (S - 1));           /* Brev(0)    */
    if (extra)
    {
        tmid_eval_brev(T1, B, S, q);
        fmpz_mul(ys + N - 2, ys + N - 2, T1);
    }
    fmpz_mul(ys + N - 1, ys + N - 1, B + 0);         /* Brev_{S-1} */

    /* ---- outputs --------------------------------------------------- */

    for (k = 1; k <= p; k++)                         /* sum -> 2k, diff -> 2k-1 */
    {
        fmpz_add(T1, ys + 2*k, ys + 2*k - 1);
        fmpz_sub(ys + 2*k - 1, ys + 2*k, ys + 2*k - 1);
        fmpz_swap(ys + 2*k, T1);
    }

    for (k = 1; k <= p; k++)
        pwk[k] = 1;
    pw = 1;                                          /* q^i */
    for (i = 0; i < w; i++)
    {
        fmpz * zi = res + i;
        if (i == 0)
        {
            fmpz_set(zi, ys + 0);
            for (k = 1; k <= p; k++)
                fmpz_add(zi, zi, ys + 2*k);
            if (extra)
                fmpz_add(zi, zi, ys + N - 2);
        }
        else
        {
            for (k = 1; k <= p; k++)
                pwk[k] *= (ulong) k;
            pw *= q;
            fmpz_zero(zi);
            for (k = 2; k <= p; k++)
                fmpz_addmul_ui(zi, ys + 2*k - (i & 1), pwk[k]);
            if (p >= 1)
                fmpz_add(zi, zi, ys + 2 - (i & 1));  /* k = 1 */
            if (extra)
                fmpz_addmul_ui(zi, ys + N - 2, pw);
        }
        if (i == w - 1)
            fmpz_add(zi, zi, ys + N - 1);
        if (tmid_G1[N] != 1)
            fmpz_divexact_ui(zi, zi, tmid_G1[N]);
        if (tmid_G2[N] != 1)
            fmpz_divexact_ui(zi, zi, tmid_G2[N]);
    }

    _fmpz_vec_clear(ys, N + 3);
}

int
_fmpz_poly_mulmid_toom_scalar_nmax(fmpz * res,
                              const fmpz * poly1, slong len1,
                              const fmpz * poly2, slong len2,
                              slong nlo, slong nhi, slong nmax)
{
    nmax = FLINT_MIN(nmax, TOOM_SCALAR_N_MAX);

    slong prodlen = len1 + len2 - 1;
    slong sa, ea, sb, eb, m1, n1, lo, hi, w, cf, cm, r;
    const fmpz *A, *B;
    slong L, S;

    /* the whole product is wanted: no trimming or windowing to do */
    if (nlo == 0 && nhi == prodlen)
        return (len1 >= len2)
             ? _fmpz_poly_mul_toom_scalar(res, poly1, len1, poly2, len2)
             : _fmpz_poly_mul_toom_scalar(res, poly2, len2, poly1, len1);

    /* trim */
    sa = FLINT_MAX(0, nlo - (len2 - 1));
    ea = FLINT_MIN(len1 - 1, nhi - 1);
    sb = FLINT_MAX(0, nlo - (len1 - 1));
    eb = FLINT_MIN(len2 - 1, nhi - 1);
    m1 = ea - sa + 1;
    n1 = eb - sb + 1;
    lo = nlo - sa - sb;
    hi = (nhi - 1) - sa - sb;
    w = hi - lo + 1;

    if (m1 >= n1)
    {
        A = poly1 + sa; L = m1;
        B = poly2 + sb; S = n1;
    }
    else
    {
        A = poly2 + sb; L = n1;
        B = poly1 + sa; S = m1;
    }

    if (S == 1)
    {
        for (r = 0; r < w; r++)
            fmpz_mul(res + r, A + lo + r, B + 0);
        return 1;
    }

    cf = m1 + n1 - 1;
    cm = w + S - 1;
    if (cf <= cm && cf <= nmax)
    {
        slong topgap = (cf - 1) - hi;
        int trunc = (cf >= 4 && FLINT_MAX(lo, topgap) > 0);
        int rv = trunc && (topgap > lo);
        slong slo = rv ? cf - 1 - hi : lo;   /* window -> [slo, slo+w) */
        fmpz * scratch;
        slong j;
        TMP_INIT;

        TMP_START;
        scratch = (fmpz *) TMP_ALLOC(cf * sizeof(fmpz));

        for (j = 0; j < cf; j++)
        {
            if (j >= slo && j < slo + w)
                scratch[j] = res[rv ? slo + w - 1 - j : j - slo];
            else
                fmpz_init(scratch + j);
        }

        if (trunc)
            toom_mul_generic(scratch, A, L, B, S, A == B && L == S,
                             rv ? topgap : lo, rv);
        else
            _fmpz_poly_mul_toom_scalar(scratch, A, L, B, S);

        for (r = 0; r < w; r++)
            res[r] = scratch[rv ? slo + w - 1 - r : slo + r];
        for (j = 0; j < cf; j++)
            if (j < slo || j >= slo + w)
                fmpz_clear(scratch + j);

        TMP_END;
        return 1;
    }

    if (cm <= nmax)
    {
        toom_mulmid_engine(res, A, L, lo, B, S, w);
        return 1;
    }

    if (cf <= nmax)
    {
        fmpz * scratch;
        slong j;
        TMP_INIT;

        TMP_START;
        scratch = (fmpz *) TMP_ALLOC(cf * sizeof(fmpz));
        for (j = 0; j < cf; j++)
        {
            if (j >= lo && j < lo + w)
                scratch[j] = res[j - lo];
            else
                fmpz_init(scratch + j);
        }
        _fmpz_poly_mul_toom_scalar(scratch, A, L, B, S);
        for (r = 0; r < w; r++)
            res[r] = scratch[lo + r];
        for (j = 0; j < cf; j++)
            if (j < lo || j >= lo + w)
                fmpz_clear(scratch + j);
        TMP_END;
        return 1;
    }

    return 0;
}

int
_fmpz_poly_mulmid_toom_scalar(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    return _fmpz_poly_mulmid_toom_scalar_nmax(res, poly1, len1, poly2, len2, nlo, nhi, TOOM_SCALAR_N_MAX);
}

int
fmpz_poly_mulmid_toom_scalar(fmpz_poly_t res, const fmpz_poly_t poly1,
                             const fmpz_poly_t poly2, slong nlo, slong nhi)
{
    slong len1 = poly1->length, len2 = poly2->length;
    slong w = nhi - nlo;
    int ok;

    if (len1 == 0 || len2 == 0 || w <= 0)
    {
        fmpz_poly_zero(res);
        return 1;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;

        fmpz_poly_init2(t, w);
        ok = _fmpz_poly_mulmid_toom_scalar(t->coeffs, poly1->coeffs, len1,
                                           poly2->coeffs, len2, nlo, nhi);
        if (ok)
        {
            _fmpz_poly_set_length(t, w);
            _fmpz_poly_normalise(t);
            fmpz_poly_swap(res, t);
        }

        fmpz_poly_clear(t);
        return ok;
    }

    fmpz_poly_fit_length(res, w);
    ok = _fmpz_poly_mulmid_toom_scalar(res->coeffs, poly1->coeffs, len1,
                                       poly2->coeffs, len2, nlo, nhi);
    if (ok)
    {
        _fmpz_poly_set_length(res, w);
        _fmpz_poly_normalise(res);
    }

    return ok;
}

