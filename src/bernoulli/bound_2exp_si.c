/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bernoulli.h"

const short bernoulli_bound_tab[256] = {
    1, -2, -4, -5, -4, -3, -1, 1, 3, 6, 10, 13,
    17, 21, 25, 30, 34, 39, 44, 49, 55, 60, 66, 71,
    77, 83, 89, 95, 102, 108, 115, 121, 128, 135, 141, 148,
    155, 162, 170, 177, 184, 192, 199, 207, 214, 222, 230, 237,
    245, 253, 261, 269, 277, 285, 294, 302, 310, 318, 327, 335,
    344, 352, 361, 370, 378, 387, 396, 405, 413, 422, 431, 440,
    449, 458, 467, 477, 486, 495, 504, 514, 523, 532, 542, 551,
    561, 570, 580, 589, 599, 608, 618, 628, 638, 647, 657, 667,
    677, 687, 697, 707, 717, 727, 737, 747, 757, 767, 777, 787,
    797, 808, 818, 828, 838, 849, 859, 870, 880, 890, 901, 911,
    922, 932, 943, 953, 964, 975, 985, 996, 1007, 1017, 1028, 1039,
    1050, 1061, 1071, 1082, 1093, 1104, 1115, 1126, 1137, 1148, 1159, 1170,
    1181, 1192, 1203, 1214, 1225, 1236, 1247, 1258, 1270, 1281, 1292, 1303,
    1315, 1326, 1337, 1349, 1360, 1371, 1383, 1394, 1405, 1417, 1428, 1440,
    1451, 1463, 1474, 1486, 1497, 1509, 1520, 1532, 1544, 1555, 1567, 1579,
    1590, 1602, 1614, 1625, 1637, 1649, 1661, 1672, 1684, 1696, 1708, 1720,
    1732, 1743, 1755, 1767, 1779, 1791, 1803, 1815, 1827, 1839, 1851, 1863,
    1875, 1887, 1899, 1911, 1923, 1935, 1948, 1960, 1972, 1984, 1996, 2008,
    2021, 2033, 2045, 2057, 2070, 2082, 2094, 2106, 2119, 2131, 2143, 2156,
    2168, 2180, 2193, 2205, 2218, 2230, 2243, 2255, 2267, 2280, 2292, 2305,
    2317, 2330, 2342, 2355, 2368, 2380, 2393, 2405, 2418, 2430, 2443, 2456,
    2468, 2481, 2494, 2506,
};

#define LOG_PREC 6

/* table of ceil(log(n,2) * M - log_prec * M) for n from M to 2M inclusive
   where M = 2^LOG_PREC */
const unsigned char log_tab[] = {
    0, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23,
    25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
    42, 43, 44, 44, 45, 46, 47, 48, 49, 50, 51, 51, 52, 53, 54, 55, 55,
    56, 57, 58, 59, 59, 60, 61, 62, 62, 63, 64, 64
};

slong
bernoulli_bound_2exp_si(ulong n)
{
    if (n % 2)
    {
        if (n == 1)
            return -WORD(1);
        else
            return LONG_MIN;
    }
    else if (n < 512)
    {
        return bernoulli_bound_tab[n / 2];
    }
    else
    {
        /* |B_n| < 4 * n! / (2*pi)^n < 4 * (n+1)^(n+1) e^(-n) / (2*pi)^n */
        mp_limb_t l, u, hi, lo;
        int b, shift;

        b = FLINT_BIT_COUNT(n + 1);
        shift = b - (LOG_PREC + 1);

        /* (n+1) * log_2(n+1) */
        u = n + 1;
        l = log_tab[((u >> shift) + 1) - (1 << LOG_PREC)];
        l += (LOG_PREC << LOG_PREC);
        umul_ppmm(hi, lo, l, u);

        if (hi || n > (UWORD(1) << (FLINT_BITS - 6)))
        {
            flint_throw(FLINT_ERROR, "bernoulli_bound_2exp_si: n too large\n");
        }

        l = (lo >> LOG_PREC) + 1;
        l += shift * u;

        /* log_2(2*pi*e) > 131 / 32 */
        return l + 2 - (131 * n) / 32;
    }
}

