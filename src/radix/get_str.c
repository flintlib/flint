/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdint.h>
#include "mpn_extras.h"
#include "radix.h"

static const uint8_t n_sizeinbase10_tab1[64] = {
    1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8,
    8, 8, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13,
    14, 14, 14, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19,
    19, 19, 20,
};

static const ulong n_sizeinbase10_tab2[FLINT_BITS] = {
    1, 1, 1, 10, 10, 10, 100, 100, 100, 1000, 1000, 1000, 1000, 10000, 10000,
    10000, 100000, 100000, 100000, 1000000, 1000000, 1000000, 1000000, 10000000,
    10000000, 10000000, 100000000, 100000000, 100000000, 1000000000,
    1000000000, 1000000000,
#if FLINT_BITS == 64
    1000000000, UWORD(10000000000), UWORD(10000000000), UWORD(10000000000),
    UWORD(100000000000), UWORD(100000000000), UWORD(100000000000),
    UWORD(1000000000000), UWORD(1000000000000), UWORD(1000000000000),
    UWORD(1000000000000), UWORD(10000000000000), UWORD(10000000000000),
    UWORD(10000000000000), UWORD(100000000000000), UWORD(100000000000000),
    UWORD(100000000000000), UWORD(1000000000000000), UWORD(1000000000000000),
    UWORD(1000000000000000), UWORD(1000000000000000), UWORD(10000000000000000),
    UWORD(10000000000000000), UWORD(10000000000000000),
    UWORD(100000000000000000), UWORD(100000000000000000),
    UWORD(100000000000000000),  UWORD(1000000000000000000),
    UWORD(1000000000000000000), UWORD(1000000000000000000),
    UWORD(1000000000000000000), UWORD(10000000000000000000),
#endif
};

static const char dec_to_str_tab[200] =
    "000102030405060708091011121314151617181920212223242526272829"
    "303132333435363738394041424344454647484950515253545556575859"
    "606162636465666768697071727374757677787980818283848586878889"
    "90919293949596979899";

FLINT_FORCE_INLINE void n_to_str_2(char * s, uint32_t x)
{
    memcpy(s, &dec_to_str_tab[2 * x], 2);
}

FLINT_FORCE_INLINE void n_to_str_3(char * s, uint32_t x)
{
    n_to_str_2(s, x / 10);
    s[2] = '0' + x % 10;
}

FLINT_FORCE_INLINE void n_to_str_4(char * s, uint32_t x)
{
    n_to_str_2(s, x / 100);
    n_to_str_2(s + 2, x % 100);
}

FLINT_FORCE_INLINE void n_to_str_8(char * s, uint32_t x)
{
    n_to_str_4(s, x / 10000);
    n_to_str_4(s + 4, x % 10000);
}

static void
n_max_decimal_limb_get_str(char * s, ulong x)
{
#if FLINT_BITS == 64
    FLINT_ASSERT(x < UWORD(10000000000000000000));

    ulong q = x / UWORD(10000000000000000);
    ulong r = x % UWORD(10000000000000000);

    n_to_str_3(s, q);
    n_to_str_8(s + 3, r / 100000000);
    n_to_str_8(s + 11, r % 100000000);
#else
    FLINT_ASSERT(x < 1000000000);

    s[0] = '0' + (x / 100000000);
    n_to_str_8(s + 1, x % 100000000);
#endif
}

static void
n_get_str_nd(char * s, ulong x, int d)
{
    int i;

    for (i = d - 1; i >= 0; i--)
    {
        s[i] = (x % 10) + '0';
        x /= 10;
    }
}

static slong
n_nonzero_sizeinbase10(ulong n)
{
    FLINT_ASSERT(n != 0);
    int b = FLINT_BIT_COUNT(n) - 1;
    return n_sizeinbase10_tab1[b] - (n < n_sizeinbase10_tab2[b]);
}

static char * _radix_decimal_get_str(char * res, nn_srcptr t, slong decimal_limbs, int negative, slong digits_per_limb)
{
    slong i;

    if (decimal_limbs == 0)
    {
        if (res == NULL)
            res = flint_malloc(2);
        res[0] = '0';
        res[1] = '\0';
    }
    else
    {
        slong nd, nleading = n_nonzero_sizeinbase10(t[decimal_limbs - 1]);

        nd = (decimal_limbs - 1) * digits_per_limb + nleading;

        if (res == NULL)
            res = flint_malloc(nd + negative + 1);

        if (negative)
            res[0] = '-';

        n_get_str_nd(res + negative, t[decimal_limbs - 1], nleading);

        if (digits_per_limb == ((FLINT_BITS == 64) ? 19 : 9))
        {
            for (i = 1; i < decimal_limbs; i++)
                n_max_decimal_limb_get_str(res + negative + nleading + (i - 1) * digits_per_limb, t[decimal_limbs - 1 - i]);
        }
        else
        {
            /* todo: fast code */
            for (i = 1; i < decimal_limbs; i++)
                n_get_str_nd(res + negative + nleading + (i - 1) * digits_per_limb, t[decimal_limbs - 1 - i], digits_per_limb);
        }
        res[negative + nd] = '\0';
    }

    return res;
}

char * radix_get_str_decimal(char * res, nn_srcptr x, slong n, int negative, const radix_t radix)
{
    if (DIGIT_RADIX(radix) == 10)
    {
        return _radix_decimal_get_str(res, x, n, negative, radix->exp);
    }
    else if (n == 1 && LIMB_RADIX(radix) < n_sizeinbase10_tab2[FLINT_BITS - 1])  /* todo: could work even for 10/20-digit input */
    {
        return _radix_decimal_get_str(res, x, 1, negative, ((FLINT_BITS == 64) ? 19 : 9));
    }
    else
    {
        /* todo: improve */
        nn_ptr t;
        slong tn;
        TMP_INIT;
        TMP_START;
        t = TMP_ALLOC(n * sizeof(ulong));
        tn = radix_get_mpn(t, x, n, radix);
        res = flint_mpn_get_str(res, 10, t, tn, negative);
        TMP_END;
        return res;
    }
}

/* descending => -(D[n-1] * b^(e*(n-1)) + ... + D[1] * b^e + D[0]) */
/* ascending  => -(D[0] + D[1] * b^e + ... + D[n-1] * b^b(e*(n-1))) */
char * radix_get_str_sum(char * res, nn_srcptr x, slong n, int negative, int ascending, const radix_t radix)
{
    if (n <= 1)
        return radix_get_str_decimal(res, x, n, negative, radix);

    /* " * b^" */
    char power_base_str[32];
    slong power_base_str_len;

    /* length of b */
    slong base_len = n_nonzero_sizeinbase10(DIGIT_RADIX(radix));

    memcpy(power_base_str, " * ", 3);
    n_get_str_nd(power_base_str + 3, DIGIT_RADIX(radix), base_len);
    power_base_str[3 + base_len] = '^';
    power_base_str_len = 3 + base_len + 1;

    slong alloc = WORD_MAX;

    if (res == NULL)
    {
        /* max length of D[i] */
        slong limb_max_len = n_nonzero_sizeinbase10(LIMB_RADIX(radix) - 1);
        /* length of largest exponent */
        slong exp_len = n_nonzero_sizeinbase10((n - 1) * radix->exp);

        /* space for 0 or -(), and null terminator */
        alloc = 4;
        /* space for limbs */
        alloc += n * limb_max_len;
        /* space for " + " and " * b^e" */
        alloc += (n - 1) * (3 + power_base_str_len + exp_len);

        res = flint_malloc(alloc);
    }

    slong rlen, i, j;

    rlen = 0;
    if (negative)
    {
        memcpy(res, "-(", 2);
        rlen = 2;
    }

    for (j = 0; j < n; j++)
    {
        i = ascending ? j : n - 1 - j;

        ulong d = x[i];

        if (d == 0)
        {
            res[rlen] = '0';
            rlen++;
        }
        else
        {
            slong l = n_nonzero_sizeinbase10(d);
            n_get_str_nd(res + rlen, d, l);
            rlen += l;
        }

        if (i != 0)
        {
            memcpy(res + rlen, power_base_str, power_base_str_len);
            rlen += power_base_str_len;

            ulong e = i * radix->exp;
            slong l = n_nonzero_sizeinbase10(e);
            n_get_str_nd(res + rlen, e, l);
            rlen += l;
        }

        if (j < n - 1)
        {
            memcpy(res + rlen, " + ", 3);
            rlen += 3;
        }
    }

    if (negative)
    {
        res[rlen] = ')';
        rlen++;
    }

    res[rlen] = '\0';

    FLINT_ASSERT(rlen + 1 <= alloc);

    return res;
}

