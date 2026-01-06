/*
    Copyright (C) 2024, 2026 Fredrik Johansson

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

FLINT_TLS_PREFIX radix_t default_decimal_ctx;
FLINT_TLS_PREFIX int default_decimal_ctx_initialized = 0;

static void
default_decimal_ctx_cleanup(void)
{
    if (default_decimal_ctx_initialized)
    {
        default_decimal_ctx_initialized = 0;
        radix_clear(default_decimal_ctx);
    }
}

radix_struct * get_default_decimal_ctx(void);

radix_struct * get_default_decimal_ctx(void)
{
    if (!default_decimal_ctx_initialized)
    {
        radix_init(default_decimal_ctx, 10, 0);
        flint_register_cleanup_function(default_decimal_ctx_cleanup);
        default_decimal_ctx_initialized = 1;
    }

    return default_decimal_ctx;
}

char * flint_mpn_get_str(char * res, int base, mp_srcptr x, mp_size_t xn, int negative)
{
    radix_struct * radix;
    nn_ptr t;
    slong alloc;
    slong decimal_limbs, i;
    slong digits_per_limb;
    TMP_INIT;

#if FLINT_HAVE_FFT_SMALL
    /* Todo: improve the radix code so that it beats GMP in this range too. */
    if ((xn > 40 && xn < 1400) || base != 10)
#endif
    {
        mpz_t tmp;
        while (xn > 0 && x[xn - 1] == 0)
            xn--;
        tmp->_mp_d = (mp_ptr) x;
        tmp->_mp_alloc = tmp->_mp_size = negative ? -xn : xn;
        if (res == NULL)
            res = flint_malloc(mpz_sizeinbase(tmp, base) + 2);
        return mpz_get_str(res, base, tmp);
    }

    TMP_START;

    radix = get_default_decimal_ctx();
    digits_per_limb = radix->exp;
    alloc = radix_set_mpn_need_alloc(xn, radix);

    t = TMP_ALLOC(alloc * sizeof(ulong));
    decimal_limbs = radix_set_mpn(t, x, xn, radix);

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
        for (i = 1; i < decimal_limbs; i++)
            n_max_decimal_limb_get_str(res + negative + nleading + (i - 1) * digits_per_limb, t[decimal_limbs - 1 - i]);

        res[negative + nd] = '\0';
    }

    TMP_END;
    return res;
}

