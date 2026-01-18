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
    slong decimal_limbs;
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
    alloc = radix_set_mpn_need_alloc(xn, radix);
    t = TMP_ALLOC(alloc * sizeof(ulong));
    decimal_limbs = radix_set_mpn(t, x, xn, radix);
    res = radix_get_str_decimal(res, t, decimal_limbs, negative, radix);
    TMP_END;
    return res;
}

