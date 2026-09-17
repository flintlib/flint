/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "mpn_mod.h"

static void
_mpn_mod_modulus_fmpz(fmpz_t p, gr_ctx_t ctx)
{
    fmpz_set_ui_array(p, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
}

/* Euler's criterion: x^((p-1)/2) is 1 for squares and -1 for nonsquares */
truth_t
mpn_mod_is_square(nn_srcptr x, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    fmpz_t p, e;
    nn_ptr t;
    truth_t res;

    if (mpn_mod_is_zero(x, ctx) == T_TRUE || mpn_mod_is_one(x, ctx) == T_TRUE)
        return T_TRUE;

    if (MPN_MOD_CTX_IS_PRIME(ctx) != T_TRUE)
        return T_UNKNOWN;

    fmpz_init(p);
    fmpz_init(e);
    _mpn_mod_modulus_fmpz(p, ctx);
    fmpz_sub_ui(e, p, 1);
    fmpz_tdiv_q_2exp(e, e, 1);

    t = flint_malloc(n * sizeof(ulong));

    if (mpn_mod_pow_fmpz(t, x, e, ctx) != GR_SUCCESS)
        res = T_UNKNOWN;
    else
        res = (mpn_mod_is_one(t, ctx) == T_TRUE) ? T_TRUE : T_FALSE;

    flint_free(t);
    fmpz_clear(p);
    fmpz_clear(e);

    return res;
}

/*
    Tonelli and Shanks. Writing p - 1 = q 2^s with q odd, the case s = 1 is
    a single exponentiation, and otherwise we walk down the 2-part of the
    group using a fixed quadratic nonresidue.
*/
int
mpn_mod_sqrt(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    fmpz_t p, q, e;
    nn_ptr scratch, c, r, t, b, z, u;
    slong s, i, j, m;
    int status = GR_SUCCESS;

    if (mpn_mod_is_zero(x, ctx) == T_TRUE || mpn_mod_is_one(x, ctx) == T_TRUE)
        return mpn_mod_set(res, x, ctx);

    if (MPN_MOD_CTX_IS_PRIME(ctx) != T_TRUE)
        return GR_UNABLE;

    if (mpn_mod_is_square(x, ctx) != T_TRUE)
    {
        mpn_mod_zero(res, ctx);
        return GR_DOMAIN;
    }

    fmpz_init(p);
    fmpz_init(q);
    fmpz_init(e);

    _mpn_mod_modulus_fmpz(p, ctx);

    /* p - 1 = q 2^s, q odd */
    fmpz_sub_ui(q, p, 1);
    s = fmpz_val2(q);
    fmpz_tdiv_q_2exp(q, q, s);

    scratch = flint_malloc(6 * n * sizeof(ulong));
    c = scratch; r = c + n; t = r + n; b = t + n; z = b + n; u = z + n;

    if (s == 1)
    {
        /* p = 3 mod 4, so the root is x^((p+1)/4) */
        fmpz_add_ui(e, p, 1);
        fmpz_tdiv_q_2exp(e, e, 2);
        status = mpn_mod_pow_fmpz(res, x, e, ctx);
        goto cleanup;
    }

    /* the smallest quadratic nonresidue; there is one well within reach */
    for (i = 2; ; i++)
    {
        status |= mpn_mod_set_ui(z, (ulong) i, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        if (mpn_mod_is_square(z, ctx) == T_FALSE)
            break;

        if (i == WORD(1) << 20)
        {
            /* only reachable if the modulus is not actually prime */
            status = GR_UNABLE;
            goto cleanup;
        }
    }

    status |= mpn_mod_pow_fmpz(c, z, q, ctx);       /* c = z^q */
    status |= mpn_mod_pow_fmpz(t, x, q, ctx);       /* t = x^q */

    fmpz_add_ui(e, q, 1);
    fmpz_tdiv_q_2exp(e, e, 1);
    status |= mpn_mod_pow_fmpz(r, x, e, ctx);       /* r = x^((q+1)/2) */

    if (status != GR_SUCCESS)
        goto cleanup;

    m = s;

    while (mpn_mod_is_one(t, ctx) != T_TRUE)
    {
        /* the least i with 0 < i < m and t^(2^i) = 1 */
        status |= mpn_mod_set(u, t, ctx);

        for (i = 1; i < m; i++)
        {
            status |= mpn_mod_sqr(u, u, ctx);

            if (mpn_mod_is_one(u, ctx) == T_TRUE)
                break;
        }

        if (status != GR_SUCCESS)
            goto cleanup;

        if (i >= m)
        {
            /* again, only reachable for a modulus that is not prime */
            status = GR_UNABLE;
            goto cleanup;
        }

        /* b = c^(2^(m - i - 1)) */
        status |= mpn_mod_set(b, c, ctx);

        for (j = 0; j < m - i - 1; j++)
            status |= mpn_mod_sqr(b, b, ctx);

        m = i;
        status |= mpn_mod_sqr(c, b, ctx);
        status |= mpn_mod_mul(t, t, c, ctx);
        status |= mpn_mod_mul(r, r, b, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;
    }

    status |= mpn_mod_set(res, r, ctx);

cleanup:
    if (status != GR_SUCCESS)
        mpn_mod_zero(res, ctx);

    flint_free(scratch);
    fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(e);

    return status;
}
