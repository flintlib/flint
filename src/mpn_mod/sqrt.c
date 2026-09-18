/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"
#include "mpn_mod.h"

static void
_mpn_mod_modulus_fmpz(fmpz_t p, gr_ctx_t ctx)
{
    fmpz_set_ui_array(p, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
}

/* A length-n mpn read as an mpz; GMP wants the top limb to be nonzero. */
static void
_mpn_mod_roinit(mpz_t z, nn_srcptr x, slong n)
{
    while (n > 0 && x[n - 1] == 0)
        n--;

    mpz_roinit_n(z, (mp_srcptr) x, n);
}

/* The Jacobi symbol (x/p); the modulus is assumed odd. */
static int
_mpn_mod_jacobi(nn_srcptr x, gr_ctx_t ctx)
{
    mpz_t xz, pz;

    _mpn_mod_roinit(xz, x, MPN_MOD_CTX_NLIMBS(ctx));
    _mpn_mod_roinit(pz, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));

    return mpz_jacobi(xz, pz);
}

/* The same for a small integer, which need not be reduced. */
static int
_mpn_mod_jacobi_ui(ulong x, gr_ctx_t ctx)
{
    mpz_t xz, pz;

    mpz_roinit_n(xz, (mp_srcptr) &x, x != 0);
    _mpn_mod_roinit(pz, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));

    return mpz_jacobi(xz, pz);
}

/*
    The Jacobi symbol, which GMP computes in quasi-linear time, rather than
    Euler's criterion: a symbol of -1 also settles the question for an odd
    modulus that is not prime, where the criterion says nothing.
*/
truth_t
mpn_mod_is_square(nn_srcptr x, gr_ctx_t ctx)
{
    if (mpn_mod_is_zero(x, ctx) == T_TRUE || mpn_mod_is_one(x, ctx) == T_TRUE)
        return T_TRUE;

    if (MPN_MOD_CTX_MODULUS(ctx)[0] & 1)
    {
        if (_mpn_mod_jacobi(x, ctx) == -1)
            return T_FALSE;

        if (MPN_MOD_CTX_IS_PRIME(ctx) == T_TRUE)
            return T_TRUE;
    }

    return T_UNKNOWN;
}

/*
    Writing p - 1 = q 2^s with q odd, s = 1 and s = 2 are closed forms and
    otherwise we walk down the 2-part of the group (Tonelli and Shanks)
    using the least quadratic nonresidue.
*/
int
mpn_mod_sqrt(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    ulong tmp[6 * MPN_MOD_MAX_LIMBS];
    nn_ptr c = tmp, r = c + n, t = r + n, b = t + n, z = b + n, u = z + n;
    fmpz_t p, q, e;
    slong s, i, j, m;
    ulong k;
    int status = GR_SUCCESS;

    if (mpn_mod_is_zero(x, ctx) == T_TRUE || mpn_mod_is_one(x, ctx) == T_TRUE)
        return mpn_mod_set(res, x, ctx);

    /* a Jacobi symbol of -1 rules out a root without knowing p to be prime */
    if (mpn_mod_is_square(x, ctx) == T_FALSE)
    {
        mpn_mod_zero(res, ctx);
        return GR_DOMAIN;
    }

    if (MPN_MOD_CTX_IS_PRIME(ctx) != T_TRUE)
        return GR_UNABLE;

    fmpz_init(p);
    fmpz_init(q);
    fmpz_init(e);

    _mpn_mod_modulus_fmpz(p, ctx);

    /* p - 1 = q 2^s, q odd */
    fmpz_sub_ui(q, p, 1);
    s = fmpz_val2(q);
    fmpz_tdiv_q_2exp(q, q, s);

    if (s == 1)
    {
        /* p = 3 mod 4, so the root is x^((p+1)/4) */
        fmpz_add_ui(e, p, 1);
        fmpz_tdiv_q_2exp(e, e, 2);
        status = mpn_mod_pow_fmpz(res, x, e, ctx);
        goto cleanup;
    }

    if (s == 2)
    {
        /*
            p = 5 mod 8. Atkin: with b = (2x)^((p-5)/8) and y = 2 x b^2,
            a root is x b (y - 1). This is one exponentiation, where the
            variant in fmpz_sqrtmod needs a second one half of the time.
        */
        fmpz_sub_ui(e, p, 5);
        fmpz_tdiv_q_2exp(e, e, 3);

        status |= mpn_mod_add(t, x, x, ctx);            /* t = 2x */
        status |= mpn_mod_pow_fmpz(b, t, e, ctx);
        status |= mpn_mod_sqr(u, b, ctx);
        status |= mpn_mod_mul(u, u, t, ctx);            /* u = y = 2 x b^2 */
        status |= mpn_mod_sub_ui(u, u, 1, ctx);
        status |= mpn_mod_mul(b, b, x, ctx);
        status |= mpn_mod_mul(res, b, u, ctx);
        goto cleanup;
    }

    /*
        x^((q-1)/2) yields both x^((q+1)/2) and x^q for a multiplication
        each, so the descent is set up with two exponentiations, not three.
    */
    fmpz_sub_ui(e, q, 1);
    fmpz_tdiv_q_2exp(e, e, 1);

    status |= mpn_mod_pow_fmpz(u, x, e, ctx);           /* u = x^((q-1)/2) */
    status |= mpn_mod_mul(r, u, x, ctx);                /* r = x^((q+1)/2) */
    status |= mpn_mod_mul(t, r, u, ctx);                /* t = x^q */

    if (status != GR_SUCCESS)
        goto cleanup;

    /*
        The least quadratic nonresidue, which is well within reach. Only
        odd candidates need testing: 2 is a residue for p = 1 mod 8, and
        then so is any even number below the first odd nonresidue.
    */
    for (k = 3; ; k += 2)
    {
        if (_mpn_mod_jacobi_ui(k, ctx) == -1)
            break;

        if (k >= (UWORD(1) << 20))
        {
            /* only reachable if the modulus is not actually prime */
            status = GR_UNABLE;
            goto cleanup;
        }
    }

    status |= mpn_mod_set_ui(z, k, ctx);
    status |= mpn_mod_pow_fmpz(c, z, q, ctx);           /* c = z^q */

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

    fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(e);

    return status;
}
