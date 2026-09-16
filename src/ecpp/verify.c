/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "ecpp.h"

/*
    Checks one step of a certificate. Returns 1 if the step is valid, i.e.
    n is prime provided q is prime; 0 otherwise.
*/
int
ecpp_verify_step(const ecpp_step_struct * s)
{
    int result = 0;
    fmpz_t t, u, k, acc;
    fmpz_mod_ctx_t ctx;
    ecpp_point_t P, R;

    if (fmpz_cmp_ui(s->n, 3) <= 0 || fmpz_is_even(s->n))
        return 0;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(k);
    fmpz_init(acc);
    fmpz_mod_ctx_init(ctx, s->n);
    ecpp_point_init(P);
    ecpp_point_init(R);

    /* q | m, m = k q with k > 1 */
    if (fmpz_cmp_ui(s->q, 1) <= 0 || !fmpz_divisible(s->m, s->q))
        goto cleanup;
    fmpz_divexact(k, s->m, s->q);
    if (fmpz_cmp_ui(k, 1) <= 0)
        goto cleanup;

    /* q > (n^{1/4} + 1)^2, i.e. (sqrt(q) - 1)^4 > n; check via floor roots */
    fmpz_sqrt(t, s->q);                     /* t = floor(sqrt(q)) */
    fmpz_sub_ui(t, t, 1);
    if (fmpz_sgn(t) <= 0)
        goto cleanup;
    fmpz_pow_ui(u, t, 4);                   /* (floor(sqrt q) - 1)^4 */
    /* need (sqrt(q) - 1)^4 > n; sqrt(q) - 1 >= floor(sqrt(q)) - 1, so this suffices */
    if (fmpz_cmp(u, s->n) <= 0)
        goto cleanup;

    /* gcd(6, n) = 1 and the curve is nonsingular: gcd(4a^3 + 27b^2, n) = 1 */
    if (fmpz_divisible_si(s->n, 3))
        goto cleanup;
    fmpz_mod_mul(t, s->a, s->a, ctx);
    fmpz_mod_mul(t, t, s->a, ctx);
    fmpz_mod_mul_ui(t, t, 4, ctx);
    fmpz_mod_mul(u, s->b, s->b, ctx);
    fmpz_mod_mul_ui(u, u, 27, ctx);
    fmpz_mod_add(t, t, u, ctx);
    fmpz_gcd(t, t, s->n);
    if (!fmpz_is_one(t))
        goto cleanup;

    /* P on the curve */
    fmpz_mod_mul(t, s->x, s->x, ctx);
    fmpz_mod_add(t, t, s->a, ctx);
    fmpz_mod_mul(t, t, s->x, ctx);
    fmpz_mod_add(t, t, s->b, ctx);
    fmpz_mod_mul(u, s->y, s->y, ctx);
    if (!fmpz_equal(t, u))
        goto cleanup;

    ecpp_point_set_affine(P, s->x, s->y);
    fmpz_one(acc);

    /* m P = O */
    ecpp_point_mul(R, P, s->m, s->a, acc, ctx);
    if (!ecpp_point_is_zero(R))
        goto cleanup;

    /* (m / q) P != O */
    ecpp_point_mul(R, P, k, s->a, acc, ctx);
    if (ecpp_point_is_zero(R))
        goto cleanup;

    /* all case distinctions were legitimate */
    fmpz_gcd(t, acc, s->n);
    if (!fmpz_is_one(t))
        goto cleanup;

    result = 1;

cleanup:
    ecpp_point_clear(P);
    ecpp_point_clear(R);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(k);
    fmpz_clear(acc);

    return result;
}

/*
    Returns 1 if the certificate proves that n is prime, 0 otherwise. The
    chain must start at n and be consecutive (n_{i+1} = q_i), and end with
    a q < 2^64 which is checked directly. For n < 2^64 the certificate is
    empty and n is checked directly.
*/
int
ecpp_verify(const ecpp_cert_t cert, const fmpz_t n)
{
    slong i;
    const fmpz * last;

    if (cert->num == 0)
        return fmpz_sgn(n) > 0 && fmpz_bits(n) <= 64 && fmpz_is_prime(n);

    if (!fmpz_equal(cert->steps[0].n, n))
        return 0;

    for (i = 0; i < cert->num; i++)
    {
        const ecpp_step_struct * s = cert->steps + i;

        if (i > 0 && !fmpz_equal(s->n, cert->steps[i - 1].q))
            return 0;

        if (!ecpp_verify_step(s))
            return 0;
    }

    last = cert->steps[cert->num - 1].q;
    return fmpz_sgn(last) > 0 && fmpz_bits(last) <= 64 && fmpz_is_prime(last);
}
