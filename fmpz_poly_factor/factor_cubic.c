/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"


/*
    z = sqrt(x) clobber x
        return absolute precision (>= 0) of z if z exists, otherwise
        return negative if x is not a square
*/
static slong binary_sqrt(fmpz_t z, fmpz_t x, slong p)
{
    slong e, new_p, n;
    fmpz_t t, tx, s;
    fmpz two = 2;
#if FLINT_WANT_ASSERT
    fmpz_t x_org;
#endif

    FLINT_ASSERT(p > 0);
    fmpz_fdiv_r_2exp(x, x, p);

    if (fmpz_is_zero(x))
    {
        fmpz_zero(z);
        return p/2;
    }

    e = fmpz_remove(x, x, &two);

    if ((e % 2) != 0)
    {
        fmpz_zero(z);
        return -WORD(1);
    }

    /* new precision after power of 2 removal */
    new_p = p - e;
    if (new_p < 4)
    {
        fmpz_one(z);
        fmpz_mul_2exp(z, z, e/2);
        return e/2 + 1;
    }

    if (fmpz_fdiv_ui(x, 8) != 1)
    {
        fmpz_zero(z);
        return -WORD(1);
    }

#if FLINT_WANT_ASSERT
    fmpz_init_set(x_org, x);
    fmpz_mul_2exp(x_org, x_org, e);
#endif

    fmpz_init(t);
    fmpz_init(tx);
    fmpz_init(s);

    /*
        n is the current precision of z = 1/sqrt(x) as measured by
            2^n divides x*z^2 - 1
    */
    n = 4;

    fmpz_fdiv_r_2exp(z, x, n + 1);
    fmpz_sub_ui(z, z, 3);
    fmpz_neg(z, z);
    fmpz_fdiv_r_2exp(z, z, n + 1);
    fmpz_fdiv_q_2exp(z, z, 1);

#if FLINT_WANT_ASSERT
    fmpz_mul(t, z, z);
    fmpz_mul(t, t, x);
    fmpz_sub_ui(t, t, 1);
    fmpz_fdiv_r_2exp(t, t, n);
    FLINT_ASSERT(fmpz_is_zero(t));
#endif

    /* iterate z <- z*(3 - x*z^2)/2 */
    while (n < new_p - 1)
    {
        n = 2*n - 2;
        fmpz_mul(t, z, z);
        fmpz_fdiv_r_2exp(tx, x, n + 1);
        fmpz_mul(s, tx, t);
        fmpz_sub_ui(s, s, 3);
        fmpz_neg(s, s);
        fmpz_fdiv_r_2exp(s, s, n + 1);
        fmpz_fdiv_q_2exp(s, s, 1);
        fmpz_mul(t, z, s);
        fmpz_fdiv_r_2exp(t, t, n);
        fmpz_swap(z, t);

    #if FLINT_WANT_ASSERT
        fmpz_mul(t, z, z);
        fmpz_mul(t, t, x);
        fmpz_sub_ui(t, t, 1);
        fmpz_fdiv_r_2exp(t, t, n);
        FLINT_ASSERT(fmpz_is_zero(t));
    #endif
    }

    fmpz_mul(t, z, x);
    fmpz_fdiv_r_2exp(t, t, new_p - 1);
    fmpz_mul_2exp(t, t, e/2);
    fmpz_swap(z, t);

#if FLINT_WANT_ASSERT
    fmpz_submul(x_org, z, z);
    fmpz_fdiv_r_2exp(t, x_org, new_p + e/2);
    FLINT_ASSERT(fmpz_is_zero(t));
    fmpz_clear(x_org);
#endif

    fmpz_clear(t);
    fmpz_clear(tx);
    fmpz_clear(s);

    return new_p - 1 + e/2;
}


static mp_limb_t fmpz_fdiv_r_2exp_flint_bits(const fmpz_t a)
{
    if (COEFF_IS_MPZ(*a))
    {
        const __mpz_struct * A = COEFF_TO_PTR(*a);
        return A->_mp_size > 0 ? A->_mp_d[0] : -A->_mp_d[0];
    }
    else
    {
        return *a;
    }
}


static void binary_cubic_lift_fac(
    fmpz_t r,
    fmpz_t s,
    const fmpz_t a,
    const fmpz_t b,
    const fmpz_t inv,
    const fmpz_t r2,            /* old r^2 */
    slong e,
    slong n,
    fmpz_t c, fmpz_t d, fmpz_t t) /* temp */
{
    fmpz_mul_2exp(c, r2, e);
    fmpz_add(c, c, a);
    fmpz_sub(c, c, s);
    fmpz_fdiv_q_2exp(c, c, n);

    fmpz_set(d, b);
    fmpz_submul(d, r, s);
    fmpz_fdiv_q_2exp(d, d, n);

    fmpz_mul(t, d, r);
    fmpz_mul_2exp(t, t, e + 1);
    fmpz_addmul(t, c, s);
    fmpz_fdiv_r_2exp(t, t, n);
    fmpz_mul(t, t, inv);
    fmpz_fdiv_r_2exp(t, t, n);
    fmpz_mul_2exp(t, t, n);
    fmpz_add(s, s, t);

    fmpz_submul(d, c, r);
    fmpz_fdiv_r_2exp(d, d, n);
    fmpz_mul(d, d, inv);
    fmpz_fdiv_r_2exp(d, d, n);
    fmpz_mul_2exp(d, d, n);
    fmpz_add(r, r, d);
}


static void binary_cubic_lift_inv(
    fmpz_t inv,
    fmpz_t r2,  /* r^2 */
    const fmpz_t r,
    const fmpz_t s,
    slong e,
    slong n,
    fmpz_t t, fmpz_t t2)  /* temp */
{
    fmpz_mul(r2, r, r);
    fmpz_mul(t, inv, inv);
    fmpz_mul_2exp(t2, r2, e + 1);
    fmpz_add(t2, t2, s);
    fmpz_mul_2exp(inv, inv, 1);
    fmpz_submul(inv, t, t2);
    fmpz_fdiv_r_2exp(inv, inv, n);
}


/*
    factor 2^e*y^3 + a*y + b = (y + r)*(2^e*y^2 - r*y + s) mod 2^p

    assuming:
        a and b are both odd and e > 0, or
        a is even and b is odd and e = 0.
*/

static slong binary_cubic_lift(
    fmpz_t r,
    fmpz_t s,
    fmpz_t inv,
    const fmpz_t a,
    const fmpz_t b,
    slong e,
    slong p)
{
    slong n;
    fmpz_t r2, c, d, t;
    mp_limb_t A, B, C, D, INV, R, R2, S, E;

    /* start with a factorization mod 2^n */
    n = 1;

    A = fmpz_fdiv_r_2exp_flint_bits(a);
    B = fmpz_fdiv_r_2exp_flint_bits(b);
    R = 1;
    S = 1;
    INV = 1;
    R2 = R*R;
    E = (e < FLINT_BITS) ? (UWORD(1) << e) : 0;

    while (n <= FLINT_BITS/2)
    {
        mp_limb_t mask = (UWORD(1) << n);
        C = (A - (S - R2*E)) >> n;
        D = (B - (R*S)) >> n;
        R += (((D - C*R)*INV) % mask) << n;
        S += (((2*E*D*R + C*S)*INV) % mask) << n;
        n *= 2;
        R2 = R*R;
        INV = 2*INV - (INV*INV)*(2*R2*E + S);
    }

    fmpz_set_ui(r, R);
    fmpz_set_ui(s, S);
    fmpz_set_ui(inv, INV);

    if (n >= p)
        return n;

    fmpz_init(t);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_init_set_ui(r2, R);
    fmpz_mul_ui(r2, r2, R);

    while (n < p)
    {
        binary_cubic_lift_fac(r, s, a, b, inv, r2, e, n, c, d, t);
        n *= 2;
        if (n < p)
            binary_cubic_lift_inv(inv, r2, r, s, e, n, t, d);
    }

#if FLINT_WANT_ASSERT
    fmpz_mul(r2, r, r);
    fmpz_mul_2exp(c, r2, e);
    fmpz_add(c, c, a);
    fmpz_sub(c, c, s);
    fmpz_fdiv_r_2exp(t, c, n);
    FLINT_ASSERT(fmpz_is_zero(t));

    fmpz_set(d, b);
    fmpz_submul(d, r, s);
    fmpz_fdiv_r_2exp(t, d, n);
    FLINT_ASSERT(fmpz_is_zero(t));
#endif

    fmpz_clear(t);
    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(r2);

    return n;
}

/*
    on input
        n is even
        0 <= s < 2^n
        0 <= r < 2^n
        a = s - r^2*t mod 2^n
        b = r*s       mod 2^n
        inv = (s + 2*r^2*2^e)^-1 mod 2^(n/2)

    on return
        0 <= s < 2^(2n)
        0 <= r < 2^(2n)
        0 <= inv < 2^n
        a = s - r^2*t mod 2^(2n)
        b = r*s       mod 2^(2n)
        inv = (s + 2*r^2*2^e)^-1 mod 2^n
*/
static slong binary_cubic_lift_continue(
    fmpz_t r,
    fmpz_t s,
    fmpz_t inv,
    const fmpz_t a,
    const fmpz_t b,
    slong e,
    slong n)
{
    fmpz_t r2, t, c, d;

    fmpz_init(r2);
    fmpz_init(t);
    fmpz_init(c);
    fmpz_init(d);

    binary_cubic_lift_inv(inv, r2, r, s, e, n, t, d);
    binary_cubic_lift_fac(r, s, a, b, inv, r2, e, n, c, d, t);

    n *= 2;

#if FLINT_WANT_ASSERT
    fmpz_mul(r2, r, r);
    fmpz_mul_2exp(c, r2, e);
    fmpz_add(c, c, a);
    fmpz_sub(c, c, s);
    fmpz_fdiv_r_2exp(t, c, n);
    FLINT_ASSERT(fmpz_is_zero(t));

    fmpz_set(d, b);
    fmpz_submul(d, r, s);
    fmpz_fdiv_r_2exp(t, d, n);
    FLINT_ASSERT(fmpz_is_zero(t));
#endif

    fmpz_clear(r2);
    fmpz_clear(t);
    fmpz_clear(c);
    fmpz_clear(d);

    return n;
}


/* return f(0)*...*f(largest_prime - 1) mod prime_product */
static mp_limb_t eval_product_mod_n(
    const fmpz_t a,
    const fmpz_t b,
    mp_limb_t prime_product,
    mp_limb_t largest_prime)
{
    nmod_t ctx;
    ulong A, B, F, G, H, P, i;

    nmod_init(&ctx, prime_product);

    A = fmpz_fdiv_ui(a, ctx.n);
    B = fmpz_fdiv_ui(b, ctx.n);
    A = nmod_add(A, nmod_add(A, A, ctx), ctx);

    /*
        f(x) = B + A*x - x^3
        g(x) = f(x + 1) - f(x)
             = A - 1 - 3*x - 3*x^2
        h(x) = g(x) - g(x + 1)
             = 6 + 6*x
    */

    F = B;
    G = nmod_sub(A, 1, ctx);
    H = 6;
    P = F;
    for (i = 1; i < largest_prime; i++)
    {
        FLINT_ASSERT(H < prime_product);
        F = nmod_add(F, G, ctx);
        G = nmod_sub(G, H, ctx);
        H += 6;
        P = nmod_mul(P, F, ctx);
    }

    return P;
}

/* signed remainder mod 2^prec */
static void _fmpz_map_from_ZZ2(fmpz_t x, slong prec)
{
    FLINT_ASSERT(prec > 0);
    fmpz_fdiv_r_2exp(x, x, prec);
    if (fmpz_bits(x) >= prec)
    {
        fmpz_neg(x, x);
        fmpz_fdiv_r_2exp(x, x, prec);
        fmpz_neg(x, x);
    }
}


/*
    write the integer roots of f(X) = X^3 - 3*a*X - b to x
    return 0: irreducible cubic
           1: (X - x[0])*(irreducible quadratic)
           2: (X - x[0])*(X - x[1])*(X - x[2])
           3: (X - x[0])*(X - x[1])^2
           4: (X - x[0])^3

    a and b are clobbered!
*/
static int _fmpz_cubic_roots(fmpz * x, fmpz_t a, fmpz_t b)
{
    slong i;
    int ret, sign_a, sign_b, sign_d;
    fmpz_t d, t1, t2, t3, t4, ta, tb, r, s, inv, z;
    slong prec;
    slong cubic_prec, sqrt_prec;
    ulong alpha, beta, alpha2, beta3;
    fmpz two = 2;

    FLINT_ASSERT(a != b);

    sign_a = fmpz_sgn(a);
    sign_b = fmpz_sgn(b);

    if (fmpz_is_pm1(b))
    {
        if (sign_a == 0)
        {
            fmpz_swap(x + 0, b);
            return 1;
        }
        else
        {
            return 0;
        }
    }

    if (fmpz_is_zero(b))
    {
        fmpz_zero(x + 0);

        if (sign_a <= 0)
            return sign_a == 0 ? 4 : 1;

        fmpz_mul_ui(a, a, 3);

        if (fmpz_is_square(a))
        {
            fmpz_sqrt(x + 1, a);
            fmpz_neg(x + 2, x + 1);
            return 2;
        }
        else
        {
            return 1;
        }
    }

    if (sign_b < 0)
        fmpz_neg(b, b);

    /* b >= 2 now, and sign_b has the sign of the original b */

    /* check irreducibility mod 2 */
    if (fmpz_is_odd(b) && fmpz_is_odd(a))
        return 0;

    /* check irreducibility mod some p */
    FLINT_ASSERT(FLINT_BITS >= 32);
    if (0 != eval_product_mod_n(a, b, UWORD(5)*7*11*13*17*19*23*29, 29))
        return 0;

    fmpz_init(d);
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(t3);
    fmpz_init(t4);
    fmpz_init(ta);
    fmpz_init(tb);
    fmpz_init(r);
    fmpz_init(s);
    fmpz_init(inv);
    fmpz_init(z);

    /* d = b^2 - 4*a^3 = discriminant/-27 */
    fmpz_mul(d, b, b);
    fmpz_mul(t1, a, a);
    fmpz_mul(t2, t1, a);
    fmpz_submul_ui(d, t2, 4);

    sign_d = fmpz_sgn(d);

    if (sign_d == 0)
    {
        FLINT_ASSERT(fmpz_divisible(b, a));
        fmpz_divexact(x + 0, b, a);
        FLINT_ASSERT(fmpz_is_even(x + 0));
        fmpz_divexact_si(x + 1, x + 0, -2);
        ret = 3;
        goto cleanup;
    }

    if (sign_d > 0)
    {
        /*
            The real root is (cbrt(4*b + 4*sqrt(d)) + cbrt(4*b - 4*sqrt(d)))/2
            which is calculated by rounding intermediates to zero.
        */

        fmpz_sqrt(t3, d);
        fmpz_add(t4, b, t3);
        fmpz_mul_ui(t4, t4, 4);
        fmpz_root(t1, t4, 3);
        fmpz_sub(t4, b, t3);
        fmpz_mul_ui(t4, t4, 4);
        fmpz_root(t2, t4, 3);
        fmpz_add(t1, t1, t2);
        fmpz_fdiv_q_2exp(x + 0, t1, 1);

        /*
            If x is the actual real root in RR, then, experimentally,
            x[0] - 1 < x < x[0] + 2, so that x[0] - 1 cannot be a root. Test
            all three numbers x[0] - 1, x[0], x[0] + 1 anyways since these
            follow from easily-proven bounds.
        */

        fmpz_mul(t2, x + 0, x + 0);
        fmpz_sub(t3, t2, a);
        fmpz_submul_ui(t2, a, 3);
        fmpz_mul(t1, t2, x + 0);
        if (fmpz_equal(t1, b))
        {
            ret = 1;
            goto cleanup;
        }

        fmpz_submul_ui(b, x + 0, 3);
        fmpz_mul_ui(t3, t3, 3);
        fmpz_add_ui(t3, t3, 1);
        fmpz_sub(t2, b, t3);
        if (fmpz_equal(t1, t2))
        {
            fmpz_add_ui(x + 0, x + 0, 1);
            ret = 1;
            goto cleanup;
        }

        fmpz_add(t2, b, t3);
        if (fmpz_equal(t1, t2))
        {
            fmpz_sub_ui(x + 0, x + 0, 1);
            ret = 1;
            goto cleanup;
        }

        ret = 0;
        goto cleanup;
    }

    /*
        Three real roots. The observations
                f(-2*sqrt(a)) < 0,   f(-sqrt(a)) > 0
                f(-b/(2*a)) > 0,     f(-b/(3*a)) < 0
                f(sqrt(a)) < 0,      f(2*sqrt(a)) > 0
        give intervals bounding the roots. They are sought in ZZ_2 here.
    */

    /* check irreducibility mod some more p */
    FLINT_ASSERT(FLINT_BITS >= 32);
    if (0 != eval_product_mod_n(a, b, UWORD(31)*37*41*43*47, 47))
    {
        ret = 0;
        goto cleanup;
    }

    /* 2*sqrt(a) is bound on absolute value of roots */
    prec = fmpz_bits(a)/2 + 3;

    fmpz_mul_si(ta, a, -3);
    fmpz_mul_si(tb, b, -1);

    FLINT_ASSERT(!fmpz_is_zero(ta));
    FLINT_ASSERT(!fmpz_is_zero(tb));

    alpha = fmpz_remove(ta, ta, &two);
    beta = fmpz_remove(tb, tb, &two);

    if (3*alpha == 2*beta)
    {
        ret = 0;
        goto cleanup;
    }
    else if (3*alpha > 2*beta)
    {
        /* only one root, which has valuation beta/3 */

        beta3 = beta/3;
        if ((beta % 3) != 0)
        {
            ret = 0;
            goto cleanup;
        }

        fmpz_mul_2exp(ta, ta, alpha - 2*beta3);
        binary_cubic_lift(r, s, inv, ta, tb, 0, prec);
        fmpz_mul_2exp(x + 0, r, beta3);
        fmpz_neg(x + 0, x + 0);
        _fmpz_map_from_ZZ2(x + 0, prec);

        fmpz_mul(t1, x + 0, x + 0);
        fmpz_submul_ui(t1, a, 3);
        fmpz_mul(t2, t1, x + 0);
        ret = fmpz_equal(t2, b) ? 1 : 0;
        goto cleanup;
    }
    else
    {
        /* there is a root with valuation beta - alpha */

        alpha2 = alpha/2;
        if ((alpha % 2) != 0)
        {
            /* there are no other roots */

            binary_cubic_lift(r, s, inv, ta, tb, 2*beta - 3*alpha, prec);
            fmpz_mul_2exp(x + 0, r, beta - alpha);
            fmpz_neg(x + 0, x + 0);
            _fmpz_map_from_ZZ2(x + 0, prec);

            fmpz_mul(t1, x + 0, x + 0);
            fmpz_submul_ui(t1, a, 3);
            fmpz_mul(t2, t1, x + 0);

            ret = fmpz_equal(t2, b) ? 1 : 0;
            goto cleanup;
        }

        cubic_prec = binary_cubic_lift(r, s, inv, ta, tb, 2*beta - 3*alpha, prec);

        fmpz_mul_2exp(x + 0, r, beta - alpha);
        fmpz_neg(x + 0, x + 0);
        _fmpz_map_from_ZZ2(x + 0, prec);

        fmpz_mul(t1, x + 0, x + 0);
        fmpz_submul_ui(t1, a, 3);
        fmpz_mul(t2, t1, x + 0);

        /* there are two roots with valuation alpha/2 */

        if (fmpz_equal(t2, b))
        {
            fmpz_mul_ui(a, a, 4);
            fmpz_submul(a, x + 0, x + 0);
            fmpz_mul_ui(a, a, 3);
            if (fmpz_is_square(a))
            {
                fmpz_sqrt(t1, a);
                fmpz_sub(x + 1, t1, x + 0);
                fmpz_add(x + 2, t1, x + 0);
                FLINT_ASSERT(fmpz_is_even(x + 1));
                FLINT_ASSERT(fmpz_is_even(x + 2));
                fmpz_divexact_ui(x + 1, x + 1, 2);
                fmpz_divexact_si(x + 2, x + 2, -2);
                ret = 2;
                goto cleanup;
            }
            else
            {
                ret = 1;
                goto cleanup;
            }
        }

        /*
            The root with valuation beta - alpha is irrational. We have already
            found it as -2^(beta-alpha)*r so factor it out and find the roots
            of the remaining quadratic with the quadratic formula:

            2^(beta-alpha-1)*r +- 2^(alpha/2)*sqrt(2^(2*beta-3*alpha-2)*r^2 - s)

            r and s are calculated to precision O(2^cubic_prec), where
            cubic_prec >= prec. This formula might lose precision in the
            relatively common cases

                (1) alpha = 0, or
                (2) alpha > 0 and 2*beta-3*alpha-2 = 0

            Since binary_cubic_lift_continue doubles cubic_prec, and binary_sqrt
            outputs a sqrt_prec with sqrt_prec >= cubic_prec/2 - 1, at most
            two jumps to try_again are required.
        */

try_again:

        FLINT_ASSERT(2*beta - 3*alpha >= 2);

        fmpz_mul(d, r, r);
        fmpz_mul_2exp(d, d, 2*beta - 3*alpha - 2);
        fmpz_sub(d, d, s);
        sqrt_prec = binary_sqrt(z, d, cubic_prec);

        if (sqrt_prec < 0)
        {
            /* roots with valuation alpha/2 are irrational */
            ret = 0;
            goto cleanup;
        }

        if (sqrt_prec + alpha2 < prec)
        {
            cubic_prec = binary_cubic_lift_continue(r, s, inv, ta, tb,
                                                 2*beta - 3*alpha, cubic_prec);
            goto try_again;
        }

        fmpz_mul_2exp(r, r, beta - alpha - 1);
        fmpz_mul_2exp(z, z, alpha2);
        fmpz_add(x + 1, r, z);
        fmpz_sub(x + 2, r, z);

        for (i = 1; i <= 2; i++)
        {
            _fmpz_map_from_ZZ2(x + i, prec);

            fmpz_mul(t1, x + i, x + i);
            fmpz_submul_ui(t1, a, 3);
            fmpz_mul(t2, t1, x + i);
            if (fmpz_equal(t2, b))
            {
                fmpz_swap(x + i, x + 0);
                ret = 1;
                goto cleanup;
            }
        }

        ret = 0;
        goto cleanup;
    }

cleanup:

    if (sign_b < 0)
    {
        fmpz_neg(x + 0, x + 0);
        fmpz_neg(x + 1, x + 1);
        fmpz_neg(x + 2, x + 2);
    }

    fmpz_clear(d);
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(t3);
    fmpz_clear(t4);
    fmpz_clear(ta);
    fmpz_clear(tb);
    fmpz_clear(r);
    fmpz_clear(s);
    fmpz_clear(inv);
    fmpz_clear(z);

    return ret;
}


/* roots of depressed cubic back to original */
static void _raise_linear_factor(
    fmpz_poly_t p,
    const fmpz_t a,
    const fmpz_t b,
    const fmpz_t root,
    fmpz_t T)
{
    FLINT_ASSERT(p->alloc >= 2);
    fmpz_mul_ui(p->coeffs + 1, a, 3);
    fmpz_sub(p->coeffs + 0, b, root);
    fmpz_gcd(T, p->coeffs + 0, p->coeffs + 1);
    fmpz_divexact(p->coeffs + 0, p->coeffs + 0, T);
    fmpz_divexact(p->coeffs + 1, p->coeffs + 1, T);
    _fmpz_poly_set_length(p, 2);
}

void _fmpz_poly_factor_cubic(fmpz_poly_factor_t fac,
                                                const fmpz_poly_t f, slong exp)
{
    fmpz_t A, B, b2, T, ac;
    fmpz r[3];
    const fmpz *a, *b, *c, *d;
    fmpz_poly_t p;

    FLINT_ASSERT(f->length == 4);
    FLINT_ASSERT(fmpz_sgn(f->coeffs + 3) > 0);

    d = f->coeffs;
    c = f->coeffs + 1;
    b = f->coeffs + 2;
    a = f->coeffs + 3;

    /* depress the cubic to X^3 - 3*A*X - B */

    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(b2);
    fmpz_init(T);
    fmpz_init(ac);
    fmpz_init(r + 0);
    fmpz_init(r + 1);
    fmpz_init(r + 2);
    fmpz_poly_init2(p, 3);

    /* A = b^2 - 3*a*c */
    fmpz_mul(ac, a, c);
    fmpz_mul(b2, b, b);
    fmpz_set(A, b2);
    fmpz_submul_ui(A, ac, 3);

    /* B = (9*a*c-2*b^2)*b - 27*a^2*d */
    fmpz_mul_ui(b2, b2, 2);
    fmpz_mul_ui(ac, ac, 9);
    fmpz_sub(B, ac, b2);
    fmpz_mul(B, B, b);
    fmpz_mul(T, a, a);
    fmpz_mul(T, T, d);
    fmpz_submul_ui(B, T, 27);

    switch (_fmpz_cubic_roots(r, A, B))
    {
        case 4:
            _raise_linear_factor(p, a, b, r + 0, T);
            fmpz_poly_factor_insert(fac, p, 3*exp);
            break;

        case 3:
            _raise_linear_factor(p, a, b, r + 0, T);
            fmpz_poly_factor_insert(fac, p, 1*exp);
            _raise_linear_factor(p, a, b, r + 1, T);
            fmpz_poly_factor_insert(fac, p, 2*exp);
            break;

        case 2:
            _raise_linear_factor(p, a, b, r + 0, T);
            fmpz_poly_factor_insert(fac, p, 1*exp);
            _raise_linear_factor(p, a, b, r + 1, T);
            fmpz_poly_factor_insert(fac, p, 1*exp);
            _raise_linear_factor(p, a, b, r + 2, T);
            fmpz_poly_factor_insert(fac, p, 1*exp);
            break;

        case 1:
            _raise_linear_factor(p, a, b, r + 0, T);
            fmpz_poly_factor_insert(fac, p, 1*exp);
            fmpz_poly_divides(p, f, p);
            fmpz_poly_factor_insert(fac, p, 1*exp);
            break;

        default:
            fmpz_poly_factor_insert(fac, f, exp);
    }

    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(b2);
    fmpz_clear(T);
    fmpz_clear(ac);
    fmpz_clear(r + 0);
    fmpz_clear(r + 1);
    fmpz_clear(r + 2);

    fmpz_poly_clear(p);
}

