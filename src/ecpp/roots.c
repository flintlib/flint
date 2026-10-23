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
#include "fmpz_mod_poly.h"
#include "ecpp.h"

/*
    Roots of polynomials of degree 2, 3 and 4 that split completely modulo
    a prime n, by radicals: a few square and cube roots instead of a
    powering of a polynomial modulo the polynomial, which costs some ten
    times more at these degrees. For degrees 3 and 4 the cube roots are
    taken in F_n when n = 1 mod 3 and in F_{n^2} when n = 2, 5 mod 9; the
    caller falls back to the generic method otherwise or on failure. n is assumed prime: for
    composite n the functions may fail, but a wrong root is harmless since
    the curve is tested afterwards.
*/

/* w^3 = z for a cube z modulo n = 1 mod 3 (Adleman-Manders-Miller for r = 3) */
static int
_cbrtmod(fmpz_t w, const fmpz_t z, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_t t, e, y, g, G, tmp, pw;
    slong s = 0, k, j, tries;
    int result = 0;

    if (fmpz_is_zero(z))
    {
        fmpz_zero(w);
        return 1;
    }

    fmpz_init(t); fmpz_init(e); fmpz_init(y); fmpz_init(g); fmpz_init(G);
    fmpz_init(tmp); fmpz_init(pw);

    /* n - 1 = 3^s t */
    fmpz_sub_ui(t, n, 1);
    while (fmpz_divisible_si(t, 3))
    {
        fmpz_divexact_si(t, t, 3);
        s++;
    }
    if (s == 0)
        goto cleanup;

    /* w = z^e with 3 e = 1 mod t; then w^3 = z y^k, y = z^t of order | 3^{s-1} */
    fmpz_set_ui(e, 3);
    fmpz_invmod(e, e, t);
    fmpz_powm(w, z, e, n);
    fmpz_powm(y, z, t, n);
    if (fmpz_is_one(y))
    {
        result = 1;
        goto cleanup;
    }
    /* k = (3 e - 1) / t */
    fmpz_mul_ui(tmp, e, 3);
    fmpz_sub_ui(tmp, tmp, 1);
    fmpz_divexact(tmp, tmp, t);
    fmpz_powm(y, y, tmp, n);        /* y^k, in the 3-Sylow subgroup, a cube there */

    /* G = g^t of order 3^s for a non-cube g */
    for (tries = 0; tries < 100; tries++)
    {
        fmpz_randm(g, state, n);
        fmpz_powm(G, g, t, n);
        fmpz_set(pw, G);
        for (j = 0; j < s - 1; j++)
            fmpz_mod_pow_ui(pw, pw, 3, ctx);
        if (!fmpz_is_one(pw) && !fmpz_is_zero(pw))
            break;
    }
    if (tries == 100)
        goto cleanup;

    /*
        Discrete logarithm of y^k = G^{3 m} to base G (Pohlig-Hellman in a
        cyclic group of order 3^s, digit by digit); then w /= G^m.
    */
    {
        fmpz_t m, Gi, target, cur;
        slong digit;
        fmpz_init(m); fmpz_init(Gi); fmpz_init(target); fmpz_init(cur);
        fmpz_zero(m);           /* accumulates 3 m */
        fmpz_set(target, y);
        fmpz_mod_inv(Gi, G, ctx);
        for (k = 0; k < s; k++)
        {
            /* digit d: (target / G^{3m})^{3^{s-1-k}} = (G^{3^{s-1}})^d */
            fmpz_powm(cur, Gi, m, n);
            fmpz_mod_mul(cur, cur, target, ctx);
            for (j = 0; j < s - 1 - k; j++)
                fmpz_mod_pow_ui(cur, cur, 3, ctx);
            fmpz_set(pw, G);
            for (j = 0; j < s - 1; j++)
                fmpz_mod_pow_ui(pw, pw, 3, ctx);    /* element of order 3 */
            for (digit = 0; digit < 3; digit++)
            {
                if (fmpz_is_one(cur))
                    break;
                fmpz_mod_mul(cur, cur, pw, ctx);    /* multiply by pw^{-1}... use pw^2 = pw^{-1} */
                fmpz_mod_mul(cur, cur, pw, ctx);
            }
            if (digit == 3)
                break;
            fmpz_set_ui(tmp, 1);
            for (j = 0; j < k; j++)
                fmpz_mul_ui(tmp, tmp, 3);
            fmpz_addmul_ui(m, tmp, digit);
        }
        if (k == s && fmpz_divisible_si(m, 3))
        {
            fmpz_divexact_si(m, m, 3);
            fmpz_powm(cur, Gi, m, n);
            fmpz_mod_mul(w, w, cur, ctx);
            fmpz_mod_pow_ui(cur, w, 3, ctx);
            result = fmpz_equal(cur, z);
        }
        fmpz_clear(m); fmpz_clear(Gi); fmpz_clear(target); fmpz_clear(cur);
    }

cleanup:
    fmpz_clear(t); fmpz_clear(e); fmpz_clear(y); fmpz_clear(g); fmpz_clear(G);
    fmpz_clear(tmp); fmpz_clear(pw);
    return result;
}

/*
    Arithmetic in F_{n^2} = F_n[s] / (s^2 + 3) for n = 2 mod 3 (where -3 is
    not a square): elements a + b s.
*/
static void
_fn2_mul(fmpz_t ra, fmpz_t rb, const fmpz_t a1, const fmpz_t b1,
                        const fmpz_t a2, const fmpz_t b2, const fmpz_mod_ctx_t ctx)
{
    fmpz_t t, u;
    fmpz_init(t); fmpz_init(u);
    fmpz_mod_mul(t, a1, a2, ctx);
    fmpz_mod_mul(u, b1, b2, ctx);
    fmpz_mod_mul_ui(u, u, 3, ctx);
    fmpz_mod_sub(t, t, u, ctx);             /* a1 a2 - 3 b1 b2 */
    fmpz_mod_mul(u, a1, b2, ctx);
    fmpz_mod_addmul(u, u, a2, b1, ctx);     /* a1 b2 + a2 b1 */
    fmpz_swap(ra, t);
    fmpz_swap(rb, u);
    fmpz_clear(t); fmpz_clear(u);
}

static void
_fn2_pow(fmpz_t ra, fmpz_t rb, const fmpz_t a, const fmpz_t b, const fmpz_t e,
                                                        const fmpz_mod_ctx_t ctx)
{
    fmpz_t xa, xb;
    slong i;
    fmpz_init(xa); fmpz_init(xb);
    fmpz_one(xa);
    fmpz_zero(xb);
    for (i = fmpz_bits(e) - 1; i >= 0; i--)
    {
        _fn2_mul(xa, xb, xa, xb, xa, xb, ctx);
        if (fmpz_tstbit(e, i))
            _fn2_mul(xa, xb, xa, xb, a, b, ctx);
    }
    fmpz_swap(ra, xa);
    fmpz_swap(rb, xb);
    fmpz_clear(xa); fmpz_clear(xb);
}

/*
    u^3 = z in F_{n^2} for n = 2, 5 mod 9 (so that 3 divides n^2 - 1 exactly
    once and a cube root of a cube z is z^e with 3 e = 1 mod (n^2 - 1) / 3).
*/
static int
_cbrt_fn2(fmpz_t ua, fmpz_t ub, const fmpz_t za, const fmpz_t zb, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_t t, e, ca, cb;
    int ok;
    fmpz_init(t); fmpz_init(e); fmpz_init(ca); fmpz_init(cb);
    fmpz_mul(t, n, n);
    fmpz_sub_ui(t, t, 1);
    fmpz_divexact_ui(t, t, 3);
    fmpz_set_ui(e, 3);
    fmpz_invmod(e, e, t);
    _fn2_pow(ua, ub, za, zb, e, ctx);
    _fn2_mul(ca, cb, ua, ub, ua, ub, ctx);
    _fn2_mul(ca, cb, ca, cb, ua, ub, ctx);
    ok = fmpz_equal(ca, za) && fmpz_equal(cb, zb);
    fmpz_clear(t); fmpz_clear(e); fmpz_clear(ca); fmpz_clear(cb);
    return ok;
}

/* a root of x^2 + c1 x + c0 */
static int
_root_quadratic(fmpz_t x, const fmpz_t c1, const fmpz_t c0, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_t d, r, h;
    int ok;

    fmpz_init(d); fmpz_init(r); fmpz_init(h);
    fmpz_mod_mul(d, c1, c1, ctx);
    fmpz_mod_mul_ui(r, c0, 4, ctx);
    fmpz_mod_sub(d, d, r, ctx);
    ok = fmpz_sqrtmod(r, d, n);
    if (ok)
    {
        fmpz_mod_sub(r, r, c1, ctx);
        fmpz_set_ui(h, 2);
        fmpz_mod_inv(h, h, ctx);
        fmpz_mod_mul(x, r, h, ctx);
    }
    fmpz_clear(d); fmpz_clear(r); fmpz_clear(h);
    return ok;
}

/* a root of x^3 + c2 x^2 + c1 x + c0 (Cardano) */
static int
_root_cubic(fmpz_t x, const fmpz_t c2, const fmpz_t c1, const fmpz_t c0,
                                    flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_t p, q, t, u, v, i3, shift;
    int ok = 0;

    fmpz_init(p); fmpz_init(q); fmpz_init(t); fmpz_init(u); fmpz_init(v);
    fmpz_init(i3); fmpz_init(shift);

    /* x = y - c2/3: y^3 + p y + q, p = c1 - c2^2/3, q = c0 - c1 c2/3 + 2 c2^3/27 */
    fmpz_set_ui(i3, 3);
    fmpz_mod_inv(i3, i3, ctx);
    fmpz_mod_mul(shift, c2, i3, ctx);           /* c2 / 3 */
    fmpz_mod_mul(t, c2, shift, ctx);            /* c2^2 / 3 */
    fmpz_mod_sub(p, c1, t, ctx);
    fmpz_mod_mul(t, shift, shift, ctx);         /* c2^2 / 9 */
    fmpz_mod_mul(t, t, shift, ctx);             /* c2^3 / 27 */
    fmpz_mod_add(q, t, t, ctx);
    fmpz_mod_mul(t, c1, shift, ctx);
    fmpz_mod_sub(q, q, t, ctx);
    fmpz_mod_add(q, q, c0, ctx);

    if (fmpz_is_zero(p))
    {
        /* y^3 = -q */
        fmpz_mod_neg(t, q, ctx);
        if (fmpz_fdiv_ui(n, 3) == 1)
            ok = _cbrtmod(u, t, state, ctx);
        else
        {
            /* cube roots are unique: (-q)^e with 3 e = 1 mod n - 1 */
            fmpz_sub_ui(v, n, 1);
            fmpz_set_ui(u, 3);
            fmpz_invmod(u, u, v);
            fmpz_powm(u, t, u, n);
            ok = 1;
        }
        if (ok)
            fmpz_mod_sub(x, u, shift, ctx);
        goto cleanup;
    }

    if (fmpz_fdiv_ui(n, 3) != 1)
    {
        /*
            n = 2 mod 3: the resolvent z = -q/2 + sqrt(disc) s / 18, with
            s = sqrt(-3), lies in F_{n^2}; for a cube root u of z the
            conjugate of u is -p/(3u), so y = u + conj(u) = 2 re(u).
        */
        fmpz_t za, zb, ua, ub, disc;
        fmpz_init(za); fmpz_init(zb); fmpz_init(ua); fmpz_init(ub); fmpz_init(disc);
        /* disc = -(4 p^3 + 27 q^2) */
        fmpz_mod_mul(t, p, p, ctx);
        fmpz_mod_mul(t, t, p, ctx);
        fmpz_mod_mul_ui(t, t, 4, ctx);
        fmpz_mod_mul(u, q, q, ctx);
        fmpz_mod_mul_ui(u, u, 27, ctx);
        fmpz_mod_add(disc, t, u, ctx);
        fmpz_mod_neg(disc, disc, ctx);
        if (fmpz_sqrtmod(disc, disc, n))
        {
            fmpz_set_ui(t, 2);
            fmpz_mod_inv(t, t, ctx);
            fmpz_mod_mul(za, q, t, ctx);
            fmpz_mod_neg(za, za, ctx);
            fmpz_set_ui(t, 18);
            fmpz_mod_inv(t, t, ctx);
            fmpz_mod_mul(zb, disc, t, ctx);
            if (_cbrt_fn2(ua, ub, za, zb, ctx))
            {
                fmpz_mod_add(x, ua, ua, ctx);
                fmpz_mod_sub(x, x, shift, ctx);
                fmpz_mod_add(t, x, c2, ctx);
                fmpz_mod_mul(t, t, x, ctx);
                fmpz_mod_add(t, t, c1, ctx);
                fmpz_mod_mul(t, t, x, ctx);
                fmpz_mod_add(t, t, c0, ctx);
                ok = fmpz_is_zero(t);
            }
        }
        fmpz_clear(za); fmpz_clear(zb); fmpz_clear(ua); fmpz_clear(ub); fmpz_clear(disc);
        goto cleanup;
    }

    /* u^3 = (-q + sqrt(q^2 + 4 p^3 / 27)) / 2, v = -p / (3 u) */
    fmpz_mod_mul(t, p, p, ctx);
    fmpz_mod_mul(t, t, p, ctx);
    fmpz_mod_mul_ui(t, t, 4, ctx);
    fmpz_mod_mul(t, t, i3, ctx);
    fmpz_mod_mul(t, t, i3, ctx);
    fmpz_mod_mul(t, t, i3, ctx);                /* 4 p^3 / 27 */
    fmpz_mod_mul(u, q, q, ctx);
    fmpz_mod_add(t, t, u, ctx);
    if (!fmpz_sqrtmod(t, t, n))
        goto cleanup;
    fmpz_mod_sub(t, t, q, ctx);
    fmpz_set_ui(u, 2);
    fmpz_mod_inv(u, u, ctx);
    fmpz_mod_mul(t, t, u, ctx);
    if (fmpz_is_zero(t))
    {
        /* the other sign */
        fmpz_mod_neg(t, q, ctx);
    }
    if (!_cbrtmod(u, t, state, ctx))
        goto cleanup;
    if (fmpz_is_zero(u))
        goto cleanup;
    fmpz_mod_mul_ui(v, u, 3, ctx);
    fmpz_mod_inv(v, v, ctx);
    fmpz_mod_mul(v, v, p, ctx);
    fmpz_mod_sub(x, u, v, ctx);                 /* u - p/(3u) */
    fmpz_mod_sub(x, x, shift, ctx);

    /* check */
    fmpz_mod_add(t, x, c2, ctx);
    fmpz_mod_mul(t, t, x, ctx);
    fmpz_mod_add(t, t, c1, ctx);
    fmpz_mod_mul(t, t, x, ctx);
    fmpz_mod_add(t, t, c0, ctx);
    ok = fmpz_is_zero(t);

cleanup:
    fmpz_clear(p); fmpz_clear(q); fmpz_clear(t); fmpz_clear(u); fmpz_clear(v);
    fmpz_clear(i3); fmpz_clear(shift);
    return ok;
}

/* a root of x^4 + c3 x^3 + c2 x^2 + c1 x + c0 (Ferrari) */
static int
_root_quartic(fmpz_t x, const fmpz_t c3, const fmpz_t c2, const fmpz_t c1,
                    const fmpz_t c0, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_t p, q, r, shift, t, u, m, s, i2, a1, a0;
    int ok = 0;

    fmpz_init(p); fmpz_init(q); fmpz_init(r); fmpz_init(shift); fmpz_init(t);
    fmpz_init(u); fmpz_init(m); fmpz_init(s); fmpz_init(i2); fmpz_init(a1); fmpz_init(a0);

    /* x = y - c3/4: y^4 + p y^2 + q y + r */
    fmpz_set_ui(i2, 2);
    fmpz_mod_inv(i2, i2, ctx);
    fmpz_mod_mul(shift, c3, i2, ctx);
    fmpz_mod_mul(shift, shift, i2, ctx);        /* c3 / 4 */
    /* p = c2 - 6 shift^2, q = c1 - 2 c2 shift + 8 shift^3,
       r = c0 - c1 shift + c2 shift^2 - 3 shift^4 */
    fmpz_mod_mul(t, shift, shift, ctx);         /* shift^2 */
    fmpz_mod_mul_ui(p, t, 6, ctx);
    fmpz_mod_sub(p, c2, p, ctx);
    fmpz_mod_mul(u, t, shift, ctx);             /* shift^3 */
    fmpz_mod_mul_ui(q, u, 8, ctx);
    fmpz_mod_mul(s, c2, shift, ctx);
    fmpz_mod_sub(q, q, s, ctx);
    fmpz_mod_sub(q, q, s, ctx);
    fmpz_mod_add(q, q, c1, ctx);
    fmpz_mod_mul(r, t, t, ctx);                 /* shift^4 */
    fmpz_mod_mul_ui(r, r, 3, ctx);
    fmpz_mod_neg(r, r, ctx);
    fmpz_mod_mul(s, c2, t, ctx);
    fmpz_mod_add(r, r, s, ctx);
    fmpz_mod_mul(s, c1, shift, ctx);
    fmpz_mod_sub(r, r, s, ctx);
    fmpz_mod_add(r, r, c0, ctx);

    if (fmpz_is_zero(q))
    {
        /* biquadratic: z = y^2 root of z^2 + p z + r */
        if (!_root_quadratic(t, p, r, ctx))
            goto cleanup;
        if (!fmpz_sqrtmod(t, t, n))
            goto cleanup;
        fmpz_mod_sub(x, t, shift, ctx);
        ok = 1;
        goto cleanup;
    }

    /*
        Resolvent: 8 m^3 + 8 p m^2 + (2 p^2 - 8 r) m - q^2 = 0, i.e.
        m^3 + p m^2 + (p^2/4 - r) m - q^2/8 = 0; then
        y^4 + p y^2 + q y + r = (y^2 + p/2 + m)^2 - (sqrt(2m) y - q / (2 sqrt(2m)))^2.
    */
    fmpz_mod_mul(a1, p, p, ctx);
    fmpz_mod_mul(a1, a1, i2, ctx);
    fmpz_mod_mul(a1, a1, i2, ctx);
    fmpz_mod_sub(a1, a1, r, ctx);
    fmpz_mod_mul(a0, q, q, ctx);
    fmpz_mod_mul(a0, a0, i2, ctx);
    fmpz_mod_mul(a0, a0, i2, ctx);
    fmpz_mod_mul(a0, a0, i2, ctx);
    fmpz_mod_neg(a0, a0, ctx);
    if (!_root_cubic(m, p, a1, a0, state, ctx))
        goto cleanup;
    if (fmpz_is_zero(m))
        goto cleanup;
    fmpz_mod_add(s, m, m, ctx);
    if (!fmpz_sqrtmod(s, s, n))                  /* s = sqrt(2m) */
        goto cleanup;
    /* y^2 - s y + (p/2 + m + q/(2s)) = 0 */
    fmpz_mod_add(t, s, s, ctx);
    fmpz_mod_inv(t, t, ctx);
    fmpz_mod_mul(t, t, q, ctx);
    fmpz_mod_mul(u, p, i2, ctx);
    fmpz_mod_add(u, u, m, ctx);
    fmpz_mod_add(u, u, t, ctx);
    fmpz_mod_neg(t, s, ctx);
    if (!_root_quadratic(x, t, u, ctx))
    {
        /* the other factor: y^2 + s y + (p/2 + m - q/(2s)) */
        fmpz_mod_sub(u, u, t, ctx);
        fmpz_mod_sub(u, u, t, ctx);
        if (!_root_quadratic(x, s, u, ctx))
            goto cleanup;
    }
    fmpz_mod_sub(x, x, shift, ctx);

    /* check */
    fmpz_mod_add(t, x, c3, ctx);
    fmpz_mod_mul(t, t, x, ctx);
    fmpz_mod_add(t, t, c2, ctx);
    fmpz_mod_mul(t, t, x, ctx);
    fmpz_mod_add(t, t, c1, ctx);
    fmpz_mod_mul(t, t, x, ctx);
    fmpz_mod_add(t, t, c0, ctx);
    ok = fmpz_is_zero(t);

cleanup:
    fmpz_clear(p); fmpz_clear(q); fmpz_clear(r); fmpz_clear(shift); fmpz_clear(t);
    fmpz_clear(u); fmpz_clear(m); fmpz_clear(s); fmpz_clear(i2); fmpz_clear(a1); fmpz_clear(a0);
    return ok;
}

/*
    A root of the monic polynomial f of degree 2, 3 or 4, assumed to split
    completely modulo the prime modulus of ctx, by radicals. Returns 1 and
    sets x on success, 0 if the degree is not handled (degrees 3 and 4
    require n = 1 mod 3) or the computation failed.
*/
int
ecpp_root_radicals(fmpz_t x, const fmpz_mod_poly_t f, flint_rand_t state,
                                                        const fmpz_mod_ctx_t ctx)
{
    slong d = fmpz_mod_poly_degree(f, ctx);
    fmpz c[4];
    slong i;
    int ok = 0;

    if (d < 2 || d > 4)
        return 0;
    if (d >= 3)
    {
        ulong n9 = fmpz_fdiv_ui(fmpz_mod_ctx_modulus(ctx), 9);
        /* n = 1 mod 3 in F_n; n = 2, 5 mod 9 through F_{n^2} */
        if (n9 % 3 != 1 && n9 != 2 && n9 != 5)
            return 0;
    }
    if (!fmpz_is_one(f->coeffs + d))
        return 0;

    for (i = 0; i < 4; i++)
        fmpz_init(c + i);
    for (i = 0; i < d; i++)
        fmpz_mod_poly_get_coeff_fmpz(c + i, f, i, ctx);

    if (d == 2)
        ok = _root_quadratic(x, c + 1, c + 0, ctx);
    else if (d == 3)
        ok = _root_cubic(x, c + 2, c + 1, c + 0, state, ctx);
    else
        ok = _root_quartic(x, c + 3, c + 2, c + 1, c + 0, state, ctx);

    for (i = 0; i < 4; i++)
        fmpz_clear(c + i);
    return ok;
}

/*
    A root of the monic polynomial f, assumed to split completely modulo
    the prime modulus of ctx: by radicals when possible, else by random
    splitting with (x + r)^((n-1)/2). Returns 1 and sets x, or 0.
*/
int
ecpp_poly_root(fmpz_t x, const fmpz_mod_poly_t f, flint_rand_t state,
                                                        const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_mod_poly_t g, base, pw, finv;
    fmpz_t e, r;
    slong tries = 0;
    int found = 0;

    if (fmpz_mod_poly_degree(f, ctx) < 1)
        return 0;
    if (fmpz_mod_poly_degree(f, ctx) == 1)
    {
        fmpz_mod_poly_get_coeff_fmpz(x, f, 0, ctx);
        fmpz_mod_neg(x, x, ctx);
        return 1;
    }
    if (fmpz_mod_poly_degree(f, ctx) <= 4 && ecpp_root_radicals(x, f, state, ctx))
        return 1;

    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(base, ctx);
    fmpz_mod_poly_init(pw, ctx);
    fmpz_mod_poly_init(finv, ctx);
    fmpz_init(e);
    fmpz_init(r);
    fmpz_mod_poly_set(g, f, ctx);
    fmpz_sub_ui(e, n, 1);
    fmpz_fdiv_q_2exp(e, e, 1);

    while (!found && fmpz_mod_poly_degree(g, ctx) > 1 && tries < 30)
    {
        slong d = fmpz_mod_poly_degree(g, ctx);

        if (d <= 4 && ecpp_root_radicals(x, g, state, ctx))
        {
            found = 1;
            break;
        }
        /* gcd(g, (x + r)^((n-1)/2) - 1) */
        fmpz_randm(r, state, n);
        fmpz_mod_poly_zero(base, ctx);
        fmpz_mod_poly_set_coeff_ui(base, 1, 1, ctx);
        fmpz_mod_poly_set_coeff_fmpz(base, 0, r, ctx);
        fmpz_mod_poly_reverse(finv, g, d + 1, ctx);
        fmpz_mod_poly_inv_series(finv, finv, d + 1, ctx);
        fmpz_mod_poly_powmod_fmpz_binexp_preinv(pw, base, e, g, finv, ctx);
        fmpz_mod_poly_sub_si(pw, pw, 1, ctx);
        fmpz_mod_poly_gcd(pw, pw, g, ctx);
        if (fmpz_mod_poly_degree(pw, ctx) >= 1 && fmpz_mod_poly_degree(pw, ctx) < d)
        {
            if (2 * fmpz_mod_poly_degree(pw, ctx) > d)
                fmpz_mod_poly_div(pw, g, pw, ctx);
            fmpz_mod_poly_make_monic(g, pw, ctx);
        }
        tries++;
    }
    if (!found && fmpz_mod_poly_degree(g, ctx) == 1)
    {
        fmpz_mod_poly_get_coeff_fmpz(x, g, 0, ctx);
        fmpz_mod_neg(x, x, ctx);
        found = 1;
    }

    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(base, ctx);
    fmpz_mod_poly_clear(pw, ctx);
    fmpz_mod_poly_clear(finv, ctx);
    fmpz_clear(e);
    fmpz_clear(r);
    return found;
}
