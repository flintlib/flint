/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"
#include "impl.h"

/*
    Binary operations where one operand is quadratic.

    A quadratic q with minimal polynomial A t^2 + B t + C is written
    as q = a + b sqrt(D) with D = B^2 - 4AC (not necessarily squarefree),
    a = -B / (2A), b = +/- 1 / (2A), where sqrt(D) denotes the principal
    square root (i sqrt(-D) if D < 0) and the sign of b is determined
    numerically.

    If both operands are quadratic, the result is computed in closed form
    in Q(sqrt(D1), sqrt(D2)) (which is a quadratic field if D1 D2 is a
    square); this requires no factoring at all.

    If only y = a + b sqrt(D) is quadratic and x has degree d >= 3 with
    minimal polynomial P, then x = L(w) where w = x op y and L is a
    linear polynomial with coefficients in Q(sqrt(D)). Writing
    G(z) = P(L(z)) = U(z) + sqrt(D) V(z), the norm H = U^2 - D V^2 is an
    annihilating polynomial of degree 2d for w. If sqrt(D) is not in Q(x),
    then P is irreducible over Q(sqrt(D)), hence so is G, so H is a power of
    an irreducible polynomial, and H is irreducible if it is squarefree.
    We certify this using a prime p (not dividing lc(P) lc(H) 2 D) such that
    P mod p is squarefree and has an irreducible factor of odd degree f
    (giving a prime of Q(x) above p with residue field F_{p^f}), D is a
    quadratic nonresidue mod p and hence in F_{p^f} (so sqrt(D) is not
    in Q(x)), and H mod p is squarefree of full degree. If no such prime is
    found quickly, we factor H.
*/

/* Write x = a + b sqrt(D) where x is quadratic. */
static void
_qqbar_get_quadratic_parts(fmpq_t a, fmpq_t b, fmpz_t D, const qqbar_t x)
{
    const fmpz * A = QQBAR_COEFFS(x) + 2;
    const fmpz * B = QQBAR_COEFFS(x) + 1;
    const fmpz * C = QQBAR_COEFFS(x);
    acb_t z;
    arb_t s;
    slong prec;
    int sign = 0;

    /* D = B^2 - 4AC */
    fmpz_mul(D, A, C);
    fmpz_mul_2exp(D, D, 2);
    fmpz_submul(D, B, B);
    fmpz_neg(D, D);

    /* a = -B / (2A) */
    fmpz_neg(fmpq_numref(a), B);
    fmpz_mul_2exp(fmpq_denref(a), A, 1);
    fmpq_canonicalise(a);

    /* sign: 2A x + B = +/- sqrt(D) */
    acb_init(z);
    arb_init(s);
    acb_set(z, QQBAR_ENCLOSURE(x));

    for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        acb_mul_fmpz(z, QQBAR_ENCLOSURE(x), A, prec);
        acb_mul_2exp_si(z, z, 1);
        acb_add_fmpz(z, z, B, prec);

        if (fmpz_sgn(D) > 0)
            arb_set(s, acb_realref(z));
        else
            arb_set(s, acb_imagref(z));

        if (arb_is_positive(s))
        {
            sign = 1;
            break;
        }
        else if (arb_is_negative(s))
        {
            sign = -1;
            break;
        }

        /* should not happen since the enclosure isolates the root */
        {
            acb_t t;
            acb_init(t);
            qqbar_enclosure_raw(t, x, 2 * prec);
            acb_mul_fmpz(z, t, A, 2 * prec);
            acb_mul_2exp_si(z, z, 1);
            acb_add_fmpz(z, z, B, 2 * prec);
            acb_clear(t);

            if (fmpz_sgn(D) > 0)
                arb_set(s, acb_realref(z));
            else
                arb_set(s, acb_imagref(z));

            if (arb_is_positive(s))
            {
                sign = 1;
                break;
            }
            else if (arb_is_negative(s))
            {
                sign = -1;
                break;
            }
        }
    }

    /* b = sign / (2A) */
    fmpz_set_si(fmpq_numref(b), sign);
    fmpz_mul_2exp(fmpq_denref(b), A, 1);
    fmpq_canonicalise(b);

    acb_clear(z);
    arb_clear(s);
}

/* Compute the enclosure for res given its (irreducible) minimal polynomial,
   where res = x op y, op = 0, 1, 2, 3. */
static void
_qqbar_set_minpoly_from_op(qqbar_t res, const fmpz_poly_t poly, const qqbar_t x, const qqbar_t y, int op)
{
    fmpz_poly_factor_t fac;

    if (fmpz_poly_degree(poly) == 1)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpz_neg(fmpq_numref(t), poly->coeffs);
        fmpz_set(fmpq_denref(t), poly->coeffs + 1);
        fmpq_canonicalise(t);
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        return;
    }

    /* mock up a factorization with a single factor */
    fac->p = (fmpz_poly_struct *) poly;
    fac->num = 1;
    fac->alloc = 1;
    fac->exp = NULL;
    *(&fac->c) = 1;

    _qqbar_binary_op_select_factor(res, fac, x, y, op);
}

/* Set poly to the primitive integer polynomial with positive leading
   coefficient proportional to the given rational polynomial. */
static void
_fmpz_poly_set_fmpq_poly_primitive(fmpz_poly_t res, const fmpq_poly_t poly)
{
    fmpq_poly_get_numerator(res, poly);
    fmpz_poly_primitive_part(res, res);
    if (res->length > 0 && fmpz_sgn(res->coeffs + res->length - 1) < 0)
        fmpz_poly_neg(res, res);
}

/* Minimal polynomial of u + v sqrt(D) (v != 0, D not a square). */
static void
_quadratic_minpoly(fmpz_poly_t res, const fmpq_t u, const fmpq_t v, const fmpz_t D)
{
    fmpq_poly_t p;
    fmpq_t t, s;

    fmpq_poly_init(p);
    fmpq_init(t);
    fmpq_init(s);

    /* z^2 - 2u z + (u^2 - v^2 D) */
    fmpq_mul(t, u, u);
    fmpq_mul(s, v, v);
    fmpq_mul_fmpz(s, s, D);
    fmpq_sub(t, t, s);
    fmpq_poly_set_coeff_fmpq(p, 0, t);
    fmpq_mul_si(t, u, -2);
    fmpq_poly_set_coeff_fmpq(p, 1, t);
    fmpq_poly_set_coeff_si(p, 2, 1);

    _fmpz_poly_set_fmpq_poly_primitive(res, p);

    fmpq_poly_clear(p);
    fmpq_clear(t);
    fmpq_clear(s);
}

/* Both x and y quadratic. */
static void
_qqbar_binary_op_quadratic_quadratic(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)
{
    fmpq_t a1, b1, a2, b2, u, v, w, t, s, r;
    fmpz_t D1, D2, D12, q;
    fmpz_poly_t poly;

    fmpq_init(a1); fmpq_init(b1); fmpq_init(a2); fmpq_init(b2);
    fmpq_init(u); fmpq_init(v); fmpq_init(w); fmpq_init(t);
    fmpq_init(s); fmpq_init(r);
    fmpz_init(D1); fmpz_init(D2); fmpz_init(D12); fmpz_init(q);
    fmpz_poly_init(poly);

    _qqbar_get_quadratic_parts(a1, b1, D1, x);
    _qqbar_get_quadratic_parts(a2, b2, D2, y);

    fmpz_mul(D12, D1, D2);

    if (fmpz_sgn(D12) > 0 && fmpz_is_square(D12))
    {
        /* Same field: sqrt(D2) = r sqrt(D1) with r = sqrt(D2/D1) > 0
           (for D1, D2 of the same sign, the principal square roots
           satisfy this). */
        fmpz_sqrt(q, D12);
        fmpz_set(fmpq_numref(r), q);
        fmpz_abs(fmpq_denref(r), D1);
        fmpq_canonicalise(r);
        fmpq_mul(b2, b2, r);

        /* x = a1 + b1 s, y = a2 + b2 s, s^2 = D1; result u + v s */
        if (op == 0)
        {
            fmpq_add(u, a1, a2);
            fmpq_add(v, b1, b2);
        }
        else if (op == 1)
        {
            fmpq_sub(u, a1, a2);
            fmpq_sub(v, b1, b2);
        }
        else
        {
            if (op == 3)
            {
                /* 1/y = (a2 - b2 s) / (a2^2 - b2^2 D1) */
                fmpq_mul(t, a2, a2);
                fmpq_mul(s, b2, b2);
                fmpq_mul_fmpz(s, s, D1);
                fmpq_sub(t, t, s);
                fmpq_div(a2, a2, t);
                fmpq_div(b2, b2, t);
                fmpq_neg(b2, b2);
            }

            /* (a1 + b1 s)(a2 + b2 s) */
            fmpq_mul(u, a1, a2);
            fmpq_mul(t, b1, b2);
            fmpq_mul_fmpz(t, t, D1);
            fmpq_add(u, u, t);
            fmpq_mul(v, a1, b2);
            fmpq_addmul(v, a2, b1);
        }

        if (fmpq_is_zero(v))
        {
            qqbar_set_fmpq(res, u);
        }
        else
        {
            _quadratic_minpoly(poly, u, v, D1);
            _qqbar_set_minpoly_from_op(res, poly, x, y, op);
        }
    }
    else
    {
        /* Q(s1, s2) with s1^2 = D1, s2^2 = D2, s3 = s1 s2.
           Result u + v s1 + w s2 + r s3. */
        if (op == 0 || op == 1)
        {
            if (op == 0)
            {
                fmpq_add(u, a1, a2);
                fmpq_set(w, b2);
            }
            else
            {
                fmpq_sub(u, a1, a2);
                fmpq_neg(w, b2);
            }
            fmpq_set(v, b1);
            fmpq_zero(r);
        }
        else
        {
            if (op == 3)
            {
                fmpq_mul(t, a2, a2);
                fmpq_mul(s, b2, b2);
                fmpq_mul_fmpz(s, s, D2);
                fmpq_sub(t, t, s);
                fmpq_div(a2, a2, t);
                fmpq_div(b2, b2, t);
                fmpq_neg(b2, b2);
            }

            /* (a1 + b1 s1)(a2 + b2 s2) */
            fmpq_mul(u, a1, a2);
            fmpq_mul(v, b1, a2);
            fmpq_mul(w, a1, b2);
            fmpq_mul(r, b1, b2);
        }

        /* The element has degree 4 unless at most one of v, w, r is
           nonzero (each nontrivial automorphism of the biquadratic field
           changes at least one term if two of them are nonzero). */
        if (fmpq_is_zero(v) + fmpq_is_zero(w) + fmpq_is_zero(r) >= 2)
        {
            if (!fmpq_is_zero(v))
                _quadratic_minpoly(poly, u, v, D1);
            else if (!fmpq_is_zero(w))
                _quadratic_minpoly(poly, u, w, D2);
            else if (!fmpq_is_zero(r))
                _quadratic_minpoly(poly, u, r, D12);

            if (fmpq_is_zero(v) && fmpq_is_zero(w) && fmpq_is_zero(r))
                qqbar_set_fmpq(res, u);
            else
                _qqbar_set_minpoly_from_op(res, poly, x, y, op);
        }
        else
        {
            /* With Z = z - u, the minimal polynomial is
               Z^4 + (2 c0 - 4 D1 v^2) Z^2 - 8 D1 D2 v w r Z
                 + c0^2 - 4 D1 D2^2 w^2 r^2,
               c0 = v^2 D1 - w^2 D2 - r^2 D1 D2. */
            fmpq_poly_t p, zu;
            fmpq_t c0, c;

            fmpq_poly_init(p);
            fmpq_poly_init(zu);
            fmpq_init(c0);
            fmpq_init(c);

            fmpq_mul(c0, v, v);
            fmpq_mul_fmpz(c0, c0, D1);
            fmpq_mul(t, w, w);
            fmpq_mul_fmpz(t, t, D2);
            fmpq_sub(c0, c0, t);
            fmpq_mul(t, r, r);
            fmpq_mul_fmpz(t, t, D12);
            fmpq_sub(c0, c0, t);

            fmpq_poly_set_coeff_si(p, 4, 1);

            fmpq_mul(t, v, v);
            fmpq_mul_fmpz(t, t, D1);
            fmpq_mul_2exp(t, t, 2);
            fmpq_mul_2exp(c, c0, 1);
            fmpq_sub(c, c, t);
            fmpq_poly_set_coeff_fmpq(p, 2, c);

            fmpq_mul(c, v, w);
            fmpq_mul(c, c, r);
            fmpq_mul_fmpz(c, c, D12);
            fmpq_mul_si(c, c, -8);
            fmpq_poly_set_coeff_fmpq(p, 1, c);

            fmpq_mul(c, w, r);
            fmpq_mul(c, c, c);
            fmpq_mul_fmpz(c, c, D12);
            fmpq_mul_fmpz(c, c, D2);
            fmpq_mul_2exp(c, c, 2);
            fmpq_mul(t, c0, c0);
            fmpq_sub(c, t, c);
            fmpq_poly_set_coeff_fmpq(p, 0, c);

            /* substitute Z = z - u */
            fmpq_poly_set_coeff_si(zu, 1, 1);
            fmpq_neg(t, u);
            fmpq_poly_set_coeff_fmpq(zu, 0, t);
            fmpq_poly_compose(p, p, zu);

            _fmpz_poly_set_fmpq_poly_primitive(poly, p);
            _qqbar_set_minpoly_from_op(res, poly, x, y, op);

            fmpq_poly_clear(p);
            fmpq_poly_clear(zu);
            fmpq_clear(c0);
            fmpq_clear(c);
        }
    }

    fmpq_clear(a1); fmpq_clear(b1); fmpq_clear(a2); fmpq_clear(b2);
    fmpq_clear(u); fmpq_clear(v); fmpq_clear(w); fmpq_clear(t);
    fmpq_clear(s); fmpq_clear(r);
    fmpz_clear(D1); fmpz_clear(D2); fmpz_clear(D12); fmpz_clear(q);
    fmpz_poly_clear(poly);
}

/* Given P and linear L = (al0 + al1 t) z + (be0 + be1 t) with t^2 = D,
   compute H = U^2 - D V^2 where P(L) = U + t V. */
static void
_qqbar_norm_compose_linear(fmpz_poly_t H, const fmpz_poly_t P,
    const fmpq_t al0, const fmpq_t al1, const fmpq_t be0, const fmpq_t be1, const fmpz_t D)
{
    fmpq_poly_t U, V, L0, L1, T, S;
    slong i, d;

    d = fmpz_poly_degree(P);

    fmpq_poly_init(U);
    fmpq_poly_init(V);
    fmpq_poly_init(L0);
    fmpq_poly_init(L1);
    fmpq_poly_init(T);
    fmpq_poly_init(S);

    fmpq_poly_set_coeff_fmpq(L0, 1, al0);
    fmpq_poly_set_coeff_fmpq(L0, 0, be0);
    fmpq_poly_set_coeff_fmpq(L1, 1, al1);
    fmpq_poly_set_coeff_fmpq(L1, 0, be1);

    /* Horner: (U + tV) <- (U + tV)(L0 + t L1) + p_i */
    fmpq_poly_set_fmpz(U, P->coeffs + d);

    for (i = d - 1; i >= 0; i--)
    {
        /* T = U L0 + D V L1, S = U L1 + V L0 */
        fmpq_poly_mul(T, V, L1);
        fmpq_poly_scalar_mul_fmpz(T, T, D);
        fmpq_poly_mul(S, U, L0);
        fmpq_poly_add(T, T, S);

        fmpq_poly_mul(S, U, L1);
        fmpq_poly_mul(V, V, L0);
        fmpq_poly_add(V, V, S);

        fmpq_poly_swap(U, T);
        fmpq_poly_add_fmpz(U, U, P->coeffs + i);
    }

    fmpq_poly_mul(U, U, U);
    fmpq_poly_mul(V, V, V);
    fmpq_poly_scalar_mul_fmpz(V, V, D);
    fmpq_poly_sub(U, U, V);

    _fmpz_poly_set_fmpq_poly_primitive(H, U);

    fmpq_poly_clear(U);
    fmpq_poly_clear(V);
    fmpq_poly_clear(L0);
    fmpq_poly_clear(L1);
    fmpq_poly_clear(T);
    fmpq_poly_clear(S);
}

/* Whether the monic polynomial f has a root mod p. */
static int
_nmod_poly_has_root(const nmod_poly_t f, nmod_poly_t t, nmod_poly_t u, nmod_poly_t finv)
{
    nmod_poly_reverse(finv, f, f->length);
    nmod_poly_inv_series(finv, finv, f->length);
    nmod_poly_zero(u);
    nmod_poly_set_coeff_ui(u, 1, 1);
    nmod_poly_powmod_ui_binexp_preinv(t, u, f->mod.n, f, finv);
    nmod_poly_sub(t, t, u);
    nmod_poly_gcd(u, t, f);
    return nmod_poly_degree(u) >= 1;
}

/* Given monic squarefree f, determine whether f has an irreducible factor
   of odd degree. */
static int
_nmod_poly_has_odd_degree_factor(const nmod_poly_t f, nmod_poly_t t, nmod_poly_t u, nmod_poly_t finv)
{
    nmod_poly_factor_t fac;
    slong * degs;
    slong i, d;
    int result = 0;

    d = nmod_poly_degree(f);

    if (d % 2 == 1 || _nmod_poly_has_root(f, t, u, finv))
        return 1;

    nmod_poly_factor_init(fac);
    degs = flint_malloc(sizeof(slong) * (d + 1));
    nmod_poly_factor_distinct_deg(fac, f, &degs);

    for (i = 0; i < fac->num && !result; i++)
        if (degs[i] % 2 == 1)
            result = 1;

    flint_free(degs);
    nmod_poly_factor_clear(fac);

    return result;
}

#define CERTIFY_TRIES 16
#define CERTIFY_GIVE_UP 2

/* Attempt to certify that H is irreducible, as described above.
   Returns 1 if certified. */
static int
_qqbar_certify_norm_irreducible(const fmpz_poly_t H, const fmpz_poly_t P, const fmpz_t D)
{
    nmod_poly_t Pp, Hp, t, u, finv;
    ulong p, Dp;
    slong tries, evidence;
    int result = 0;

    nmod_poly_init(Pp, 2);
    nmod_poly_init(Hp, 2);
    nmod_poly_init(t, 2);
    nmod_poly_init(u, 2);
    nmod_poly_init(finv, 2);

    /* start with small primes (> 2 deg(H)) for speed */
    p = FLINT_MAX(1000, 2 * fmpz_poly_degree(H));
    evidence = 0;

    for (tries = 0; tries < CERTIFY_TRIES && !result && evidence < CERTIFY_GIVE_UP; tries++)
    {
        nmod_t mod;

        p = n_nextprime(p, 1);
        nmod_init(&mod, p);

        Dp = fmpz_fdiv_ui(D, p);
        if (Dp == 0)
            continue;

        nmod_poly_set_mod(Pp, mod);
        nmod_poly_set_mod(Hp, mod);
        nmod_poly_set_mod(t, mod);
        nmod_poly_set_mod(u, mod);
        nmod_poly_set_mod(finv, mod);

        fmpz_poly_get_nmod_poly(Pp, P);
        if (nmod_poly_degree(Pp) != fmpz_poly_degree(P))
            continue;

        nmod_poly_make_monic(Pp, Pp);

        if (n_jacobi_unsigned(Dp, p) == 1)
        {
            /* If P has a root mod p and D is a residue, this is (weak)
               evidence that sqrt(D) is in Q(x) (in which case D is a
               residue for all p where P has a root); we give up after a
               few such primes. */
            if (_nmod_poly_has_root(Pp, t, u, finv))
                evidence++;
            continue;
        }

        /* P squarefree mod p (p unramified in Q(x), not dividing the index) */
        nmod_poly_derivative(t, Pp);
        nmod_poly_gcd(u, Pp, t);
        if (!nmod_poly_is_one(u))
            continue;

        /* P must have an irreducible factor of odd degree f mod p; the
           nonresidue D is then not a square in F_{p^f} */
        if (!_nmod_poly_has_odd_degree_factor(Pp, t, u, finv))
            continue;

        /* H squarefree mod p with full degree */
        fmpz_poly_get_nmod_poly(Hp, H);
        if (nmod_poly_degree(Hp) != fmpz_poly_degree(H))
            continue;
        nmod_poly_derivative(t, Hp);
        nmod_poly_gcd(u, Hp, t);
        if (!nmod_poly_is_one(u))
            break;   /* if H mod p is not squarefree, H is likely not squarefree */

        result = 1;
    }

    nmod_poly_clear(Pp);
    nmod_poly_clear(Hp);
    nmod_poly_clear(t);
    nmod_poly_clear(u);
    nmod_poly_clear(finv);

    return result;
}

/* y quadratic, x of degree >= 3; op = 0, 1, 2, 3 */
static void
_qqbar_binary_op_general_quadratic(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)
{
    fmpq_t a, b, al0, al1, be0, be1, N;
    fmpz_t D;
    fmpz_poly_t H;
    fmpz_poly_factor_t fac;

    fmpq_init(a); fmpq_init(b);
    fmpq_init(al0); fmpq_init(al1); fmpq_init(be0); fmpq_init(be1);
    fmpq_init(N);
    fmpz_init(D);
    fmpz_poly_init(H);
    fmpz_poly_factor_init(fac);

    _qqbar_get_quadratic_parts(a, b, D, y);

    /* x = L(w) where w = x op y */
    if (op == 0)
    {
        /* x = w - y */
        fmpq_one(al0);
        fmpq_neg(be0, a);
        fmpq_neg(be1, b);
    }
    else if (op == 1)
    {
        /* x = w + y */
        fmpq_one(al0);
        fmpq_set(be0, a);
        fmpq_set(be1, b);
    }
    else if (op == 2)
    {
        /* x = w / y = w (a - b t) / (a^2 - b^2 D) */
        fmpq_mul(N, a, a);
        fmpq_mul(al1, b, b);
        fmpq_mul_fmpz(al1, al1, D);
        fmpq_sub(N, N, al1);
        fmpq_div(al0, a, N);
        fmpq_div(al1, b, N);
        fmpq_neg(al1, al1);
    }
    else
    {
        /* x = w y */
        fmpq_set(al0, a);
        fmpq_set(al1, b);
    }

    _qqbar_norm_compose_linear(H, QQBAR_POLY(x), al0, al1, be0, be1, D);

    if (_qqbar_certify_norm_irreducible(H, QQBAR_POLY(x), D))
    {
        _qqbar_set_minpoly_from_op(res, H, x, y, op);
    }
    else
    {
        fmpz_poly_factor(fac, H);
        _qqbar_binary_op_select_factor(res, fac, x, y, op);
    }

    fmpq_clear(a); fmpq_clear(b);
    fmpq_clear(al0); fmpq_clear(al1); fmpq_clear(be0); fmpq_clear(be1);
    fmpq_clear(N);
    fmpz_clear(D);
    fmpz_poly_clear(H);
    fmpz_poly_factor_clear(fac);
}

/* Assumes that neither x nor y is rational and that at least one of
   them is quadratic. */
void
_qqbar_binary_op_quadratic(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)
{
    slong dx, dy;

    dx = qqbar_degree(x);
    dy = qqbar_degree(y);

    if (dx == 2 && dy == 2)
    {
        _qqbar_binary_op_quadratic_quadratic(res, x, y, op);
    }
    else if (dy == 2)
    {
        _qqbar_binary_op_general_quadratic(res, x, y, op);
    }
    else
    {
        /* x quadratic */
        if (op == 0 || op == 2)
        {
            _qqbar_binary_op_general_quadratic(res, y, x, op);
        }
        else if (op == 1)
        {
            /* x - y = -(y - x) */
            _qqbar_binary_op_general_quadratic(res, y, x, 1);
            qqbar_neg(res, res);
        }
        else
        {
            /* x / y = 1 / (y / x) */
            _qqbar_binary_op_general_quadratic(res, y, x, 3);
            qqbar_inv(res, res);
        }
    }
}
