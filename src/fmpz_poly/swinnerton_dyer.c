/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*
    Exact computation of the Swinnerton-Dyer polynomial

        S_n(x) = prod_{e in {+1,-1}^n} (x - e_1 sqrt(p_1) - ... - e_n sqrt(p_n))

    by iterated resultants, without any floating-point arithmetic.

    Let Q_0(x) = x and Q_i(x) = Q_{i-1}(x - sqrt(p_i)) * Q_{i-1}(x + sqrt(p_i)),
    so that Q_n = S_n.  For i >= 1 the roots of Q_i are symmetric under
    negation, so Q_i(x) = G_i(x^2) and we work with G_i, of half the length.

    With u = x^2 and tau = 2 sqrt(p) x, so that tau^2 = 4 p u,

        Q(x + sqrt(p)) Q(x - sqrt(p)) = G(u + p + tau) G(u + p - tau).

    Writing G(u + p + tau) = A(u) + tau B(u) in Z[u][tau]/(tau^2 - 4 p u),
    the conjugation tau -> -tau gives the second factor, hence

        G_new(u) = A(u)^2 - 4 p u B(u)^2.

    The shift by p + tau is done by divide and conquer, splitting
    G = G_lo + u^h G_hi:

        G(u + p + tau) = G_lo(u + p + tau) + (u + p + tau)^h G_hi(u + p + tau).

    Since u + p + tau = (sqrt(u) + sqrt(p))^2, the power (u + p + tau)^h
    = P + tau R follows from the binomial theorem:

        P = sum_{k even} C(2h,k)     p^(k/2)     u^(h - k/2)
        R = sum_{k odd}  C(2h,k) / 2 p^((k-1)/2) u^(h - (k+1)/2)

    where C(2h,k) is even for odd k, so R has integer coefficients.
*/

/* Below this length, Horner's rule in Z[u][tau]/(tau^2 - 4pu) is used. */
#define SHIFT_CUTOFF 32

/* A + tau B = G(u + p + tau), G = sum_{j<len} c[j] u^j, by Horner's rule. */
static void
_shift_horner(fmpz_poly_t A, fmpz_poly_t B, const fmpz * c, slong len,
    const fmpz_t p, const fmpz_t p4)
{
    fmpz_poly_t T, U;
    fmpz_t t;
    slong j;

    fmpz_poly_init(T);
    fmpz_poly_init(U);
    fmpz_init(t);

    fmpz_poly_set_fmpz(A, c + len - 1);
    fmpz_poly_zero(B);

    for (j = len - 2; j >= 0; j--)
    {
        /* (A + tau B)(u + p + tau) + c_j
              = ((u + p) A + 4pu B + c_j) + tau (A + (u + p) B) */
        fmpz_poly_shift_left(T, B, 1);
        fmpz_poly_shift_left(U, A, 1);
        fmpz_poly_scalar_addmul_fmpz(U, A, p);
        fmpz_poly_scalar_addmul_fmpz(U, T, p4);

        fmpz_poly_get_coeff_fmpz(t, U, 0);
        fmpz_add(t, t, c + j);
        fmpz_poly_set_coeff_fmpz(U, 0, t);

        fmpz_poly_scalar_addmul_fmpz(T, B, p);
        fmpz_poly_add(B, T, A);
        fmpz_poly_swap(A, U);
    }

    fmpz_poly_clear(T);
    fmpz_poly_clear(U);
    fmpz_clear(t);
}

/* P + tau R = (u + p + tau)^h mod (tau^2 - 4pu), for h >= 1. */
static void
_binomial_pair(fmpz_poly_t P, fmpz_poly_t R, slong h, const fmpz_t p)
{
    fmpz_t b, pw;       /* b = C(2h,k),  pw = p^floor(k/2) */
    slong k;

    fmpz_poly_fit_length(P, h + 1);
    fmpz_poly_fit_length(R, h);

    fmpz_init_set_ui(b, 1);
    fmpz_init_set_ui(pw, 1);

    fmpz_one(P->coeffs + h);

    for (k = 1; k <= 2 * h; k++)
    {
        fmpz_mul_ui(b, b, 2 * h - k + 1);
        fmpz_divexact_ui(b, b, k);

        if (k % 2 == 0)
        {
            fmpz_mul(pw, pw, p);
            fmpz_mul(P->coeffs + (h - k / 2), b, pw);
        }
        else
        {
            fmpz_mul(R->coeffs + (h - (k + 1) / 2), b, pw);
            fmpz_tdiv_q_2exp(R->coeffs + (h - (k + 1) / 2),
                             R->coeffs + (h - (k + 1) / 2), 1);
        }
    }

    _fmpz_poly_set_length(P, h + 1);    /* leading coefficient 1 */
    _fmpz_poly_set_length(R, h);        /* leading coefficient h */

    fmpz_clear(b);
    fmpz_clear(pw);
}

/* A + tau B = G(u + p + tau) mod (tau^2 - 4pu), G = sum_{j<len} c[j] u^j. */
static void
_shift(fmpz_poly_t A, fmpz_poly_t B, const fmpz * c, slong len,
    const fmpz_t p, const fmpz_t p4)
{
    fmpz_poly_t A1, B1, P, R, X, Y, Z;
    slong h;

    if (len <= SHIFT_CUTOFF)
    {
        _shift_horner(A, B, c, len, p, p4);
        return;
    }

    h = len / 2;

    fmpz_poly_init(A1);
    fmpz_poly_init(B1);
    fmpz_poly_init(P);
    fmpz_poly_init(R);
    fmpz_poly_init(X);
    fmpz_poly_init(Y);
    fmpz_poly_init(Z);

    _shift(A, B, c, h, p, p4);                  /* G_lo(u + p + tau) */
    _shift(A1, B1, c + h, len - h, p, p4);      /* G_hi(u + p + tau) */
    _binomial_pair(P, R, h, p);                 /* (u + p + tau)^h   */

    /* (P + tau R)(A1 + tau B1), with three multiplications (Karatsuba) */
    fmpz_poly_mul(X, P, A1);
    fmpz_poly_mul(Y, R, B1);
    fmpz_poly_add(P, P, R);
    fmpz_poly_add(A1, A1, B1);
    fmpz_poly_mul(Z, P, A1);
    fmpz_poly_sub(Z, Z, X);
    fmpz_poly_sub(Z, Z, Y);                     /* Z = P B1 + R A1     */
    fmpz_poly_shift_left(Y, Y, 1);
    fmpz_poly_scalar_addmul_fmpz(X, Y, p4);     /* X = P A1 + 4pu R B1 */

    fmpz_poly_add(A, A, X);
    fmpz_poly_add(B, B, Z);

    fmpz_poly_clear(A1);
    fmpz_poly_clear(B1);
    fmpz_poly_clear(P);
    fmpz_poly_clear(R);
    fmpz_poly_clear(X);
    fmpz_poly_clear(Y);
    fmpz_poly_clear(Z);
}

void
_fmpz_poly_swinnerton_dyer(fmpz * T, ulong n)
{
    slong j, N = (WORD(1) << n);

    _fmpz_vec_zero(T, N + 1);

    if (n == 0)
    {
        fmpz_one(T + 1);
    }
    else
    {
        fmpz_poly_t G, A, B, S;
        fmpz_t p, p4;
        ulong i, prime;

        fmpz_poly_init(G);
        fmpz_poly_init(A);
        fmpz_poly_init(B);
        fmpz_poly_init(S);
        fmpz_init(p);
        fmpz_init(p4);

        /* G_1 = u - 2 */
        fmpz_poly_set_coeff_si(G, 0, -2);
        fmpz_poly_set_coeff_si(G, 1, 1);

        prime = 2;
        for (i = 1; i < n; i++)
        {
            prime = n_nextprime(prime, 0);
            fmpz_set_ui(p, prime);
            fmpz_mul_2exp(p4, p, 2);

            _shift(A, B, G->coeffs, G->length, p, p4);

            fmpz_poly_sqr(G, A);
            fmpz_poly_sqr(S, B);
            fmpz_poly_shift_left(S, S, 1);
            fmpz_poly_scalar_submul_fmpz(G, S, p4);  /* G = A^2 - 4pu B^2 */
        }

        for (j = 0; j < G->length; j++)
            fmpz_swap(T + 2 * j, G->coeffs + j);

        fmpz_poly_clear(G);
        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
        fmpz_poly_clear(S);
        fmpz_clear(p);
        fmpz_clear(p4);
    }
}

void
fmpz_poly_swinnerton_dyer(fmpz_poly_t poly, ulong n)
{
    slong N = (WORD(1) << n);
    fmpz_poly_fit_length(poly, N + 1);
    _fmpz_poly_swinnerton_dyer(poly->coeffs, n);
    _fmpz_poly_set_length(poly, N + 1);
}
