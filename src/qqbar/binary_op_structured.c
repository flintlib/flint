/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_vec.h"
#include "fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "fmpq.h"
#include "arb.h"
#include "acb.h"
#include "qqbar.h"
#include "impl.h"

/*
    Binary operations on elements of a common cyclotomic field.

    If x and y lie in Q(zeta_n), we write x = S(zeta_n), y = T(zeta_n)
    with S, T in Q[t] of degree < phi(n) (found using LLL, verified
    rigorously), and compute x op y by linear algebra in Q(zeta_n).
    Cyclotomic numbers usually have very small coefficients in the power
    basis of zeta_n, even when the power basis of x itself is badly
    conditioned (e.g. x = 1 + 2 zeta_n).

    To find n, we use that if Q(x) is contained in Q(zeta_n), then every
    prime p = 1 mod n splits completely in Q(x), so the minimal polynomial
    of x splits into distinct linear factors mod p (for p not dividing the
    leading coefficient or index). We enumerate n (n != 2 mod 4) with
    deg(x) | phi(n), deg(y) | phi(n), in order of increasing phi(n), and
    take the first n passing this test for a few primes. Since the minimal
    polynomial of a random algebraic number of degree d splits completely
    mod p with probability about 1/d!, this is a cheap filter.
*/

#define CYCLO_PHI_MAX 128
#define CYCLO_TEST_PRIMES 4

/* 1: all roots distinct and in F_p; 0: not split; -1: inconclusive */
static int
_nmod_poly_splits(const fmpz_poly_t P, ulong p)
{
    nmod_poly_t f, finv, t, u;
    nmod_t mod;
    slong d = fmpz_poly_degree(P);
    int result;

    nmod_init(&mod, p);
    nmod_poly_init_mod(f, mod);
    nmod_poly_init_mod(finv, mod);
    nmod_poly_init_mod(t, mod);
    nmod_poly_init_mod(u, mod);

    fmpz_poly_get_nmod_poly(f, P);

    if (nmod_poly_degree(f) != d)
    {
        result = -1;
    }
    else
    {
        nmod_poly_make_monic(f, f);
        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);
        nmod_poly_zero(u);
        nmod_poly_set_coeff_ui(u, 1, 1);
        nmod_poly_powmod_ui_binexp_preinv(t, u, p, f, finv);
        nmod_poly_sub(t, t, u);
        nmod_poly_gcd(u, t, f);

        if (nmod_poly_degree(u) == d)
        {
            result = 1;
        }
        else
        {
            /* inconclusive if f is not squarefree */
            nmod_poly_derivative(t, f);
            nmod_poly_gcd(u, t, f);
            result = nmod_poly_is_one(u) ? 0 : -1;
        }
    }

    nmod_poly_clear(f);
    nmod_poly_clear(finv);
    nmod_poly_clear(t);
    nmod_poly_clear(u);

    return result;
}

/* Test whether P and Q plausibly have roots in Q(zeta_n). */
static int
_qqbar_cyclotomic_split_test(const fmpz_poly_t P, const fmpz_poly_t Q, ulong n)
{
    ulong p, k;
    slong passed, tries;
    int r;

    passed = 0;

    for (k = 1 + 1000 / n, tries = 0; passed < CYCLO_TEST_PRIMES && tries < 200; k++)
    {
        p = k * n + 1;

        if (!n_is_prime(p))
            continue;

        tries++;

        r = _nmod_poly_splits(P, p);
        if (r == 0)
            return 0;
        if (r == -1)
            continue;

        if (Q != NULL)
        {
            r = _nmod_poly_splits(Q, p);
            if (r == 0)
                return 0;
            if (r == -1)
                continue;
        }

        passed++;
    }

    return passed == CYCLO_TEST_PRIMES;
}

/* Returns n such that x, y plausibly lie in Q(zeta_n), or 0. */
static ulong
_qqbar_find_common_cyclotomic_field(const qqbar_t x, const qqbar_t y)
{
    slong dx, dy, phi, i, j;
    ulong n, nmax;
    int same;
    ulong * tab;
    ulong result = 0;

    dx = qqbar_degree(x);
    dy = qqbar_degree(y);
    same = fmpz_poly_equal(QQBAR_POLY(x), QQBAR_POLY(y));

    if (FLINT_MAX(dx, dy) > CYCLO_PHI_MAX)
        return 0;

    /* phi(n) >= n / (e^gamma log log n + 3 / log log n) > 128 for n > 640 */
    nmax = 5 * CYCLO_PHI_MAX;

    /* sieve for phi */
    tab = flint_malloc(sizeof(ulong) * (nmax + 1));
    for (i = 0; i <= (slong) nmax; i++)
        tab[i] = i;
    for (i = 2; i <= (slong) nmax; i++)
        if (tab[i] == (ulong) i)
            for (j = i; j <= (slong) nmax; j += i)
                tab[j] = (tab[j] / i) * (i - 1);

    for (phi = FLINT_MAX(dx, dy); phi <= CYCLO_PHI_MAX && result == 0; phi++)
    {
        if (phi % dx != 0 || phi % dy != 0)
            continue;

        for (n = 3; n <= nmax && result == 0; n++)
        {
            if (n % 4 == 2 || tab[n] != (ulong) phi)
                continue;

            if (_qqbar_cyclotomic_split_test(QQBAR_POLY(x), same ? NULL : QQBAR_POLY(y), n))
                result = n;
        }
    }

    flint_free(tab);
    return result;
}

/* Express x in Q(zeta_n) (real = 0) or Q(zeta_n + zeta_n^-1) (real = 1),
   as a polynomial S in gen = zeta_n or gen = 2 cos(2 pi / n), using LLL at
   the given precision. For the real case, the lattice is built from the
   basis 1, 2 cos(2 pi k / n), 1 <= k < m, which gives small coefficients
   for cyclotomic numbers (unlike the power basis of 2 cos(2 pi / n)). */
static int
_qqbar_cyclotomic_express(fmpq_poly_t S, const qqbar_t gen, ulong n, int real, const qqbar_t x, slong prec)
{
    slong m, k;
    acb_ptr vec;
    fmpz * rel;
    fmpq_t q;
    int found = 0;

    m = qqbar_degree(gen);

    vec = _acb_vec_init(m + 1);
    rel = _fmpz_vec_init(m + 1);

    fmpq_init(q);

    for (k = 0; k < m; k++)
    {
        fmpz_set_ui(fmpq_numref(q), 2 * k);
        fmpz_set_ui(fmpq_denref(q), n);
        fmpq_canonicalise(q);

        if (real)
        {
            if (k == 0)
                acb_one(vec);
            else
            {
                arb_cos_pi_fmpq(acb_realref(vec + k), q, prec);
                arb_mul_2exp_si(acb_realref(vec + k), acb_realref(vec + k), 1);
            }
        }
        else
        {
            arb_sin_cos_pi_fmpq(acb_imagref(vec + k), acb_realref(vec + k), q, prec);
        }
    }

    qqbar_enclosure_raw(vec + m, x, prec);
    if (real)
        arb_zero(acb_imagref(vec + m));

    if (_qqbar_acb_lindep(rel, vec, m + 1, 1, prec) && !fmpz_is_zero(rel + m))
    {
        fmpq_poly_zero(S);

        if (real)
        {
            fmpz_poly_t C0, C1, C2, T;
            fmpz_poly_init(C0);
            fmpz_poly_init(C1);
            fmpz_poly_init(C2);
            fmpz_poly_init(T);

            /* C_0 = 2, C_1 = t, C_{k+1} = t C_k - C_{k-1} */
            fmpz_poly_set_si(C0, 2);
            fmpz_poly_set_coeff_si(C1, 1, 1);
            fmpz_poly_set_fmpz(T, rel);

            for (k = 1; k < m; k++)
            {
                fmpz_poly_scalar_addmul_fmpz(T, C1, rel + k);
                fmpz_poly_shift_left(C2, C1, 1);
                fmpz_poly_sub(C2, C2, C0);
                fmpz_poly_swap(C0, C1);
                fmpz_poly_swap(C1, C2);
            }

            fmpq_poly_set_fmpz_poly(S, T);

            fmpz_poly_clear(C0);
            fmpz_poly_clear(C1);
            fmpz_poly_clear(C2);
            fmpz_poly_clear(T);
        }
        else
        {
            fmpz_poly_t T;
            T->coeffs = rel;
            T->length = m;
            T->alloc = m;
            _fmpz_poly_normalise(T);
            fmpq_poly_set_fmpz_poly(S, T);
        }

        fmpq_poly_scalar_div_fmpz(S, S, rel + m);
        fmpq_poly_neg(S, S);

        found = qqbar_equal_fmpq_poly_val(x, S, gen);
    }

    _acb_vec_clear(vec, m + 1);
    _fmpz_vec_clear(rel, m + 1);
    fmpq_clear(q);

    return found;
}

/* Compute h(gen) where h = A op B (polynomials in gen); the op codes
   are as for _qqbar_binary_op_subfield. */
static void
_qqbar_evaluate_op(qqbar_t res, fmpq_poly_t A, const fmpq_poly_t B, const qqbar_t gen, int op)
{
    fmpq_poly_t F, G, U, V;

    fmpq_poly_init(F);
    fmpq_poly_init(G);
    fmpq_poly_init(U);
    fmpq_poly_init(V);

    fmpq_poly_set_fmpz_poly(F, QQBAR_POLY(gen));

    if (op == 0)
        fmpq_poly_add(A, A, B);
    else if (op == 1)
        fmpq_poly_sub(A, A, B);
    else if (op == 2)
    {
        fmpq_poly_mul(A, A, B);
        fmpq_poly_rem(A, A, F);
    }
    else if (op == 4 || op == 6)
    {
        /* (x - y)^2 or -(x - y)^2 / 4 */
        fmpq_poly_sub(A, A, B);
        fmpq_poly_mul(A, A, A);
        fmpq_poly_rem(A, A, F);
        if (op == 6)
            fmpq_poly_scalar_div_si(A, A, -4);
    }
    else
    {
        fmpq_poly_xgcd(G, U, V, B, F);
        fmpq_poly_mul(A, A, U);
        fmpq_poly_rem(A, A, F);
    }

    qqbar_evaluate_fmpq_poly(res, A, gen);

    fmpq_poly_clear(F);
    fmpq_poly_clear(G);
    fmpq_poly_clear(U);
    fmpq_poly_clear(V);
}

/* Numerical value of x op y. */
static void
_qqbar_op_enclosure(acb_t w, const qqbar_t x, const qqbar_t y, int op, slong prec)
{
    acb_t a, b;
    acb_init(a);
    acb_init(b);

    qqbar_enclosure_raw(a, x, prec);
    qqbar_enclosure_raw(b, y, prec);

    if (op == 0)
        acb_add(w, a, b, prec);
    else if (op == 1)
        acb_sub(w, a, b, prec);
    else if (op == 2)
        acb_mul(w, a, b, prec);
    else
        acb_div(w, a, b, prec);

    acb_clear(a);
    acb_clear(b);
}

/*
    Attempt to find the minimal polynomial of w = x op y by guessing an
    annihilating polynomial g of degree m using LLL, and proving that
    g(w) = 0 without computing a resultant.

    Let f be the minimal polynomial of w (primitive in Z[z]) and suppose
    that g(w) != 0. Then g(v) != 0 for all conjugates v of w, and
    Res(f, g) = lc(f)^m prod_{f(v) = 0} g(v) is a nonzero integer.
    Since f divides the integer polynomial
    H = L prod_{a,b} (z - (a op b)) with L = lc(P)^dy lc(Q')^dx (where
    Q' = Q, or the reversal of Q for division), lc(f) divides L, and
    deg(f) <= N = dx dy. If R bounds the absolute values of all roots of
    H, we have |g(v)| <= |g|_1 max(1, R)^m for all conjugates v, so
    |g(w)| >= |L|^(-m) (|g|_1 max(1, R)^m)^(-(N-1)). Hence if a numerical
    evaluation shows that |g(w)| is smaller than this bound, we must have
    g(w) = 0. The minimal polynomial is then the irreducible factor of g
    which vanishes at w.
*/
static int
_qqbar_guess_op(qqbar_t res, const qqbar_t x, const qqbar_t y, int op, slong m, slong prec)
{
    acb_ptr vec;
    acb_t w;
    fmpz * rel;
    fmpz_poly_t g, Q2;
    fmpz_t RP, RQ, R;
    slong i, dx, dy, N, B, lR, lg, llc;
    int success = 0;

    dx = qqbar_degree(x);
    dy = qqbar_degree(y);
    N = dx * dy;

    vec = _acb_vec_init(m + 1);
    rel = _fmpz_vec_init(m + 1);
    acb_init(w);
    fmpz_poly_init(g);

    _qqbar_op_enclosure(w, x, y, op, prec);
    _acb_vec_set_powers(vec, w, m + 1, prec);

    if (_qqbar_acb_lindep(rel, vec, m + 1, 1, prec))
    {
        for (i = 0; i <= m; i++)
            fmpz_poly_set_coeff_fmpz(g, i, rel + i);

        fmpz_poly_primitive_part(g, g);
        if (g->length > 0 && fmpz_sgn(fmpz_poly_lead(g)) < 0)
            fmpz_poly_neg(g, g);
    }

    if (fmpz_poly_degree(g) >= 1)
    {
        acb_t gw;
        mag_t mm;

        acb_init(gw);
        mag_init(mm);
        fmpz_poly_init(Q2);
        fmpz_init(RP);
        fmpz_init(RQ);
        fmpz_init(R);

        m = fmpz_poly_degree(g);

        fmpz_poly_set(Q2, QQBAR_POLY(y));
        if (op == 3)
            fmpz_poly_reverse(Q2, Q2, dy + 1);

        fmpz_poly_bound_roots(RP, QQBAR_POLY(x));
        fmpz_poly_bound_roots(RQ, Q2);

        if (op == 2 || op == 3)
            fmpz_mul(R, RP, RQ);
        else
            fmpz_add(R, RP, RQ);

        lR = fmpz_bits(R);
        llc = dy * fmpz_bits(QQBAR_COEFFS(x) + dx) + dx * fmpz_bits(Q2->coeffs + dy);
        lg = FLINT_ABS(fmpz_poly_max_bits(g)) + FLINT_BIT_COUNT(m + 1);

        /* log2 of the bound (upper estimate, all terms rounded up) */
        B = m * llc + (N - 1) * (lg + m * lR) + 1;

        /* cheap sanity check at a higher precision before the
           expensive high-precision evaluation: a spurious relation
           found by LLL will typically only vanish to about prec bits */
        _qqbar_op_enclosure(w, x, y, op, 3 * prec);
        arb_fmpz_poly_evaluate_acb(gw, g, w, 3 * prec);

        if (acb_contains_zero(gw))
        {
            slong vprec = B + 64 + m * (lR + 2);

            _qqbar_op_enclosure(w, x, y, op, vprec);
            arb_fmpz_poly_evaluate_acb(gw, g, w, vprec);
            acb_get_mag(mm, gw);

            if (mag_cmp_2exp_si(mm, -B) < 0)
            {
                /* g(w) = 0 is proved */
                fmpz_poly_factor_t fac;
                fmpz_poly_factor_init(fac);
                fmpz_poly_factor(fac, g);
                _qqbar_binary_op_select_factor(res, fac, x, y, op);
                fmpz_poly_factor_clear(fac);
                success = 1;
            }
        }

        acb_clear(gw);
        mag_clear(mm);
        fmpz_poly_clear(Q2);
        fmpz_clear(RP);
        fmpz_clear(RQ);
        fmpz_clear(R);
    }

    _acb_vec_clear(vec, m + 1);
    _fmpz_vec_clear(rel, m + 1);
    acb_clear(w);
    fmpz_poly_clear(g);

    return success;
}

/* LLL cost units: roughly 1.5e-8 seconds. */
static double
_lll_cost(slong dim, slong prec)
{
    /* the second term accounts for numerical evaluation overhead */
    return (pow(dim, 2.4) + 20.0 * dim) * prec;
}

/* Precision cap for a stream (relations with more than 1024 bits per
   coefficient are unlikely to be found within the budget anyway, and this
   prevents cheap low-dimensional streams from running to huge precision). */
#define STREAM_MAX_PREC(dim) (1024 * (dim))

#define STREAM_SUBFIELD 0
#define STREAM_CYCLOTOMIC 1
#define STREAM_GUESS 2
#define MAX_STREAMS 32
#define GUESS_MAX_DEGREE 64

/*
    Try to compute x op y without factoring a composed polynomial, using
    one of the following methods, each of which requires an LLL search that
    may fail:

    * Subfield: express y in Q(x) or x in Q(y) and use linear algebra in
      that field (see binary_op_subfield.c).
    * Cyclotomic: express x and y in Q(zeta_n) or its real subfield and use
      linear algebra there.
    * Guess: find the minimal polynomial of the result directly and verify
      it rigorously (see _qqbar_guess_op); we try degrees m = d, d/2, d/4
      where d = max(dx, dy) (if the subfield test suggests that one field
      contains the other, but there is no common cyclotomic field).

    We maintain a list of "streams" of LLL attempts (each with a lattice
    dimension and a doubling sequence of precisions) and always perform the
    cheapest next attempt, subject to a total budget proportional to the
    estimated cost of the fallback algorithm.
*/
int
_qqbar_binary_op_structured(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)
{
    slong dx, dy, d, N, bits, phi, i, j, best, num;
    int swapped, sub_ok, cyc_ok, real, success, have_x;
    double budget, spent, cost, best_cost;
    ulong n;
    qqbar_t gen;
    fmpq_poly_t S, T;
    int stype[MAX_STREAMS];
    slong sdim[MAX_STREAMS], sprec[MAX_STREAMS];
    int sactive[MAX_STREAMS];

    dx = qqbar_degree(x);
    dy = qqbar_degree(y);
    d = FLINT_MAX(dx, dy);

    if (dx < 2 || dy < 2)
        return 0;

    if (fmpz_poly_equal(QQBAR_POLY(x), QQBAR_POLY(y)))
        N = dx * (dx - 1) / 2;
    else
        N = dx * dy;

    if (N < 36)
        return 0;

    bits = dy * FLINT_ABS(fmpz_poly_max_bits(QQBAR_POLY(x)))
         + dx * FLINT_ABS(fmpz_poly_max_bits(QQBAR_POLY(y))) + dx * dy;

    /* 60% of the estimated fallback cost 1e-7 N B seconds, in units of
       1.5e-8 seconds (see _lll_cost) */
    budget = 4.0 * (double) N * (double) bits;
    spent = 0;

    /* Subfield method */
    swapped = -1;
    if (dx == dy)
        swapped = FLINT_ABS(fmpz_poly_max_bits(QQBAR_POLY(y))) <
                  FLINT_ABS(fmpz_poly_max_bits(QQBAR_POLY(x)));
    else if (dx % dy == 0)
        swapped = 0;
    else if (dy % dx == 0)
        swapped = 1;

    sub_ok = 0;
    if (swapped != -1)
    {
        const qqbar_struct * gen1 = swapped ? y : x;
        const qqbar_struct * other1 = swapped ? x : y;
        sub_ok = _qqbar_subfield_modular_test(QQBAR_POLY(gen1), QQBAR_POLY(other1),
            fmpz_poly_equal(QQBAR_POLY(gen1), QQBAR_POLY(other1)));
    }

    /* Cyclotomic method */
    n = 0;
    cyc_ok = 0;
    real = 0;
    phi = 0;

    if (dx * dy >= 64)
    {
        n = _qqbar_find_common_cyclotomic_field(x, y);

        if (n != 0)
        {
            cyc_ok = 1;
            real = (qqbar_sgn_im(x) == 0) && (qqbar_sgn_im(y) == 0);
            phi = n_euler_phi(n);
        }
    }

    if (!sub_ok && !cyc_ok)
        return 0;

    /* Without cyclotomic structure, relations with small coefficients
       are less likely (e.g. towers of nested square roots, where y is often
       in Q(x) but only with huge coefficients), so we use a smaller
       budget. */
    if (!cyc_ok)
        budget *= 0.5;

    num = 0;

    if (sub_ok)
    {
        stype[num] = STREAM_SUBFIELD;
        sdim[num] = d + 1;
        num++;
    }

    if (cyc_ok)
    {
        stype[num] = STREAM_CYCLOTOMIC;
        sdim[num] = (real ? phi / 2 : phi) + 1;
        num++;
    }

    /* Guess streams: only when there is no cyclotomic structure (the
       cyclotomic method is more reliable), for degrees m = d, d/2, d/4
       (typical for towers of quadratic extensions, where y is in Q(x) and
       the result often generates a proper subfield). */
    if (op <= 3 && !cyc_ok)
    {
        for (j = 1; j <= 4; j *= 2)
        {
            if (d % j == 0 && d / j < N && d / j <= GUESS_MAX_DEGREE && d / j >= 2)
            {
                stype[num] = STREAM_GUESS;
                sdim[num] = d / j + 1;
                num++;
            }
        }
    }

    for (i = 0; i < num; i++)
    {
        sprec[i] = 64;
        sactive[i] = 1;
    }

    qqbar_init(gen);
    fmpq_poly_init(S);
    fmpq_poly_init(T);

    if (cyc_ok)
    {
        if (real)
        {
            qqbar_cos_pi(gen, 2, n);
            qqbar_mul_2exp_si(gen, gen, 1);
        }
        else
        {
            qqbar_root_of_unity(gen, 1, n);
        }
    }

    have_x = 0;
    success = 0;

    while (!success)
    {
        best = -1;
        best_cost = 0;

        for (i = 0; i < num; i++)
        {
            if (sactive[i] && sprec[i] > STREAM_MAX_PREC(sdim[i]))
                sactive[i] = 0;

            if (!sactive[i])
                continue;

            cost = _lll_cost(sdim[i], sprec[i]);
            if (stype[i] == STREAM_CYCLOTOMIC && !have_x)
                cost *= 2;

            if (best == -1 || cost < best_cost)
            {
                best = i;
                best_cost = cost;
            }
        }

        if (best == -1)
            break;

        if (spent + best_cost > budget)
            break;

        spent += best_cost;

        if (stype[best] == STREAM_SUBFIELD)
        {
            if (_qqbar_binary_op_subfield(res, x, y, op, swapped, sprec[best], sprec[best], 0) == 1)
                success = 1;
        }
        else if (stype[best] == STREAM_CYCLOTOMIC)
        {
            if (!have_x)
                have_x = _qqbar_cyclotomic_express(S, gen, n, real, x, sprec[best]);

            if (have_x && _qqbar_cyclotomic_express(T, gen, n, real, y, sprec[best]))
            {
                _qqbar_evaluate_op(res, S, T, gen, op);
                success = 1;
            }
        }
        else
        {
            success = _qqbar_guess_op(res, x, y, op, sdim[best] - 1, sprec[best]);
        }

        sprec[best] *= 2;
    }

    qqbar_clear(gen);
    fmpq_poly_clear(S);
    fmpq_poly_clear(T);

    return success;
}
