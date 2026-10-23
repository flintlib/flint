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
#include "fmpq_poly.h"
#include "acb.h"
#include "qqbar.h"
#include "impl.h"

/*
    Binary operations x op y where one operand lies in the field generated
    by the other.

    If y lies in K = Q(x), say y = R(x) with R in Q[t], deg(R) < deg(x),
    then x op y = h(x) where h = t op R(t) mod minpoly(x), and the
    minimal polynomial of h(x) can be computed by linear algebra in K
    (see _qqbar_evaluate_fmpq_poly) instead of by factoring a composed
    polynomial of degree deg(x) deg(y). The latter is very expensive when
    the Galois group is small (e.g. for elements of cyclotomic fields),
    since the composed polynomial then splits into many factors both over Z
    and modulo every prime, making Zassenhaus/van Hoeij recombination
    costly.

    Finding R is done heuristically using LLL, and the result is verified
    rigorously (as in qqbar_express_in_field). Since a failed LLL search is
    not free, we first apply a cheap modular test giving a necessary
    condition for y to lie in Q(x).

    Let p be a prime not dividing the leading coefficients and
    discriminants of P = minpoly(x) and Q = minpoly(y). Then p is
    unramified in both fields and (Dedekind) the degrees of the irreducible
    factors of P mod p are the residue degrees f of the primes of Q(x)
    above p. If y = R(x), then reducing modulo such a prime shows that
    Q has a root in F_{p^f}. If moreover P = Q and y != x, the images of x
    and y in F_{p^f} are distinct roots of Q (since Q is squarefree mod p),
    so Q must have at least two roots in F_{p^f}. The number of roots
    of Q in F_{p^f} is the sum of deg(g) over irreducible factors g of
    Q mod p with deg(g) | f.
*/

/* Number of primes to use in the first (cheap, residue degree 1 only) and
   second (all residue degrees) phase of the modular test. In the first
   phase, we use a number of primes proportional to the degree: when
   Q(x) is Galois of degree d, only a fraction 1/d of all primes split
   completely in Q(x), and if Q(x, y) has degree 2d, a fraction 1/(2d) of
   all primes certify that y is not in Q(x) (e.g. sums of square roots). */
#define SUBFIELD_TEST_PRIMES_1(d) FLINT_MAX(32, FLINT_MIN(4 * (d), 400))
#define SUBFIELD_TEST_PRIMES_2 8

/* Number of distinct roots of the monic polynomial f over F_p. */
static slong
_qqbar_nmod_poly_num_roots(const nmod_poly_t f, nmod_poly_t t, nmod_poly_t u)
{
    nmod_poly_t finv;
    slong r;

    if (nmod_poly_degree(f) == 1)
        return 1;

    nmod_poly_init_mod(finv, f->mod);
    nmod_poly_reverse(finv, f, f->length);
    nmod_poly_inv_series(finv, finv, f->length);

    /* t = x^p - x mod f */
    nmod_poly_zero(u);
    nmod_poly_set_coeff_ui(u, 1, 1);
    nmod_poly_powmod_ui_binexp_preinv(t, u, f->mod.n, f, finv);
    nmod_poly_sub(t, t, u);
    nmod_poly_gcd(u, t, f);

    r = nmod_poly_degree(u);

    nmod_poly_clear(finv);
    return r;
}

/* Given monic squarefree poly, write cnt[g] = number of irreducible factors
   of degree g, for 0 <= g <= deg. */
static void
_qqbar_nmod_poly_factor_degree_counts(slong * cnt, const nmod_poly_t poly)
{
    nmod_poly_factor_t fac;
    slong * degs;
    slong i, n;

    n = nmod_poly_degree(poly);

    for (i = 0; i <= n; i++)
        cnt[i] = 0;

    nmod_poly_factor_init(fac);
    degs = flint_malloc(sizeof(slong) * (n + 1));

    nmod_poly_factor_distinct_deg(fac, poly, &degs);

    for (i = 0; i < fac->num; i++)
        cnt[degs[i]] += nmod_poly_degree(fac->p + i) / degs[i];

    flint_free(degs);
    nmod_poly_factor_clear(fac);
}

static int
_qqbar_nmod_poly_set_fmpz_poly_monic(nmod_poly_t res, const fmpz_poly_t poly, nmod_poly_t tmp1, nmod_poly_t tmp2, int check_squarefree)
{
    fmpz_poly_get_nmod_poly(res, poly);

    if (nmod_poly_degree(res) != fmpz_poly_degree(poly))
        return 0;

    nmod_poly_make_monic(res, res);

    if (!check_squarefree)
        return 1;

    nmod_poly_derivative(tmp1, res);
    nmod_poly_gcd(tmp2, res, tmp1);

    return nmod_poly_is_one(tmp2);
}

/* Returns 0 if y = R(x) (with y != x if P = Q) is impossible, where
   P = minpoly(x), Q = minpoly(y); returns 1 if the test is inconclusive.

   In the first phase, we only use the condition for residue degree f = 1,
   i.e. if P has a root mod p, then Q must have a root (two roots if P = Q)
   mod p. This only requires one modular exponentiation with P per prime.
   In the second phase, we use distinct-degree factorization of P and Q
   to test the condition for all residue degrees f. This is needed, for
   instance, when Q(x) is a Galois field of large degree: then P almost
   never has a root mod p, but it splits into factors of equal degree f,
   and Q must then have a factor of degree dividing f. */
int
_qqbar_subfield_modular_test(const fmpz_poly_t P, const fmpz_poly_t Q, int same)
{
    nmod_poly_t Pp, Qp, t, u;
    slong tries, passed, rP, rQ, dP, dQ, f, g, N, num_primes, informative = 0;
    slong * cntP, * cntQ;
    ulong p;
    int result = 1, phase;

    dP = fmpz_poly_degree(P);
    dQ = fmpz_poly_degree(Q);

    cntP = flint_malloc(sizeof(slong) * (dP + 1));
    cntQ = flint_malloc(sizeof(slong) * (dQ + 1));

    p = UWORD(1) << 12;
    nmod_poly_init(Pp, 2);
    nmod_poly_init(Qp, 2);
    nmod_poly_init(t, 2);
    nmod_poly_init(u, 2);

    for (phase = 1; phase <= 2 && result; phase++)
    {
        /* The second phase is relatively expensive; skip it if the
           first phase already found enough primes of residue degree 1. */
        if (phase == 2 && informative >= SUBFIELD_TEST_PRIMES_2)
            break;

        num_primes = (phase == 1) ? SUBFIELD_TEST_PRIMES_1(dP) : SUBFIELD_TEST_PRIMES_2;

        for (tries = passed = 0; passed < num_primes && tries < 4 * num_primes + 100 && result; tries++)
        {
            nmod_t mod;

            p = n_nextprime(p, 1);
            nmod_init(&mod, p);

            nmod_poly_set_mod(Pp, mod);
            nmod_poly_set_mod(Qp, mod);
            nmod_poly_set_mod(t, mod);
            nmod_poly_set_mod(u, mod);

            /* In the first phase, we skip the squarefreeness check for
               speed. The test may then give a wrong answer if p happens to
               be ramified or divide the index of Z[x] in the ring of
               integers, but this is harmless: we only lose the chance to
               use the fast algorithm. */
            if (!_qqbar_nmod_poly_set_fmpz_poly_monic(Pp, P, t, u, phase == 2))
                continue;

            if (!same && !_qqbar_nmod_poly_set_fmpz_poly_monic(Qp, Q, t, u, phase == 2))
                continue;

            if (phase == 1)
            {
                rP = _qqbar_nmod_poly_num_roots(Pp, t, u);

                if (rP >= 1)
                    informative++;

                if (same)
                {
                    if (rP == 1)
                        result = 0;
                }
                else if (rP >= 1)
                {
                    rQ = _qqbar_nmod_poly_num_roots(Qp, t, u);
                    if (rQ == 0)
                        result = 0;
                }
            }
            else
            {
                _qqbar_nmod_poly_factor_degree_counts(cntP, Pp);

                if (same)
                {
                    for (g = 0; g <= dQ; g++)
                        cntQ[g] = cntP[g];
                }
                else
                {
                    _qqbar_nmod_poly_factor_degree_counts(cntQ, Qp);
                }

                /* For each residue degree f, count roots of Q in F_{p^f}. */
                for (f = 1; f <= dP && result; f++)
                {
                    if (cntP[f] == 0)
                        continue;

                    N = 0;
                    for (g = 1; g <= FLINT_MIN(f, dQ); g++)
                        if (f % g == 0)
                            N += g * cntQ[g];

                    if (N < 1 + (same != 0))
                        result = 0;
                }
            }

            passed++;
        }
    }

    nmod_poly_clear(Pp);
    nmod_poly_clear(Qp);
    nmod_poly_clear(t);
    nmod_poly_clear(u);
    flint_free(cntP);
    flint_free(cntQ);

    return result;
}

/* Attempt to find R with y = R(x), searching up to the given precision.
   This is similar to qqbar_express_in_field, but only uses the raw
   enclosures (qqbar_get_acb may call qqbar_sub for exactness detection,
   which could lead to infinite recursion here). */
int
_qqbar_express_in_field_search(fmpq_poly_t R, const qqbar_t x, const qqbar_t y, slong min_prec, slong max_prec)
{
    slong d, prec;
    acb_ptr vec;
    acb_t z;
    int found = 0;

    d = qqbar_degree(x);

    acb_init(z);
    vec = _acb_vec_init(d + 1);

    for (prec = min_prec; prec <= max_prec && !found; prec *= 2)
    {
        qqbar_enclosure_raw(z, x, prec);
        _acb_vec_set_powers(vec, z, d, prec);
        qqbar_enclosure_raw(vec + d, y, prec);

        fmpq_poly_fit_length(R, d + 1);

        if (_qqbar_acb_lindep(R->coeffs, vec, d + 1, 1, prec) &&
            !fmpz_is_zero(R->coeffs + d))
        {
            fmpz_neg(R->den, R->coeffs + d);
            _fmpq_poly_set_length(R, d);
            _fmpq_poly_normalise(R);
            fmpq_poly_canonicalise(R);

            found = qqbar_equal_fmpq_poly_val(y, R, x);
        }
        else
        {
            fmpq_poly_zero(R);
        }
    }

    acb_clear(z);
    _acb_vec_clear(vec, d + 1);

    return found;
}

/* Attempts to compute x op y (op = 0, 1, 2, 3 for x + y, x - y, x y, x / y;
   op = 4 for (x - y)^2 and op = 6 for -(x - y)^2 / 4) by expressing
   y as an element of Q(x) (swapped = 0) or x as an element of Q(y)
   (swapped = 1). The LLL search uses precisions min_prec, 2 min_prec, ...
   up to max_prec. Returns 1 if
   successful, 0 if the LLL search failed, and -1 if the modular test
   (performed only if test is set) proves that the method does
   not apply. */
int
_qqbar_binary_op_subfield(qqbar_t res, const qqbar_t x, const qqbar_t y, int op, int swapped, slong min_prec, slong max_prec, int test)
{
    const qqbar_struct * gen = swapped ? y : x;
    const qqbar_struct * other = swapped ? x : y;
    fmpq_poly_t R, A, B, P, G, S;
    int success = 0;

    /* The fallback is cheap for small problems. */
    if (qqbar_degree(x) * qqbar_degree(y) < 36 ||
        (fmpz_poly_equal(QQBAR_POLY(x), QQBAR_POLY(y)) &&
            qqbar_degree(x) * (qqbar_degree(x) - 1) / 2 < 36))
        return 0;

    if (test && !_qqbar_subfield_modular_test(QQBAR_POLY(gen), QQBAR_POLY(other),
            fmpz_poly_equal(QQBAR_POLY(gen), QQBAR_POLY(other))))
        return -1;

    fmpq_poly_init(R);

    if (_qqbar_express_in_field_search(R, gen, other, min_prec, max_prec))
    {
        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(P);

        fmpq_poly_set_fmpz_poly(P, QQBAR_POLY(gen));

        /* A = x, B = y as polynomials in gen */
        if (swapped)
        {
            fmpq_poly_swap(A, R);
            fmpq_poly_set_coeff_ui(B, 1, 1);
        }
        else
        {
            fmpq_poly_set_coeff_ui(A, 1, 1);
            fmpq_poly_swap(B, R);
        }

        if (op == 0)
        {
            fmpq_poly_add(A, A, B);
        }
        else if (op == 1)
        {
            fmpq_poly_sub(A, A, B);
        }
        else if (op == 2)
        {
            fmpq_poly_mul(A, A, B);
            fmpq_poly_rem(A, A, P);
        }
        else if (op == 4 || op == 6)
        {
            /* (x - y)^2 or -(x - y)^2 / 4 */
            fmpq_poly_sub(A, A, B);
            fmpq_poly_mul(A, A, A);
            fmpq_poly_rem(A, A, P);
            if (op == 6)
            {
                fmpq_poly_scalar_div_si(A, A, -4);
            }
        }
        else
        {
            fmpq_poly_init(G);
            fmpq_poly_init(S);
            /* B is nonzero mod P since y != 0 */
            fmpq_poly_xgcd(G, S, R, B, P);
            fmpq_poly_mul(A, A, S);
            fmpq_poly_rem(A, A, P);
            fmpq_poly_clear(G);
            fmpq_poly_clear(S);
        }

        qqbar_evaluate_fmpq_poly(res, A, gen);
        success = 1;

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(P);
    }

    fmpq_poly_clear(R);

    return success;
}
