/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

static ulong
_capelli_smallest_prime_factor(ulong d)
{
    ulong t;
    if (d % 2 == 0)
        return 2;
    for (t = 3; t * t <= d; t += 2)
        if (d % t == 0)
            return t;
    return d;
}

/*
    Given E = theta^((q^f - 1)/p) reduced modulo Pf, a squarefree product of
    degree-f irreducibles mod q in whose residue fields theta is a unit,
    return 1 iff some component certifies theta a non-p-th power there.
    For p = 2 the component values are +-1, so a witness is a component
    where E = -1, detected by gcd(E + 1, Pf); for odd p any component value
    different from 1 is a witness, i.e. E - 1 nonzero modulo Pf.
*/
static int
_capelli_euler_witness(const nmod_poly_t E, const nmod_poly_t Pf, ulong p,
                       nmod_poly_t tmp)
{
    nmod_poly_zero(tmp);
    nmod_poly_set_coeff_ui(tmp, 0, 1);

    if (p == 2)
    {
        nmod_poly_add(tmp, E, tmp);
        if (nmod_poly_is_zero(tmp))
            return 1;
        nmod_poly_gcd(tmp, tmp, Pf);
        return nmod_poly_degree(tmp) >= 1;
    }
    else
    {
        nmod_poly_sub(tmp, E, tmp);
        return !nmod_poly_is_zero(tmp);
    }
}

static void
_capelli_preinv(nmod_poly_t finv, const nmod_poly_t f)
{
    nmod_poly_reverse(finv, f, f->length);
    nmod_poly_inv_series(finv, finv, f->length);
}

/*
    Given T irreducible over Q and p prime, attempt to certify that T(x^p)
    is irreducible over Q.  T is assumed primitive with positive leading
    coefficient (as produced by fmpz_poly_factor); the degree-1 exact test
    relies on this to read theta = -T(0)/T(1) in lowest terms, and the
    assumption spares callers with low degree and huge coefficients a
    redundant integer gcd.  The non-underscore chaining version normalises
    arbitrary input first.

    TODO (sqrt(theta) reconstruction): make the p = 2 certificate two-sided
    for deg T >= 3.  When the witness search only ever sees "theta is a
    square" evidence, attempt to construct gamma in K = Q[y]/(T) with
    gamma^2 = theta: at a prime q minimising the number r of local factors,
    compute componentwise square roots in the residue fields, Hensel-lift
    to q^k with k from a Mignotte bound on the factor g in
    T(x^2) = +- g(x) g(-x), resolve the 2^(r-1) componentwise sign
    ambiguity (cheap pruning by rational reconstruction of a single
    coefficient for moderate r; fields with abelian splitting behaviour
    have no low-r primes and would need recursive multiquadratic descent),
    then verify gamma^2 = theta exactly.  Success yields the factorisation
    with both halves irreducible by construction, so certification would
    terminate with an answer in both directions instead of falling back
    after an exhausted budget.

    By Capelli's theorem (prime exponent case), T(x^p) is irreducible over Q
    iff x^p - theta is irreducible over K = Q[y]/(T(y)) iff theta (the image
    of y) is not a p-th power in K.

    We search for a rational prime q with q not dividing lc(T)*T(0) and
    T squarefree mod q.  Then K (x) Q_q is an unramified etale algebra whose
    components have residue fields F_q[y]/(T_i(y)) for the irreducible
    factors T_i of T mod q, and theta is a unit in every component.  A global
    p-th power theta = gamma^p forces gamma to be a unit as well, so theta
    reduces to a p-th power in every residue field.  Hence a single component
    F_{q^f} with p | q^f - 1 and theta-bar^((q^f-1)/p) != 1 proves that theta
    is not a p-th power in K, certifying irreducibility of T(x^p).

    The search performs the cheap distinct-degree phase of factorisation
    mod q and one Euler-criterion power per degree class (no equal-degree
    splitting is needed: for p = 2 the witness condition on the product P_f
    of all degree-f components is deg gcd(E + 1, P_f) > 0, and for odd p it
    is (E - 1) != 0 mod P_f, where E = theta-bar^((q^f-1)/p) mod P_f).

    Returns 1 if a certificate was found (T(x^p) is then provably
    irreducible), and 0 if the search budget was exhausted, in which case
    nothing is claimed.

    Note on the budget: for typical Galois groups a witness prime has
    positive density (often ~ 1/2) among primes with a small-degree
    component, so the first few candidates succeed.  The hard case is a
    multiquadratic-type (elementary abelian) splitting field, e.g. the
    Swinnerton-Dyer polynomials, where the witness density is 2^-n =
    1/(2 deg T): there, all Frobenius elements are exposed only through
    rational primes, so roughly 2 deg T candidates are required.  The budget
    is chosen linear in deg T to cover this case; the per-candidate cost is
    a handful of modular power/gcd operations, so even an exhausted search
    (which happens when T(x^p) is in fact reducible) costs a small fraction
    of the factorisation it precedes.
*/
int
_fmpz_poly_factor_inflation_is_irreducible_capelli(const fmpz_poly_t T, ulong p)
{
    slong m = fmpz_poly_degree(T);
    slong candA = 0, rawB = 0, budgetB, powbits = 0;
    int found = 0, phaseB_ok;
    n_primes_t piter;

    if (m < 1)
        return 0;

    if (m == 1)
    {
        /* K = Q and theta = -T(0)/T(1) in lowest terms (T is primitive).
           By Capelli (p prime), x^p - theta is irreducible over Q iff theta
           is not a p-th power in Q; for p = 2 a negative theta is never a
           square.  This exact O(size) test replaces the prime scan, and is
           effectively two-sided: on a p-th power the caller's fallback
           performs the (genuine) split at trivial cost. */
        fmpz_t u, r, w;
        int upow, vpow;

        if (fmpz_is_zero(T->coeffs + 0))
            return 0;

        fmpz_init(u);
        fmpz_init(r);
        fmpz_init(w);

        /* theta = u / v with u = -T(0), v = T(1) > 0; primitivity of T
           makes this lowest terms, which the p-th power criterion below
           requires (u and v must be separately p-th powers) */
        fmpz_neg(u, T->coeffs + 0);

        if (p == 2 && fmpz_sgn(u) < 0)
            found = 1;                      /* theta < 0: nonsquare */
        else
        {
            fmpz_abs(w, u);
            upow = fmpz_root(r, w, (slong) p);
            vpow = 1;
            if (!fmpz_is_one(T->coeffs + 1))
                vpow = fmpz_root(r, T->coeffs + 1, (slong) p);
            found = !(upow && vpow);
        }

        fmpz_clear(u);
        fmpz_clear(r);
        fmpz_clear(w);

        return found;
    }

    if (m == 2 && p == 2)
    {
        /* K is a quadratic field, where theta in K^2 has an exact integer
           criterion: with T = a y^2 + b y + c irreducible (a > 0), theta
           is a square in K iff a c is a perfect square, say s^2 with
           s >= 0, and a (2 e s - b) is a perfect square for one of
           e = +-1; equivalently T(x^2) admits the norm-form factorisation
           a (x^2 + beta x + nu)(x^2 - beta x + nu), which is the only way
           it can factor when T is irreducible.  Like the degree-1 case
           this is O(size) and two-sided.  (Degree-2 T is always Galois,
           so the mixed-degree evidence bail below can never fire for it;
           this test removes the scan for the ubiquitous quadratic-field
           inputs altogether.) */
        fmpz_t s2, w;
        int e;

        fmpz_init(s2);
        fmpz_init(w);

        fmpz_mul(w, T->coeffs + 2, T->coeffs + 0);
        if (fmpz_sgn(w) < 0 || !fmpz_is_square(w))
            found = 1;
        else
        {
            found = 1;
            fmpz_sqrt(s2, w);
            for (e = 0; e < 2; e++)
            {
                fmpz_mul_2exp(w, s2, 1);
                if (e)
                    fmpz_neg(w, w);
                fmpz_sub(w, w, T->coeffs + 1);
                fmpz_mul(w, w, T->coeffs + 2);
                if (fmpz_sgn(w) >= 0 && fmpz_is_square(w))
                    found = 0;              /* theta = gamma^2: reducible */
            }
        }

        fmpz_clear(s2);
        fmpz_clear(w);

        return found;
    }

    /* Phase A is cheap at every scale and, on inputs arising in practice
       (e.g. qqbar arithmetic), succeeds often enough that it is worth
       attempting unconditionally.  Phase B, in contrast, is only worth its
       possible exhaustion when the fallback would be genuinely expensive:
       below (inflated) degree ~64 even the hard abelian-type inputs
       factor in a few milliseconds, which is what an exhausted phase-B
       scan would itself cost. */
    phaseB_ok = (p * (ulong) m > 64);

    /* Phase B (see below) can succeed only for inputs whose splitting field
       is an abelian group of exponent p; a reducible T(x^p) with such local
       behaviour makes phase B pure waste.  The common source of that trap
       is cyclotomic-type inputs, which are self-reciprocal with constant
       term of absolute value 1, whereas the hard Swinnerton-Dyer-type
       inputs have huge constant terms (norms of sums of radicals); so we
       exclude plausibly-cyclotomic inputs from phase B up front.  Search
       policy only; certificates are unaffected. */
    if (fmpz_is_pm1(T->coeffs + 0))
    {
        fmpz_poly_t rev;
        fmpz_poly_init(rev);
        fmpz_poly_reverse(rev, T, T->length);
        if (fmpz_poly_equal(rev, T))
            phaseB_ok = 0;
        else
        {
            fmpz_poly_neg(rev, rev);
            if (fmpz_poly_equal(rev, T))
                phaseB_ok = 0;
        }
        fmpz_poly_clear(rev);
    }

    /* Phase B budget: the hard irreducible inputs (elementary-abelian-type
       splitting fields, e.g. Swinnerton-Dyer) have witness density
       1/(p * deg T) among rational primes, so the budget must be a small
       multiple of deg T.  Exhausting it (which happens when T(x^p) is in
       fact reducible) costs one Frobenius power and one gcd mod q per
       prime, a small fraction of the factorisation that follows. */
    budgetB = FLINT_MIN(512 + 32 * m, WORD(24576));

    n_primes_init(piter);

    while (!found && rawB <= budgetB && (candA < 16 || phaseB_ok))
    {
        ulong q = n_primes_next(piter);
        nmod_poly_t Tq, dT, sf, R, X, xp, Pf, E, tmp, Fi;
        slong f, fmax;

        if (q == 2)
            continue;
        if (fmpz_fdiv_ui(T->coeffs + m, q) == 0)
            continue;
        if (fmpz_fdiv_ui(T->coeffs + 0, q) == 0)
            continue;
        if (candA >= 16 && p != 2 && (q - 1) % p != 0)
            continue;   /* phase B tests degree-1 residue fields: need p | q - 1 */

        nmod_poly_init(Tq, q);
        nmod_poly_init(dT, q);
        nmod_poly_init(sf, q);
        nmod_poly_init(R, q);
        nmod_poly_init(X, q);
        nmod_poly_init(xp, q);
        nmod_poly_init(Pf, q);
        nmod_poly_init(E, q);
        nmod_poly_init(tmp, q);
        nmod_poly_init(Fi, q);

        fmpz_poly_get_nmod_poly(Tq, T);
        nmod_poly_set_coeff_ui(xp, 1, 1);            /* xp = x */

        if (candA < 16)
        {
            /* Phase A (generic inputs): distinct-degree scan with an Euler
               test per degree class.  When T(x^p) is irreducible and the
               Galois group of its splitting field is not an abelian
               p-group, a witness has positive density (typically about
               1/2) among these candidates, so a handful suffice. */
            nmod_poly_derivative(dT, Tq);
            nmod_poly_gcd(sf, Tq, dT);

            if (nmod_poly_degree(sf) == 0)
            {
                slong mleqp = 0, dfseen = 0, ctested = 0;

                candA++;
                fmax = FLINT_MAX(WORD(6), FLINT_MIN(m, (slong) p));

                nmod_poly_set(R, Tq);
                nmod_poly_rem(X, xp, R);
                _capelli_preinv(Fi, R);

                for (f = 1; f <= fmax && nmod_poly_degree(R) >= 1 && !found; f++)
                {
                    nmod_poly_powmod_ui_binexp_preinv(X, X, q, R, Fi);
                    nmod_poly_sub(tmp, X, xp);
                    nmod_poly_gcd(Pf, tmp, R);

                    if (nmod_poly_degree(Pf) >= 1)
                    {
                        fmpz_t e;

                        dfseen++;
                        if (f <= (slong) p)
                            mleqp += nmod_poly_degree(Pf);

                        fmpz_init_set_ui(e, q);
                        fmpz_pow_ui(e, e, (ulong) f);
                        fmpz_sub_ui(e, e, 1);

                        if (fmpz_fdiv_ui(e, p) == 0)
                        {
                            fmpz_divexact_ui(e, e, p);
                            ctested += nmod_poly_degree(Pf) / f;
                            nmod_poly_rem(tmp, xp, Pf);
                            _capelli_preinv(Fi, Pf);   /* recomputed below if R shrinks */
                            nmod_poly_powmod_fmpz_binexp_preinv(E, tmp, e, Pf, Fi);

                            found = _capelli_euler_witness(E, Pf, p, tmp);
                        }

                        fmpz_clear(e);

                        if (!found && f < fmax)
                        {
                            nmod_poly_divrem(R, tmp, R, Pf);   /* exact */
                            nmod_poly_rem(X, X, R);
                            _capelli_preinv(Fi, R);
                        }
                    }
                }

                /* Phase B assumes every Frobenius element has order <= p
                   (abelian exponent-p splitting field); any mass in local
                   degrees > p at any good prime disproves that. */
                if (!found && mleqp != m)
                    phaseB_ok = 0;

                /* Sequential evidence bail: at a prime whose factor degrees
                   are mixed, T is provably not Galois and the residue
                   components behave nearly independently, so under
                   irreducibility of T(x^p) each tested component is a
                   non-p-th-power with probability about 1 - 1/p.  Some 12
                   consecutive p-th-power components with no witness make
                   irreducibility very unlikely; stop scanning and let the
                   fallback do the (probable) real factorisation.  Primes
                   with all component degrees equal are Galois-consistent
                   and carry essentially no such evidence -- for abelian
                   inputs (e.g. Swinnerton-Dyer) "all components are powers"
                   has probability close to 1 under irreducibility -- so
                   they are deliberately not counted; this keeps the hard
                   abelian certificates unaffected. */
                if (!found && dfseen >= 2)
                {
                    powbits += ctested;
                    if (powbits >= 12)
                    {
                        candA = 16;
                        phaseB_ok = 0;
                    }
                }
            }
        }
        else
        {
            /* Phase B (specialist): essentially only inputs whose splitting
               field has an abelian p-group Galois group survive phase A
               with T(x^p) irreducible.  For those, every completion of the
               splitting field has residue degree <= p over F_q, so a
               witness can only sit in a degree-1 residue field (a prime of
               K of residue degree p splits, rather than staying inert, in
               K(theta^(1/p))).  Hence: one Frobenius power, the product of
               the linear factors, and an Euler test only when that product
               is nonempty, which for the hard inputs is rare -- e.g. for
               Swinnerton-Dyer inputs a fraction 2^(1-n) of primes.  The
               squarefreeness check (needed for the soundness of the
               certificate, not for the search) is deferred until a
               potential witness is on the table. */
            rawB++;

            nmod_poly_rem(X, xp, Tq);   /* deg(T) = 1 edge case */
            _capelli_preinv(Fi, Tq);
            nmod_poly_powmod_ui_binexp_preinv(X, X, q, Tq, Fi);
            nmod_poly_sub(tmp, X, xp);
            nmod_poly_gcd(Pf, tmp, Tq);

            if (nmod_poly_degree(Pf) >= 1)
            {
                nmod_poly_derivative(dT, Tq);
                nmod_poly_gcd(sf, Tq, dT);

                if (nmod_poly_degree(sf) == 0)
                {
                    nmod_poly_rem(tmp, xp, Pf);
                    _capelli_preinv(Fi, Pf);
                    nmod_poly_powmod_ui_binexp_preinv(E, tmp, (q - 1) / p, Pf, Fi);

                    found = _capelli_euler_witness(E, Pf, p, tmp);
                }
            }
        }

        nmod_poly_clear(Tq);
        nmod_poly_clear(dT);
        nmod_poly_clear(sf);
        nmod_poly_clear(R);
        nmod_poly_clear(X);
        nmod_poly_clear(xp);
        nmod_poly_clear(Pf);
        nmod_poly_clear(E);
        nmod_poly_clear(tmp);
        nmod_poly_clear(Fi);
    }

    n_primes_clear(piter);

    return found;
}


/*
    Given T irreducible over Q and d >= 1, return 1 if T(x^d) is certified
    irreducible over Q, and 0 if no certificate was found (in which case
    nothing is claimed: T(x^d) may or may not be irreducible).  Chains the
    prime-step certificate over the prime factorisation of d, via
    T(x^d) = U(x^(d/p)) with U = T(x^p).
*/
int
fmpz_poly_factor_inflation_is_irreducible_capelli(const fmpz_poly_t T, ulong d)
{
    fmpz_poly_t U;
    int res = 1;

    if (d == 0 || fmpz_poly_degree(T) < 1)
        return 0;

    fmpz_poly_init(U);

    /* the underscore prime-step assumes primitive input with positive
       leading coefficient; normalise once (inflation preserves both) */
    {
        fmpz_t cont;
        fmpz_init(cont);
        fmpz_poly_content(cont, T);
        if (fmpz_sgn(fmpz_poly_lead(T)) < 0)
            fmpz_neg(cont, cont);
        fmpz_poly_scalar_divexact_fmpz(U, T, cont);
        fmpz_clear(cont);
    }

    while (d > 1 && res)
    {
        ulong p = _capelli_smallest_prime_factor(d);

        if (_fmpz_poly_factor_inflation_is_irreducible_capelli(U, p))
        {
            fmpz_poly_inflate(U, U, p);
            d /= p;
        }
        else
            res = 0;
    }

    fmpz_poly_clear(U);
    return res;
}
