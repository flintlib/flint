/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    L(s, chi) for integer s. Let chi be primitive modulo q > 1 with parity
    d (chi(-1) = (-1)^d).

    * s = 1 - n <= 0: L(1-n, chi) = -B_{n,chi} / n, which is zero unless
      n = d mod 2 (generalized Bernoulli numbers).

    * s = n >= 1 with n = d mod 2 ("same parity"):

        L(n, chi) = (-1)^(1+(n-d)/2) tau(chi) / (2 i^d) (2 pi/q)^n B_{n,conj(chi)} / n!

    * s = 1, d = 0: L(1, chi) = -(tau(chi)/q) sum_a conj(chi(a)) log|2 sin(pi a/q)|.

    * s >= 2 with n != d mod 2 ("opposite parity"): no closed form is known;
      we use the Euler product for large s, Borwein's algorithm with binary
      splitting (acb_dirichlet_l_ui_borwein_bsplit) at high precision, and
      the approximate functional equation or the Hurwitz zeta function
      otherwise.

    Imprimitive characters are reduced to the primitive character inducing
    them, and q = 1 is the Riemann zeta function.
*/

#include <math.h>
#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "bernoulli.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"

/* res = sum_c v[c] exp(2 pi i c / m) (sign = +1) or exp(-2 pi i c / m)
   (sign = -1); v has length m */
static void
eval_classes(acb_t res, arb_srcptr v, ulong m, int sign, slong prec)
{
    if (m == 1)
    {
        arb_set(acb_realref(res), v);
        arb_zero(acb_imagref(res));
    }
    else if (m == 2)
    {
        arb_sub(acb_realref(res), v, v + 1, prec);
        arb_zero(acb_imagref(res));
    }
    else
    {
        acb_ptr w;
        acb_t z;
        ulong c;

        w = _acb_vec_init(m);
        acb_init(z);

        for (c = 0; c < m; c++)
            arb_set(acb_realref(w + c), v + c);

        acb_unit_root(z, m, prec);
        if (sign < 0)
            acb_conj(z, z);

        _acb_poly_evaluate(res, w, m, z, prec);

        _acb_vec_clear(w, m);
        acb_clear(z);
    }
}

/*
    Generalized Bernoulli number B_{n,psi} with psi = chi (conj = 0) or
    conj(chi) (conj = 1), for chi primitive modulo q >= 3 and n = parity
    mod 2:

        B_{n,psi} = q^(n-1) sum_{a=1}^{q} psi(a) B_n(a/q)
                  = 2 q^(n-1) sum_{a < q/2} psi(a) B_n(a/q)

    (the terms a and q - a are equal). With q^(n-1) B_n(a/q) =
    (1/(q D)) sum_k C_k a^(n-k), C_k = D binom(n,k) B_k q^k integers, the
    values are evaluated exactly by Horner's rule and summed exactly over
    the classes of residues with the same character value; only the final
    combination of at most ord(chi) class sums is numerical.
*/
static void
bernoulli_char(acb_t res, ulong n, const dirichlet_group_t G, const dirichlet_char_t chi, int conj, slong prec)
{
    ulong q = G->q, a, m, c;
    ulong * cls;
    slong k, wp;
    fmpq * B;
    fmpz * C;
    fmpz * U;
    fmpz_t D, t, v;
    fmpq_t r;
    arb_ptr V;

    m = dirichlet_order_char(G, chi);

    B = _fmpq_vec_init(n + 1);
    C = _fmpz_vec_init(n + 1);
    U = _fmpz_vec_init(m);
    fmpz_init(D);
    fmpz_init(t);
    fmpz_init(v);
    fmpq_init(r);

    bernoulli_fmpq_vec_no_cache(B, 0, n + 1);

    /* C_k = binom(n,k) B_k q^k, scaled by the common denominator D */
    fmpz_one(D);
    for (k = 0; k <= (slong) n; k++)
        if (!fmpq_is_zero(B + k))
            fmpz_lcm(D, D, fmpq_denref(B + k));

    for (k = 0; k <= (slong) n; k++)
    {
        if (fmpq_is_zero(B + k))
            continue;

        fmpz_bin_uiui(t, n, k);
        fmpz_mul(t, t, fmpq_numref(B + k));
        fmpz_divexact(v, D, fmpq_denref(B + k));
        fmpz_mul(t, t, v);
        fmpz_ui_pow_ui(v, q, k);
        fmpz_mul(C + k, t, v);
    }

    /* class indices of chi(a) for a < q/2 */
    cls = flint_malloc(sizeof(ulong) * (q / 2 + 1));
    dirichlet_chi_vec_order(cls, G, chi, m, q / 2 + 1);

    for (a = 1; 2 * a < q; a++)
    {
        c = cls[a];

        if (c == DIRICHLET_CHI_NULL)
            continue;

        if (conj && c != 0)
            c = m - c;

        /* v = sum_k C_k a^(n-k) */
        fmpz_set(v, C + 0);
        for (k = 1; k <= (slong) n; k++)
        {
            fmpz_mul_ui(v, v, a);
            fmpz_add(v, v, C + k);
        }

        fmpz_add(U + c, U + c, v);
    }

    flint_free(cls);

    /* res = 2 / (q D) sum_c U_c exp(2 pi i c / m) */
    fmpz_mul_ui(D, D, q);

    if (m <= 2)
    {
        if (m == 2)
            fmpz_sub(U, U, U + 1);
        fmpz_mul_2exp(U, U, 1);
        fmpq_set_fmpz_frac(r, U, D);
        arb_set_fmpq(acb_realref(res), r, prec);
        arb_zero(acb_imagref(res));
    }
    else
    {
        /* the class sums have about the same size as the result times
           the number of classes; allow for cancellation */
        wp = prec + 2 * FLINT_BIT_COUNT(m) + 10;

        V = _arb_vec_init(m);
        for (c = 0; c < m; c++)
            arb_set_round_fmpz(V + c, U + c, wp);

        eval_classes(res, V, m, 1, wp);

        acb_mul_2exp_si(res, res, 1);
        acb_div_fmpz(res, res, D, prec);

        _arb_vec_clear(V, m);
    }

    _fmpq_vec_clear(B, n + 1);
    _fmpz_vec_clear(C, n + 1);
    _fmpz_vec_clear(U, m);
    fmpz_clear(D);
    fmpz_clear(t);
    fmpz_clear(v);
    fmpq_clear(r);
}

/* L(s, chi) by the closed forms for s <= 0 and for s >= 1 of the same
   parity as chi; chi primitive modulo q >= 3. */
static void
l_si_bernoulli(acb_t res, slong s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    int d = dirichlet_parity_char(G, chi);
    ulong q = G->q;
    slong wp;

    if (s <= 0)
    {
        ulong n = 1 - s;

        if ((n % 2) != (ulong) d)
        {
            acb_zero(res);
            return;
        }

        /* L(1-n, chi) = -B_{n,chi} / n */
        bernoulli_char(res, n, G, chi, 0, prec + 10);
        acb_div_si(res, res, -(slong) n, prec);
    }
    else
    {
        ulong n = s;
        acb_t tau;
        arb_t t;

        acb_init(tau);
        arb_init(t);

        /* the Gauss sum loses a few bits */
        wp = prec + 2 * FLINT_BIT_COUNT(q) + 20;

        bernoulli_char(res, n, G, chi, 1, wp);

        acb_dirichlet_gauss_sum(tau, G, chi, wp);
        acb_mul(res, res, tau, wp);

        /* (2 pi / q)^n / n! */
        arb_const_pi(t, wp);
        arb_mul_2exp_si(t, t, 1);
        arb_div_ui(t, t, q, wp);
        arb_pow_ui(t, t, n, wp);
        acb_mul_arb(res, res, t, wp);
        arb_fac_ui(t, n, wp);
        acb_div_arb(res, res, t, wp);

        /* (-1)^(1+(n-d)/2) / (2 i^d) */
        acb_mul_2exp_si(res, res, -1);
        if (d)
            acb_div_onei(res, res);
        if ((1 + (n - d) / 2) % 2)
            acb_neg(res, res);

        if (dirichlet_char_is_real(G, chi))
            arb_zero(acb_imagref(res));

        acb_set_round(res, res, prec);

        acb_clear(tau);
        arb_clear(t);
    }
}

/*
    L(1, chi) for chi even and primitive modulo q >= 3:

        L(1, chi) = -(2 tau(chi) / q) sum_{a < q/2} conj(chi(a)) log(2 sin(pi a / q)).

    The sines are the imaginary parts of the powers of exp(pi i / q), and
    the logarithm is taken of the product of the sines in each class of
    residues with the same character value.
*/
static void
l_one_even(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    ulong q = G->q, a, m, c;
    ulong * cls;
    slong wp;
    arb_ptr P;
    acb_t w, z, tau;
    arb_t t;

    m = dirichlet_order_char(G, chi);

    /* guard bits for the accumulated error in the powers and for
       cancellation between the classes */
    wp = prec + 3 * FLINT_BIT_COUNT(q) + 20;

    P = _arb_vec_init(m);
    acb_init(w);
    acb_init(z);
    acb_init(tau);
    arb_init(t);

    for (c = 0; c < m; c++)
        arb_one(P + c);

    acb_unit_root(w, 2 * q, wp);
    acb_one(z);

    cls = flint_malloc(sizeof(ulong) * (q / 2 + 1));
    dirichlet_chi_vec_order(cls, G, chi, m, q / 2 + 1);

    for (a = 1; 2 * a < q; a++)
    {
        acb_mul(z, z, w, wp);

        c = cls[a];

        if (c == DIRICHLET_CHI_NULL)
            continue;

        /* 2 sin(pi a / q) */
        arb_mul_2exp_si(t, acb_imagref(z), 1);
        arb_mul(P + c, P + c, t, wp);
    }

    flint_free(cls);

    if (m == 2)
    {
        /* one logarithm: log(P_0 / P_1) */
        arb_div(acb_realref(res), P + 0, P + 1, wp);
        arb_log(acb_realref(res), acb_realref(res), wp);
        arb_zero(acb_imagref(res));
    }
    else
    {
        for (c = 0; c < m; c++)
            arb_log(P + c, P + c, wp);

        eval_classes(res, P, m, -1, wp);
    }

    acb_dirichlet_gauss_sum(tau, G, chi, wp);
    acb_mul(res, res, tau, wp);
    acb_mul_2exp_si(res, res, 1);
    acb_div_ui(res, res, q, wp);
    acb_neg(res, res);

    if (dirichlet_char_is_real(G, chi))
        arb_zero(acb_imagref(res));

    acb_set_round(res, res, prec);

    _arb_vec_clear(P, m);
    acb_clear(w);
    acb_clear(z);
    acb_clear(tau);
    arb_clear(t);
}

/* whether acb_dirichlet_l_euler_product reaches the full accuracy for
   integer s >= 2: it needs the primes up to about 2^(prec/s) and stops at
   100 + prec^(3/2) */
static int
euler_product_ok(slong s, slong prec)
{
    double wp = prec + 10;
    return s >= 2 && (wp / s) * 0.6931471805599453 < log(0.5 * (100 + wp * sqrt(wp)));
}

/* L(s, chi) for chi primitive modulo q >= 3 */
static void
l_si_primitive(acb_t res, slong s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    ulong q = G->q;
    int d = dirichlet_parity_char(G, chi);

    if (s <= 0)
    {
        ulong n = 1 - s;

        if ((n % 2) != (ulong) d)
        {
            acb_zero(res);
        }
        else if (n < 30 || prec >= 0.84 * pow(n, 1.4))   /* measured cutoff */
        {
            l_si_bernoulli(res, s, G, chi, prec);
        }
        else
        {
            /* functional equation and the Euler product (which lose some
               bits in the gamma function and powers) */
            acb_t t;
            acb_init(t);
            acb_set_si(t, s);
            _acb_dirichlet_l_generic(res, t, G, chi, prec + 2 * FLINT_BIT_COUNT(n) + 10);
            acb_set_round(res, res, prec);
            acb_clear(t);
        }
    }
    else if (s % 2 == d)
    {
        /* closed form, unless s is so large that the Euler product
           converges immediately (measured cutoff) */
        if (euler_product_ok(s, prec) &&
            prec < 3.0 * pow(s, 1.25) * (1.0 + 0.1 * log(q) * 1.4426950408889634))
        {
            acb_t t;
            acb_init(t);
            acb_set_si(t, s);
            acb_dirichlet_l_euler_product(res, t, G, chi, prec + 10);
            acb_set_round(res, res, prec);
            acb_clear(t);
        }
        else
        {
            l_si_bernoulli(res, s, G, chi, prec);
        }
    }
    else if (s == 1)
    {
        l_one_even(res, G, chi, prec);
    }
    else
    {
        double sd = s, qd = q;
        acb_t t;

        acb_init(t);
        acb_set_si(t, s);

        if (prec <= 10 * sd && euler_product_ok(s, prec))
        {
            acb_dirichlet_l_euler_product(res, t, G, chi, prec + 10);
        }
        else if (q == 4 && s == 2)
        {
            arb_const_catalan(acb_realref(res), prec);
            arb_zero(acb_imagref(res));
        }
        /* timings at s = 2: Borwein's algorithm is the fastest at all
           precisions for q <= 101, from about 150 bits for q = 331 and
           from about 1150 bits for q = 1009 (the approximate functional
           equation being faster below) */
        else if (prec >= FLINT_MAX(1.2 * sd * sd, 0.0006 * qd * qd * sd))
        {
            acb_dirichlet_l_ui_borwein_bsplit(res, s, G, chi, prec);
        }
        else if (q >= 30 && (prec >= 64 || q >= 500) &&
            (sd <= 10 + ((q >= 100) ? 0.05 : 0.02) * prec || q >= 500))
        {
            fmpq_t u;
            fmpq_init(u);
            fmpq_set_si(u, s, 1);
            acb_dirichlet_l_fmpq_afe(res, u, G, chi, prec);
            fmpq_clear(u);
        }
        else
        {
            _acb_dirichlet_l_generic(res, t, G, chi, prec + 10);
        }

        acb_set_round(res, res, prec);
        acb_clear(t);
    }
}

void
acb_dirichlet_l_si(acb_t res, slong s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    ulong q, cond;

    q = (G == NULL) ? 1 : G->q;

    /* avoid overflow in 1 - s */
    if (s < -(WORD(1) << (FLINT_BITS - 4)) || s > (WORD(1) << (FLINT_BITS - 4)))
    {
        acb_t t;
        acb_init(t);
        acb_set_si(t, s);
        _acb_dirichlet_l_generic(res, t, G, chi, prec);
        acb_clear(t);
        return;
    }

    cond = (q == 1) ? 1 : dirichlet_conductor_char(G, chi);

    if (cond == 1)
    {
        /* zeta(s) times the Euler factors of the principal character */
        if (s == 1)
        {
            acb_indeterminate(res);
            return;
        }

        if (s >= 0)
        {
            arb_zeta_ui(acb_realref(res), s, prec + (q > 1) * 10);
            arb_zero(acb_imagref(res));
        }
        else
        {
            acb_set_si(res, s);
            acb_dirichlet_zeta(res, res, prec + (q > 1) * 10);
        }
    }
    else if (cond == q)
    {
        l_si_primitive(res, s, G, chi, prec);
        return;
    }
    else
    {
        dirichlet_group_t G0;
        dirichlet_char_t chi0;
        slong wp;

        dirichlet_subgroup_init(G0, G, cond);
        dirichlet_char_init(chi0, G0);
        dirichlet_char_lower(chi0, G0, chi, G);

        /* the Euler factors can cancel when s = 0 */
        wp = prec + 10 + ((s == 0) ? FLINT_BIT_COUNT(G->expo) : 0);

        l_si_primitive(res, s, G0, chi0, wp);

        if (q > 1 && !acb_is_zero(res))
        {
            /* multiply by prod_{p | q} (1 - chi0(p) p^(-s)) */
            acb_t f, c;
            slong k;
            n_factor_t fac;

            acb_init(f);
            acb_init(c);

            n_factor_init(&fac);
            n_factor(&fac, q, 1);

            for (k = 0; k < fac.num; k++)
            {
                ulong p = fac.p[k];

                if (cond % p == 0)
                    continue;

                acb_dirichlet_chi(c, G0, chi0, p, wp);

                if (s >= 0)
                {
                    arb_ui_pow_ui(acb_realref(f), p, s, wp);
                    arb_zero(acb_imagref(f));
                    acb_div(c, c, f, wp);
                }
                else
                {
                    arb_ui_pow_ui(acb_realref(f), p, -s, wp);
                    arb_zero(acb_imagref(f));
                    acb_mul(c, c, f, wp);
                }

                acb_sub_ui(c, c, 1, wp);
                acb_neg(c, c);
                acb_mul(res, res, c, wp);
            }

            if (dirichlet_char_is_real(G, chi))
                arb_zero(acb_imagref(res));

            acb_clear(f);
            acb_clear(c);
        }

        acb_set_round(res, res, prec);

        dirichlet_char_clear(chi0);
        dirichlet_group_clear(G0);
        return;
    }

    if (q > 1 && !acb_is_zero(res) && acb_is_finite(res))
    {
        /* principal character: multiply by prod_{p | q} (1 - p^(-s)) */
        slong k, wp = prec + 10;
        arb_t f;
        n_factor_t fac;

        arb_init(f);

        n_factor_init(&fac);
        n_factor(&fac, q, 1);

        for (k = 0; k < fac.num; k++)
        {
            ulong p = fac.p[k];

            if (s >= 0)
            {
                arb_ui_pow_ui(f, p, s, wp);
                arb_inv(f, f, wp);
            }
            else
            {
                arb_ui_pow_ui(f, p, -s, wp);
            }

            arb_sub_ui(f, f, 1, wp);
            arb_neg(f, f);
            acb_mul_arb(res, res, f, wp);
        }

        acb_set_round(res, res, prec);
        arb_clear(f);
    }
}
