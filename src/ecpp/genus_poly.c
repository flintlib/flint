/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "arb.h"
#include "acb.h"
#include "arb_poly.h"
#include "acb_modular.h"
#include "ecpp.h"

/*
    The factor of the Hilbert class polynomial over the genus field.

    Let D = p_1^* ... p_g^* be a fundamental discriminant (with q0 in
    {-4, 8, -8} counted among the p_i^* when D is even), K = Q(sqrt D), H
    its Hilbert class field and L = K(sqrt(p_1^*), ..., sqrt(p_g^*)) the
    genus field, the fixed field of the subgroup Cl^2 of the class group
    Cl = Gal(H/K); [L : K] = 2^{g-1}. Accordingly

        H_D = prod_{gamma in Cl / Cl^2} F_gamma,
        F_gamma = prod_{[a] in the coset gamma} (X - j(a)),

    where the F_gamma are the conjugates over K of F_0, the factor of the
    principal genus, whose coefficients lie in L. By Artin reciprocity the
    class [a] acts on sqrt(p_i^*) by the genus character chi_i([a]) =
    (p_i^* / a). Writing a coefficient of F_0 as sum_S kappa_S b_S with
    b_S = prod_{i in S} sqrt(p_i^*) over the subsets S of {1, ..., g-1}
    and kappa_S in K, the coefficients of the conjugates are
    sum_S chi_gamma(S) kappa_S b_S, so that kappa_S b_S is recovered by a
    character sum over the genera. The kappa_S = u + v sqrt(D) have 2-power
    denominators.

    Modulo a prime n with all chi_i(n) = 1, the square roots of the p_i^*
    exist and F_0 maps to a factor of H_D mod n of degree
    o = h / 2^{g-1}, whose roots are j-invariants of curves with CM by
    the maximal order of K just as well as those of H_D. Finding a root of
    F_0 instead of H_D costs a fraction 1 / 2^{g-1} of the work.

    Input: D, the pstar[0..g-1] (the p_i^* and possibly q0), sqrts[i] a
    square root of pstar[i] modulo n. On success (return 1) F is monic of
    degree o with F | H_D mod n; returns 0 if the coefficients could not be
    identified (the caller then falls back to H_D).
*/

/* genus of the form (a, b, c): bit i set if chi_i = -1 */
static slong
_form_genus(slong a, slong c, const slong * pstar, slong g)
{
    slong i, code = 0;
    fmpz_t t, u;
    fmpz_init(t);
    fmpz_init(u);
    for (i = 0; i < g; i++)
    {
        slong p = FLINT_ABS(pstar[i]);
        slong x = (a % (p == 4 || p == 8 ? 2 : p) != 0) ? a : c;
        fmpz_set_si(t, pstar[i]);
        fmpz_set_si(u, x);
        if (fmpz_kronecker(t, u) == -1)
            code |= WORD(1) << i;
    }
    fmpz_clear(t);
    fmpz_clear(u);
    return code;
}

int
ecpp_class_poly_genus(fmpz_mod_poly_t F, slong D, const slong * pstar, slong g,
                                const fmpz * sqrts, const fmpz_mod_ctx_t ctx)
{
    slong i, k, a, b, c, ac, nforms, ngenus, o, prec, S, gamma, den;
    slong * forms;          /* a, b (negative: b and -b), c, genus */
    slong * genus_of;       /* genus code -> index 0 .. 2^{g-1} - 1, or -1 */
    slong * genus_code;
    arb_poly_struct * Fg;
    arb_t sqrtD;
    acb_t z, bS, sqrtDnum, kappa;
    arb_t t, u, v;
    fmpz_t U, V, inv2, tmp, tmp2;
    double lgh, lgmax;
    int success = 1;

    if (g < 2)
        return 0;

    /* reduced forms (Cohen 5.3.5), as in acb_modular_hilbert_class_poly */
    nforms = 0;
    forms = NULL;
    {
        slong alloc = 0;
        b = D & 1;
        do
        {
            ac = (b * b - D) / 4;
            a = FLINT_MAX(b, 1);
            do
            {
                if (ac % a == 0 && n_gcd(n_gcd(a, b), ac / a) == 1)
                {
                    c = ac / a;
                    if (nforms >= alloc)
                    {
                        alloc = FLINT_MAX(8, 2 * alloc);
                        forms = flint_realloc(forms, 4 * alloc * sizeof(slong));
                    }
                    forms[4 * nforms + 0] = a;
                    forms[4 * nforms + 1] = (a == b || a * a == ac || b == 0) ? b : -b;
                    forms[4 * nforms + 2] = c;
                    forms[4 * nforms + 3] = _form_genus(a, c, pstar, g);
                    nforms++;
                }
                a++;
            }
            while (a * a <= ac);
            b += 2;
        }
        while (3 * b * b <= -D);
    }

    /* the genera that occur (2^{g-1} of them) */
    ngenus = WORD(1) << (g - 1);
    genus_of = flint_malloc((WORD(1) << g) * sizeof(slong));
    genus_code = flint_malloc(ngenus * sizeof(slong));
    for (i = 0; i < (WORD(1) << g); i++)
        genus_of[i] = -1;
    k = 0;
    for (i = 0; i < nforms; i++)
    {
        slong code = forms[4 * i + 3];
        if (genus_of[code] == -1)
        {
            if (k == ngenus)
            {
                flint_free(forms);
                flint_free(genus_of);
                flint_free(genus_code);
                return 0;   /* should not happen */
            }
            genus_code[k] = code;
            genus_of[code] = k++;
        }
    }
    if (k != ngenus || genus_of[0] != 0)
    {
        flint_free(forms);
        flint_free(genus_of);
        flint_free(genus_code);
        return 0;
    }

    /* precision from the largest factor height (heuristic as for H_D) */
    lgmax = 0.0;
    for (gamma = 0; gamma < ngenus; gamma++)
    {
        lgh = 0.0;
        for (i = 0; i < nforms; i++)
            if (genus_of[forms[4 * i + 3]] == gamma)
                lgh += (forms[4 * i + 1] < 0 ? 2.0 : 1.0) / forms[4 * i];
        lgmax = FLINT_MAX(lgmax, lgh);
    }
    prec = 3.141593 * sqrt((double) -D) * lgmax * 1.442696;
    prec = prec * 1.01 + 40 + 2 * g;

    Fg = flint_malloc(ngenus * sizeof(arb_poly_struct));
    for (gamma = 0; gamma < ngenus; gamma++)
        arb_poly_init(Fg + gamma);
    arb_init(sqrtD);
    acb_init(z); acb_init(bS); acb_init(sqrtDnum); acb_init(kappa);
    arb_init(t); arb_init(u); arb_init(v);
    fmpz_init(U); fmpz_init(V); fmpz_init(inv2); fmpz_init(tmp); fmpz_init(tmp2);

    for (;;)
    {
        arb_poly_t lin;
        arb_poly_init(lin);

        arb_set_si(sqrtD, -D);
        arb_sqrt(sqrtD, sqrtD, prec);
        for (gamma = 0; gamma < ngenus; gamma++)
            arb_poly_one(Fg + gamma);

        for (i = 0; i < nforms; i++)
        {
            a = forms[4 * i];
            b = forms[4 * i + 1];
            gamma = genus_of[forms[4 * i + 3]];
            arb_set_si(acb_realref(z), -FLINT_ABS(b));
            arb_set(acb_imagref(z), sqrtD);
            acb_div_si(z, z, 2 * a, prec);
            acb_modular_j(z, z, prec);
            if (b < 0)
            {
                /* (x^2 - 2 re(j) x + |j|^2) */
                arb_poly_fit_length(lin, 3);
                arb_mul(lin->coeffs, acb_realref(z), acb_realref(z), prec);
                arb_addmul(lin->coeffs, acb_imagref(z), acb_imagref(z), prec);
                arb_mul_2exp_si(lin->coeffs + 1, acb_realref(z), 1);
                arb_neg(lin->coeffs + 1, lin->coeffs + 1);
                arb_one(lin->coeffs + 2);
                _arb_poly_set_length(lin, 3);
            }
            else
            {
                arb_poly_fit_length(lin, 2);
                arb_neg(lin->coeffs, acb_realref(z));
                arb_one(lin->coeffs + 1);
                _arb_poly_set_length(lin, 2);
            }
            arb_poly_mul(Fg + gamma, Fg + gamma, lin, prec);
        }
        arb_poly_clear(lin);

        o = arb_poly_degree(Fg + 0);
        for (gamma = 1; gamma < ngenus && success; gamma++)
            if (arb_poly_degree(Fg + gamma) != o)
                success = 0;
        if (!success)
            break;

        /* numerical square roots: sqrt(p^*) real or i sqrt(-p^*); sqrt(D) = product */
        acb_one(sqrtDnum);
        for (i = 0; i < g; i++)
        {
            arb_set_si(t, FLINT_ABS(pstar[i]));
            arb_sqrt(t, t, prec);
            if (pstar[i] > 0)
                acb_mul_arb(sqrtDnum, sqrtDnum, t, prec);
            else
            {
                acb_mul_arb(sqrtDnum, sqrtDnum, t, prec);
                acb_mul_onei(sqrtDnum, sqrtDnum);
            }
        }

        /* modular square root of D consistent with the numerical one */
        fmpz_one(tmp2);
        for (i = 0; i < g; i++)
            fmpz_mod_mul(tmp2, tmp2, sqrts + i, ctx);

        fmpz_set_ui(inv2, UWORD(1) << (g + 1));
        fmpz_mod_inv(inv2, inv2, ctx);

        fmpz_mod_poly_zero(F, ctx);
        fmpz_mod_poly_set_coeff_ui(F, o, 1, ctx);

        for (k = 0; k < o && success; k++)
        {
            fmpz_t coeff;
            fmpz_init(coeff);

            for (S = 0; S < ngenus && success; S++)
            {
                /* t_S = 2^{-(g-1)} sum_gamma chi_gamma(S) c_{gamma,k} */
                arb_zero(t);
                for (gamma = 0; gamma < ngenus; gamma++)
                {
                    slong code = genus_code[gamma], bitsSg = code & S, par = 0;
                    while (bitsSg)
                    {
                        par ^= 1;
                        bitsSg &= bitsSg - 1;
                    }
                    if (par)
                        arb_sub(t, t, arb_poly_get_coeff_ptr(Fg + gamma, k), prec);
                    else
                        arb_add(t, t, arb_poly_get_coeff_ptr(Fg + gamma, k), prec);
                }
                arb_mul_2exp_si(t, t, -(g - 1));

                /* b_S numerically and modulo n; den = prod_{i in S} |p_i^*| */
                acb_one(bS);
                fmpz_one(tmp);
                den = 1;
                for (i = 0; i < g - 1; i++)
                {
                    if (S & (WORD(1) << i))
                    {
                        arb_set_si(u, FLINT_ABS(pstar[i]));
                        arb_sqrt(u, u, prec);
                        acb_mul_arb(bS, bS, u, prec);
                        if (pstar[i] < 0)
                            acb_mul_onei(bS, bS);
                        fmpz_mod_mul(tmp, tmp, sqrts + i, ctx);
                        den *= FLINT_ABS(pstar[i]);
                    }
                }

                /* kappa_S = t_S / b_S = (U + V sqrt(D)) / (2^{g+1} den): dividing
                   an algebraic integer by b_S brings in the primes of S */
                acb_set_arb(kappa, t);
                acb_div(kappa, kappa, bS, prec);
                arb_set(u, acb_realref(kappa));
                arb_div(v, acb_imagref(kappa), acb_imagref(sqrtDnum), prec);
                arb_mul_2exp_si(u, u, g + 1);
                arb_mul_2exp_si(v, v, g + 1);
                arb_mul_si(u, u, den, prec);
                arb_mul_si(v, v, den, prec);
                if (!arb_get_unique_fmpz(U, u) || !arb_get_unique_fmpz(V, v))
                {
                    success = 0;
                    break;
                }

                /* coefficient += b_S (U + V sqrt(D)) / (2^{g+1} den) */
                fmpz_mod_set_fmpz(U, U, ctx);
                fmpz_mod_set_fmpz(V, V, ctx);
                fmpz_mod_mul(V, V, tmp2, ctx);
                fmpz_mod_add(U, U, V, ctx);
                fmpz_mod_mul(U, U, inv2, ctx);
                fmpz_set_si(V, den);
                fmpz_mod_inv(V, V, ctx);
                fmpz_mod_mul(U, U, V, ctx);
                fmpz_mod_mul(U, U, tmp, ctx);
                fmpz_mod_add(coeff, coeff, U, ctx);
            }

            fmpz_mod_poly_set_coeff_fmpz(F, k, coeff, ctx);
            fmpz_clear(coeff);
        }

        if (success)
            break;

        /* not enough precision: retry, at most a few times */
        if (prec > 20 * (3.141593 * sqrt((double) -D) * lgmax * 1.442696 + 100))
            break;
        success = 1;
        prec = prec * 1.5 + 50;
    }

    for (gamma = 0; gamma < ngenus; gamma++)
        arb_poly_clear(Fg + gamma);
    flint_free(Fg);
    arb_clear(sqrtD);
    acb_clear(z); acb_clear(bS); acb_clear(sqrtDnum); acb_clear(kappa);
    arb_clear(t); arb_clear(u); arb_clear(v);
    fmpz_clear(U); fmpz_clear(V); fmpz_clear(inv2); fmpz_clear(tmp); fmpz_clear(tmp2);
    flint_free(forms);
    flint_free(genus_of);
    flint_free(genus_code);

    return success;
}
