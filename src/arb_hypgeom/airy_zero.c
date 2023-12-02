/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
#else
# include <math.h>
#endif

/*
https://dlmf.nist.gov/9.9

a_k ~ -T(3/8 pi (4k-1))
a'_k ~ -U(3/8 pi (4k-3))
b_k ~ -T(3/8 pi (4k-3))
b'_k ~ -U(3/8 pi (4k-1))

For a_k and b_k, the u^8 and u^10 truncations are known to give lower
bounds. [G. Pittaluga and L. Sacripante (1991) Inequalities for the
zeros of the Airy functions. SIAM J. Math. Anal. 22 (1), pp. 260–267.]

We don't have proofs for a'_k and b'_k. However, in that case, we can just
do a single interval Newton step to verify that we have isolated a
zero (the enclosure must be for the correct zero due to sandwiching).
*/

#define AI 0
#define BI 1
#define AI_PRIME 2
#define BI_PRIME 3

static const double initial[4][10] = {{
 -658118728906175.0,-1150655474581104.0,-1553899449042978.0,-1910288501594969.0,
 -2236074816421182.0,-2539650438812533.0,-2826057838960988.0,
 -3098624122012011.0,-3359689702679955.0,-3610979637739094.0},{
 -330370902027041.0,-920730911234245.0,-1359731821477101.0,-1736658984124319.0,
 -2076373934490092.0,-2390271103799312.0,-2684763040788193.0,
 -2963907159065113.0,-3230475233555475.0,-3486466475611047.0},{
 -286764727967452.0,-914286338795679.0,-1356737313209586.0,-1734816794389239.0,
 -2075083421171399.0,-2389296605766914.0,-2683990299959380.0,
 -2963272965051282.0,-3229941298662311.0,-3486008018531685.0},{
 -645827356227815.0,-1146491233835383.0,-1551601459626981.0,-1908764696253222.0,
 -2234961611612173.0,-2538787015856429.0,-2825360342097020.0,
 -3098043823061022.0,-3359196018589429.0,-3610552233837226.0,
}};

void
_arb_hypgeom_airy_zero(arb_t res, const fmpz_t n, int which, slong prec)
{
    slong asymp_accuracy, wp;

    if (fmpz_cmp_ui(n, 10) <= 0)
    {
        if (fmpz_sgn(n) <= 0)
        {
            flint_throw(FLINT_ERROR, "Airy zero only defined for index >= 1\n");
        }

        /* The asymptotic expansions work well except when n == 1, so
           use precomputed starting intervals (also for the first
           few larger n as a small optimization). */
        arf_set_d(arb_midref(res), ldexp(initial[which][fmpz_get_ui(n)-1], -48));
        mag_set_d(arb_radref(res), ldexp(1.0, -48));
        asymp_accuracy = 48;
    }
    else
    {
        arb_t z, u, u2, u4, s;
        fmpz_t c;

        arb_init(z);
        arb_init(u);
        arb_init(u2);
        arb_init(u4);
        arb_init(s);
        fmpz_init(c);

        if (which == AI || which == BI_PRIME)
            asymp_accuracy = 13 + 10 * (fmpz_bits(n) - 1);
        else
        {
            fmpz_sub_ui(c, n, 1);
            asymp_accuracy = 13 + 10 * (fmpz_bits(c) - 1);
        }

        wp = asymp_accuracy + 8;
        /* Reduce precision since we may not need to do any Newton steps. */
        if (which == AI || which == BI)
            wp = FLINT_MIN(wp, prec + 8);

        arb_const_pi(z, wp);
        fmpz_mul_2exp(c, n, 2);
        if (which == AI || which == BI_PRIME)
            fmpz_sub_ui(c, c, 1);
        else
            fmpz_sub_ui(c, c, 3);
        fmpz_mul_ui(c, c, 3);
        arb_mul_fmpz(z, z, c, wp);
        arb_mul_2exp_si(z, z, -3);

        arb_inv(u, z, wp);
        arb_mul(u2, u, u, wp);
        arb_mul(u4, u2, u2, wp);

        if (which == AI || which == BI)
        {
            /* u^8 truncation gives lower bound */
            arb_mul_si(s, u4, -108056875, wp);
            arb_addmul_si(s, u2, 6478500, wp);
            arb_add_si(s, s, -967680, wp);
            arb_mul(s, s, u4, wp);
            arb_addmul_si(s, u2, 725760, wp);
            arb_div_ui(s, s, 6967296, wp);

            /* u^10 term gives upper bound */
            arb_mul(u4, u4, u4, 10);
            arb_mul(u4, u4, u2, 10);
            arb_mul_ui(u4, u4, 486, 10);
        }
        else
        {
            /* u^8 truncation gives upper bound */
            arb_mul_si(s, u4, 18683371, wp);
            arb_addmul_si(s, u2, -1087338, wp);
            arb_add_si(s, s, 151200, wp);
            arb_mul(s, s, u4, wp);
            arb_addmul_si(s, u2, -181440, wp);
            arb_div_ui(s, s, 1244160, wp);

            /* u^10 term gives lower bound */
            arb_mul(u4, u4, u4, 10);
            arb_mul(u4, u4, u2, 10);
            arb_mul_ui(u4, u4, 477, 10);
            arb_neg(u4, u4);
        }

        arb_mul_2exp_si(u4, u4, -1);
        arb_add(s, s, u4, wp);
        arb_add_error(s, u4);

        arb_add_ui(s, s, 1, wp);
        arb_root_ui(z, z, 3, wp);
        arb_mul(z, z, z, wp);
        arb_mul(res, z, s, wp);

        arb_neg(res, res);

        arb_clear(z);
        arb_clear(u);
        arb_clear(u2);
        arb_clear(u4);
        arb_clear(s);
        fmpz_clear(c);
    }

    /* Do interval Newton steps for refinement. Important: for the
       primed zeros, we need to do at least one interval Newton step to
       validate the initial (tentative) inclusion. */
    if (asymp_accuracy < prec || (which == AI_PRIME || which == BI_PRIME))
    {
        arb_t f, fprime, root;
        mag_t C, r;
        slong * steps;
        slong step, extraprec;

        arb_init(f);
        arb_init(fprime);
        arb_init(root);
        mag_init(C);
        mag_init(r);
        steps = flint_malloc(sizeof(slong) * FLINT_BITS);

        extraprec = 0.25 * fmpz_bits(n) + 8;
        wp = asymp_accuracy + extraprec;

        /* C = |f''| or |f'''| on the initial interval given by res */
        /* f''(x) = xf(x) */
        /* f'''(x) = xf'(x) + f(x) */
        if (which == AI || which == AI_PRIME)
            arb_hypgeom_airy(f, fprime, NULL, NULL, res, wp);
        else
            arb_hypgeom_airy(NULL, NULL, f, fprime, res, wp);

        if (which == AI || which == BI)
            arb_mul(f, f, res, wp);
        else
            arb_addmul(f, fprime, res, wp);

        arb_get_mag(C, f);

        step = 0;
        steps[step] = prec;

        while (steps[step] / 2 > asymp_accuracy - extraprec)
        {
            steps[step + 1] = steps[step] / 2;
            step++;
        }

        arb_set(root, res);

        for ( ; step >= 0; step--)
        {
            wp = steps[step] + extraprec;
            wp = FLINT_MAX(wp, arb_rel_accuracy_bits(root) + extraprec);

            /* store radius, set root to the midpoint */
            mag_set(r, arb_radref(root));
            mag_zero(arb_radref(root));

            if (which == AI || which == AI_PRIME)
                arb_hypgeom_airy(f, fprime, NULL, NULL, root, wp);
            else
                arb_hypgeom_airy(NULL, NULL, f, fprime, root, wp);

            /* f, f' = f', xf */
            if (which == AI_PRIME || which == BI_PRIME)
            {
                arb_mul(f, f, root, wp);
                arb_swap(f, fprime);
            }

            /* f'([m+/-r]) = f'(m) +/- f''([m +/- r]) * r */
            mag_mul(r, C, r);
            arb_add_error_mag(fprime, r);
            arb_div(f, f, fprime, wp);
            arb_sub(root, root, f, wp);

            /* Verify inclusion so that C is still valid, and for the
               primed zeros also to make sure that the initial
               intervals really were correct. */
            if (!arb_contains(res, root))
            {
                flint_printf("unexpected: no containment of Airy zero\n");
                arb_indeterminate(root);
                break;
            }
        }

        arb_set(res, root);

        arb_clear(f);
        arb_clear(fprime);
        arb_clear(root);
        mag_clear(C);
        mag_clear(r);
        flint_free(steps);
    }

    arb_set_round(res, res, prec);
}

void
arb_hypgeom_airy_zero(arb_t ai, arb_t aip, arb_t bi, arb_t bip, const fmpz_t n, slong prec)
{
    if (ai != NULL)
        _arb_hypgeom_airy_zero(ai, n, AI, prec);
    if (aip != NULL)
        _arb_hypgeom_airy_zero(aip, n, AI_PRIME, prec);
    if (bi != NULL)
        _arb_hypgeom_airy_zero(bi, n, BI, prec);
    if (bip != NULL)
        _arb_hypgeom_airy_zero(bip, n, BI_PRIME, prec);
}
