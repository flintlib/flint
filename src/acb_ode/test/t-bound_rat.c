#include "test_helpers.h"
#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "acb_types.h"
#include "arb.h"
#include "mag.h"

/* based on a draft generated using Claude Opus 4.8 */

static void
ref_values(arb_ptr res,
           acb_poly_struct * indder, acb_poly_struct * numder,
           slong n, slong mult, slong ord, slong prec)
{
    acb_t nn, val, fact;
    acb_poly_t P, G, R;
    arb_t s, ab;

    acb_init(nn);
    acb_init(val);
    acb_init(fact);
    acb_poly_init(P);
    acb_poly_init(G);
    acb_poly_init(R);
    arb_init(s);
    arb_init(ab);

    acb_set_si(nn, n);

    acb_one(fact);
    for (slong i = 0; i < ord; i++)
    {
        if (i > 0)
            acb_mul_si(fact, fact, i, prec);
        acb_poly_evaluate(val, numder + i, nn, prec);
        acb_div(val, val, fact, prec);
        acb_poly_set_coeff_acb(P, i, val);
    }

    acb_one(fact);
    for (slong t = 1; t <= mult; t++)
        acb_mul_si(fact, fact, t, prec);
    for (slong j = 0; j < ord; j++)
    {
        if (j > 0)
            acb_mul_si(fact, fact, mult + j, prec);
        acb_poly_evaluate(val, indder + (mult + j), nn, prec);
        acb_div(val, val, fact, prec);
        acb_poly_set_coeff_acb(G, j, val);
    }

    acb_poly_div_series(R, P, G, ord, prec);

    arb_zero(res + 0);
    for (slong t = 1; t < ord; t++)
    {
        acb_poly_get_coeff_acb(val, R, t - 1);
        acb_abs(ab, val, prec);
        arb_mul_si(ab, ab, n, prec);
        arb_add(res + t, res + t - 1, ab, prec);
    }

    arb_clear(ab);
    arb_clear(s);
    acb_poly_clear(R);
    acb_poly_clear(G);
    acb_poly_clear(P);
    acb_clear(fact);
    acb_clear(val);
    acb_clear(nn);
}


TEST_FUNCTION_START(acb_ode_bound_rat, state)
{
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        // ulong seed1, seed2;
        // flint_rand_get_seed(&seed1, &seed2, state);
        // flint_printf("flint_rand_set_seed(state, 0x%lx, 0x%lx);\n",
        //              seed1, seed2);

        arb_t b, ref_past, ref_asy;
        acb_t junk;
        acb_poly_t ind;
        acb_ode_exponents_t expos;
        acb_ode_ind_lbound_t ind_lbound;
        acb_ode_stairs_t stairs;

        acb_init(junk);
        acb_ode_exponents_init(expos);
        acb_poly_init(ind);
        arb_init(b);
        acb_ode_ind_lbound_init(ind_lbound);
        acb_ode_stairs_init(stairs);
        arb_init(ref_past);
        arb_init(ref_asy);

        acb_randtest(junk, state, 20, 20);

        slong prec = 2 + n_randint(state, 40);

        /* denominator */

        slong elen = 1 + n_randint(state, 8);
        slong disp = n_randint(state, 50);
        acb_ode_exponents_randtest(expos, state, elen, disp, prec, 4);

        slong grp = n_randint(state, expos->ngroups);
        acb_ptr leader0 = expos->groups[grp].leader;

        for (slong g = 0; g < expos->ngroups; g++)
        {
            if (g != grp)
                acb_sub(expos->groups[g].leader, expos->groups[g].leader,
                        leader0, prec);
        }
        acb_zero(leader0);

        acb_ode_group_struct * group = expos->groups + grp;
        slong shiftmax = group->shifts[group->nshifts - 1].n;

        acb_ode_indicial_polynomial_from_exponents(ind, expos, prec);

        /* numerators (we only check the numerator of index i0) */

        slong nnum = n_randint(state, 16) ? 1 : n_randint(state, 4);
        slong i0 = n_randint(state, nnum);
        acb_poly_struct * num = _acb_poly_vec_init(nnum);
        for (slong i = 0; i < nnum; i++)
            acb_poly_randtest(num + i, state, n_randint(state, ind->length),
                              prec, 4);

        /* some data used in the computation of reference values */

        slong ordmax = acb_ode_group_nlogs(group, WORD_MAX) + 1
                       + n_randint(state, 3);

        slong multmax = 0;
        for (slong s = 0; s < group->nshifts; s++)
            multmax = FLINT_MAX(multmax, group->shifts[s].mult);

        acb_poly_struct * indder = _acb_poly_vec_init(ordmax + multmax);
        acb_poly_set(indder + 0, ind);
        for (slong k = 0; k + 1 < ordmax + multmax; k++)
            acb_poly_derivative(indder + k + 1, indder + k, prec);

        acb_poly_struct * numder = _acb_poly_vec_init(ordmax);
        if (nnum > 0)
            acb_poly_set(numder + 0, num + i0);
        for (slong k = 0; k + 1 < ordmax; k++)
            acb_poly_derivative(numder + k + 1, numder + k, prec);

        acb_ode_ind_lbound_precompute(ind_lbound, expos, grp, prec);
        acb_ode_stairs_precompute(stairs, num, nnum, ind, group, ind_lbound,
                                  prec);

        arb_ptr ref_n = _arb_vec_init(ordmax + 1);
        arb_ptr ref_max = _arb_vec_init(ordmax + 1);
        mag_ptr res = _mag_vec_init(nnum);

        for (slong n = 1000 + n_randint(state, 10);
             n >= 0;
             n = n > 2 * shiftmax + 1 ? n / 2 : n - 1)
        {
            if (n_randint(state, 8))
                continue;

            /* reference values */

            slong mult = acb_ode_group_multiplicity(group, n);
            ref_values(ref_n, indder, numder, n, mult, ordmax, prec + 10);

            if (!mult)
                for (slong ord = 0; ord < ordmax; ord++)
                    arb_max(ref_max + ord, ref_max + ord, ref_n + ord, prec);

            /* test acb_ode_bound_rat_ref_vec */

            for (slong ord = 0; ord < ordmax; ord++)
            {
                acb_ode_bound_rat_ref_vec(res, num, nnum, ind, n, mult, ord,
                                          prec);

                if (nnum == 0)
                    continue;

                arb_zero(b);
                mag_set(arb_radref(b), res + i0);
                if (!arb_overlaps(b, ref_n + ord))
                {
                    flint_printf("FAIL (bound_rat_ref_vec)\n\n");
                    flint_printf("expos="); acb_ode_exponents_println(expos);
                    flint_printf("ind = %{acb_poly}\n", ind);
                    flint_printf("num = %{acb_poly}\n", num + i0);
                    flint_printf("n = %wd, ord = %wd\n", n, ord);
                    flint_printf("res = %{mag}\n", res + i0);
                    flint_printf("ref(n) = %{arb}\n", ref_n + ord);
                    flint_abort();
                }
            }

            /* test acb_ode_bound_rat_ordinary_vec */

            if (mult)
                continue;

            for (slong ord = 0; ord < ordmax; ord++)
            {
                acb_ode_bound_rat_ordinary_vec(res, num, nnum, ind, ind_lbound,
                                               n, ord, prec);

                if (nnum == 0)
                    continue;

                arb_zero(b);
                mag_set(arb_radref(b), res + i0);
                if (!arb_overlaps(b, ref_max + ord))
                {
                    flint_printf("FAIL (bound_rat_ordinary_vec)\n\n");
                    flint_printf("expos="); acb_ode_exponents_println(expos);
                    flint_printf("ind = %{acb_poly}\n", ind);
                    flint_printf("num = %{acb_poly}\n", num + i0);
                    flint_printf("n = %wd, ord = %wd\n", n, ord);
                    flint_printf("res = %{mag}\n", res + i0);
                    flint_printf("max ref(n) = %{arb}\n", ref_max + ord);
                    flint_abort();
                }
            }
        }

        _arb_vec_zero(ref_max, ordmax);

        for (slong n = 1000 + n_randint(state, 10);
             n >= 0;
             n = n > 2 * shiftmax + 1 ? n / 2 : n - 1)
        {
            if (n_randint(state, 3))
                continue;

            /* reference values */

            slong mult = acb_ode_group_multiplicity(group, n);

            ref_values(ref_n, indder, numder, n, mult, ordmax, prec + 10);

            slong tau = acb_ode_group_nlogs(group, n);
            arb_max(ref_asy, ref_asy, ref_n + tau, prec);

            if (mult)
            {
                _arb_vec_zero(ref_max, ordmax);
                arb_set(ref_past, ref_asy);
            }
            else
            {
                for (slong ord = 0; ord < ordmax; ord++)
                    arb_max(ref_max + ord, ref_max + ord, ref_n + ord, prec);
            }

            /* test acb_ode_bound_rat_vec */

            for (slong ord = 0; ord < ordmax; ord++)
            {
                acb_swap(leader0, junk);
                acb_ode_bound_rat_vec(res, num, nnum, ind, group, ind_lbound,
                                      stairs, n, ord, prec);
                acb_swap(leader0, junk);

                if (nnum == 0)
                    continue;

                mag_set(arb_radref(b), res + i0);
                if (!(arb_overlaps(b, ref_max + ord)
                      && arb_overlaps(b, ref_past)))
                {
                    flint_printf("FAIL (bound_rat_vec)\n\n");
                    flint_printf("expos="); acb_ode_exponents_println(expos);
                    flint_printf("ind = %{acb_poly}\n", ind);
                    flint_printf("num = %{acb_poly}\n", num + i0);
                    flint_printf("n = %wd, ord = %wd\n", n, ord);
                    flint_printf("res = %{mag}\n", res + i0);
                    flint_printf("max ref(n) up to next shift = %{arb}\n",
                                 ref_max + ord);
                    flint_printf("max ref(n) past next shift = %{arb}\n",
                                 ref_past);
                    flint_abort();
                }
            }
        }

        _mag_vec_clear(res, nnum);
        _acb_poly_vec_clear(numder, ordmax);
        _acb_poly_vec_clear(indder, ordmax + multmax);
        _arb_vec_clear(ref_n, ordmax + 1);
        _arb_vec_clear(ref_max, ordmax + 1);
        arb_clear(ref_past);
        arb_clear(ref_asy);
        acb_ode_stairs_clear(stairs);
        acb_ode_ind_lbound_clear(ind_lbound);
        _acb_poly_vec_clear(num, nnum);
        acb_poly_clear(ind);
        acb_ode_exponents_clear(expos);
        acb_clear(junk);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
