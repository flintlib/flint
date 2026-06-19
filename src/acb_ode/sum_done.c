#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_poly.h"
#include "acb_ode.h"
#include "arb_poly.h"

static int
sol_has_ini_past(const acb_ode_sol_t sol, slong n, const acb_ode_group_t group)
{
    for (slong s = 0; s < group->nshifts; s++)
        if (group->shifts[s].n >= n)
            for (slong k = 0; k < group->shifts[s].mult; k++)
                if (!acb_is_zero(acb_mat_entry(sol->extini, s, k)))
                    return 1;
    return 0;
}

static void
_mag_vec_inf(mag_ptr mag, slong len)
{
    for (slong i = 0; i < len; i++)
        mag_inf(mag + i);
}

void
acb_ode_cvest_println(const acb_ode_cvest_t cvest)
{
    flint_printf("neglogterm=%f cvg_rate=%f terms_wanted=%f accuracy=%ld "
                 "loss_rate=%f terms_full_prec=%f prec_wanted=%ld\n",
                 cvest->neglogterm, cvest->cvg_rate, cvest->terms_wanted,
                 cvest->accuracy, cvest->loss_rate, cvest->terms_full_prec,
                 cvest->prec_wanted);
}


/* todo: support the case when only part of the residual is known */
static void
bound_normalized_residuals(arb_poly_struct * nres_maj,
                           acb_ode_sum_struct * sum,
                           slong prec)
{
    acb_t cst;
    arf_t tmparf;
    acb_poly_t ind_nd;
    TMP_INIT;
    TMP_START;

    acb_init(cst);
    arf_init(tmparf);

    acb_poly_struct * dop_lc = sum->dop + sum->dop_len - 1;

    acb_poly_init(ind_nd);
    slong nlogs0 = acb_ode_sum_max_nlogs(sum);
    acb_ptr rhs_term = TMP_ALLOC(sizeof(acb_struct) * nlogs0);
    acb_ptr nres_term = _acb_vec_init(nlogs0);

    for (slong m = 0; m < sum->nsols; m++)
    {
        arb_poly_fit_length(nres_maj + m, sum->dop_clen);
        _arb_poly_set_length(nres_maj + m, sum->dop_clen);
    }

    for (slong d = 0, s = 0; d < sum->dop_clen - 1; d++)
    {
        slong n = sum->n;

        /* Multiplicity of λ+n+d as an indicial root */
        slong mult = 0;
        for (; s < sum->dop_len - 1; s++) {
            if (sum->group->shifts[s].n == n)
            {
                mult = sum->group->shifts[s].mult;
                break;
            }
        }

        /* ind_nd/cst = jet of monic indicial polynomial at λ+n+d */
        // XXX change acb_ode_sum_max_nlogs to take k (=n+d here) or k - n instead of k - n0
        slong nlogs_d = acb_ode_sum_max_nlogs_xn(sum, n - sum->n0 + d);

        acb_ode_poly_taylor_shift_aps_trunc(ind_nd, sum->ind, sum->group->leader,
                                            n, nlogs_d + mult, prec);
        acb_set_round(cst, dop_lc->coeffs, prec);
        acb_neg(cst, cst);

        acb_poly_fit_length(ind_nd, nlogs_d + mult);
        for (slong m = 0; m < sum->nsols; m++)
        {
            // XXX skip solutions that have converged?

            // XXX tighten?
            slong nlogs_dm = FLINT_MIN(sum->sol[m].nlogs, nlogs_d);
            // flint_printf("n=%ld sol[%ld].nlogs=%ld nlogs_d=%ld\n", sum->n, m, sum->sol[m].nlogs, nlogs_d);

            /* nres_term = coeff of x^{λ+n+d} in normalized residual of sol m */
            for (slong k = 0; k < nlogs_dm; k++)
                rhs_term[k] = (sum->sol[m].series + k)->coeffs[n + d - sum->n0];
            _acb_ode_poly_negdivrevhigh(nres_term, ind_nd->coeffs + mult,
                                        cst, rhs_term, nlogs_dm, prec);

            /* nres_maj[m] = common majorant of the coefficients of log(x)^k/k!
             * in the normalized residual of solution m */
            for (slong k = 0; k < nlogs_dm; k++)
            {
                arf_ptr c = arb_midref((nres_maj + m)->coeffs + d);
                acb_get_abs_ubound_arf(tmparf, nres_term + k, prec);
                arf_max(c, c, tmparf);
            }
        }
    }

    for (slong m = 0; m < sum->nsols; m++)
        _arb_poly_normalise(nres_maj + m);

    _acb_vec_clear(nres_term, nlogs0);
    acb_poly_clear(ind_nd);
    acb_clear(cst);
    arf_clear(tmparf);
    TMP_END;
}


static void
update_tail_bounds(acb_ode_sum_struct * sum,
                   acb_ode_bound_t bound, acb_ode_group_bound_t gbound,
                   slong prec)
{
    arb_t rad;
    arb_poly_t jet_tb;
    TMP_INIT;
    TMP_START;

    arb_init(rad);
    arb_poly_init(jet_tb);
    arb_poly_struct * nres_maj = TMP_ALLOC(sizeof(arb_poly_t) * sum->nsols);
    for (slong m = 0; m < sum->nsols; m++)
        arb_poly_init(nres_maj + m);

    /* ore_algebra uses a bound on nlogs:
    slong nlogs = acb_ode_group_nlogs(sum->group, sum->n);
    */

    // TODO Only call precompute_integrand if the current bounds are not tight
    // enough or nlogs has changed. (We can also use the larger nlogs from the
    // start, but this leads to worse bounds before nlogs increases.)
    acb_ode_bound_precompute_integrand(sum->itg_pol, sum->itg_num, sum->group,
                                       bound, gbound, sum->n,
                                       acb_ode_sum_max_nlogs(sum), MAG_BITS);

    bound_normalized_residuals(nres_maj, sum, prec);

    /* ore_algebra uses a single nres_maj for all solutions:

    for (slong m = 1; m < sum->nsols; m++)
    {
        arb_poly_fit_length(nres_maj, nres_maj[m].length);
        for (slong i = 0; i < nres_maj[m].length; i++)
            arb_max(nres_maj->coeffs + i, nres_maj->coeffs + i,
                    nres_maj[m].coeffs + i, prec);
        nres_maj->length = FLINT_MAX(nres_maj->length, nres_maj[m].length);
    }
    _arb_poly_normalise(nres_maj);
    */

    for (slong m = 0; m < sum->nsols; m++)
    {
        slong i;

        /* tail_bound_jet only yields a bound for n past the nonzero initial
           values of the solution */
        if (sol_has_ini_past(sum->sol + m, sum->n, sum->group))
        {
            _mag_vec_inf((sum->sol + m)->tb, sum->nder);
            continue;
        }

        arf_set_mag(arb_midref(rad), sum->mag);
        // flint_printf("n=%ld m=%ld nres_maj=%{arb_poly}\n", sum->n, m, nres_maj + m);
        acb_ode_tail_bound_jet_precomp(jet_tb, bound, sum->n, sum->itg_pol,
                                       sum->itg_num, nres_maj + m, rad,
                                       sum->nder, MAG_BITS);
        for (i = 0; i < jet_tb->length; i++)
            arb_get_mag((sum->sol + m)->tb + i, jet_tb->coeffs + i);
        for (; i < sum->nder; i++)
            mag_zero((sum->sol + m)->tb + i);
    }

    arb_poly_clear(jet_tb);
    arb_clear(rad);
    for (slong m = 0; m < sum->nsols; m++)
        arb_poly_clear(nres_maj + m);
    TMP_END;
}

/* stride = terms since last check */
int
acb_ode_sum_done(acb_ode_sum_struct * sum, slong stride,
                 /* in the future sum_done may refine bound, gbound */
                 acb_ode_bound_t bound, acb_ode_group_bound_t gbound,
                 slong prec)
{
    mag_t eps, est, sum_mag, sum_rad;
    mag_init(eps);
    mag_init(est);
    mag_init(sum_mag);
    mag_init(sum_rad);

    int stop = 1;
    int have_tail_bounds = 0;

    slong deg = sum->dop_clen - 1;
    FLINT_ASSERT(sum->n >= sum->n0 + deg);  /* todo: raise this limitation */

    /* detect solutions that have converged or will not improve */

    for (slong m = 0; m < sum->nsols; m++)
    {
        acb_ode_sol_struct * sol = sum->sol + m;

        /* todo: non-rigorous check with less than deg terms?
           (maybe using known terms from the previous residuals in addition to
           the available coefficients) */
        slong off = sum->n - deg - sum->n0;
        slong accuracy = acb_ode_sol_estimate_terms(est, sol, off, deg,
                                                    sum->magpow);
        slong sum_accuracy = acb_ode_sol_estimate_sums(sum_mag, sum_rad, sol);
        acb_ode_cvest_update(sol->cvest, sol->cvest,
                             est, accuracy, stride, prec, sum->wp);

        /* Attempt to detect cases where it is better to bail out even though,
           technically, we may still be able to improve the results (typically
           entire functions with a large hump to climb). XXX Not sure if this
           belongs here; maybe find a better criterion. */
        if (mag_cmp_2exp_si(sum_rad, 4 * prec) >= 0
            && (sum_accuracy < 0 || sum->n > prec * FLINT_BIT_COUNT(prec)))
        {
            mag_inf(eps);
        }
        /* sum_accuracy may be <= 0 and still have a chance to improve,
           typically when the terms are still increasing.

           If accuracy <= 0, the approximation of the sum is no longer going to
           get better, but the error bound still might (e.g., if the evaluation
           point contains zero), roughly until the radius of the partial sum
           exceeds the tail bound. */
        else if (accuracy >= 0 || sum_accuracy >= 0)
        {
            /* stopping condition: ~prec correct bits OR working precision
               insufficient to do better */
            mag_mul_2exp_si(eps, sum_mag, -prec);
            mag_max(eps, eps, sum_rad);
            mag_mul_2exp_si(eps, eps, -8);
        }
        else
        {
            /* todo: The point about the bound having a chance to improve still
               holds. A better heuristic (maybe based on the rate at which the
               tail bounds are improving?) should improve evaluations on wide
               intervals (typically containing 0). */
            mag_inf(eps);
        }

        // flint_printf("n=%ld ", sum->n); acb_ode_cvest_println(sol->cvest);
        // if (sum->n % (10 * stride) == 0)
        //     flint_printf("n=%ld prec=%ld m=%ld est=%{mag} accuracy=%ld sum_mag=%{mag} sum_rad=%{mag} sum_accuracy=%ld eps=%{mag}\n", sum->n, prec, m, est, accuracy, sum_mag, sum_rad, sum_accuracy, eps);

        /* quick non-rigorous check */
        if (mag_cmp(est, eps) > 0 && accuracy >= 0)
        {
            stop = 0;
            // flint_printf("cont\n");
            continue;
        }

        if (sum->flags & ACB_ODE_APPROX)
        {
            sol->done = 1;
            // flint_printf("done (approx)\n");
            continue;
        }

        /* rigorous check */
        if (!have_tail_bounds)
        {
            // XXX wasteful when only some sols have converged
            update_tail_bounds(sum, bound, gbound, MAG_BITS);
            have_tail_bounds = 1;
        }

        slong max_shift = sum->group->shifts[sum->group->nshifts - 1].n;

        // XXX also check derivatives?
        if (sum->nder > 0 && mag_cmp(sol->tb + 0, eps) <= 0)
        {
            // flint_printf("done\n");
            sol->done = 1;  // XXX currently unused
        }
        else if (sum_accuracy == -ARF_PREC_EXACT)
        {
            sol->done = 1;  // XXX currently unused
        }
        /* Crude heuristic to detect situations where the tail bound might as
           well be infinite. XXX Should somehow take into account the rate at
           which the bounds improve, and maybe also indicial roots from other
           groups. */
        else if (sum->n > max_shift + sum->dop_clen
                 && mag_get_d_log2_approx(sol->tb + 0)
                        > mag_get_d_log2_approx(eps) + prec)
        {
            sol->done = 1;  // XXX currently unused
        }
        else
        {
            // flint_printf("cont\n");
            stop = 0;
        }

        if (stop || sum->n % (10 * stride) == 0)
            flint_printf("n=%ld m=%ld tb=%{mag*}\n", sum->n, m, sol->tb, sol->nder);
    }

    mag_clear(eps);
    mag_clear(est);
    mag_clear(sum_mag);
    mag_clear(sum_rad);

    return stop;
}
