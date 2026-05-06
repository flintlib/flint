#include "arf.h"
#include "acb.h"
#include "acb_ode.h"

void
acb_ode_ind_lbound_init(acb_ode_ind_lbound_t ind_lbound)
{
    ind_lbound->length = 0;
    ind_lbound->r = NULL;
}

void
acb_ode_ind_lbound_clear(acb_ode_ind_lbound_t ind_lbound)
{
    for (slong i = 0; i < ind_lbound->length; i++)
    {
        acb_clear(ind_lbound->r[i].root);
        mag_clear(ind_lbound->r[i].global_lbound);
    }
    flint_free(ind_lbound->r);
}

static void
find_min(slong * n_min, mag_ptr lb, const acb_ode_group_t group,
         const acb_t rt, slong mult, slong prec)
{
    acb_t tmpacb;
    arb_t crit_n;
    arf_t x0, x1;
    mag_t tmpmag;

    acb_init(tmpacb);
    arb_init(crit_n);
    arf_init(x0);
    arf_init(x1);
    mag_init(tmpmag);

    arb_sqr(crit_n, acb_realref(rt), prec);
    arb_addmul(crit_n, acb_imagref(rt), acb_imagref(rt), prec);
    arb_div(crit_n, crit_n, acb_realref(rt), prec);

    if (mag_cmp_2exp_si(arb_radref(crit_n), 3) > 0)
    {
        * n_min = WORD_MAX;
        /* Lower bound valid for all real n > 0 */
        arb_get_mag_lower(lb, acb_imagref(rt));
        acb_get_mag(tmpmag, rt);
        mag_div(lb, lb, tmpmag);
    }
    else
    {
        arb_get_interval_arf(x0, x1, crit_n, prec);
        slong n0 = arf_get_si(x0, ARF_RND_DOWN);
        slong n1 = arf_get_si(x1, ARF_RND_UP);

        /* Extend the interval so that the endpoints are ordinary indices */
        for (slong s = 0; s < group->nshifts;)
        {
            slong exn = group->shifts[s].n;
            if (exn == n0)
                n0--, s = FLINT_MAX(0, s-1);
            else if (exn == n1)
                n1++, s++;
            else
                s++;
        }

        n0 = FLINT_MAX(n0, 1);
        * n_min = n1;

        mag_one(lb);
        for (slong n = n0, s = 0; n <= n1; n++)
        {
            /* Skip exceptional indices (which would lead to lb = 0 and will be
             * handled separately) */
            while (s + 1 < group->nshifts && group->shifts[s].n < n)
                s++;
            if (n == group->shifts[s].n)
                continue;

            /* Evaluate |1-rt/n| */
            acb_div_si(tmpacb, rt, n, prec);
            acb_sub_si(tmpacb, tmpacb, 1, prec);
            acb_get_mag_lower(tmpmag, tmpacb);
            mag_min(lb, lb, tmpmag);
        }
        mag_pow_ui_lower(lb, lb, mult);
    }

    acb_clear(tmpacb);
    arb_clear(crit_n);
    arf_clear(x0);
    arf_clear(x1);
    mag_clear(tmpmag);
}

void
acb_ode_ind_lbound_precompute(acb_ode_ind_lbound_t ind_lbound,
                              const acb_ode_exponents_t expos,
                              slong grp,  // group index in expos
                              slong prec)
{
    acb_t leader, rt;
    /* flint_printf("== acb_ode_ind_lbound_precompute ==\n"); */

    acb_init(leader);
    acb_init(rt);

    slong length = 0;
    for (slong g = 0; g < expos->ngroups; g++)
        length += expos->groups[g].nshifts;

    if (ind_lbound->r != NULL)
        acb_ode_ind_lbound_clear(ind_lbound);
    ind_lbound->r = flint_malloc(sizeof(*(ind_lbound->r)) * length);

    slong i = 0;
    for (slong g = 0; g < expos->ngroups; g++)
    {
        if (g == grp)
            acb_zero(leader);
        else
            acb_sub(leader, expos->groups[g].leader, expos->groups[grp].leader, prec);

        for (slong s = 0; s < expos->groups[g].nshifts; s++)
        {
            /* A root of the shifted indicial polynomial */
            acb_add_si(rt, leader, expos->groups[g].shifts[s].n, prec);

            /* flint_printf("g=%ld s=%ld rt=%{acb}\n", g, s, rt); */

            /* When Re(α) ≤ 0, the sequence |1-α/n| decreases to 1. */
            if (arb_is_nonpositive(acb_realref(rt)))
                continue;

            /* Otherwise, it first decreases to its minimum (which may be 0
             * if α is an integer), then increases to 1. We compute the
             * minimum and a value of n after which the sequence is
             * nondecreasing. */
            acb_init(ind_lbound->r[i].root);
            mag_init(ind_lbound->r[i].global_lbound);

            acb_set(ind_lbound->r[i].root, rt);
            ind_lbound->r[i].mult = expos->groups[g].shifts[s].mult;

            find_min(&ind_lbound->r[i].n_min, ind_lbound->r[i].global_lbound,
                     expos->groups + g, rt, ind_lbound->r[i].mult, prec);

            /* flint_printf("mult=%ld n_min=%ld lbound=%{mag}\n",
                         ind_lbound->r[i].mult, ind_lbound->r[i].n_min,
                         ind_lbound->r[i].global_lbound + i); */

            i++;
        }
    }
    ind_lbound->length = i;

    acb_clear(rt);
    acb_clear(leader);
}
