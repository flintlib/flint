#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_poly.h"
#include "acb_ode.h"
#include "arb_poly.h"
#include "arf.h"
#include "fmpz_vec.h"


void
acb_ode_sum_init(acb_ode_sum_t sum, slong dop_len, slong npts,
                 slong nsols, slong nder)
{
    slong dop_order = dop_len - 1;

    sum->dop = _acb_poly_vec_init(dop_len);
    sum->dop_len = dop_len;

    acb_ode_group_init(sum->group, dop_order);
    for (slong s = 0; s < dop_order; s++)
        sum->group->shifts[s].n = -1;
    acb_poly_init(sum->ind);

    sum->sol = flint_malloc(nsols * sizeof(acb_ode_sol_struct));
    /* using dop_order as a bound
     * - for max possible log prec (will need updating to support inhomogeneous
     *   equations),
     * - for number of initial value positions */
    for (slong m = 0; m < nsols; m++)
        acb_ode_sol_init(sum->sol + m, dop_order, dop_order, npts, nder);
    sum->nsols = nsols;

    sum->pts = _acb_vec_init(npts);
    sum->pows = _acb_vec_init(npts * nder);
    sum->shifted_sums = flint_malloc(npts);
    sum->npts = npts;
    sum->nder = nder;

    mag_init(sum->mag);
    mag_init(sum->magpow);
    mag_init(sum->cvrad);
    arb_poly_init(sum->itg_pol);
    arb_poly_init(sum->itg_num);

    sum->binom_n = _fmpz_vec_init(nder);

    sum->flags = 0;
    sum->have_precomputed = 0;
}

void
acb_ode_sum_clear(acb_ode_sum_t sum)
{
    _acb_poly_vec_clear(sum->dop, sum->dop_len);

    acb_ode_group_clear(sum->group);
    acb_poly_clear(sum->ind);

    for (slong m = 0; m < sum->nsols; m++)
        acb_ode_sol_clear(sum->sol + m);
    flint_free(sum->sol);

    _acb_vec_clear(sum->pts, sum->npts);
    _acb_vec_clear(sum->pows, sum->npts * sum->nder);
    flint_free(sum->shifted_sums);

    mag_clear(sum->mag);
    mag_clear(sum->magpow);
    mag_clear(sum->cvrad);
    arb_poly_clear(sum->itg_pol);
    arb_poly_clear(sum->itg_num);

    _fmpz_vec_clear(sum->binom_n, sum->nder);
}

/* Operator */

void
acb_ode_sum_set_diffop(acb_ode_sum_t sum, const acb_poly_t dop, slong dop_len,
                       const mag_t cvrad)
{
    FLINT_ASSERT(dop_len > 0);
    _acb_poly_vec_set(sum->dop, dop, dop_len);
    mag_set(sum->cvrad, cvrad);
}

/* Solution group */

void
acb_ode_sum_set_ordinary(acb_ode_sum_t sum)
{
    acb_zero(sum->group->leader);
    for (slong n = 0; n < sum->dop_len - 1; n++)
    {
        sum->group->shifts[n].n = n;
        sum->group->shifts[n].mult = 1;
    }
}

void
acb_ode_sum_set_group(acb_ode_sum_t sum, const acb_ode_group_t group)
{
    acb_ode_group_set(sum->group, group);
}

/* Initial values */

void
acb_ode_sum_set_ini_echelon(acb_ode_sum_t sum)
{
    for (int m = 0; m < sum->nsols; m++)
        acb_ode_sol_unit_ini(sum->sol + m, m, sum->group->shifts);
}


void
acb_ode_sum_set_ini_highest(acb_ode_sum_t sum)
{
    for (int m = 0; m < sum->nsols; m++)
    {
        acb_mat_zero(sum->sol[m].extini);
        acb_one(acb_mat_entry(sum->sol[m].extini, m,
                              sum->group->shifts[m].mult - 1));
    }
}

/* Evaluation points */

void
acb_ode_sum_set_points(acb_ode_sum_t sum, acb_srcptr pts, slong npts)
{
    _acb_vec_set(sum->pts, pts, npts);
}

void
acb_ode_sum_precompute(acb_ode_sum_t sum)
{
    if (sum->have_precomputed)
        return;

    slong dop_len = sum->dop_len;

    arf_t abspt;
    arf_init(abspt);

    sum->dop_clen = _acb_poly_vec_length(sum->dop, dop_len);

    /* XXX _gr_ore_poly_indicial_polynomial? */
    acb_poly_fit_length(sum->ind, dop_len);
    _acb_poly_set_length(sum->ind, dop_len);
    for (slong i = 0; i < dop_len; i++)
        acb_poly_get_coeff_acb((sum->ind)->coeffs + i, sum->dop + i, 0);
    _acb_poly_normalise(sum->ind);

    for (slong j = 0; j < sum->npts; j++)  /* XXX move? */
    {
        /* We can trade some full-width muls for muls by small integers by
           multiplying everyone by ξ^n, but this optimization is only
           valid when ξ does not contain 0. */
        sum->shifted_sums[j] = !acb_contains_zero(sum->pts + j);
    }

    _acb_vec_get_mag(sum->mag, sum->pts, sum->npts);
    FLINT_ASSERT(mag_cmp(sum->mag, sum->cvrad) <= 0);

    sum->have_precomputed = 1;

    arf_clear(abspt);
}
