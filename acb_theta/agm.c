
#include "acb_theta.h"

static void
agm_get_conv_rates(arf_t c, arf_t r, acb_srcptr v, slong n, slong prec)
{
    arb_t eps;
    slong lowprec = ACB_THETA_AGM_LOWPREC;

    arb_init(eps);

    acb_theta_agm_rel_dist(eps, v, n, lowprec, prec);
    arb_get_ubound_arf(r, eps, lowprec);
    acb_theta_agm_conv_rate(c, r, r, lowprec);

    arb_clear(eps);
}

void
acb_theta_agm(acb_t res, acb_srcptr a, acb_srcptr roots, slong nb_bad,
              slong g, slong prec)
{
    acb_ptr v;
    acb_t scal;
    arf_t c, r;
    arf_t err;
    slong nb_good;
    slong n = 1 << g;
    slong k;

    v = _acb_vec_init(n);
    acb_init(scal);
    arf_init(c);
    arf_init(r);
    arf_init(err);

    _acb_vec_set(v, a, n);

    for (k = 0; k < nb_bad; k++)
    {
        acb_theta_agm_step_bad(v, v, roots + k * n, g, prec);
    }

    /* Get convergence rate */
    agm_get_conv_rates(c, r, v, n, prec);
    nb_good = acb_theta_agm_nb_good_steps(c, r, prec);

    /* Perform half the steps */
    acb_set(scal, &v[0]);
    _acb_vec_scalar_div(v, v, n, scal, prec);

    for (k = 0; k < nb_good / 2; k++)
    {
        acb_theta_agm_step_good(v, v, g, prec);
    }

    /* Readjust convergence rate */
    agm_get_conv_rates(c, r, v, n, prec);
    nb_good = acb_theta_agm_nb_good_steps(c, r, prec);

    /* Perform remaining steps */
    for (k = 0; k < nb_good - 1; k++)
    {
        acb_theta_agm_step_good(v, v, g, prec);
    }

    if (nb_good > 0)
        acb_theta_agm_step_last(res, v, g, prec);
    else
        acb_set(res, &v[0]);

    /* Rescale, add relative error */
    acb_mul(res, res, scal, prec);
    arf_one(err);
    arf_mul_2exp_si(err, err, -prec);
    acb_one(scal);
    acb_add_error_arf(scal, err);
    acb_mul(res, res, scal, prec);

    _acb_vec_clear(v, n);
    acb_clear(scal);
    arf_clear(c);
    arf_clear(r);
    arf_clear(err);
}
