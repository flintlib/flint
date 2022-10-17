
#include "acb_theta.h"

static void
agm_ext_get_conv_rate(arf_t c1, arf_t c2, arf_t r, acb_srcptr a,
        slong n, slong prec)
{    
    arb_t eps;
    slong lowprec = ACB_THETA_AGM_LOWPREC;

    arb_init(eps);
    
    acb_theta_agm_max_abs(eps, a, 2*n, lowprec);
    arb_get_ubound_arf(c1, eps, lowprec);
    acb_theta_agm_min_abs(eps, a, 2*n, lowprec);
    arb_get_lbound_arf(c2, eps, lowprec);
    acb_theta_agm_rel_dist(eps, a+n, n, lowprec, prec);
    arb_get_ubound_arf(r, eps, lowprec);

    acb_theta_agm_ext_conv_rate(c1, c2, r, r, c2, c1, lowprec);    

    arb_clear(eps);
}

void
acb_theta_agm_ext(acb_t r, acb_t s, acb_srcptr a, acb_srcptr roots,
        slong nb_bad, slong g, slong prec)
{  
    acb_ptr v;
    acb_t scal;
    acb_t rel;
    arf_t c1, c2, u;
    arf_t err;
    fmpz_t exp;
    slong n = 1<<g;
    slong nb1, nb2;
    slong k;

    v = _acb_vec_init(2*n);
    acb_init(scal);
    acb_init(rel);
    arf_init(c1);
    arf_init(c2);
    arf_init(u);
    arf_init(err);
    fmpz_init(exp);

    _acb_vec_set(v, a, 2*n);
  
    for (k = 0; k < nb_bad; k++)
    {
	acb_theta_agm_ext_step_bad(v, v, roots + k*2*n, g, prec);
    }
    
    /* Get convergence rate */
    agm_ext_get_conv_rate(c1, c2, u, v, n, prec);
    nb1 = acb_theta_agm_nb_good_steps(c1, u, prec);
    
    /* Perform half the steps */
    acb_set(scal, &v[n]);
    _acb_vec_scalar_div(v, v, 2*n, scal, prec);
    
    for (k = 0; k < nb1/2; k++)
    {
	acb_theta_agm_ext_step_good(v, v, g, prec);
    }
    
    /* Readjust convergence rate */
    agm_ext_get_conv_rate(c1, c2, u, v, n, prec);
    nb2 = acb_theta_agm_nb_good_steps(c1, u, prec);
    nb2 = FLINT_MAX(1, nb2);

    /* Perform remaining steps, at least 1 */
    for (k = 0; k < nb2; k++)
    {
        acb_theta_agm_ext_step_good(v, v, g, prec);
    }
    
    /* Set regular Borchardt */
    acb_set(s, &v[n]);
    arf_one(err);
    arf_mul_2exp_si(err, err, -prec);
    acb_one(rel);
    acb_add_error_arf(rel, err);    
    acb_mul(s, s, rel, prec);
    
    /* Set extended Borchardt */
    acb_theta_agm_ext_step_last(r, s, v, g, prec);
    acb_theta_agm_ext_rel_err(err, c2, u, nb2 + 1, prec);
    acb_one(rel);
    acb_add_error_arf(rel, err);
    acb_mul(r, r, rel, prec);

    /* Renormalize */
    acb_mul(s, s, scal, prec);
    fmpz_one(exp);
    fmpz_mul_2exp(exp, exp, nb_bad + nb1/2 + nb2 + 1);
    acb_pow_fmpz(r, r, exp, prec);

    _acb_vec_clear(v, 2*n);
    acb_clear(scal);
    acb_clear(rel);
    arf_clear(c1);
    arf_clear(c2);
    arf_clear(u);
    arf_clear(err);
    fmpz_clear(exp);
}
