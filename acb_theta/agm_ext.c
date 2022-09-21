
#include "acb_theta.h"

void
acb_theta_agm_ext(acb_t r, acb_t s, acb_srcptr a, acb_srcptr roots,
        slong nb_bad, slong nb_good, slong g, slong prec)
{  
    acb_ptr v;
    acb_t scal;
    fmpz_t exp;
    arb_t abs;
    arf_t err;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong n = 1<<g;
    slong k;

    v = _acb_vec_init(2*n);
    acb_init(scal);
    fmpz_init(exp);
    arb_init(abs);
    arf_init(err);

    _acb_vec_set(v, a, 2*n);
  
    for (k = 0; k < nb_bad; k++)
    {
	acb_theta_agm_ext_step_bad(v, v, all_roots + k*2*n, g, prec);
    }

    acb_set(scal, &v[0]);
    _acb_vec_scalar_div(v, v, 2*n, scal, prec);
    
    for (k = 0; k < nb_good; k++)
    {
	acb_theta_agm_ext_step_good(v, v, g, prec);
    }

    /* Set regular Borchardt */
    acb_set(s, &v[n]);
    arf_one(err);
    arf_mul_2exp_si(err, err, -prec);
    acb_add_error_arf(s, err);
    acb_mul(s, &v[0], scal, prec);

    /* Set extended Borchardt */
    acb_div(r, &v[0], &v[n], prec);
    acb_abs(abs, r);
    arb_get_ubound_arf(err, arf, lowprec);
    arf_mul_2exp_si(err, err, -prec);
    acb_add_error_arf(r, err);
    
    fmpz_one(exp);
    fmpz_mul_2exp(exp, exp, nb_good + nb_bad);
    acb_pow_fmpz(r, r, exp, prec);
        
    flint_printf("(agm_ext) Reached agms\n");
    acb_printd(r, 10); flint_printf("\n");
    acb_printd(s, 10); flint_printf("\n");

    _acb_vec_clear(v, 2*n);
    acb_clear(scal);
    fmpz_clear(exp);
    arb_clear(abs);
    arf_clear(err);
}
