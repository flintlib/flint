
#include "acb_theta.h"

void
acb_theta_agm_ext(acb_t r, acb_t s, acb_srcptr a, acb_srcptr roots,
        const arf_t c, const arf_t e, slong nb_bad, slong nb_good, slong g,
        slong prec)
{  
    acb_ptr v;
    acb_t scal;
    acb_t rel;
    fmpz_t exp;
    arf_t err;
    slong n = 1<<g;
    slong k;

    v = _acb_vec_init(2*n);
    acb_init(scal);
    acb_init(rel);
    fmpz_init(exp);
    arf_init(err);

    _acb_vec_set(v, a, 2*n);
  
    for (k = 0; k < nb_bad; k++)
    {
	acb_theta_agm_ext_step_bad(v, v, roots + k*2*n, g, prec);
    }

    acb_set(scal, &v[n]);
    _acb_vec_scalar_div(v, v, 2*n, scal, prec);

    /* Make good steps, at least 1 */
    nb_good = FLINT_MAX(1, nb_good);    
    for (k = 0; k < nb_good; k++)
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
    acb_theta_agm_ext_rel_err(err, c, e, nb_good+1, prec);
    acb_one(rel);
    acb_add_error_arf(rel, err);
    acb_mul(r, r, rel, prec);

    /* Renormalize */
    acb_mul(s, s, scal, prec);
    fmpz_one(exp);
    fmpz_mul_2exp(exp, exp, nb_bad + nb_good);
    acb_pow_fmpz(r, r, exp, prec);
        
    flint_printf("(agm_ext) Reached agms\n");
    acb_printd(r, 10); flint_printf("\n");
    acb_printd(s, 10); flint_printf("\n");

    _acb_vec_clear(v, 2*n);
    acb_clear(scal);
    acb_clear(rel);
    fmpz_clear(exp);
    arf_clear(err);
}
