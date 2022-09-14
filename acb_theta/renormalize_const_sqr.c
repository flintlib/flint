
#include "acb_theta.h"

void acb_theta_renormalize_const_sqr(acb_t scal, acb_srcptr th2,
        const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_bad = acb_theta_agm_nb_bad_steps(tau, lowprec);
    slong nb_good;
    acb_ptr roots;
    arf_t rel_err;
    slong n = 1<<g;
    slong k;

    roots = _acb_vec_init(n * nb_bad);
    arf_init(rel_err);

    nb_good = acb_theta_agm_nb_good_steps(rel_err, g, prec);
    for (k = 0; k < nb_bad; k++)
    {
        acb_theta_naive_const(&roots[k*n], tau, lowprec);
    }
    if (nb_bad > 0) /* Renormalize lowprec square roots */
    {
        acb_sqrt(scal, &th2[0], 2*lowprec);
        acb_div(scal, &roots[0], scal, lowprec);
        _acb_vec_scalar_mul(roots, roots, n*nb_bad, scal, lowprec);
    }

    acb_theta_agm(scal, th2, roots, rel_err, nb_bad, nb_good, g, prec);
    acb_inv(scal, scal, prec);

    _acb_vec_clear(roots, n * nb_bad);
    arf_clear(rel_err);    
}
