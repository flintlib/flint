
#include "acb_theta.h"

void acb_theta_renormalize_sqr(acb_t scal_z, acb_t scal_0, acb_srcptr th2_z,
        acb_srcptr th2_0, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_bad = acb_theta_agm_ext_nb_bad_steps(z, tau, lowprec);
    slong nb_good;
    acb_mat_t w;
    acb_ptr th2;
    acb_ptr z_0;
    acb_ptr roots;
    acb_t scal;
    arf_t rel_err;
    slong n = 1<<g;
    slong k;

    acb_mat_init(w, g, g);
    th2 = _acb_vec_init(2*n);
    z_0 = _acb_vec_init(2*g);
    roots = _acb_vec_init(2*n*nb_bad);
    acb_init(scal);
    arf_init(rel_err);

    acb_mat_set(w, tau);
    _acb_vec_set(th2, th2_z, n);
    _acb_vec_set(th2+n, th2_0, n);
    
    nb_good = acb_theta_agm_ext_nb_good_steps(rel_err, g, prec);
    for (k = 0; k < nb_bad; k++)
    {
        acb_theta_naive(&roots[2*k*n], z_0, 2, w, lowprec);
        acb_mat_scalar_mul_2exp_si(w, w, 1);
    }
    if (nb_bad > 0) /* Renormalize lowprec square roots */
    {
        acb_sqrt(scal, &th2[n], 2*lowprec);
        acb_div(scal, scal, &roots[n], lowprec);
        _acb_vec_scalar_mul(roots, roots, 2*n*nb_bad, scal, lowprec);

        acb_sqrt(scal, &th2[0], 2*lowprec);
        acb_div(scal, scal, &roots[n], lowprec);
        acb_div(scal, scal, &roots[0], lowprec);
        for (k = 0; k < nb_bad; k++)
        {
            _acb_vec_scalar_mul(&roots[2*n*k], &roots[2*n*k], n, scal, lowprec);
            acb_sqrt(scal, scal, lowprec);
        }
    }
    
    acb_theta_agm_ext(scal_z, scal_0, th2, roots, rel_err, nb_bad, nb_good, g,
            prec);
    acb_mul(scal_z, scal_z, scal_0, prec);
    acb_inv(scal_z, scal_z, prec);
    acb_inv(scal_0, scal_0, prec);

    acb_mat_clear(w);
    _acb_vec_clear(th2, 2*n);
    _acb_vec_clear(z_0, 2*g);
    _acb_vec_clear(roots, 2*n*nb_bad);
    acb_clear(scal);
    arf_clear(rel_err);
}
