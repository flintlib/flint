
#include "acb_theta.h"

void acb_theta_renormalize_sqr(acb_t scal_z, acb_t scal_0, acb_srcptr th2_z,
        acb_srcptr th2_0, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_bad = 1 + acb_theta_agm_ext_nb_bad_steps(z, tau, lowprec);
    slong nb_good;
    acb_ptr a;
    acb_ptr th2;
    acb_ptr roots;
    acb_t scal;
    arf_t c, c_ext, e;
    slong n = 1<<g;
    slong k;

    a = _acb_vec_init(2*n);
    th2 = _acb_vec_init(2*n);
    roots = _acb_vec_init(2*n*nb_bad);
    acb_init(scal);
    arf_init(c);
    arf_init(c_ext);
    arf_init(e);

    _acb_vec_set(th2, th2_z, n);
    _acb_vec_set(th2+n, th2_0, n);

    /* Compute lowprec square roots */
    acb_theta_agm_ext_roots(roots, z, tau, nb_bad, prec);
    
    /* Renormalize lowprec square roots */
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

    /* Set convergence rate */
    for (k = 0; k < 2*n; k++)
    {
        acb_sqr(&a[k], &roots[2*(nb_bad-1) + k], lowprec);
    }
    acb_theta_agm_ext_conv_rate(c_ext, c, e, a, g, lowprec);
    nb_good = acb_theta_agm_nb_good_steps(c, e, prec);
    
    acb_theta_agm_ext(scal_z, scal_0, th2, roots, c_ext, e, nb_bad, nb_good, g,
            prec);
    acb_mul(scal_z, scal_z, scal_0, prec);
    acb_inv(scal_z, scal_z, prec);
    acb_inv(scal_0, scal_0, prec);

    _acb_vec_clear(a, 2*n);
    _acb_vec_clear(th2, 2*n);
    _acb_vec_clear(roots, 2*n*nb_bad);
    acb_clear(scal);
    arf_clear(c);
    arf_clear(c_ext);
    arf_clear(e);
}
