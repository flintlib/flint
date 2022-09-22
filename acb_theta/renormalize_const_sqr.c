
#include "acb_theta.h"

void acb_theta_renormalize_const_sqr(acb_t scal, acb_srcptr th2,
        const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_bad = 1 + acb_theta_agm_nb_bad_steps(tau, lowprec);
    slong nb_good;
    acb_mat_t w;
    acb_ptr roots;
    acb_ptr a;
    arf_t c, e;
    slong n = 1<<g;
    slong k;

    acb_mat_init(w, g, g);
    roots = _acb_vec_init(n * nb_bad);
    a = _acb_vec_init(n);
    arf_init(c);
    arf_init(e);

    acb_mat_set(w, tau);

    /* Compute lowprec roots */
    for (k = 0; k < nb_bad; k++)
    {
        acb_theta_naive_const(&roots[k*n], w, lowprec);
        acb_mat_scalar_mul_2exp_si(w, w, 1);
    }
    
    /* Renormalize lowprec square roots */
    acb_sqrt(scal, &th2[0], 2*lowprec);
    acb_div(scal, scal, &roots[0], lowprec);
    _acb_vec_scalar_mul(roots, roots, n*nb_bad, scal, lowprec);
    
    /* Set convergence rate */
    for (k = 0; k < n; k++)
    {
        acb_sqr(&a[k], &roots[(nb_bad-1)*n + k], lowprec);
    }
    acb_theta_agm_conv_rate(c, e, a, g, lowprec);    
    nb_good = acb_theta_agm_nb_good_steps(c, e, prec);

    acb_theta_agm(scal, th2, roots, nb_bad, nb_good, g, prec);
    acb_inv(scal, scal, prec);

    acb_mat_clear(w);
    _acb_vec_clear(roots, n * nb_bad);
    _acb_vec_clear(a, n);
    arf_clear(c);
    arf_clear(e);
}
