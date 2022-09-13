
#include "acb_theta.h"

static void
worker_dim0(acb_ptr th, const acb_t term, slong* coords, slong g,
        ulong ab, slong ord, slong prec, slong fullprec)
{
    acb_t x;
    slong sgn;

    acb_init(x);
  
    sgn = acb_theta_dot(ab, coords, g) % 4;
  
    acb_set(x, term);
    if (sgn == 1) acb_mul_onei(x, x);
    else if (sgn == 2) acb_neg(x, x);
    else if (sgn == 3) acb_div_onei(x, x);
    
    acb_add(th, th, x, fullprec);    
    acb_clear(x);
}

void
acb_theta_naive_ind(acb_t th, ulong ab, acb_srcptr z, const acb_mat_t tau,
        slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    arf_t eps;
    acb_t c;
    acb_ptr new_z;
    int all = 0;
    slong ord = 0;
    slong k = 0;
    slong nb = 1;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, 1, g);
    arf_init(eps);
    acb_init(c);
    new_z = _acb_vec_init(g);

    acb_theta_naive_ellipsoid(E, eps, c, new_z, ab, all, ord, z, 1, tau, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, new_z, tau, E, prec);
    
    acb_theta_naive_worker(th, nb, c, eps, E, D, k, ab, ord, prec,
            worker_dim0);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    arf_clear(eps);
    acb_clear(c);
    _acb_vec_clear(new_z, g);
}
