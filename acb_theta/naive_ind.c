
#include "acb_theta.h"

static void worker_dim0(acb_ptr th, const acb_t term, slong* coords, slong g,
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
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    arf_t epsilon;
    int all = 0;
    int unif = 0;
    slong ord = 0;
    acb_ptr exp_z;
    acb_mat_t lin_powers;
    acb_t cofactor;
    slong fullprec;
    slong g = acb_mat_nrows(tau);
    slong k;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, g);
    arf_init(epsilon);
    exp_z = _acb_vec_init(g);
    acb_mat_init(lin_powers, g, g);
    acb_init(cofactor);

    acb_theta_naive_ellipsoid(E, epsilon, ab, all, unif, ord, z, tau, prec);  
    fullprec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, tau, E, fullprec);

    acb_mat_set(lin_powers, acb_theta_precomp_exp_mat(D));
    for (k = 0; k < g; k++)
    {
        acb_mul_2exp_si(&exp_z[k], &z[k], 1);
        acb_exp_pi_i(&exp_z[k], &exp_z[k], prec);
    }
    acb_one(cofactor);

    acb_zero(th);
    acb_theta_naive_worker_rec(th, lin_powers, E, D, z, cofactor, ab, ord,
            fullprec, fullprec, worker_dim0);
    acb_add_error_arf(th, epsilon);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    arf_clear(epsilon);
    _acb_vec_clear(exp_z, g);
    acb_mat_clear(lin_powers);
    acb_clear(cofactor);  
}
