
#include "acb_theta.h"

static void
worker_dim0(acb_ptr th, const acb_t term, slong * coords, slong g,
            ulong ab, slong ord, slong prec, slong fullprec)
{
    acb_t x;
    slong sgn;
    ulong b;
    slong n = 1 << g;

    acb_init(x);

    for (b = 0; b < n; b++)
    {
        sgn = acb_theta_dot(b, coords, g) % 4;

        acb_set(x, term);
        if (sgn == 1)
            acb_mul_onei(x, x);
        else if (sgn == 2)
            acb_neg(x, x);
        else if (sgn == 3)
            acb_div_onei(x, x);

        acb_add(&th[b], &th[b], x, fullprec);
    }

    /*
       flint_printf("(naive) Coords");
       for (b = 0; b < g; b++) flint_printf(" %wd", coords[b]);
       flint_printf(": "); acb_printd(term, 10); flint_printf("\n");
     */

    acb_clear(x);
}

void
acb_theta_naive(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau,
                slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    arf_struct *eps;
    acb_ptr c;
    acb_ptr new_z;
    int all = 0;
    slong ord = 0;
    ulong ab = 0;
    slong nb = 1 << g;
    slong k;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, nb_z, g);
    eps = flint_malloc(nb_z * sizeof(arf_struct));
    for (k = 0; k < nb_z; k++)
        arf_init(&eps[k]);
    c = _acb_vec_init(nb_z);
    new_z = _acb_vec_init(nb_z * g);

    acb_theta_naive_ellipsoid(E, eps, c, new_z, ab, all, ord, z, nb_z,
                              tau, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, new_z, tau, E, prec);

    for (k = 0; k < nb_z; k++)
    {
        acb_theta_naive_worker(&th[k * nb], nb, &c[k], &eps[k], E, D, k, ab,
                               ord, prec, worker_dim0);
    }

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    for (k = 0; k < nb_z; k++)
        arf_clear(&eps[k]);
    flint_free(eps);
    _acb_vec_clear(c, nb_z);
    _acb_vec_clear(new_z, nb_z * g);
}
