
#include "acb_theta.h"

static int
acb_theta_eld_contains_rec(const acb_theta_eld_t E, slong* pt)
{
    slong d = acb_theta_eld_dim(E);
    slong step = acb_theta_eld_step(E);
    slong c = pt[d-1];
    slong k;

    if (c < acb_theta_eld_min(E)
            || c > acb_theta_eld_max(E)
            || (acb_theta_eld_max(E) - c) % acb_theta_eld_step(E) != 0)
    {
        return 0;
    }
    else if (d == 1)
    {
        return 1;
    }
    else if (c >= acb_theta_eld_mid(E))
    {
        k = (c - acb_theta_eld_mid(E))/step;
        return acb_theta_eld_contains_rec(acb_theta_eld_rchild(E, k), pt);
    }
    else
    {
        k = (acb_theta_eld_mid(E) - step - c)/step;
        return acb_theta_eld_contains_rec(acb_theta_eld_lchild(E, k), pt);
    }
}

int
acb_theta_eld_contains(const acb_theta_eld_t E, slong* pt)
{
    slong g = acb_theta_eld_ambient_dim(E);
    slong d = acb_theta_eld_dim(E);
    slong k;

    if (acb_theta_eld_nb_pts(E) == 0) return 0;

    for (k = d; k < g; k++)
    {
        if (pt[k] != acb_theta_eld_coord(E, k)) return 0;
    }

    return acb_theta_eld_contains_rec(E, pt);
}
