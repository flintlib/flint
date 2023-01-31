
#include "acb_theta.h"

slong
acb_theta_naive_fullprec(const acb_theta_eld_t E, slong prec)
{
    return prec + ceil(ACB_THETA_NAIVE_FULLPREC_ADDLOG
                       * n_flog(1 + acb_theta_eld_nb_pts(E), 2));
}
