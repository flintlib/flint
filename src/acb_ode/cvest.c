#include <math.h>

#include "acb_ode.h"
#include "arf.h"

void
acb_ode_cvest_init(acb_ode_cvest_t cvest)
{
    cvest->neglogterm = 0.;
    cvest->cvg_rate = 0.;
    cvest->terms_wanted = INFINITY;
    cvest->accuracy = WORD_MAX; // (?)
    cvest->loss_rate = 1.; // (?)
    cvest->terms_full_prec = INFINITY;
    cvest->prec_wanted = ARF_PREC_EXACT;
}

void
acb_ode_cvest_clear(acb_ode_cvest_t cvest)
{
    return;
}

void
acb_ode_cvest_update(acb_ode_cvest_t cvest,  /* aliasing allowed */
                     const acb_ode_cvest_t old,
                     mag_t est, slong accuracy,
                     slong stride,
                     slong prec, slong work_prec)
{
    FLINT_ASSERT(stride > 0);

    acb_ode_cvest_struct _old = * old;

    /* estimated absolute accuracy of the sum in bits based on the magnitude of
       the last few terms */
    cvest->neglogterm = - mag_get_d_log2_approx(est);

    /* running estimate of convergence rate in bits/term,
       exponential discounting */
    double cvg_rate_stride = (cvest->neglogterm - _old.neglogterm)/stride;
    cvest->cvg_rate = _old.cvg_rate/2 + cvg_rate_stride/2;

    /* number of terms still needed for full accuracy at the current rate of
       convergence */
    cvest->terms_wanted = cvest->cvg_rate > 0.
                        ? (prec - cvest->neglogterm)/cvest->cvg_rate
                        : INFINITY;

    /* relative accuracy of the last few terms */
    cvest->accuracy = accuracy;

    /* running estimate of interval growth (lost bits/term),
       exponential discounting
       (work_prec may have decreased since _old.accuracy was computed) */
    double loss_stride = FLINT_MIN(_old.accuracy, work_prec) - cvest->accuracy;
    double loss_rate_stride = FLINT_MAX(0., loss_stride)/stride;
    cvest->loss_rate = _old.loss_rate/2 + loss_rate_stride/2;

    /* estimated number of terms that can be computed before interval radii
       exceed the target error */
    double abs_acc = cvest->neglogterm + cvest->accuracy;
    // todo: double abs_tol = prec - log(sum_est); ?
    cvest->terms_full_prec = (abs_acc - prec)/cvest->loss_rate;

    /* suggested working precision */
    if (cvest->terms_full_prec <= cvest->terms_wanted
        || cvest->loss_rate >= _old.loss_rate * 1.0125)
        cvest->prec_wanted = ARF_PREC_EXACT;
    else if (cvest->terms_wanted > 0)
        cvest->prec_wanted = prec - 0.875 * cvest->neglogterm
                           + cvest->terms_wanted * cvest->loss_rate;
    else
        cvest->prec_wanted = prec + 4;
}
