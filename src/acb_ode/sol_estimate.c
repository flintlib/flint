#include "acb.h"
#include "acb_ode.h"
#include "arf.h"
#include "mag.h"

slong
acb_ode_sol_estimate_terms(mag_t est,
                           const acb_ode_sol_t sol, slong off, slong len,
                           const mag_t radpow)
{
    mag_t m;

    mag_init(m);

    slong accuracy = ARF_PREC_EXACT;
    mag_zero(est);

    for (slong k = 0; k < sol->nlogs; k++)
    {
        for (slong i = 0; i < len; i++)
        {
            acb_ptr c = (sol->series + k)->coeffs + off + i;

            acb_get_mag(m, c);
            mag_add(est, est, m);

            accuracy = FLINT_MIN(accuracy, acb_rel_accuracy_bits(c));
        }
    }

    mag_mul(est, est, radpow);

    mag_clear(m);
    return accuracy;
}

slong
acb_ode_sol_estimate_sums(mag_t mag, mag_t rad, const acb_ode_sol_t sol)
{
    mag_t tmp;

    mag_init(tmp);

    mag_zero(mag);
    mag_zero(rad);
    slong accuracy = ARF_PREC_EXACT;

    for (slong p = 0; p < sol->npts; p++)
        for (slong k = 0; k < sol->nlogs; k++)
            for (slong i = 0; i < sol->nder; i++)
            {
                acb_ptr c = acb_ode_sol_sum_ptr(sol, p, k, i);

                acb_get_mag(tmp, c);
                mag_max(mag, mag, tmp);

                mag_max(rad, rad, &c->real.rad);
                mag_max(rad, rad, &c->imag.rad);

                accuracy = FLINT_MIN(accuracy, acb_rel_accuracy_bits(c));
            }

    mag_clear(tmp);

    return accuracy;
}
