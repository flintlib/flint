
#include "acb_theta.h"

void
acb_theta_agm_ext_rel_err(arf_t err, const arf_t c2, const arf_t r,
                          slong nb_good, slong prec)
{
    fmpz_t exp;
    arb_t x, y, z;

    fmpz_init(exp);
    arb_init(x);
    arb_init(y);
    arb_init(z);

    if (nb_good < 1)
    {
        flint_printf
            ("acb_theta_agm_ext_rel_err: Error (need at least 1 good step)\n");
        fflush(stdout);
        flint_abort();
    }
    arb_set_arf(x, r);
    arb_mul_2exp_si(x, x, 2);
    arb_sub_si(x, x, 1, prec);
    if (!arb_is_negative(x))
    {
        flint_printf("acb_theta_agm_ext_rel_err: Error (decay too slow)\n");
        arf_printd(r, 10);
        flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_one(exp);
    fmpz_mul_2exp(exp, exp, FLINT_MAX(nb_good - 1, 0));
    arb_set_arf(x, r);
    arb_pow_fmpz(x, x, exp, prec);
    arb_mul_arf(z, x, c2, prec);

    arb_mul_si(y, x, -2, prec);
    arb_add_si(y, y, 1, prec);
    arb_div(x, x, y, prec);

    arb_add_si(y, z, 2, prec);
    arb_sub_si(z, z, 1, prec);
    arb_mul_si(z, z, -2, prec);
    arb_div(y, y, z, prec);
    arb_mul(x, x, y, prec);

    arb_mul_arf(x, x, c2, prec);

    arb_expm1(x, x, prec);
    arb_get_ubound_arf(err, x, prec);

    fmpz_clear(exp);
    arb_clear(x);
    arb_clear(y);
    arb_clear(z);
}
