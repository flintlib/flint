#include "acb_types.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "acb.h"

void
acb_ode_randtest_acb(acb_poly_struct * dop, acb_struct * lcroots,
                     acb_ode_exponents_t expos,
                     flint_rand_t state,
                     slong disp, slong lcdeg, slong clen, slong len, slong prec)
{
    if (len <= 0)
        return;

    acb_poly_t tmp;
    acb_poly_init(tmp);

    /* start with an Euler equation */

    /* todo: occasionally generate irregular singular operators (but with the
       correct exponents!) */
    slong expos_len = len - 1;

    acb_ode_exponents_randtest(expos, state, expos_len, disp, prec, 3);

    acb_ode_indicial_polynomial_from_exponents(tmp, expos, prec);
    for (slong i = 0; i < tmp->length; i++)
        acb_poly_set_acb(dop + i, tmp->coeffs + i);

    // flint_printf("expos_len=%wd ind=%{acb_poly} expos=", expos_len, tmp);
    acb_ode_exponents_println(expos);

    /* adjust the leading coefficient */

    for (slong i = 0; i < lcdeg; i++)
    {
        if (i > 0 && n_randint(state, 8))
            acb_set(lcroots + i, lcroots + i - 1);
        acb_randtest_precise(lcroots + i, state, prec, 5);
    }
    acb_poly_product_roots(tmp, lcroots, lcdeg, prec);

    acb_ptr cst = tmp->coeffs + 0;
    for (slong i = 0; i < len; i++)
        acb_poly_scalar_mul(dop + i, dop + i, cst, prec);
    acb_poly_swap(dop + len - 1, tmp);

    /* perturb without changing the indicial polynomial */

    if (clen > 1 && n_randint(state, 16))
    {
        for (slong i = 0; i < len - 1; i++)
        {
            acb_poly_randtest(tmp, state, clen - 1, prec, 3);
            acb_poly_shift_left(tmp, tmp, 1);
            acb_poly_add(dop + i, dop + i, tmp, prec);
        }
    }

    acb_poly_clear(tmp);
}


/* todo: gr version, NOT returning the roots and exponents */
