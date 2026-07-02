#include "acb_poly.h"

void
acb_ode_poly_taylor_shift_aps_trunc(acb_poly_t g,
                                    const acb_poly_t f, acb_srcptr a,
                                    slong n, slong len, slong prec)
{
    acb_poly_fit_length(g, len);

    if (len == 1 && acb_is_zero(a))  /* covers ordinary points */
    {
        acb_ptr c = g->coeffs;
        acb_zero(c);

        for (slong i = acb_poly_degree(f); i >= 0; i--)
        {
            acb_mul_si(c, c, n, prec);
            acb_add(c, c, f->coeffs + i, prec);
        }

        _acb_poly_set_length(g, 1);
        _acb_poly_normalise(g);
    }
    else
    {
        acb_t b;
        acb_init(b);

        acb_add_si(b, a, n, prec);
        acb_poly_taylor_shift(g, f, b, prec);  /* wasteful */
        acb_poly_truncate(g, len);

        acb_clear(b);
    }
}

