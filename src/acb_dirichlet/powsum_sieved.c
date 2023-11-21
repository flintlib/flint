/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "acb_poly.h"

#define POWER(_k) (powers + (((_k)-1)/2) * (len))
#define DIVISOR(_k) (divisors[((_k)-1)/2])

void
acb_dirichlet_powsum_sieved(acb_ptr z, const acb_t s, ulong n, slong len, slong prec)
{
    slong * divisors;
    slong powers_alloc;
    slong i, j, k, ibound, power_of_two, horner_point;
    ulong kprev;
    int critical_line, integer;

    acb_ptr powers;
    acb_ptr t, u, x;
    acb_ptr p1, p2;
    arb_t logk, v, w;

    if (n <= 1)
    {
        acb_set_ui(z, n);
        _acb_vec_zero(z + 1, len - 1);
        return;
    }

    critical_line = arb_is_exact(acb_realref(s)) &&
        (arf_cmp_2exp_si(arb_midref(acb_realref(s)), -1) == 0);

    integer = arb_is_zero(acb_imagref(s)) && arb_is_int(acb_realref(s));

    divisors = flint_calloc(n / 2 + 1, sizeof(slong));
    powers_alloc = (n / 6 + 1) * len;
    powers = _acb_vec_init(powers_alloc);

    ibound = n_sqrt(n);
    for (i = 3; i <= ibound; i += 2)
        if (DIVISOR(i) == 0)
            for (j = i * i; j <= n; j += 2 * i)
                DIVISOR(j) = i;

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);
    x = _acb_vec_init(len);
    arb_init(logk);
    arb_init(v);
    arb_init(w);

    power_of_two = 1;
    while (power_of_two * 2 <= n)
        power_of_two *= 2;
    horner_point = n / power_of_two;

    _acb_vec_zero(z, len);

    kprev = 1;
    acb_dirichlet_powsum_term(x, logk, &kprev, s, 2,
        integer, critical_line, len, prec);

    for (k = 1; k <= n; k += 2)
    {
        /* t = k^(-s) */
        if (DIVISOR(k) == 0)
        {
            acb_dirichlet_powsum_term(t, logk, &kprev, s, k,
                integer, critical_line, len, prec);
        }
        else
        {
            p1 = POWER(DIVISOR(k));
            p2 = POWER(k / DIVISOR(k));

            if (len == 1)
                acb_mul(t, p1, p2, prec);
            else
                _acb_poly_mullow(t, p1, len, p2, len, len, prec);
        }

        if (k * 3 <= n)
            _acb_vec_set(POWER(k), t, len);

        _acb_vec_add(u, u, t, len, prec);

        while (k == horner_point && power_of_two != 1)
        {
            _acb_poly_mullow(t, z, len, x, len, len, prec);
            _acb_vec_add(z, t, u, len, prec);

            power_of_two /= 2;
            horner_point = n / power_of_two;
            horner_point -= (horner_point % 2 == 0);
        }
    }

    _acb_poly_mullow(t, z, len, x, len, len, prec);
    _acb_vec_add(z, t, u, len, prec);

    flint_free(divisors);
    _acb_vec_clear(powers, powers_alloc);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
    _acb_vec_clear(x, len);
    arb_clear(logk);
    arb_clear(v);
    arb_clear(w);
}

