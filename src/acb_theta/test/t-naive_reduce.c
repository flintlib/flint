/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("naive_reduce....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: special values of z */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 5);
        slong nb_z = n_randint(state, 10);
        slong bits = n_randint(state, 5);
        slong prec = 100 + n_randint(state, 200);
        acb_mat_t tau;
        arb_mat_t Y, cho;
        acb_ptr z, new_z, c;
        arb_ptr v, offset;
        arb_t pi;
        acb_t t, x;
        slong *n, *zero;
        slong err_exp = - 10 - n_randint(state, 20);
        slong k, j;
        int res;

        acb_mat_init(tau, g, g);
        arb_mat_init(Y, g, g);
        arb_mat_init(cho, g, g);
        z = _acb_vec_init(g * nb_z);
        new_z = _acb_vec_init(g * nb_z);
        c = _acb_vec_init(nb_z);
        v = _arb_vec_init(g * nb_z);
        offset = _arb_vec_init(g);
        arb_init(pi);
        acb_init(t);
        acb_init(x);
        n = flint_malloc(g * nb_z * sizeof(slong));
        zero = flint_malloc(g * sizeof(slong));

        /* Set tau, cho, Y */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_get_imag(cho, tau);
        arb_const_pi(pi, prec);
        arb_mat_scalar_mul_arb(cho, cho, pi, prec);
        arb_mat_cho(cho, cho, prec);
        arb_mat_transpose(cho, cho);
        acb_mat_get_imag(Y, tau);

        /* Test: if z are real, new_z = z, c = 1 and offset = 0 */
        for (k = 0; k < g * nb_z; k++)
        {
            arb_randtest_precise(acb_realref(&z[k]), state, prec, bits);
        }
        acb_theta_naive_reduce(offset, new_z, c, z, nb_z, tau, cho, prec);

        res = 1;
        for (k = 0; k < nb_z; k++)
        {
            res = res && acb_is_one(&c[k]);            
        }

        if (!_arb_vec_is_zero(offset, g)
            || !res
            || !_acb_vec_equal(new_z, z, g * nb_z))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        /* Test: if im(z) = - Y . (even integral vector n) + small error,
           then terms for 2 * n and 0 correspond */
        for (j = 0; j < g; j++)
        {
            zero[j] = 0;
        }
        for (k = 0; k < nb_z; k++)
        {
            for (j = k * g; j < (k + 1) * g; j++)
            {
                n[j] = 2 * n_randint(state, 10);
                arb_set_si(&v[j], n[j]);
            }
            arb_mat_vector_mul_col(v + k * g, Y, v + k * g, prec);
            for (j = k * g; j < (k + 1) * g; j++)
            {
                arb_urandom(acb_imagref(&z[j]), state, prec);
                arb_mul_2exp_si(acb_imagref(&z[j]), acb_imagref(&z[j]), err_exp);
                arb_sub(acb_imagref(&z[j]), acb_imagref(&z[j]), &v[j], prec);
                n[j] *= 2;
            }
        }
        acb_theta_naive_reduce(offset, new_z, c, z, nb_z, tau, cho, prec);

        for (k = 0; k < nb_z; k++)
        {
            acb_theta_naive_term(x, z + k * g, tau, n + k * g, prec);
            acb_theta_naive_term(t, new_z + k * g, tau, zero, prec);
            acb_mul(t, t, &c[k], prec);
            
            if (!acb_overlaps(x, t))
            {
                flint_printf("FAIL (value, k = %wd)\n", k);
                flint_printf("tau:\n");
                acb_mat_printd(tau, 10);
                flint_printf("z:\n");
                _acb_vec_printd(z + k * g, g, 10);
                flint_printf("\nValues:\n");
                acb_printd(x, 10);
                flint_printf("\n");
                acb_printd(t, 10);
                flint_printf("\n");
                acb_printd(&c[k], 10);
                flint_printf("\nNew z:\n");
                _acb_vec_printd(new_z + k * g, g, 10);
                flint_printf("\n");
                flint_abort();
            }
        }
        
        arb_mat_inv(cho, cho, prec);
        arb_mat_vector_mul_col(offset, cho, offset, prec);
        for (k = 0; k < g; k++)
        {
            arb_mul_2exp_si(&offset[k], &offset[k], - err_exp - 2);
            arb_sub_si(&offset[k], &offset[k], 1, prec);
            if (!arb_is_negative(&offset[k]))
            {
                flint_printf("FAIL (offset)\n");
                arb_printd(&offset[k], 10);
                flint_printf("\n");
                flint_abort();
            }
        }
                
        acb_mat_clear(tau);
        arb_mat_clear(Y);
        arb_mat_clear(cho);
        _acb_vec_clear(z, g * nb_z);
        _acb_vec_clear(new_z, g * nb_z);
        _acb_vec_clear(c, nb_z);
        _arb_vec_clear(v, g * nb_z);
        _arb_vec_clear(offset, g);
        arb_clear(pi);
        acb_clear(t);
        acb_clear(x);
        flint_free(n);
        flint_free(zero);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

