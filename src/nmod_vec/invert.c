/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"

void _nmod_vec_invert_naive(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    for (ulong k = 0; k < len; k++)
        res[k] = nmod_inv(vec[k], mod);
}

void _nmod_vec_invert_generic(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    if (len > 1)
    {
        ulong k;

        /* fill tmp with tmp[k] = vec[0] * ... * vec[k] */
        nn_ptr tmp = _nmod_vec_init(len);
        tmp[0] = vec[0];
        for (k = 1; k < len; k++)
            tmp[k] = nmod_mul(vec[k], tmp[k-1], mod);

        /* current inverse is (vec[0] * ... * vec[len-1])^{-1} */
        k = len - 1;
        ulong inv = nmod_inv(tmp[k], mod);

        /* at the end of iteration k, 
         * res[k] = (vec[0] * ... * vec[k-1]) * (vec[0] * ... * vec[k])^{-1}
         *        = vec[k]^{-1}
         * inv = (vec[0] * ... * vec[k-1])^{-1}
         **/
        for (; k > 0; k--)
        {
            res[k] = nmod_mul(inv, tmp[k-1], mod);
            inv = nmod_mul(inv, vec[k], mod);
        }
        res[0] = inv;
        _nmod_vec_clear(tmp);
    }
    else if (len == 1)
    {
        res[0] = nmod_inv(vec[0], mod);
    }
}

void _nmod_vec_invert_shoup(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    if (len > 1)
    {
        ulong k;

        /* fill tmp with 
         *    tmp[3*k+0] = precomp quo for mul by vec[k]
         *    tmp[3*k+1] = vec[0] * ... * vec[k]
         *    tmp[3*k+2] = precomp quo for mul by tmp[3*k+1]
         * temporary rem = precomp rem for mul by vec[k]
         */
        ulong rem;
        nn_ptr tmp = _nmod_vec_init(3*len);

        tmp[0] = n_mulmod_precomp_shoup(vec[0], mod.n);
        tmp[1] = vec[0];
        tmp[2] = tmp[0];
        for (k = 1; k < len; k++)
        {
            n_mulmod_precomp_shoup_quo_rem(tmp+3*k, &rem, vec[k], mod.n);
            n_mulmod_and_precomp_shoup(tmp+3*k+1, tmp+3*k+2, vec[k], tmp[3*k-2], tmp[3*k], rem, tmp[3*k-1], mod.n);
        }

        /* current inverse is (vec[0] * ... * vec[len-1])^{-1} */
        ulong inv = nmod_inv(tmp[3*len-2], mod);

        /* at the end of iteration k, 
         * res[k] = (vec[0] * ... * vec[k-1]) * (vec[0] * ... * vec[k])^{-1}
         *        = vec[k]^{-1}
         * inv = (vec[0] * ... * vec[k-1])^{-1}
         **/
        for (k = len-1; k > 0; k--)
        {
            res[k] = n_mulmod_shoup(tmp[3*k-2], inv, tmp[3*k-1], mod.n);
            inv = n_mulmod_shoup(vec[k], inv, tmp[3*k], mod.n);
        }
        res[0] = inv;
        _nmod_vec_clear(tmp);
    }
    else if (len == 1)
    {
        res[0] = nmod_inv(vec[0], mod);
    }
}

void _nmod_vec_invert(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    if (NMOD_CAN_USE_SHOUP(mod))
        _nmod_vec_invert_shoup(res, vec, len, mod);
    else
        _nmod_vec_invert_generic(res, vec, len, mod);
}
