/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <gmp.h>
#include "longlong.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_extras.h"
#include "mpn_mod.h"
#include "gr.h"

static int _gr_lucas_chain(gr_ptr Vm, gr_ptr Vm1, gr_srcptr A, const fmpz_t m, gr_ctx_t ctx)
{
    gr_ptr t;
    slong i, B, mn;
    int status = GR_SUCCESS;
    nn_srcptr md;
    ulong mtmp;
    slong i_limb, i_bit;
    int FLINT_SET_BUT_UNUSED(msgn);

    FMPZ_GET_MPN_READONLY(msgn, mn, md, mtmp, *m);
    B = (mn - 1) * FLINT_BITS + FLINT_BIT_COUNT(md[mn - 1]);

    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);
    gr_method_binary_op_ui sub_ui = GR_BINARY_OP_UI(ctx, SUB_UI);

    GR_TMP_INIT(t, ctx);

    status |= gr_set_ui(Vm, 2, ctx);
    status |= gr_set(Vm1, A, ctx);

    for (i = B - 1; i >= 0; i--)
    {
       i_limb = i / FLINT_BITS;
       i_bit = i % FLINT_BITS;

       if (md[i_limb] & (UWORD(1) << i_bit)) /* 1 in binary repn */
       {
          status |= mul(t, Vm, Vm1, ctx);
          status |= sub(Vm, t, A, ctx);
          status |= mul(t, Vm1, Vm1, ctx);
          status |= sub_ui(Vm1, t, 2, ctx);
       }
       else /* 0 in binary repn */
       {
          status |= mul(t, Vm, Vm1, ctx);
          status |= sub(Vm1, t, A, ctx);
          status |= mul(t, Vm, Vm, ctx);
          status |= sub_ui(Vm, t, 2, ctx);
       }
    }

    GR_TMP_CLEAR(t, ctx);
    return status;
}

void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A,
                                         const fmpz_t m, const fmpz_t n)
{
    gr_ctx_t ctx;
    gr_ptr rVm, rVm1, rA;
    int status = GR_SUCCESS;

    if (fmpz_size(n) == 1)
        gr_ctx_init_nmod(ctx, fmpz_get_ui(n));
    else if (gr_ctx_init_mpn_mod(ctx, n) != GR_SUCCESS)
        gr_ctx_init_fmpz_mod(ctx, n);

    GR_TMP_INIT3(rVm, rVm1, rA, ctx);

    status |= gr_set_fmpz(rA, A, ctx);
    status |= _gr_lucas_chain(rVm, rVm1, rA, m, ctx);
    status |= gr_get_fmpz(Vm, rVm, ctx);
    status |= gr_get_fmpz(Vm1, rVm1, ctx);
    GR_MUST_SUCCEED(status);

    GR_TMP_CLEAR3(rVm, rVm1, rA, ctx);
    gr_ctx_clear(ctx);
    return;
}

void fmpz_lucas_chain_full(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t B,
                                         const fmpz_t m, const fmpz_t n)
{
    fmpz_t t, Q;
    slong i, bits = fmpz_sizeinbase(m, 2);

    fmpz_init(t);
    fmpz_init(Q);
    fmpz_set_ui(Q, 1);
    fmpz_set_ui(Vm, 2);
    fmpz_set(Vm1, A);

    for (i = bits - 1; i >= 0; i--)
    {
       if (fmpz_tstbit(m, i)) /* 1 in binary repn */
       {
          fmpz_mul(t, Vm1, Vm);
          fmpz_submul(t, Q, A);
          fmpz_mod(Vm, t, n);

          fmpz_mul(Vm1, Vm1, Vm1);
          fmpz_mul_ui(t, Q, 2);
          fmpz_mul(t, t, B);
          fmpz_sub(Vm1, Vm1, t);
          fmpz_mod(Vm1, Vm1, n);

          fmpz_mul(Q, Q, Q);
          fmpz_mul(Q, Q, B);
          fmpz_mod(Q, Q, n);
       } else /* 0 in binary repn */
       {
          fmpz_mul(t, Vm, Vm1);
          fmpz_submul(t, Q, A);
          fmpz_mod(Vm1, t, n);

          fmpz_mul(t, Vm, Vm);
          fmpz_submul_ui(t, Q, 2);
          fmpz_mod(Vm, t, n);

          fmpz_mul(Q, Q, Q);
          fmpz_mod(Q, Q, n);
       }
    }

    fmpz_clear(Q);
    fmpz_clear(t);
}

/* Compute U_{2m}, U_{2m + 1} given U_m, U_{m + 1} */
void fmpz_lucas_chain_double(fmpz_t U2m, fmpz_t U2m1, const fmpz_t Um,
                            const fmpz_t Um1, const fmpz_t A, const fmpz_t B,
                            const fmpz_t n)
{
   fmpz_t t, t2;

   fmpz_init(t);
   fmpz_init(t2);

   fmpz_mul_2exp(t, Um1, 1); /* U_m(2U_{m+1) - AU_m) */
   fmpz_submul(t, Um, A);
   fmpz_mul(t, t, Um);

   fmpz_mul(U2m1, Um1, Um1); /* U_{m+1}^2 - BU_m^2 */
   fmpz_mul(t2, Um, Um);
   fmpz_submul(U2m1, t2, B);
   fmpz_mod(U2m1, U2m1, n);

   fmpz_mod(U2m, t, n);

   fmpz_clear(t);
   fmpz_clear(t2);
}

/*
   Compute U_{m + n}, U_{m + n + 1} given U_m, U_{m + 1} and
   U_n, U_{n + 1}
*/
void fmpz_lucas_chain_add(fmpz_t Umn, fmpz_t Umn1, const fmpz_t Um,
                            const fmpz_t Um1, const fmpz_t Un,
                            const fmpz_t Un1, const fmpz_t A, const fmpz_t B,
                            const fmpz_t n)
{
   fmpz_t t, t2;

   fmpz_init(t);
   fmpz_init(t2);

   fmpz_mul(t, Un, A); /* U_nU_{m + 1} - BU_m(AU_n - U_{n + 1})/B */
   fmpz_sub(t, Un1, t);
   fmpz_mul(t, t, Um);
   fmpz_addmul(t, Un, Um1);

   fmpz_mul(Umn1, Un1, Um1); /* U_{n + 1}U_{m + 1} - BU_mU_n */
   fmpz_mul(t2, Um, Un);
   fmpz_submul(Umn1, t2, B);
   fmpz_mod(Umn1, Umn1, n);

   fmpz_mod(Umn, t, n);

   fmpz_clear(t);
   fmpz_clear(t2);
}

/* Compute U_{km}, U_{km + 1} from U_m, U_{m + 1}, k > 0 */
void fmpz_lucas_chain_mul(fmpz_t Ukm, fmpz_t Ukm1,
                        const fmpz_t Um, const fmpz_t Um1,
                         const fmpz_t A, const fmpz_t B, const fmpz_t k,
                         const fmpz_t n)
{
   slong i = 0, b = fmpz_sizeinbase(k, 2);
   fmpz_t t, t1;

   fmpz_init(t);
   fmpz_init(t1);

   fmpz_set(Ukm, Um);
   fmpz_set(Ukm1, Um1);

   while (!fmpz_tstbit(k, i))
   {
      fmpz_lucas_chain_double(Ukm, Ukm1, Ukm, Ukm1, A, B, n);
      i++;
   }

   i++;

   if (i < b)
   {
      fmpz_set(t, Ukm);
      fmpz_set(t1, Ukm1);
   }

   while (i < b)
   {
      fmpz_lucas_chain_double(t, t1, t, t1, A, B, n);
      if (fmpz_tstbit(k, i))
         fmpz_lucas_chain_add(Ukm, Ukm1, Ukm, Ukm1, t, t1, A, B, n);
      i++;
   }

   fmpz_clear(t);
   fmpz_clear(t1);
}

/* Compute U_m, U_{m + 1} from V_m, V_{m + 1} */
void fmpz_lucas_chain_VtoU(fmpz_t Um, fmpz_t Um1,
        const fmpz_t Vm, const fmpz_t Vm1,
        const fmpz_t A, const fmpz_t FLINT_UNUSED(B), const fmpz_t Dinv,
        const fmpz_t n)
{
   fmpz_t t;

   fmpz_init(t);

   fmpz_mul_2exp(t, Vm1, 1); /* (2V_{m + 1} - AV_m) / D */
   fmpz_submul(t, Vm, A);
   fmpz_mul(t, t, Dinv);

   fmpz_set(Um1, Vm);
   fmpz_mod(Um, t, n);

   fmpz_addmul(Um1, Um, A); /* (V_m + AU_m) / 2 */
   if (!fmpz_is_even(Um1))
      fmpz_add(Um1, Um1, n);
   fmpz_tdiv_q_2exp(Um1, Um1, 1);

   fmpz_mod(Um1, Um1, n);

   fmpz_clear(t);
}
