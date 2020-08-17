/*
    Copyright 2009 Jason Moxham
    Copyright 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

int main(void)
{
    ulong xn, yn, b, zn, c, dn;
    gmp_randstate_t rands;
    int k, cc;
    mp_limb_t xp[10000], dp[10000], qp[10000], yp[10000];
    mp_limb_t rp[10000], zp[10000], tp[10000], tb;
    int result = 1;
    
    FLINT_TEST_INIT(state);

    gmp_randinit_default(rands);    
 
    flint_printf("mulmod_2expp1_basecase....");
    fflush(stdout);

    b = 1;
    tb = 1;
    tb <<= b;
    
    for ( ; b < 600; b++, tb*=2) {
       xn = BITS_TO_LIMBS(b);
       k = xn*GMP_NUMB_BITS - b;
       
       if (tb == 0 || tb > GMP_NUMB_MASK) 
          tb=1;
       
       mpn_zero(dp, xn);
       mpn_com_n(dp, dp, xn);
       dp[xn-1] &= GMP_NUMB_MASK >> k; /* dp is 2^b-1 */
       dn = xn;
       dp[xn] = mpn_add_1(dp, dp, xn, 2);
       
       if (dp[dn] != 0)
          dn++; /* dp is 2^b+1 */ 
       
       for (c = 0; c < 20; c++) {
          mpn_random2(xp, xn);
          mpn_random2(yp, xn);
          xp[xn-1] &= GMP_NUMB_MASK >> k;
          yp[xn-1] &= GMP_NUMB_MASK >> k;
          mpn_mul_n(zp, xp, yp, xn);
          zn = xn*2;
          MPN_NORM(zp, zn);
          
          if (zn >= dn)
             mpn_tdiv_qr(qp, rp, 0, zp, zn, dp, dn);
          else
             mpn_copyi(rp, zp, dn);
        
          cc = tp[xn] = flint_mpn_mulmod_2expp1_basecase(tp, xp, yp, 0, b, qp);
        
          if (cc != 0 && dn == xn)
             tp[xn-1] |= tb;
        
          result = (mpn_cmp(tp, rp, dn) == 0);
          if (!result) {
             flint_printf("FAIL:\n");
             flint_printf("b = %wd\n", b);
             abort();
          }        
       }
    }

    b = 1;
    tb = 1;
    tb <<= b;

    for ( ; b < 600; b++, tb*=2) {
       xn = BITS_TO_LIMBS(b);
       k = xn*GMP_NUMB_BITS - b;
      
       if (tb == 0 || tb > GMP_NUMB_MASK)
          tb = 1;
      
       mpn_zero(dp, xn);
       mpn_com_n(dp, dp, xn);
       dp[xn-1] &= GMP_NUMB_MASK >> k; /* dp is 2^b-1 */
       dn = xn;
       dp[xn] = mpn_add_1(dp, dp, xn, 2);
       
       if (dp[dn] != 0)
          dn++; /* dp is 2^b+1 */
     
       for (c = 0; c < 20; c++) {
          mpn_random2(xp, xn);
          mpn_zero(yp, xn); /* set yp to 2^b */
          xp[xn-1] &= GMP_NUMB_MASK >> k;
          yp[xn-1] &= GMP_NUMB_MASK >> k;
          yn = xn;
          
          if (tb == 1)
             yn++;
          
          yp[yn-1] = tb;
         
          mpn_mul(zp, yp, yn, xp, xn);
          zn = xn*2;
          MPN_NORM(zp, zn);
          mpn_zero(yp, xn); /* set yp to 2^b */
           
          if (zn >= dn)
             mpn_tdiv_qr(qp, rp, 0, zp, zn, dp, dn);
          else
             mpn_copyi(rp, zp, dn);
         
          cc = tp[xn] = flint_mpn_mulmod_2expp1_basecase(tp, xp, yp, 1, b, qp);    

          if (cc != 0 && dn == xn)
             tp[xn-1] |= tb; 

          result = (mpn_cmp(tp, rp, dn) == 0);
          if (!result)
          {
             flint_printf("FAIL:\n");
             flint_printf("b = %wd\n", b);
             abort();
          }        
       }
    }

    b = 1;
    tb = 1;
    tb <<= b;

    for ( ; b < 600; b++, tb*=2) {
       xn = BITS_TO_LIMBS(b);
       k = xn*GMP_NUMB_BITS - b;
      
       if (tb == 0 || tb > GMP_NUMB_MASK)
          tb = 1;
    
       mpn_zero(dp, xn);
       mpn_com_n(dp, dp, xn);
       dp[xn-1] &= GMP_NUMB_MASK >> k; /* dp is 2^b-1 */
       dn = xn;
       dp[xn] = mpn_add_1(dp, dp, xn, 2);
      
       if (dp[dn] != 0)
          dn++; /* dp is 2^b+1 */

       for (c = 0; c < 20; c++) {
          mpn_random2(xp, xn);
          mpn_zero(yp, xn); /* set yp to 2^b */
          xp[xn-1] &= GMP_NUMB_MASK >> k;
          yp[xn-1] &= GMP_NUMB_MASK >> k;
          yn = xn;
          
          if (tb == 1)
             yn++;
          
          yp[yn-1] = tb;
         
          mpn_mul(zp, yp, yn, xp, xn);
          zn = xn*2;
          MPN_NORM(zp, zn);
          mpn_zero(yp, xn); /* set yp to 2^b */
        
          if (zn >= dn)
             mpn_tdiv_qr(qp, rp, 0, zp, zn, dp, dn);
          else
             mpn_copyi(rp, zp, dn);
         
          cc = tp[xn] = flint_mpn_mulmod_2expp1_basecase(tp, yp, xp, 2, b, qp);
          
          if (cc != 0 && dn == xn)
             tp[xn-1] |= tb;
         
          result = (mpn_cmp(tp, rp, dn) == 0);
          if (!result) {
             flint_printf("FAIL\n");
             flint_printf("b = %wd\n", b);
             abort();
          }        
       }
    }

    rp[0] = 1;
    mpn_zero(rp + 1, 1000);
    b = 1;
    tb = 1;
    tb <<= b; 

    for ( ; b < 600; b++, tb*=2) {
       xn = BITS_TO_LIMBS(b);
       k = xn*GMP_NUMB_BITS - b;
      
       if (tb == 0 || tb > GMP_NUMB_MASK)
          tb = 1;
    
       mpn_zero(dp, xn);
       mpn_com_n(dp, dp, xn);
       dp[xn-1] &= GMP_NUMB_MASK >> k; /* dp is 2^b-1 */
       dn = xn;
       dp[xn] = mpn_add_1(dp, dp, xn, 2);
       
       if (dp[dn] != 0)
          dn++; /* dp is 2^b+1 */ 
     
       for (c = 0; c < 1; c++) {
          mpn_zero(xp, xn);
          mpn_zero(yp, xn); /* set xp, yp to 2^b */
          xp[xn-1] &= GMP_NUMB_MASK >> k;
          yp[xn-1] &= GMP_NUMB_MASK >> k;
          cc = tp[xn] = flint_mpn_mulmod_2expp1_basecase(tp, yp, xp, 3, b, qp);
         
          if (cc != 0 && dn == xn)
             tp[xn-1] |= tb;
       
          result = (mpn_cmp(tp, rp, dn) == 0);
          if (!result) {
             flint_printf("FAIL\n");
             flint_printf("b = %wd\n", b);
             abort();
          }        
       }
    }

    gmp_randclear(rands);
    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
