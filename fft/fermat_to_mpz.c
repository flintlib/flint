/* 

Copyright 2009, 2011 William Hart. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY William Hart ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL William Hart OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of William Hart.

*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"

void fermat_to_mpz(mpz_t m, mp_limb_t * i, mp_size_t limbs)
{
   mp_limb_signed_t hi;
   
   mpz_realloc(m, limbs + 1);
   flint_mpn_copyi(m->_mp_d, i, limbs + 1);
   hi = i[limbs];
   if (hi < 0L)
   {
      mpn_neg_n(m->_mp_d, m->_mp_d, limbs + 1);
      m->_mp_size = limbs + 1;
      while ((m->_mp_size) && (!m->_mp_d[m->_mp_size - 1])) 
         m->_mp_size--;
      m->_mp_size = -m->_mp_size;
   } else
   {
      m->_mp_size = limbs + 1;
      while ((m->_mp_size) && (!m->_mp_d[m->_mp_size - 1])) 
         m->_mp_size--;
   }
}
