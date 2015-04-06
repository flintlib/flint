#include<stdio.h>
#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#include "flint.h"
#include "ulong_extras.h"

int n_is_perfect_power(ulong * root, ulong n)
{	
    static unsigned char primes[18] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
	
    static int mod32[32] = {262143, 262143, 0, 262142, 1, 262142, 0, 262142, 2, 262143, 0,
    							262142, 0, 262142, 0, 262142, 1, 262143, 0, 262142, 0, 262142, 0,
    							262142, 2, 262143, 0, 262142, 0, 262142, 0, 262142 };
    							
    static int mod27[27] = {262143, 262143, 262140, 0, 262141, 262140, 0, 262141, 262142,
    							1, 262143, 262140, 0, 262141, 262140, 0, 262141, 262142, 0, 262143,
    							262140, 0, 262141, 262140, 0, 262141, 262142 };
   
							
    /* if we can't find any other answer we just return that n = n ^ 1 */
    
    int t = mod32[n%32];
    if (!t) return 0;

    t &= mod27[n%27];
    if (!t) return 0;
    
    ulong rem;
    
    int i, pivot = 1, power = 1;
    for (i = 0; i < 18; i++, pivot = pivot << 1)
    { 
    	if (t & pivot ){
    		n_rootrem(&rem, n, primes[i]);
    		while (rem == 0 ){
    			n = pow((double) n, 1.0 / ((double) primes[i]) );
    			power *= primes[i];
			n_rootrem(&rem, n, primes[i]);
    		}	
    	}
    }
    
    root = n;
    if( power != 0) return power;
    return 0;
}
