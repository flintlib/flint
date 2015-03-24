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
   
							
    *root = n;
    /* if we can't find any other answer we just return that n = n ^ 1 */
    
    int t;
    t = mod32[n%32];
    if (!t) return 1;

    t &= mod27[n%27];
    if (!t) return 1;
    
    ulong rem;
    
    int i = 17;
    while(i + 1){
	printf("%d %d %d %d\n",t,1 << i, n, primes[i]);
    	if(t & (1 << i) ){
    		n_rootrem(&rem, n, primes[i]);
    		if(!(rem) ){
    			*root = pow((double) n, 1.0 / ((double) primes[i]) );
    			return primes[i];
    		}	
    	}
    	--i;
    }
    
    return 1;
}
