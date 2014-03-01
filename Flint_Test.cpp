#include "flint.h"
#include "nmod_mat.h"
#include<stdio.h>
#include<time.h>
#include<unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
using namespace std;
#define MMSI 2

size_t aflint_sprintf(char * s, const char * str, ...)
{
   va_list ap;
   size_t len = strlen(str);
   char * str2 = (char *)flint_malloc(len + 1);
   int w1 = 0, w2 = 0;
   void * w3;
   double d;
   ulong wu;
   slong w;
   int args, floating;
   size_t ret;

   /* deal with first substring */
   size_t n = strcspn(str, "%");
   strncpy(str2, str, n);
   str2[n] = '\0';
   ret = sprintf(s, "%s", str2);
   len -= n;
   str += n;

   va_start(ap, str);
   while (len) /* deal with fmt spec prefixed strings */
   {
      n = strcspn(str + 2, "%") + 2; /* be sure to skip a %% */
      strncpy(str2, str, n);
      str2[n] = '\0';
      switch (str[1])
      {
      case 'w':
         if (str[2] == 'x')
         {
            wu = (ulong) va_arg(ap, ulong);
            ret += sprintf(s + ret, WORD_FMT "x", wu);
            ret += sprintf(s + ret, "%s", str2 + 3);
         } else if (str[2] == 'u')
         {
            wu = (ulong) va_arg(ap, ulong);
            ret += sprintf(s + ret, WORD_FMT "u", wu);
            ret += sprintf(s + ret, "%s", str2 + 3);
         } else if (str[2] == 'd')
         {
            w = (slong) va_arg(ap, slong);
            ret += sprintf(s + ret, WORD_FMT "d", w);
            ret += sprintf(s + ret, "%s", str2 + 3);
         } else
         {
            w = (slong) va_arg(ap, slong);
            ret += sprintf(s + ret, WORD_FMT "d", w);
            ret += sprintf(s + ret, "%s", str2 + 2);
         }
         break;
         case '%':
            ret +=sprintf(s+ret,"%s",str2+1);
            break;
      default: /* pass to sprintf */
         args = parse_fmt(&floating, str2);
         if (args)
         {
            if (args == 3)
               w1 = va_arg(ap, int);
            if (args >= 2)
               w2 = va_arg(ap, int);
            if (floating)
            {
               d = va_arg(ap, double);
               if (args == 2)
                  ret += sprintf(s + ret, str2, w2, d);
               else if (args == 3)
                  ret += sprintf(s + ret, str2, w1, w2, d);
               else
                  ret += sprintf(s + ret, str2, d);
            } else
            {
               w3 = va_arg(ap, void *);
               if (args == 2)
                  ret += sprintf(s + ret, str2, w2, w3);
               else if (args == 3)
                  ret += sprintf(s + ret, str2, w1, w2, w3);
               else
                  ret += sprintf(s + ret, str2, w3);
            }
         } else ret += sprintf(s + ret,"%s",str2); /* zero args */
      }

      len -= n;
      str += n;
   }

   va_end(ap);
   flint_free(str2);

   return ret;
}

int main()
{
    nmod_mat_t A,B,C;
    mp_limb_t mod=1000000007;
    nmod_mat_init(A,MMSI,MMSI,mod);
    nmod_mat_init(B,MMSI,MMSI,mod);
    nmod_mat_init(C,MMSI,MMSI,mod);
    flint_rand_t st;
    flint_randinit(st);
    nmod_mat_randfull(A,st);
    nmod_mat_randfull(B,st);
    flint_randclear(st);

    //nmod_mat_print_pretty(A);flint_printf("\n\n");
    //nmod_mat_print_pretty(B);flint_printf("\n\n");
    nmod_mat_mul(B,A,A);
    nmod_mat_mul(C,B,A);
    nmod_mat_print_pretty(C);flint_printf("\n\n");
    //printf("%f\n",mul_time);
    nmod_mat_pow(C,A,3);
    nmod_mat_print_pretty(C);flint_printf("\n\n");
    nmod_mat_entry(C,0,0)=((nmod_mat_entry(A,0,0)*nmod_mat_entry(A,0,0))%1000000007) + ((nmod_mat_entry(A,0,1)*nmod_mat_entry(A,1,0))%1000000007);
    nmod_mat_entry(C,0,0)=nmod_mat_entry(C,0,0)%1000000007;
    nmod_mat_entry(C,0,1)=((nmod_mat_entry(A,0,0)*nmod_mat_entry(A,0,1))%1000000007) + ((nmod_mat_entry(A,0,1)*nmod_mat_entry(A,1,1))%1000000007);
    nmod_mat_entry(C,0,1)=nmod_mat_entry(C,0,1)%1000000007;
    nmod_mat_entry(C,1,0)=((nmod_mat_entry(A,1,0)*nmod_mat_entry(A,0,0))%1000000007) + ((nmod_mat_entry(A,1,1)*nmod_mat_entry(A,1,0))%1000000007);
    nmod_mat_entry(C,1,0)=nmod_mat_entry(C,1,0)%1000000007;
    nmod_mat_entry(C,1,1)=((nmod_mat_entry(A,1,0)*nmod_mat_entry(A,0,1))%1000000007) + ((nmod_mat_entry(A,1,1)*nmod_mat_entry(A,1,1))%1000000007);
    nmod_mat_entry(C,1,1)=nmod_mat_entry(C,1,1)%1000000007;
    nmod_mat_print_pretty(C);flint_printf("\n\n");
    return 0;
}
