/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "flint.h"

int flint_fscanf(FILE * f, const char * str, ...)
{
   va_list ap;
   size_t len = strlen(str);
   char * str2 = flint_malloc(len + 1);
   int * w1 = NULL, * w2 = NULL;
   void * w3;
   double * d;
   ulong * wu;
   slong * w;
   int args, floating;
   int ret;

   /* deal with first substring */
   size_t n = strcspn(str, "%");
   strncpy(str2, str, n);
   str2[n] = '\0';
   ret = 0;
   if (!fread(str2, 1, n, f) && n > 0)
      goto cleanup;
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
            wu = (ulong *) va_arg(ap, ulong *);
            ret += fscanf(f, WORD_FMT "x", wu);
            if (!fread(str2 + 3, 1, n - 3, f) && n > 3)
               goto cleanup;
         } else if (str[2] == 'u')
         {
            wu = (ulong *) va_arg(ap, ulong *);
            ret += fscanf(f, WORD_FMT "u", wu);
            if (!fread(str2 + 3, 1, n - 3, f) && n > 3)
               goto cleanup;
         } else if (str[2] == 'd')
         {
            w = (slong *) va_arg(ap, slong *);
            ret += fscanf(f, WORD_FMT "d", w);
            if (!fread(str2 + 3, 1, n - 3, f) && n > 3)
               goto cleanup;
         } else
         {
            w = (slong *) va_arg(ap, slong *);
            ret += fscanf(f, WORD_FMT "d", w);
            if (!fread(str2 + 2, 1, n - 2, f) && n > 2)
               goto cleanup;
         }
         break;
      default: /* pass to printf */
         args = parse_fmt(&floating, str2);
         if (args) 
         {
            if (args == 3)
               w1 = va_arg(ap, int *);
            if (args >= 2)
               w2 = va_arg(ap, int *);
            if (floating)
            {
               d = va_arg(ap, double *);
               if (args == 2)
                  ret += fscanf(f, str2, w2, d);
               else if (args == 3)
                  ret += fscanf(f, str2, w1, w2, d);
               else
                  ret += fscanf(f, str2, d);
            } else
            {
               w3 = va_arg(ap, void *);
               if (args == 2)
                  ret += fscanf(f, str2, w2, w3);
               else if (args == 3)
                  ret += fscanf(f, str2, w1, w2, w3);
               else
                  ret += fscanf(f, str2, w3);
            }
         } else 
         {
            if (!fread(str2, 1, n, f) && n > 0) /* zero args */
               goto cleanup;
         }
               
      }

      len -= n;
      str += n;
   }

   va_end(ap);

cleanup:
   flint_free(str2);

   return ret;
}
