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

int flint_fprintf(FILE * f, const char * str, ...)
{
   va_list ap;
   size_t len = strlen(str);
   char * str2 = flint_malloc(len + 1);
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
   ret = fprintf(f, "%s", str2);
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
            ret += fprintf(f, WORD_FMT "x", wu);
            ret += fprintf(f, "%s", str2 + 3);
         } else if (str[2] == 'u')
         {
            wu = (ulong) va_arg(ap, ulong);
            ret += fprintf(f, WORD_FMT "u", wu);
            ret += fprintf(f, "%s", str2 + 3);
         } else if (str[2] == 'd')
         {
            w = (slong) va_arg(ap, slong);
            ret += fprintf(f, WORD_FMT "d", w);
            ret += fprintf(f, "%s", str2 + 3);
         } else
         {
            w = (slong) va_arg(ap, slong);
            ret += fprintf(f, WORD_FMT "d", w);
            ret += fprintf(f, "%s", str2 + 2);
         }
         break;
      default: /* pass to fprintf */
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
                  ret += fprintf(f, str2, w2, d);
               else if (args == 3)
                  ret += fprintf(f, str2, w1, w2, d);
               else
                  ret += fprintf(f, str2, d);
            } else
            {
               w3 = va_arg(ap, void *);
               if (args == 2)
                  ret += fprintf(f, str2, w2, w3);
               else if (args == 3)
                  ret += fprintf(f, str2, w1, w2, w3);
               else
                  ret += fprintf(f, str2, w3);
            }
         } else ret += fprintf(f, "%s", str2); /* zero args */
      }

      len -= n;
      str += n;
   }

   va_end(ap);
   flint_free(str2);

   return (int) ret;
}
