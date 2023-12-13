/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "flint.h"

/* return number of arguments called for by a specific format specifier */
int parse_fmt(int * floating, const char * fmt)
{
   int args = 1;

   fmt++; /* skip % */

   if (fmt[0] == '%')
      return 0; /* special case, print % */

   if (fmt[0] == ' ' || fmt[0] == '+' || fmt[0] == '-')
      fmt++; /* skip flag */

   if (fmt[0] == '*')
   {
      args++;
      fmt++; /* skip * */
   } else
      while (isdigit((unsigned char) fmt[0]))
         fmt++; /* skip width */

   if (fmt[0] == '.')
   {
      fmt++; /* skip . */
      if (fmt[0] == '*')
      {
         args++;
         fmt++; /* skip * */
      } else
         while (isdigit((unsigned char) fmt[0]))
            fmt++; /* skip precision */
   }

   if (fmt[0] == 'h' || fmt[0] == 'l' || fmt[0] == 'L')
      fmt++; /* skip length */

   if (fmt[0] == 'e' || fmt[0] == 'E' || fmt[0] == 'f' || fmt[0] == 'g' || fmt[0] == 'G')
      (*floating) = 1;
   else
      (*floating) = 0;

   return args;
}

FLINT_WARN_UNUSED int flint_scanf(const char * str, ...)
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
   if (!fread(str2, 1, n, stdin) && n > 0)
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
            ret += scanf(WORD_FMT "x", wu);
            if (!fread(str2 + 3, 1, n - 3, stdin) && n > 3)
               goto cleanup;
         } else if (str[2] == 'u')
         {
            wu = (ulong *) va_arg(ap, ulong *);
            ret += scanf(WORD_FMT "u", wu);
            if (!fread(str2 + 3, 1, n - 3, stdin) && n > 3)
               goto cleanup;
         } else if (str[2] == 'd')
         {
            w = (slong *) va_arg(ap, slong *);
            ret += scanf(WORD_FMT "d", w);
            if (!fread(str2 + 3, 1, n - 3, stdin) && n > 3)
               goto cleanup;
         } else
         {
            w = (slong *) va_arg(ap, slong *);
            ret += scanf(WORD_FMT "d", w);
            if (!fread(str2 + 2, 1, n - 2, stdin) && n > 2)
               goto cleanup;
         }
         break;
      default: /* pass to scanf */
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
                  ret += scanf(str2, w2, d);
               else if (args == 3)
                  ret += scanf(str2, w1, w2, d);
               else
                  ret += scanf(str2, d);
            } else
            {
               w3 = va_arg(ap, void *);
               if (args == 2)
                  ret += scanf(str2, w2, w3);
               else if (args == 3)
                  ret += scanf(str2, w1, w2, w3);
               else
                  ret += scanf(str2, w3);
            }
         } else
         {
            if (!fread(str2, 1, n, stdin) && n > 0) /* zero args */
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
