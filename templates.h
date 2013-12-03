/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen
 
******************************************************************************/

#ifndef TEMPLATES_H_
#define TEMPLATES_H_

#define CAT(X,Y) X##_##Y
#define CAT3(X,Y,Z) X##_##Y##_##Z
#define CAT4(X,Y,Z,A) X##_##Y##_##Z##_##A
#define _CAT(X,Y) _##X##_##Y
#define _CAT3(X,Y,Z) _##X##_##Y##_##Z
#define _CAT4(X,Y,Z,A) _##X##_##Y##_##Z##_##A
#define __CAT(X,Y) __##X##_##Y
#define __CAT3(X,Y,Z) __##X##_##Y##_##Z
#define __CAT4(X,Y,Z,A) __##X##_##Y##_##Z##_##A

#define TEMPLATE(X,Y) CAT(X,Y)
#define TEMPLATE3(X,Y,Z) CAT3(X,Y,Z)
#define TEMPLATE4(X,Y,Z,A) CAT4(X,Y,Z,A)
#define _TEMPLATE(X,Y) _CAT(X,Y)
#define _TEMPLATE3(X,Y,Z) _CAT3(X,Y,Z)
#define _TEMPLATE4(X,Y,Z,A) _CAT4(X,Y,Z,A)
#define __TEMPLATE(X,Y) __CAT(X,Y)
#define __TEMPLATE3(X,Y,Z) __CAT3(X,Y,Z)
#define __TEMPLATE4(X,Y,Z,A) __CAT4(X,Y,Z,A)

#define TEMPLATE_STR(x)   #x
#define TEMPLATE_PRINT(x) flint_printf("%s", TEMPLATE_STR(x));
#define TEMPLATE_PRINTF(s, x) flint_printf(s, TEMPLATE_STR(x));

#endif
