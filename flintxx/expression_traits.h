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

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

// This file contains helpers recognising expression templates

#ifndef CXX_EXPRESSION_TRAITS_H
#define CXX_EXPRESSION_TRAITS_H

#include "mp.h"
#include "traits.h"

namespace flint {
namespace operations {
// These are the operation tags the expression class creates directly.

// unary operations
struct immediate { };
struct negate { };
struct complement { };

// binary operations
struct plus { };
struct minus { };
struct times { };
struct divided_by { };
struct modulo { };
struct shift { }; // left
struct binary_and { };
struct binary_or { };
struct binary_xor { };
} // operations

namespace traits {
template<class T, class Enable = void>
struct is_expression : mp::false_ { };
template<class T>
struct is_expression<T, typename T::IS_EXPRESSION_MARKER> : mp::true_ { };

template<class T>
struct _is_immediate_expr
    : _is_convertible<
          typename basetype<T>::type::operation_t,
          operations::immediate
        >
{ };

// Compute if T is an expression, with operation "immediate"
template<class T, class Enable = void>
struct is_immediate_expr : _is_immediate_expr<T> { };
template<class T>
struct is_immediate_expr<T,
  typename mp::enable_if<mp::not_<is_expression<T> > >::type>
    : mp::false_ { };

// Compute if T is an immediate expression, *or not an expression at all*
template<class T>
struct is_immediate
    : mp::or_<mp::not_<is_expression<T> >, is_immediate_expr<T> > { };

// Compute if T is a non-immediate expression
template<class T>
struct is_lazy_expr
    : mp::and_<is_expression<T>, mp::not_<is_immediate_expr<T> > > { };

// Compute if Expr is an expression with prescribed evaluated type "T"
template<class Expr, class T, class Enable = void>
struct is_T_expr
    : mp::equal_types<typename Expr::evaluated_t, T> { };
template<class Expr, class T>
struct is_T_expr<Expr, T,
    typename mp::disable_if<traits::is_expression<Expr> >::type>
    : false_ { };

// Decide if an expressing yielding From can be directly evaluated into To.
// To be further specialised!
template<class To, class From, class Enable = void>
struct can_evaluate_into : mp::false_ { };
template<class T>
struct can_evaluate_into<T, T> : mp::true_ { };

// Decide if we should use temporary merging
template<class Expr>
struct use_temporary_merging : mp::true_ { };
} // traits
} // flint

#endif
