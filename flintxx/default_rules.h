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

// This file contains default rule implementations

#ifndef CXX_DEFAULT_RULES_H
#define CXX_DEFAULT_RULES_H

#include "mp.h"
#include "expression.h" // because we want to reuse binary_op_helper etc
#include "expression_traits.h"
#include "evaluation_tools.h"

namespace flint {
namespace rules {
// Composite binary operators
// These rules implement binary operators by implementing both arguments
// separately, then performing the operation on the evaluated types by
// instantiating the appropriate rule again.
//
// Hence to evaluate expressions like a + (b + c), it suffices to write
// rules for composition of two immediates.

namespace rdetail {
template<class Op, class Data1, class Data2, class Enable = void>
struct binary_should_enable { };

template<class Op, class Data1, class Data2>
struct binary_should_enable<Op, Data1, Data2,
    typename mp::enable_if<mp::or_<
        traits::is_lazy_expr<Data1>,
        traits::is_lazy_expr<Data2>
      > >::type>
{
    typedef tools::evaluate_2<Data1, Data2> ev2_t;
    typedef typename ev2_t::return1_t return1_t;
    typedef typename ev2_t::return2_t return2_t;
    typedef typename detail::binary_op_helper<return1_t, Op, return2_t>::cond cond;
    typedef mp::enable_if<cond> enable;
};
}

template<bool result_is_temporary, class Op, class Data1, class Data2>
struct evaluation<
    Op, tuple<Data1, tuple<Data2, empty_tuple> >, result_is_temporary, 0,
    typename rdetail::binary_should_enable<Op, Data1, Data2>::enable::type>
{
    typedef tools::evaluate_2<Data1, Data2> ev2_t;
    typedef typename ev2_t::return1_t return1_t;
    typedef typename ev2_t::return2_t return2_t;
    typedef typename ev2_t::temporaries_t temporaries_t;

    typedef detail::binary_op_helper<return1_t, Op, return2_t> binop_helper;
    typedef typename binop_helper::return_t::evaluated_t return_t;
    typedef typename mp::make_tuple<Data1, Data2>::type data_t;

    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output)
    {
        ev2_t ev2(temps, input.first(), input.second());
        *output = binop_helper::make(ev2.get1(), ev2.get2());
    }
};

template<class Op, class Data, bool result_is_temporary>
struct evaluation<Op, tuple<Data, empty_tuple>, result_is_temporary, 0,
    typename mp::enable_if<traits::is_lazy_expr<Data> >::type>
{
    typedef typename mp::find_evaluation<
        typename Data::operation_t, typename Data::data_t,
        true
      >::type ev_t;
    typedef typename ev_t::return_t interm_t;
    typedef typename ev_t::temporaries_t interm_temps_t;

    typedef detail::unary_op_helper<Op, interm_t> unop_helper;
    typedef typename unop_helper::return_t::evaluated_t return_t;
    typedef mp::merge_tuple<
        typename mp::make_tuple<interm_t*>::type,
        interm_temps_t> merger;
    typedef typename merger::type mixed_temporaries_t;

    template<class Ignored>
    struct do_merging : traits::use_temporary_merging<interm_t> { };
    typedef typename mp::if_<do_merging<void>,
            mixed_temporaries_t, interm_temps_t>::type temporaries_t;

    typedef typename mp::make_tuple<Data>::type data_t;

    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output,
        typename mp::enable_if<do_merging<Return> >::type* = 0)
    {
        interm_t* interm = merger::get_first(temps).head;
        ev_t::doit(input.head._data(), merger::get_second(temps), interm);
        *output = unop_helper::make(*interm);
    }

    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output,
        typename mp::disable_if<do_merging<Return> >::type* = 0)
    {
        interm_t interm = rules::instantiate_temporaries<
            typename data_t::head_t, interm_t>::get(input.head);
        ev_t::doit(input.head._data(), temps, &interm);
        *output = unop_helper::make(interm);
    }
};

// Automatically invoke binary_expression or commutative_binary_expression
namespace rdetail {
template<class Expr1, class Op, class Expr2, class Enable = void>
struct inverted_binary_expression
{
  typedef commutative_binary_expression<Expr2, Op, Expr1> wrapped_t;
  typedef typename wrapped_t::return_t return_t;
  template<class Return>
  static void doit(Return& to, const Expr1& e1, const Expr2& e2)
  {
    return wrapped_t::doit(to, e2, e1);
  }
};

template<template<class E1, class O, class E2, class En> class BE,
    class Data1, class Op, class Data2>
struct binary_expr_helper
{
    typedef typename traits::basetype<Data1>::type data1_t;
    typedef typename traits::basetype<Data2>::type data2_t;
    typedef BE<data1_t, Op, data2_t, void> wrapped_t;
    typedef typename wrapped_t::return_t return_t;
    typedef empty_tuple temporaries_t;
    typedef typename mp::make_tuple<Data1, Data2>::type data_t;
    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output)
    {
        wrapped_t::doit(*output, input.first(), input.second());
    }
};
} // rdetail

template<bool result_is_temporary, class Op, class Data1, class Data2>
struct evaluation<
    Op, tuple<Data1, tuple<Data2, empty_tuple> >, result_is_temporary, 0,
    typename mp::enable_if<
        mp::and_<
            traits::is_immediate<typename traits::basetype<Data1>::type>,
            mp::and_<
                traits::is_immediate<typename traits::basetype<Data2>::type>,
                mp::or_<
                    traits::is_implemented<binary_expression<
                        typename traits::basetype<Data1>::type,
                        Op,
                        typename traits::basetype<Data2>::type
                      > >,
                    mp::or_<
                        traits::is_implemented<commutative_binary_expression<
                            typename traits::basetype<Data1>::type,
                            Op,
                            typename traits::basetype<Data2>::type
                          > >,
                        traits::is_implemented<commutative_binary_expression<
                            typename traits::basetype<Data2>::type,
                            Op,
                            typename traits::basetype<Data1>::type
                          > >
                      >
                  >
              >
          >
      >::type>
    : mp::if_<
            traits::is_implemented<binary_expression<
                typename traits::basetype<Data1>::type,
                Op,
                typename traits::basetype<Data2>::type
              > >,
            rdetail::binary_expr_helper<binary_expression, Data1, Op, Data2>,
            typename mp::if_<
                traits::is_implemented<commutative_binary_expression<
                    typename traits::basetype<Data1>::type,
                    Op,
                    typename traits::basetype<Data2>::type
                  > >,
                rdetail::binary_expr_helper<
                    commutative_binary_expression, Data1, Op, Data2>,
                rdetail::binary_expr_helper<
                    rdetail::inverted_binary_expression, Data1, Op, Data2>
              >::type
          >::type
{ };


// Automatically invoke unary_expression
template<bool result_is_temporary, class Op, class Data>
struct evaluation<Op, tuple<Data, empty_tuple>, result_is_temporary, 0,
    typename mp::enable_if<
        traits::is_implemented<
            unary_expression<Op, typename traits::basetype<Data>::type> > >::type>
{
    typedef unary_expression<Op, typename traits::basetype<Data>::type> wrapped_t;
    typedef typename wrapped_t::return_t return_t;
    typedef empty_tuple temporaries_t;
    typedef typename mp::make_tuple<Data>::type data_t;
    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output)
    {
        wrapped_t::doit(*output, input.head);
    }
};

// Instantiating temporaries

namespace rdetail {
template<class T>
struct evaluated_type_pred
{
    template<class Expr>
    struct type : mp::equal_types<typename Expr::evaluated_t, T> { };
};
}

template<class Expr, class T>
struct use_default_temporary_instantiation : mp::true_ { };

template<class Expr, class T, class Enable>
struct instantiate_temporaries
{
    static T get(const Expr& e)
    {
        return T();
    }
};

template<class Expr, class T>
struct instantiate_temporaries<Expr, T, typename mp::enable_if<mp::and_<
    use_default_temporary_instantiation<Expr, T>,
    tools::has_subexpr<rdetail::evaluated_type_pred<T>, Expr> > >::type>
{
    static T get(const Expr& e)
    {
        return tools::find_subexpr<rdetail::evaluated_type_pred<T> >(e)
            .create_temporary();
    }
};

} // rules
} // flint


#endif
