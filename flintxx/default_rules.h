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

template<class Op, class Args>
struct can_evaluate_tuple : traits::is_implemented<
    typename mp::find_evaluation<Op,
        typename tools::evaluated_args_tuple<Args>::type, true>::type> { };
template<class Op, class Args, class Enable = void>
struct should_evaluate_tuple : can_evaluate_tuple<Op, Args> { };
template<class Op, class Args>
struct should_evaluate_tuple<Op, Args,
    typename mp::disable_if<tools::count_nonimm<Args> >::type> : mp::false_ { };

template<class Op, class Data1, class Data2>
struct binary_should_enable
{
    typedef mp::enable_if<should_evaluate_tuple<
              Op, typename mp::make_tuple<Data1, Data2>::type> > enable;
};
}

template<bool result_is_temporary, class Op, class Data>
struct evaluation<Op, Data, result_is_temporary, 0,
    typename mp::enable_if<rdetail::should_evaluate_tuple<Op, Data> >::type>
{
    typedef tools::evaluate_n<Data> evn_t;
    typedef typename evn_t::evtup_t evtup_t;
    typedef typename evn_t::temporaries_t temporaries_t;
    typedef typename mp::find_evaluation<
              Op, evtup_t, result_is_temporary>::type rule_t;
    typedef typename rule_t::return_t return_t;

    template<class Return>
    static void doit(const Data& input, temporaries_t temps, Return* output)
    {
        evn_t ev(input, temps);
        rule_t::doit(ev.gettuple(), empty_tuple(), output);
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


// Automatically invoke threeary_expression
template<bool result_is_temporary, class Op,
    class Data1, class Data2, class Data3>
struct evaluation<Op, tuple<Data1, tuple<Data2, tuple<Data3, empty_tuple> > >,
    result_is_temporary, 0,
    typename mp::enable_if<
        traits::is_implemented<
            threeary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type> > >::type>
{
    typedef threeary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type> wrapped_t;
    typedef typename wrapped_t::return_t return_t;
    typedef empty_tuple temporaries_t;
    typedef typename mp::make_tuple<Data1, Data2, Data3>::type data_t;
    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output)
    {
        wrapped_t::doit(*output, mp::tuple_get<data_t, 0>::get(input),
                mp::tuple_get<data_t, 1>::get(input),
                mp::tuple_get<data_t, 2>::get(input));
    }
};
// Automatically invoke fourary_expression
template<bool result_is_temporary, class Op,
    class Data1, class Data2, class Data3, class Data4>
struct evaluation<Op, tuple<Data1, tuple<Data2, tuple<Data3, tuple<Data4, empty_tuple> > > >,
    result_is_temporary, 0,
    typename mp::enable_if<
        traits::is_implemented<
            fourary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type,
                typename traits::basetype<Data4>::type> > >::type>
{
    typedef fourary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type,
                typename traits::basetype<Data4>::type> wrapped_t;
    typedef typename wrapped_t::return_t return_t;
    typedef empty_tuple temporaries_t;
    typedef typename mp::make_tuple<Data1, Data2, Data3, Data4>::type data_t;
    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output)
    {
        wrapped_t::doit(*output, mp::tuple_get<data_t, 0>::get(input),
                mp::tuple_get<data_t, 1>::get(input),
                mp::tuple_get<data_t, 2>::get(input),
                mp::tuple_get<data_t, 3>::get(input));
    }
};
// Automatically invoke fiveary_expression
template<bool result_is_temporary, class Op,
    class Data1, class Data2, class Data3, class Data4, class Data5>
struct evaluation<Op, tuple<Data1, tuple<Data2, tuple<Data3, tuple<Data4, tuple<Data5, empty_tuple> > > > >,
    result_is_temporary, 0,
    typename mp::enable_if<
        traits::is_implemented<
            fiveary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type,
                typename traits::basetype<Data4>::type,
                typename traits::basetype<Data5>::type> > >::type>
{
    typedef fiveary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type,
                typename traits::basetype<Data4>::type,
                typename traits::basetype<Data5>::type> wrapped_t;
    typedef typename wrapped_t::return_t return_t;
    typedef empty_tuple temporaries_t;
    typedef typename mp::make_tuple<Data1, Data2, Data3, Data4, Data5>::type data_t;
    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output)
    {
        wrapped_t::doit(*output, mp::tuple_get<data_t, 0>::get(input),
                mp::tuple_get<data_t, 1>::get(input),
                mp::tuple_get<data_t, 2>::get(input),
                mp::tuple_get<data_t, 3>::get(input),
                mp::tuple_get<data_t, 4>::get(input));
    }
};
// Automatically invoke sixary_expression
template<bool result_is_temporary, class Op,
    class Data1, class Data2, class Data3, class Data4, class Data5, class Data6>
struct evaluation<Op, tuple<Data1, tuple<Data2, tuple<Data3, tuple<Data4, tuple<Data5, tuple<Data6, empty_tuple> > > > > >,
    result_is_temporary, 0,
    typename mp::enable_if<
        traits::is_implemented<
            sixary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type,
                typename traits::basetype<Data4>::type,
                typename traits::basetype<Data5>::type,
                typename traits::basetype<Data6>::type> > >::type>
{
    typedef sixary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type,
                typename traits::basetype<Data4>::type,
                typename traits::basetype<Data5>::type,
                typename traits::basetype<Data6>::type> wrapped_t;
    typedef typename wrapped_t::return_t return_t;
    typedef empty_tuple temporaries_t;
    typedef typename mp::make_tuple<Data1, Data2, Data3, Data4, Data5, Data6>::type data_t;
    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output)
    {
        wrapped_t::doit(*output, mp::tuple_get<data_t, 0>::get(input),
                mp::tuple_get<data_t, 1>::get(input),
                mp::tuple_get<data_t, 2>::get(input),
                mp::tuple_get<data_t, 3>::get(input),
                mp::tuple_get<data_t, 4>::get(input),
                mp::tuple_get<data_t, 5>::get(input));
    }
};
// Automatically invoke sevenary_expression
template<bool result_is_temporary, class Op,
    class Data1, class Data2, class Data3, class Data4, class Data5, class Data6, class Data7>
struct evaluation<Op, tuple<Data1, tuple<Data2, tuple<Data3, tuple<Data4, tuple<Data5, tuple<Data6, tuple<Data7, empty_tuple> > > > > > >,
    result_is_temporary, 0,
    typename mp::enable_if<
        traits::is_implemented<
            sevenary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type,
                typename traits::basetype<Data4>::type,
                typename traits::basetype<Data5>::type,
                typename traits::basetype<Data6>::type,
                typename traits::basetype<Data7>::type> > >::type>
{
    typedef sevenary_expression<Op,
                typename traits::basetype<Data1>::type,
                typename traits::basetype<Data2>::type,
                typename traits::basetype<Data3>::type,
                typename traits::basetype<Data4>::type,
                typename traits::basetype<Data5>::type,
                typename traits::basetype<Data6>::type,
                typename traits::basetype<Data7>::type> wrapped_t;
    typedef typename wrapped_t::return_t return_t;
    typedef empty_tuple temporaries_t;
    typedef typename mp::make_tuple<Data1, Data2, Data3, Data4, Data5, Data6, Data7>::type data_t;
    template<class Return>
    static void doit(const data_t& input, temporaries_t temps, Return* output)
    {
        wrapped_t::doit(*output, mp::tuple_get<data_t, 0>::get(input),
                mp::tuple_get<data_t, 1>::get(input),
                mp::tuple_get<data_t, 2>::get(input),
                mp::tuple_get<data_t, 3>::get(input),
                mp::tuple_get<data_t, 4>::get(input),
                mp::tuple_get<data_t, 5>::get(input),
                mp::tuple_get<data_t, 6>::get(input));
    }
};

// Instantiating temporaries

namespace rdetail {
template<class T>
struct evaluated_type_pred
{
    template<class Expr>
    struct type : mp::equal_types<
              typename tools::evaluation_helper<Expr>::type, T> { };
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
    traits::is_expression<T>,
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
