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

#ifndef FMPZ_MATXX_H
#define FMPZ_MATXX_H FMPZ_MATXX_H

#include <algorithm> // max

#include "fmpz_mat.h"

#include "fmpzxx.h"

namespace flint {
// TODO move to stdmath
FLINT_DEFINE_UNOP(transpose)
FLINT_DEFINE_UNOP(trace)

namespace detail {
template<class Operation>
struct fmpz_matxx_outsize;
template<class Mat>
struct fmpz_matxx_traits
{
    static slong rows(const Mat& m)
    {
        return fmpz_matxx_outsize<typename Mat::operation_t>::rows(m);
    }
    static slong cols(const Mat& m)
    {
        return fmpz_matxx_outsize<typename Mat::operation_t>::cols(m);
    }
};
} // detail

template<class Operation, class Data>
class fmpz_matxx_expression
    : public expression<derived_wrapper<fmpz_matxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpz_matxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpz_matxx_traits<fmpz_matxx_expression> traits_t;

    FLINTXX_DEFINE_BASICS(fmpz_matxx_expression)
    FLINTXX_DEFINE_CTORS(fmpz_matxx_expression)
    FLINTXX_DEFINE_C_REF(fmpz_matxx_expression, fmpz_mat_struct, _mat)

public:
    // TODO
    fmpzxx_ref at(slong i, slong j)
        {return fmpzxx_ref::make(fmpz_mat_entry(_mat(), i, j));}

    slong rows() const {return traits_t::rows(*this);}
    slong cols() const {return traits_t::cols(*this);}

    evaluated_t create_temporary() const
    {
        return evaluated_t(rows(), cols());
    }
};

namespace detail {
struct fmpz_mat_data;
} // detail

typedef fmpz_matxx_expression<operations::immediate, detail::fmpz_mat_data> fmpz_matxx;
typedef fmpz_matxx_expression<operations::immediate,
            flint_classes::ref_data<fmpz_matxx, fmpz_mat_struct> > fmpz_matxx_ref;
typedef fmpz_matxx_expression<operations::immediate, flint_classes::srcref_data<
    fmpz_matxx, fmpz_matxx_ref, fmpz_mat_struct> > fmpz_matxx_srcref;

namespace detail {
struct fmpz_mat_data
{
    typedef fmpz_mat_t& data_ref_t;
    typedef const fmpz_mat_t& data_srcref_t;

    fmpz_mat_t inner;

    fmpz_mat_data(slong m, slong n)
    {
        fmpz_mat_init(inner, m, n);
    }

    fmpz_mat_data(const fmpz_mat_data& o)
    {
        fmpz_mat_init_set(inner, o.inner);
    }

    ~fmpz_mat_data() {fmpz_mat_clear(inner);}
};
} // detail

#define FMPZ_MATXX_COND_S FLINTXX_COND_S(fmpz_matxx)
#define FMPZ_MATXX_COND_T FLINTXX_COND_T(fmpz_matxx)

namespace traits {
template<> struct use_temporary_merging<fmpz_matxx> : mp::false_ { };

template<class T> struct is_fmpz_matxx : mp::or_<
     traits::is_T_expr<T, fmpz_matxx>,
     flint_classes::is_source<fmpz_matxx, T> > { };
} // traits
namespace mp {
template<class T1, class T2 = void, class T3 = void, class T4 = void>
struct all_fmpz_matxx : mp::and_<all_fmpz_matxx<T1>, all_fmpz_matxx<T2, T3, T4> > { };
template<class T>
struct all_fmpz_matxx<T, void, void, void> : traits::is_fmpz_matxx<T> { };

template<class Out, class T1, class T2 = void, class T3 = void, class T4 = void>
struct enable_all_fmpz_matxx
    : mp::enable_if<all_fmpz_matxx<T1, T2, T3, T4>, Out> { };
} // mp

namespace detail {
template<class Operation>
struct fmpz_matxx_outsize_generic
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head.cols();
    }

    template<class Data1, class Data2>
    static slong rows(const fmpz_matxx_expression<
            Operation, tuple<Data1, tuple<Data2, empty_tuple> > >& m,
            typename mp::enable_if<traits::is_fmpz_matxx<Data2> >::type* = 0)
    {
        return m._data().tail.head.rows();
    }
    template<class Data1, class Data2>
    static slong cols(const fmpz_matxx_expression<
            Operation, tuple<Data1, tuple<Data2, empty_tuple> > >& m,
            typename mp::enable_if<traits::is_fmpz_matxx<Data2> >::type* = 0)
    {
        return m._data().tail.head.cols();
    }
};
template<class Operation>
struct fmpz_matxx_outsize : fmpz_matxx_outsize_generic<Operation> { };
template<>
struct fmpz_matxx_outsize<operations::immediate>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return fmpz_mat_nrows(m._mat());
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return fmpz_mat_ncols(m._mat());
    }
};

template<class Mat>
struct both_mat : mp::all_fmpz_matxx<
    typename traits::basetype<typename Mat::data_t::head_t>::type,
    typename traits::basetype<typename Mat::data_t::tail_t::head_t>::type> { };
template<>
struct fmpz_matxx_outsize<operations::times>
{
    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::enable_if<both_mat<Mat> >::type* = 0)
    {
        return m._data().head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::enable_if<both_mat<Mat> >::type* = 0)
    {
        return m._data().tail.head.cols();
    }

    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::disable_if<both_mat<Mat> >::type* = 0)
    {
        return fmpz_matxx_outsize_generic<operations::times>::rows(m);
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::disable_if<both_mat<Mat> >::type* = 0)
    {
        return fmpz_matxx_outsize_generic<operations::times>::cols(m);
    }
};

template<>
struct fmpz_matxx_outsize<operations::transpose_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head.cols();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head.rows();
    }
};
}

namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpz_mat_mul(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZXX_COND_S,
        fmpz_mat_scalar_mul_fmpz(to._mat(), e1._mat(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpz_mat_add(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(transpose_op, fmpz_matxx, FMPZ_MATXX_COND_S,
        fmpz_mat_transpose(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(trace_op, fmpzxx, FMPZ_MATXX_COND_S,
        fmpz_mat_trace(to._fmpz(), from._mat()))
} // rules
} // flint

#endif
