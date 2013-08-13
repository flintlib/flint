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

// Common code shared among matrix classes

#ifndef FLINTXX_MATRIX_H
#define FLINTXX_MATRIX_H

#include "flint_classes.h"
#include "mp.h"
#include "rules.h"
#include "traits.h"
#include "tuple.h"

namespace flint {
FLINT_DEFINE_BINOP(mat_solve)
FLINT_DEFINE_THREEARY(mat_at)
FLINT_DEFINE_UNOP(transpose)
FLINT_DEFINE_UNOP(trace)
FLINT_DEFINE_UNOP(det)
FLINT_DEFINE_UNOP(charpoly)
FLINT_DEFINE_UNOP(nullspace)

template<class T> struct matrix_traits : rules::UNIMPLEMENTED { };
// override this for the non-reference type of your choice
// {
//     template<class M> static slong rows(const M&);
//     template<class M> static slong cols(const M&);
//     template<class M> static ??? at(const M&, slong, slong);
//     template<class M> static ??? at(M&, slong, slong);
// };
namespace traits {
template <class T, class Enable = void> struct is_mat : is_implemented<
    matrix_traits<typename flint_classes::to_nonref<T>::type> > { };
template<class T> struct is_mat<T,
    typename mp::disable_if<flint_classes::is_flint_class<T> >::type>
        : mp::false_ { };
} // traits
namespace matrices {
namespace mdetail {
template<class Data> struct second_is_mat_data : mp::false_ { };
template<class Data1, class Data2>
struct second_is_mat_data<tuple<Data1, tuple<Data2, empty_tuple> > >
    : traits::is_mat<typename traits::basetype<Data2>::type> { };
template<class Expr> struct second_is_mat
    : second_is_mat_data<typename Expr::data_t> { };

template<class Mat>
struct both_mat : mp::and_<
    traits::is_mat<
        typename traits::basetype<typename Mat::data_t::head_t>::type>,
    traits::is_mat<
        typename traits::basetype<typename Mat::data_t::tail_t::head_t>::type>
  > { };

template<class Expr> struct immediate_traits
    : matrix_traits<typename flint_classes::to_nonref<Expr>::type> { };
} // mdetail

template<class Operation>
struct outsize_generic
{
    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::disable_if<mdetail::second_is_mat<Mat> >::type* = 0)
    {
        return m._data().head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::disable_if<mdetail::second_is_mat<Mat> >::type* = 0)
    {
        return m._data().head.cols();
    }

    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::enable_if<mdetail::second_is_mat<Mat> >::type* = 0)
    {
        return m._data().tail.head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::enable_if<mdetail::second_is_mat<Mat> >::type* = 0)
    {
        return m._data().tail.head.cols();
    }
};

template<class Operation>
struct outsize : outsize_generic<Operation> { };

template<>
struct outsize<operations::immediate>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return mdetail::immediate_traits<Mat>::rows(m);
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return mdetail::immediate_traits<Mat>::cols(m);
    }
};

template<>
struct outsize<operations::times>
{
    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::enable_if<mdetail::both_mat<Mat> >::type* = 0)
    {
        return m._data().head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::enable_if<mdetail::both_mat<Mat> >::type* = 0)
    {
        return m._data().tail.head.cols();
    }

    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::disable_if<mdetail::both_mat<Mat> >::type* = 0)
    {
        return outsize_generic<operations::times>::rows(m);
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::disable_if<mdetail::both_mat<Mat> >::type* = 0)
    {
        return outsize_generic<operations::times>::cols(m);
    }
};

template<>
struct outsize<operations::transpose_op>
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

template<>
struct outsize<operations::nullspace_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head.cols();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head.cols();
    }
};

template<unsigned n>
struct outsize<operations::ltuple_get_op<n> >
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return outsize<
            typename Mat::data_t::head_t::operation_t>::rows(m._data().head);
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return outsize<
            typename Mat::data_t::head_t::operation_t>::cols(m._data().head);
    }
};

// This is not actually a matrix expression, but called by the above ...
template<>
struct outsize<operations::mat_solve_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().second().rows();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().second().cols();
    }
};

namespace mdetail {
struct base_traits
{
    template<class M>
    static slong rows(const M& m)
    {
        return matrices::outsize<typename M::operation_t>::rows(m);
    }
    template<class M>
    static slong cols(const M& m)
    {
        return matrices::outsize<typename M::operation_t>::cols(m);
    }
};
} // mdetail

template<class Mat>
struct generic_traits : mdetail::base_traits
{
    template<class T, class U>
    struct at
    {
        typedef FLINT_THREEARY_ENABLE_RETTYPE(mat_at, Mat, T, U)
            entry_ref_t;
        typedef entry_ref_t entry_srcref_t;

        static entry_srcref_t get(const Mat& m, T i, U j)
        {
            return mat_at(m, i, j);
        }
    };
};

template<class Mat, class Srcref>
struct generic_traits_srcref : generic_traits<Mat>
{
    template<class T, class U>
    struct at
    {
        typedef Srcref entry_ref_t;
        typedef Srcref entry_srcref_t;

        template<class M>
        static Srcref get(const M& m, T i, U j)
        {
            return Srcref::make(mdetail::immediate_traits<M>::at(m, i, j));
        }
    };
};

template<class Mat, class Ref, class Srcref>
struct generic_traits_ref : generic_traits_srcref<Mat, Srcref>
{
    template<class T, class U>
    struct at
        : generic_traits_srcref<Mat, Srcref>::template at<T, U>
    {
        typedef Ref entry_ref_t;

        template<class M>
        static Ref get(M& m, T i, U j)
        {
            return Ref::make(mdetail::immediate_traits<M>::at(m, i, j));
        }
    };
};
} // matrices
} // flint

// Define rows(), cols(), create_temporary() and at() methods.
// For this to work, Traits must behave like the above generic traits.
// Also your matrix class must have a static create_temporary_rowscols function.
// See fmpz_mat for an example.
#define FLINTXX_DEFINE_MATRIX_METHODS(Traits) \
template<class T, class U> \
typename Traits::template at<T, U>::entry_ref_t at(T i, U j) \
    {return Traits::template at<T, U>::get(*this, i, j);} \
template<class T, class U> \
typename Traits::template at<T, U>::entry_ref_t at(T i, U j) const \
    {return Traits::template at<T, U>::get(*this, i, j);} \
\
slong rows() const {return Traits::rows(*this);} \
slong cols() const {return Traits::cols(*this);} \
evaluated_t create_temporary() const \
{ \
    return create_temporary_rowscols(rows(), cols()); \
}

// Disable temporary merging. Requires create_temporary_rowscols.
#define FLINTXX_DEFINE_TEMPORARY_RULES(Matrix) \
template<class Expr> \
struct use_default_temporary_instantiation<Expr, Matrix> : mp::false_ { }; \
template<class Expr> \
struct instantiate_temporaries<Expr, Matrix, \
    typename mp::disable_if<flint_classes::is_Base<Matrix, Expr> >::type> \
{ \
    /* The only case where this should ever happen is if Expr is an ltuple */ \
    /* TODO static assert this */ \
    static Matrix get(const Expr& e) \
    { \
        typedef typename Expr::operation_t op_t; \
        return Matrix::create_temporary_rowscols( \
                matrices::outsize<op_t>::rows(e), \
                matrices::outsize<op_t>::cols(e)); \
    } \
}; \
template<class Expr> \
struct instantiate_temporaries<Expr, Matrix, \
    typename mp::enable_if<flint_classes::is_Base<Matrix, Expr> >::type> \
{ \
    static Matrix get(const Expr& e) \
    { \
        return e.create_temporary(); \
    } \
};

#endif
