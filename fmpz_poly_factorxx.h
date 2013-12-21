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

#ifndef FMPZ_POLY_FACTORXX_H
#define FMPZ_POLY_FACTORXX_H


#include "fmpz_poly.h"
#include "nmod_polyxx.h"

namespace flint {
class fmpz_poly_factorxx
{
private:
    fmpz_poly_factor_t inner;

public:
    fmpz_poly_factorxx() {fmpz_poly_factor_init(inner);}
    explicit fmpz_poly_factorxx(slong alloc)
        {fmpz_poly_factor_init2(inner, alloc);}
    ~fmpz_poly_factorxx() {fmpz_poly_factor_clear(inner);}

    fmpz_poly_factorxx(const fmpz_poly_factorxx& o)
    {
        fmpz_poly_factor_init2(inner, o.size());
        fmpz_poly_factor_set(inner, o.inner);
    }

    bool operator==(const fmpz_poly_factorxx& o)
    {
        if(o.content() != content() || o.size() != size())
            return false;
        for(ulong i = 0;i < size();++i)
            if(p(i) != o.p(i) || exp(i) != o.exp(i))
            return false;
        return true;
    }

    fmpz_poly_factorxx& operator=(const fmpz_poly_factorxx& o)
    {
        fmpz_poly_factor_set(inner, o.inner);
        return *this;
    }

    ulong size() const {return inner->num;}
    slong exp(slong i) const {return inner->exp[i];}
    slong& exp(slong i) {return inner->exp[i];}
    fmpz_polyxx_srcref p(slong i) const
        {return fmpz_polyxx_srcref::make(inner->p + i);}
    fmpz_polyxx_ref p(slong i) {return fmpz_polyxx_ref::make(inner->p + i);}
    fmpzxx_srcref content() const
        {return fmpzxx_srcref::make(& inner->c);}
    fmpzxx_ref content() {return fmpzxx_ref::make(& inner->c);}

    fmpz_poly_factor_t& _data() {return inner;}
    const fmpz_poly_factor_t& _data() const {return inner;}

    void realloc(slong a) {fmpz_poly_factor_realloc(inner, a);}
    void fit_length(slong a) {fmpz_poly_factor_fit_length(inner, a);}

    void print() const {fmpz_poly_factor_print(inner);}

    template<class Fmpz_poly>
    void insert(const Fmpz_poly& p, slong e,
            typename mp::enable_if<traits::is_fmpz_polyxx<Fmpz_poly> >::type* = 0)
        {fmpz_poly_factor_insert(_data(), p.evaluate()._poly(), e);}

    void concat(const fmpz_poly_factorxx& o)
        {fmpz_poly_factor_concat(_data(), o._data());}

    template<class Fmpz_poly>
    void set_factor_squarefree(const Fmpz_poly& p,
            typename mp::enable_if<traits::is_fmpz_polyxx<Fmpz_poly> >::type* = 0)
        {fmpz_poly_factor_squarefree(_data(), p.evaluate()._poly());}

    template<class Fmpz_poly, class Fmpz>
    void set_factor_zassenhaus_recombination(
            const fmpz_poly_factorxx& lifted_fac, const Fmpz_poly& F,
            const Fmpz& P, slong exp,
            typename mp::enable_if<traits::is_fmpz_polyxx<Fmpz_poly> >::type* = 0,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        fmpz_poly_factor_zassenhaus_recombination(_data(), lifted_fac._data(),
                F.evaluate()._poly(), P.evaluate()._fmpz(), exp);
    }

    template<class Fmpz_poly>
    void set_factor_zassenhaus(const Fmpz_poly& p,
            typename mp::enable_if<traits::is_fmpz_polyxx<Fmpz_poly> >::type* = 0)
        {fmpz_poly_factor_zassenhaus(_data(), p.evaluate()._poly());}

    template<class Fmpz_poly>
    void set_hensel_lift_once(const Fmpz_poly& f,
            const nmod_poly_factorxx& local_fac, slong N,
            typename mp::enable_if<traits::is_fmpz_polyxx<Fmpz_poly> >::type* = 0)
    {
        fmpz_poly_hensel_lift_once(_data(), f.evaluate()._poly(),
                local_fac._data(), N);
    }
};

template<class Fmpz_poly>
inline fmpz_poly_factorxx factor_squarefree(const Fmpz_poly& p,
        typename mp::enable_if<traits::is_fmpz_polyxx<Fmpz_poly> >* = 0)
{
    fmpz_poly_factorxx res;
    res.set_factor_squarefree(p);
    return res;
}
template<class Fmpz_poly>
inline fmpz_poly_factorxx factor_zassenhaus(const Fmpz_poly& p,
        typename mp::enable_if<traits::is_fmpz_polyxx<Fmpz_poly> >* = 0)
{
    fmpz_poly_factorxx res;
    res.set_factor_zassenhaus(p);
    return res;
}

template<class Fmpz_poly>
inline fmpz_poly_factorxx hensel_lift_once(const Fmpz_poly& f,
        const nmod_poly_factorxx& local_fac, slong N,
        typename mp::enable_if<traits::is_fmpz_polyxx<Fmpz_poly> >* = 0)
{
    fmpz_poly_factorxx res;
    res.set_hensel_lift_once(f, local_fac, N);
    return res;
}

inline void print(const fmpz_poly_factorxx& f)
{
    f.print();
}
} // flint

#endif
