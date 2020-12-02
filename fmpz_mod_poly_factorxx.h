/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_POLY_FACTORXX_H
#define FMPZ_MOD_POLY_FACTORXX_H


#include "fmpz_mod_poly.h"

namespace flint {
class fmpz_mod_poly_factorxx
{
private:
    fmpz_mod_poly_factor_t inner;
    fmpz_modxx_ctx_srcref ctx;

public:

    fmpz_mod_poly_factorxx(fmpz_modxx_ctx_srcref c)
        : ctx(c)
    {
        fmpz_mod_poly_factor_init(inner, ctx._ctx());
    }

    ~fmpz_mod_poly_factorxx()
    {
        fmpz_mod_poly_factor_clear(inner, ctx._ctx());
    }

    fmpz_mod_poly_factorxx(const fmpz_mod_poly_factorxx& o)
        : ctx(o.ctx)
    {
        fmpz_mod_poly_factor_init(inner, ctx._ctx());
        fmpz_mod_poly_factor_set(inner, o.inner, ctx._ctx());
    }

    bool operator==(const fmpz_mod_poly_factorxx& o)
    {
        if(o.size() != size())
            return false;
        for(slong i = 0;i < size();++i)
            if(p(i) != o.p(i) || exp(i) != o.exp(i))
            return false;
        return true;
    }

    fmpz_mod_poly_factorxx& operator=(const fmpz_mod_poly_factorxx& o)
    {
        fmpz_mod_poly_factor_set(inner, o.inner, ctx._ctx());
        return *this;
    }

    slong size() const {return inner->num;}
    slong exp(slong i) const {return inner->exp[i];}
    slong& exp(slong i) {return inner->exp[i];}
    fmpz_mod_polyxx p(slong i) const
    {
        fmpz_mod_polyxx p(ctx);
        fmpz_mod_poly_set(p._poly(), inner->poly + i, ctx._ctx());
        return p;
    }

    fmpz_mod_poly_factor_t& _data() {return inner;}
    const fmpz_mod_poly_factor_t& _data() const {return inner;}

    void realloc(slong a) {fmpz_mod_poly_factor_realloc(inner, a, ctx._ctx());}
    void fit_length(slong a) {fmpz_mod_poly_factor_fit_length(inner, a, ctx._ctx());}

    void print() const {fmpz_mod_poly_factor_print(inner, ctx._ctx());}

    template<class Fmpz_mod_poly>
    void insert(const Fmpz_mod_poly& p, slong e,
            typename mp::enable_if<traits::is_fmpz_mod_polyxx<Fmpz_mod_poly> >::type* = 0)
        {fmpz_mod_poly_factor_insert(_data(), p.evaluate()._poly(), e, ctx._ctx());}

    void concat(const fmpz_mod_poly_factorxx& o)
        {fmpz_mod_poly_factor_concat(_data(), o._data(), ctx._ctx());}

    void pow(slong exp) {fmpz_mod_poly_factor_pow(_data(), exp, ctx._ctx());}

#define FMPZ_MOD_POLY_FACTORXX_DEFINE_SET_FACTOR(name) \
    template<class Fmpz_mod_poly> \
    void set_##name(const Fmpz_mod_poly& p, \
            typename mp::enable_if<traits::is_fmpz_mod_polyxx<Fmpz_mod_poly> >::type* = 0) \
        {fmpz_mod_poly_##name(_data(), p.evaluate()._poly(), ctx._ctx());}

    FMPZ_MOD_POLY_FACTORXX_DEFINE_SET_FACTOR(factor)
    FMPZ_MOD_POLY_FACTORXX_DEFINE_SET_FACTOR(factor_squarefree)
    FMPZ_MOD_POLY_FACTORXX_DEFINE_SET_FACTOR(factor_cantor_zassenhaus)
    FMPZ_MOD_POLY_FACTORXX_DEFINE_SET_FACTOR(factor_berlekamp)
    FMPZ_MOD_POLY_FACTORXX_DEFINE_SET_FACTOR(factor_kaltofen_shoup)

    template<class Fmpz_mod_poly>
    bool set_factor_equal_deg_probab(frandxx& state, const Fmpz_mod_poly& p, slong d,
            typename mp::enable_if<traits::is_fmpz_mod_polyxx<Fmpz_mod_poly> >::type* = 0)
    {
        return fmpz_mod_poly_factor_equal_deg_prob(_data(), state._data(),
                p.evaluate()._poly(), d);
    }
    template<class Fmpz_mod_poly>
    void set_factor_equal_deg(const Fmpz_mod_poly& p, slong d,
            typename mp::enable_if<traits::is_fmpz_mod_polyxx<Fmpz_mod_poly> >::type* = 0)
    {
        fmpz_mod_poly_factor_equal_deg(_data(), p.evaluate()._poly(), d, ctx._ctx());
    }

    template<class Fmpz_mod_poly>
    void set_factor_distinct_deg(const Fmpz_mod_poly& p, std::vector<slong>& degs,
            typename mp::enable_if<traits::is_fmpz_mod_polyxx<Fmpz_mod_poly> >::type* = 0)
    {
        slong* dgs = &degs.front();
        fmpz_mod_poly_factor_distinct_deg(_data(), p.evaluate()._poly(), &dgs, ctx._ctx());
    }
};

#define FMPZ_MOD_POLY_FACTORXX_DEFINE_FACTOR(name) \
template<class Fmpz_mod_poly> \
fmpz_mod_poly_factorxx name(const Fmpz_mod_poly& p, \
        typename mp::enable_if<traits::is_fmpz_mod_polyxx<Fmpz_mod_poly> >::type* = 0) \
{ \
    fmpz_mod_poly_factorxx res(p.evaluate().get_ctx()); \
    res.set_##name(p); \
    return res; \
}
FMPZ_MOD_POLY_FACTORXX_DEFINE_FACTOR(factor)
FMPZ_MOD_POLY_FACTORXX_DEFINE_FACTOR(factor_squarefree)
FMPZ_MOD_POLY_FACTORXX_DEFINE_FACTOR(factor_cantor_zassenhaus)
FMPZ_MOD_POLY_FACTORXX_DEFINE_FACTOR(factor_berlekamp)
FMPZ_MOD_POLY_FACTORXX_DEFINE_FACTOR(factor_kaltofen_shoup)

// TODO do we want global versions of factor_distinct_deg etc?

inline void print(const fmpz_mod_poly_factorxx& f)
{
    f.print();
}
} // flint

#endif
