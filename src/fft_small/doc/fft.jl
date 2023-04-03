# This file contains an illustration of David Harvey's truncated fft fitted
# with reduced twiddle factor loads. It also contains ifft_trunc_formula which
# prints formulas that are needed in the inverse truncated transform.

using Nemo

L = 2
Qw,w = PolynomialRing(QQ, "x")
K,w = NumberField(w^(2^(L-1))+1, "w")
R,X = PolynomialRing(K, vcat([Symbol(:x,i) for i in 0:2^L-1]))

function revbits(a::Int, l::Int)
    @assert 0 <= a < 2^l
    return evalpoly(2, reverse(digits(a, base=2, pad=l)))
end

# on input:
#  x is assumed to have length at least 2^k
#  x[1+itrunc] ... x[1+2^k-1] are not read and assumed to be zero
# on output (checked by an assert):
#  x[1+otrunc] ... x[1+2^k-1] are not defined
#  x[1+i] is evalpoly(w^revbits(j*2^k+i,L), x[1:itrunc])
# all accesses to x are offset in the real x by I and strided by S
function tfft!(x::Vector, I::Int, S::Int, k::Int, j::Int, itrunc::Int, otrunc::Int)
    @assert k >= 0
    @assert 0 <= itrunc <= 2^k
    @assert 0 <= otrunc <= 2^k
    if otrunc < 1
        return
    elseif itrunc < 1
        for i in 0:otrunc-1
            x[1+I+S*i] = zero(R)
        end
        return
    end

    # for answer check
    inx = elem_type(R)[x[1+I+S*i] for i in 0:itrunc-1]

    if k < 1
        return
    elseif k == 1
        ww = w^revbits(2*j+0,L)
        x0 = 0 < itrunc ? x[1+I+0*S] : zero(R)
        x1 = 1 < itrunc ? x[1+I+1*S] : zero(R)
        x[1+I+0*S] = 0 < otrunc ? x0 + ww*x1 : zero(R)
        x[1+I+1*S] = 1 < otrunc ? x0 - ww*x1 : zero(R)
    else
        # any k1 + k2 = k with k1, k2 >= 1 will do
        k1 = fld(k, 2)
        k2 = k - k1

        n1, n2 = divrem(otrunc, 2^k2)
        z1, z2 = divrem(itrunc, 2^k2)
        n1p = n1 + (n2 != 0)
        z2p = min(2^k2, itrunc)

        # columns   !!! j is not changing in this loop !!!
        for a in 0:z2p-1
            tfft!(x, I + a*S, S*2^k2, k1, j, z1 + (a < z2), n1p)
        end

        # rows
        for b in 0:n1-1
            tfft!(x, I + S*b*2^k2, S, k2, j*2^k1 + b, z2p, 2^k2)
        end
        if n2 > 0
            tfft!(x, I + S*n1*2^k2, S, k2, j*2^k1 + n1, z2p, n2)
        end
    end

    # check answer
    for i in 0:otrunc-1
        if x[1+I+S*i] != evalpoly(R(w^revbits(j*2^k+i,L)), inx)
            error("error at i = $i with k = $k")
        end
    end
end

#=
L = 2^k
suppose a[z] = ... = a[L-1] = 0
n <= z
1 <= n + f <= L
IFFT(L, zeta, z, n, f; ah[0], ..., ah[n-1], L*a[n], ..., L*a[z-1])
returns with (L*a[0], ..., L*a[n-1])        if f = 0
             (L*a[0], ..., L*a[n-1], ah[n]) if f = 1

    ifft_trunc_formula generates the input -> output map as a n+f by z matrix.
    Special case for L = 2:

        Harvey                                          !!! Here !!!
        ah[0] =       a[0] + a[1]                      ah[0] = a[0] + w*a[1]
        ah[1] = zeta*(a[0] - a[1])                     ah[1] = a[0] - w*a[1]

    if n = 2,
        2*a[0] = ah[0] + zeta^-1*ah[1]                2*a[0] =  ah[0] + ah[1]
        2*a[1] = ah[0] - zeta^-1*ah[1]                2*a[1] = (ah[0] - ah[1])*w^-1

    if n = 1, f = 1, z = 2
        ah[0], 2*a[1]
        2*a[0] = 2*ah[0] - 2*a[1]                     2*a[0] = 2*ah[0] - w*2*a[1]
         ah[1] = zeta*(ah[0] - 2*a[1])                 ah[1] =   ah[0] - w*2*a[1]

    if n = 1, f = 1, z = 1
        2*a[0] = 2*ah[0]
         ah[1] = zeta*(ah[0])

    if n = 0, f = 1, z = 2
        ah[0] = (2*a[0] + 2*a[1])/2

    if n = 0, f = 1, z = 1
        ah[0] = (2*a[0])/2
=#
function itfft!(x::Vector, I::Int, S::Int, k::Int, j::Int, z::Int, n::Int, f::Bool)
    @assert n <= z
    @assert 1 <= z <= 2^k
    @assert 1 <= n + f <= 2^k
    if k < 1
        return
    elseif k == 1
        ww = w^revbits(2*j+0,L)
        if n == 2
            u = x[1+I+0*S]
            v = x[1+I+1*S]
            x[1+I+0*S] =  u + v
            x[1+I+1*S] = (u - v)*ww^-1
        elseif n == 1
            u = x[1+I+0*S]
            v = z == 2 ? x[1+I+1*S] : zero(R)
            x[1+I+0*S] = 2*u - ww*v
            if f
                x[1+I+1*S] = u - ww*v
            end
        elseif n == 0
            u = x[1+I+0*S]
            v = z == 2 ? x[1+I+1*S] : zero(R)
            x[1+I+0*S] = (u + ww*v)*inv(K(2))
        else
            error("case n=$n, f=$f, z=$z not implemented")
        end
    else
        # any k1 + k2 = k with k1, k2 >= 1 will do
        k1 = fld(k, 2)
        k2 = k - k1

        n1, n2 = divrem(n, 2^k2)
        z1, z2 = divrem(z, 2^k2)
        fp = n2 + f > 0
        z2p = min(2^k2, z)
        m = min(n2, z2); mp = max(n2, z2)

        # complete rows
        for b in 0:n1-1
            itfft!(x, I + S*b*2^k2, S, k2, j*2^k1 + b, 2^k2, 2^k2, false)
        end

        # rightmost columns !!! j is not changing in this loop !!!
        for a in n2:z2p-1
            itfft!(x, I + S*a, S*2^k2, k1, j, z1 + (a < mp), n1, fp)
        end

        # last partial row
        if fp
            itfft!(x, I + S*n1*2^k2, S, k2, j*2^k1 + n1, z2p, n2, f)
        end

        # leftmost columns  !!! j is not changing in this loop !!!
        for a in 0:n2-1
            itfft!(x, I + S*a, S*2^k2, k1, j, z1 + (a < m), n1 + 1, false)
        end
    end
end

println("-------- fft ----------")
for itrunc in 1:2^L
    @show itrunc
    for otrunc in 1:2^L
        x = copy(X)
        tfft!(x, 0, 1, L, 0, itrunc, otrunc)
        for i in 0:otrunc-1
            @assert x[1+i] == evalpoly(R(w^revbits(i,L)), X[1:itrunc])
        end
    end
end

println("-------- ifft ---------")
for trunc in 1:2^L
    @show trunc
    x = copy(X)
    tfft!(x, 0, 1, L, 0, trunc, trunc)
    itfft!(x, 0, 1, L, 0, trunc, trunc, false)
#show(stdout, "text/plain", x); println();
    for i in 0:trunc-1
        @assert x[1+i] == 2^L*X[1+i]
    end
end

function ifft_trunc_formula(k::Int, n::Int, z::Int, f::Bool)
    @assert n <= z
    @assert 1 <= z <= 2^k
    @assert 1 <= n + f <= 2^k
    l = 2^k

    local Qr,r = PolynomialRing(QQ, "r")
    local K,r = NumberField(r^(2^(k-1)) + 1, "r")
    local Rwr,w = PolynomialRing(K, "w")
    r = Rwr(r)
    local F = FractionField(Rwr)
    r = F(r)
    w = F(w)

    M = zero_matrix(F, l, l)
    for i in 0:n-1, j in 0:l-1
        M[1+i,1+j] = (r^revbits(i,k)*w)^j
    end
    for i in n:l-1
        M[1+i,1+i] = l
    end

    N = zero_matrix(F, n+f, l)
    for i in 0:n-1
        N[1+i,1+i] = l
    end
    if f
        for j in 0:l-1
            N[1+n,1+j] = (r^revbits(n,k)*w)^j
        end
    end

    println("\nk = $k, n = $n, z = $z, f = $f")
    show(stdout, "text/plain", (N*inv(M))[:, 1:z])
    println()
    nothing
end

println("\n-------- radix 2 --------------------")
ifft_trunc_formula(1, 2,2,false)
ifft_trunc_formula(1, 1,2,true)
ifft_trunc_formula(1, 1,2,false)
ifft_trunc_formula(1, 1,1,true)
ifft_trunc_formula(1, 1,1,false)

println("\n-------- radix 4 (r^2 = -1) ---------")
ifft_trunc_formula(2, 4,4,false)
ifft_trunc_formula(2, 1,1,false)
ifft_trunc_formula(2, 1,1,true)
ifft_trunc_formula(2, 2,4,false)
ifft_trunc_formula(2, 2,2,false)
ifft_trunc_formula(2, 3,3,false)
ifft_trunc_formula(2, 3,3,true)
ifft_trunc_formula(2, 2,2,true)
ifft_trunc_formula(2, 1,4,false)
ifft_trunc_formula(2, 3,4,false)
ifft_trunc_formula(2, 0,4,true)
ifft_trunc_formula(2, 2,4,true)
ifft_trunc_formula(2, 3,4,true)
ifft_trunc_formula(2, 1,4,true)

nothing
