#!/usr/bin/env python3
"""
Generates fixed-limb-size dot product kernels and classical polynomial
multiplication loops for FLINT's mpn_extras module.

    python3 dev/gen_mpn_dot_rev.py kernels     > src/mpn_extras/dot_rev_kernels.c
    python3 dev/gen_mpn_dot_rev.py loops       > src/mpn_extras/poly_mulmid_kernels.c
    python3 dev/gen_mpn_dot_rev.py prototypes  (table declarations for mpn_extras.h)

Dot product kernels
-------------------

    _mpn_dot_{n1}x{n2}_{s}(res, a, b, len)                 forward: a[k] b[k]
    _mpn_dot_rev_{n1}x{n2}_{s}(res, a, b, len)             reversed: a[k] b[len-1-k]
    _mpn_dot_strided_{n1}x{n2}_{s}(res, a, astride, b, bstride, len)

(and _signed versions for two's complement inputs) compute
{res, s} = sum_{k=0}^{len-1} a_k * b_k mod 2^(64 s), where a and b are
arrays of n1-limb resp. n2-limb integers (contiguous, or with the given
limb strides, which may be negative), for n1 >= n2,
s in {n1+n2-1, n1+n2, n1+n2+1}. Dedicated forward and reversed kernels
are generated for n1 <= 2; for larger sizes the strided kernels are used
through inline wrappers, the stride overhead being negligible there.

Each partial product A_p * B_q (limb offset o = p + q) is accumulated in a
3-limb carry chain associated with the offset o (the third limb absorbs
carries); the two topmost offsets share the last chain since carries out
of limb s - 1 are discarded anyway. Using one short chain per offset keeps
the number of live registers small (spills are the main bottleneck on
x86_64) while splitting the carry chains for instruction-level parallelism.
The products are emitted interleaved between chains so that dependent
additions are not adjacent. (Alternatives that were tried and found to be
slower: one chain per diagonal p - q, which needs fewer additions but more
registers; and splitting the low limb of a chain between products.)

For signed (two's complement) inputs, the top x top product is computed
with smul_ppmm and sign-extended into the chain pads, while the mixed
top x low products are corrected by subtracting (sign mask & low limb) of
the other operand at the appropriate offset, directly into the main
chains (such chains may become negative and are sign-extended in the
final combination), or, when s = n1 + n2 - 1, via a separate accumulator.

Classical multiplication loops
------------------------------

    _mpn_poly_mulmid_{n1}x{n2}_{s}[_signed](res, f, flen, g, glen, nlo, nhi)
    _mpn_poly_sqrmid_{n}x{n}_{s}[_signed](res, f, flen, nlo, nhi)

compute coefficients [nlo, nhi) of the product f * g (or f^2) with the
dot product body inlined in the coefficient loop, so that no function
calls are made and, for squaring, the doubling of the cross terms and
the addition of the square term happen in registers. These are
generated for n1 <= 3 only; larger sizes go through the kernels or
through flint_mpn_mul. For n1 <= 2, each term's full product is computed
with a minimal number of additions and added to a single accumulator,
which was measured to be faster than the chain-per-offset scheme.

    _mpn_poly_mulmid_{n1}x{n2}_{s}_signmag, _mpn_poly_sqrmid_{n}x{n}_{s}_signmag

are versions for sign-magnitude coefficients, stored as a sign limb
(0 or 1) followed by the magnitude, so that no sign bit is needed in the
magnitude; each (unsigned) product is negated conditionally on the
product sign before being accumulated. The output is two's complement.
"""

import re, sys

MAXN = 4        # strided kernels are generated for n1, n2 <= MAXN
ASMKN = 4       # kernels for n1 >= ASMKN use assembly products (see gen_kernel)
DEDN = 2        # dedicated forward and reversed kernels for n1 <= DEDN
MAXS = 8        # longest carry chains supported by longlong.h
LOOPN = 3       # inlined loops are generated for n1 <= LOOPN
ACCN = 2        # loops for n1 <= ACCN accumulate full products (see use_acc)
ASMN = (4,)     # "_short" loops (assembly products, register accumulation) for these n1
                # (for 3 limbs, the column loops with memory operands win at all lengths)
SIGNMAGN = 4    # sign-magnitude loops for n1 <= SIGNMAGN


def addmacro(L):
    return "add_" + "s" * L + "a" * (2 * L)


def submacro(L):
    return "sub_" + "d" * L + "m" * L + "s" * L


def emit_add(out, dst, src, indent):
    L = len(dst)
    if L == 1:
        out.append(f"{indent}{dst[0]} += {src[0]};")
    else:
        d = ", ".join(reversed(dst))
        s = ", ".join(reversed(src))
        out.append(f"{indent}{addmacro(L)}({d}, {d}, {s});")


def emit_sub(out, dst, src, indent):
    L = len(dst)
    if L == 1:
        out.append(f"{indent}{dst[0]} -= {src[0]};")
    else:
        d = ", ".join(reversed(dst))
        s = ", ".join(reversed(src))
        out.append(f"{indent}{submacro(L)}({d}, {d}, {s});")



def prune_unused(text):
    """Drop declared-but-unused variables from a generated function,
    together with pure assignments (loads, mask computations) whose
    results are never read: small kernel shapes use only a subset of the
    full register set, and the assembly-product kernels read the
    operands from memory so that the register loads only feed the sign
    corrections. Iterates so that removing an assignment can free its
    inputs in turn."""
    lines = text.split("\n")
    while True:
        decls = []
        for l in lines:
            m = re.match(r"\s*(?:ulong|slong) ([A-Za-z0-9_, ]+);$", l)
            if m:
                decls += m.group(1).split(", ")
        # pure single-assignment statements per variable
        assigns = {}
        for i, l in enumerate(lines):
            m = re.match(r"\s*([A-Za-z0-9_]+) = ([^=;]|[!=<>]=)*;$", l)
            if m and m.group(1) in decls:
                assigns.setdefault(m.group(1), []).append(i)
        # a umul_ppmm whose high output is never read becomes a plain
        # low multiply, giving good codegen regardless of how umul_ppmm
        # itself is implemented
        for i, l in enumerate(lines):
            m = re.match(r"(\s*)umul_ppmm\(([A-Za-z0-9_]+), ([A-Za-z0-9_]+), ([^,]+), ([^;)]+)\);(.*)$", l)
            if m and m.group(2) in decls:
                v1 = m.group(2)
                uses = 0
                for j, l2 in enumerate(lines):
                    if j == i or re.match(r"\s*(?:ulong|slong) [A-Za-z]", l2):
                        continue
                    uses += len(re.findall(r"\b%s\b" % v1, l2))
                if uses == 0:
                    lines[i] = f"{m.group(1)}{m.group(3)} = ({m.group(4)}) * ({m.group(5)});{m.group(6)}"
        drop_lines = set()
        drop_vars = set()
        for v in decls:
            uses = 0
            for i, l in enumerate(lines):
                if re.match(r"\s*(?:ulong|slong) [A-Za-z]", l):
                    continue
                n = len(re.findall(r"\b%s\b" % v, l))
                if i in assigns.get(v, []):
                    n -= 1
                uses += n
            if uses == 0:
                drop_vars.add(v)
                drop_lines.update(assigns.get(v, []))
        if not drop_vars:
            break
        out = []
        for i, l in enumerate(lines):
            if i in drop_lines:
                continue
            m = re.match(r"(\s*)(ulong|slong) ([A-Za-z0-9_, ]+);$", l)
            if m:
                keep = [v for v in m.group(3).split(", ") if v not in drop_vars]
                if not keep:
                    continue
                l = f"{m.group(1)}{m.group(2)} {', '.join(keep)};"
            out.append(l)
        lines = out
    return "\n".join(lines)

class Kernel:
    """Code fragments for one (n1, n2, s, signed) dot product kernel."""

    def __init__(self, n1, n2, s, signed):
        assert n1 >= n2 >= 1 and n1 + n2 - 1 <= s <= n1 + n2 + 1 and s <= MAXS
        self.n1, self.n2, self.s, self.signed = n1, n2, s, signed

        # one 3-limb chain per offset o <= top; the offsets above top share
        # the last chain
        top = max(s - 3, 0)
        self.chains = []
        for o in range(0, top + 1):
            prods = [(p, q, p + q) for p in range(n1) for q in range(n2)
                     if (p + q == o if o < top else p + q >= top)]
            prods.sort(key=lambda t: (t[2], t[0]))
            if prods:
                self.chains.append(dict(lo=o, L=min(o + 3, s) - o, prods=prods, negative=False))

        # sign corrections: subtract (mask_a & B_q) at offset n1 + q for
        # q < n2 - 1, and (mask_b & A_p) at offset n2 + p for p < n1 - 1
        self.ca_m = max(0, min(n2 - 1, s - n1)) if signed else 0
        self.cb_m = max(0, min(n1 - 1, s - n2)) if signed else 0
        self.inline_sign = signed and s >= n1 + n2
        # operands referenced in memory (folded into the multiply
        # instructions) instead of being loaded into registers, which
        # relieves register pressure in the wider unsigned column kernels
        self.memops = False
        if self.inline_sign:
            self.ca_L = self.cb_L = 0
        else:
            self.ca_L = min(self.ca_m + 1, s - n1) if self.ca_m > 0 else 0
            self.cb_L = min(self.cb_m + 1, s - n2) if self.cb_m > 0 else 0

        # which chains receive inline subtractions (they may go negative)
        if self.inline_sign:
            for q in range(self.ca_m):
                self.chains[self.chain_for_offset(n1 + q)]["negative"] = True
            for p in range(self.cb_m):
                self.chains[self.chain_for_offset(n2 + p)]["negative"] = True

    def chain_for_offset(self, o):
        for ci, c in enumerate(self.chains):
            if c["lo"] == o:
                return ci
        return len(self.chains) - 1

    def declarations(self, asm=False, loop=False, signmag=False):
        n1, n2, s = self.n1, self.n2, self.s
        decl = []
        decl.append(", ".join(f"A{p}" for p in range(n1)) + ", " + ", ".join(f"B{q}" for q in range(n2)))
        decl.append("p0, p1")
        if loop and (self.use_acc() or signmag):
            decl.append(", ".join(f"P{j}" for j in range(min(n1 + n2, s))) + ", t0, t1")
        for ci, c in enumerate(self.chains):
            decl.append(", ".join(f"c{ci}_{j}" for j in range(c["L"])))
        if self.ca_L > 0:
            decl.append(", ".join(f"ca{j}" for j in range(self.ca_L)))
        if self.cb_L > 0:
            decl.append(", ".join(f"cb{j}" for j in range(self.cb_L)))
        if self.ca_m > 0 or (asm and self.signed):
            decl.append("ma")
        if self.cb_m > 0 or (asm and self.signed):
            decl.append("mb")
        if self.signed:
            decl.append("sx")
        if signmag:
            decl.append("ms, cnt")
        return [f"ulong {d};" for d in decl]

    def init(self, indent):
        out = []
        for ci, c in enumerate(self.chains):
            out.append(indent + " = ".join(f"c{ci}_{j}" for j in range(c["L"])) + " = 0;")
        if self.ca_L > 0:
            out.append(indent + " = ".join(f"ca{j}" for j in range(self.ca_L)) + " = 0;")
        if self.cb_L > 0:
            out.append(indent + " = ".join(f"cb{j}" for j in range(self.cb_L)) + " = 0;")
        return out

    def loads(self, indent, aptr="a", bptr="b", asm=False, signmag=False):
        """Load A and B from the given pointers (bptr may equal aptr). With
        signmag, each coefficient is a sign limb followed by the magnitude,
        and the product sign mask is loaded into ms."""
        out = []
        off = 1 if signmag else 0
        if signmag:
            if bptr == aptr:
                out.append(f"{indent}ms = 0;")
            else:
                out.append(f"{indent}ms = -({aptr}[0] ^ {bptr}[0]);")
        for p in range(self.n1):
            out.append(f"{indent}A{p} = {aptr}[{p + off}];")
        for q in range(self.n2):
            if bptr == aptr:
                out.append(f"{indent}B{q} = A{q};")
            else:
                out.append(f"{indent}B{q} = {bptr}[{q + off}];")
        if self.ca_m > 0 or (asm and self.signed):
            out.append(f"{indent}ma = FLINT_SIGN_EXT(A{self.n1 - 1});")
        if self.cb_m > 0 or (asm and self.signed):
            out.append(f"{indent}mb = FLINT_SIGN_EXT(B{self.n2 - 1});")
        return out

    def body(self, indent, square=False):
        """Accumulate the products of the loaded A and B. With square=True
        (requires n1 == n2 and B == A), each product A_p A_q with p < q is
        computed once and added twice."""
        n1, n2, s, signed = self.n1, self.n2, self.s, self.signed
        out = []

        # interleave the products of the different chains
        items = []
        for ci, c in enumerate(self.chains):
            for idx, (p, q, o) in enumerate(c["prods"]):
                if square and p > q:
                    continue
                items.append((idx, ci, p, q, o))
        items.sort(key=lambda t: (t[0], t[1]))

        for (idx, ci, p, q, o) in items:
            c = self.chains[ci]
            lo, L = c["lo"], c["L"]
            r = o - lo
            is_toptop = signed and p == n1 - 1 and q == n2 - 1
            dst = [f"c{ci}_{j}" for j in range(r, L)]
            n = len(dst)
            Ap = f"a[{p}]" if self.memops else f"A{p}"
            Bq = f"b[{q}]" if self.memops else f"B{q}"
            if o == s - 1:
                out.append(f"{indent}p0 = {Ap} * {Bq};")
                src = ["p0"]
            else:
                mul = "smul_ppmm" if is_toptop else "umul_ppmm"
                out.append(f"{indent}{mul}(p1, p0, {Ap}, {Bq});")
                if n == 2:
                    src = ["p0", "p1"]
                elif is_toptop:
                    out.append(f"{indent}sx = FLINT_SIGN_EXT(p1);")
                    src = ["p0", "p1"] + ["sx"] * (n - 2)
                else:
                    src = ["p0", "p1"] + ["UWORD(0)"] * (n - 2)
            emit_add(out, dst, src, indent)
            if square and p < q:
                emit_add(out, dst, src, indent)

        # sign corrections (with square=True, the two symmetric sets of
        # corrections coincide, so each is applied twice)
        reps = 2 if square else 1
        if self.inline_sign:
            for q in range(self.ca_m):
                o = n1 + q
                ci = self.chain_for_offset(o)
                c = self.chains[ci]
                dst = [f"c{ci}_{j}" for j in range(o - c["lo"], c["L"])]
                for _ in range(reps):
                    emit_sub(out, dst, [f"(ma & B{q})"] + ["UWORD(0)"] * (len(dst) - 1), indent)
            if not square:
                for p in range(self.cb_m):
                    o = n2 + p
                    ci = self.chain_for_offset(o)
                    c = self.chains[ci]
                    dst = [f"c{ci}_{j}" for j in range(o - c["lo"], c["L"])]
                    emit_sub(out, dst, [f"(mb & A{p})"] + ["UWORD(0)"] * (len(dst) - 1), indent)
        if self.ca_L > 0:
            limbs = [f"(ma & B{q})" for q in range(self.ca_m)] + ["UWORD(0)"] * (self.ca_L - self.ca_m)
            for _ in range(reps):
                emit_add(out, [f"ca{j}" for j in range(self.ca_L)], limbs, indent)
        if self.cb_L > 0 and not square:
            limbs = [f"(mb & A{p})" for p in range(self.cb_m)] + ["UWORD(0)"] * (self.cb_L - self.cb_m)
            emit_add(out, [f"cb{j}" for j in range(self.cb_L)], limbs, indent)
        return out

    def combine(self, indent, r="r", accumulate=False):
        """Set (or add to) r0..r_{s-1} the value of the chains."""
        n1, n2, s = self.n1, self.n2, self.s
        out = []

        def pad(ci, c):
            if c["negative"] and c["lo"] + c["L"] < s:
                out.append(f"{indent}sx = FLINT_SIGN_EXT(c{ci}_{c['L'] - 1});")
                return "sx"
            return "UWORD(0)"

        first = self.chains[0]
        assert first["lo"] == 0
        ext = pad(0, first)
        if accumulate:
            limbs = [f"c0_{j}" if j < first["L"] else ext for j in range(s)]
            emit_add(out, [f"{r}{j}" for j in range(s)], limbs, indent)
        else:
            for j in range(s):
                out.append(f"{indent}{r}{j} = {f'c0_{j}' if j < first['L'] else ext};")
        for ci, c in enumerate(self.chains[1:], start=1):
            lo, L = c["lo"], c["L"]
            ext = pad(ci, c)
            limbs = [f"c{ci}_{j}" if j < L else ext for j in range(s - lo)]
            emit_add(out, [f"{r}{j}" for j in range(lo, s)], limbs, indent)
        if self.ca_L > 0:
            limbs = [f"ca{j}" if j < self.ca_L else "UWORD(0)" for j in range(s - n1)]
            emit_sub(out, [f"{r}{j}" for j in range(n1, s)], limbs, indent)
        if self.cb_L > 0:
            limbs = [f"cb{j}" if j < self.cb_L else "UWORD(0)" for j in range(s - n2)]
            emit_sub(out, [f"{r}{j}" for j in range(n2, s)], limbs, indent)
        return out

    def name(self, kind, direction="rev"):
        sfx = "_signed" if self.signed else ""
        if kind == "kernel":
            mid = {"fwd": "", "rev": "_rev", "strided": "_strided"}[direction]
            return f"_flint_mpn_dot{mid}_{self.n1}x{self.n2}_{self.s}{sfx}"
        if kind == "mul":
            return f"_flint_mpn_poly_mulmid_{self.n1}x{self.n2}_{self.s}{sfx}"
        if kind == "sqr":
            return f"_flint_mpn_poly_sqrmid_{self.n1}x{self.n2}_{self.s}{sfx}"
        if kind == "mul_short":
            return f"_flint_mpn_poly_mulmid_{self.n1}x{self.n2}_{self.s}{sfx}_short"
        if kind == "sqr_short":
            return f"_flint_mpn_poly_sqrmid_{self.n1}x{self.n2}_{self.s}{sfx}_short"
        if kind == "mul_signmag":
            return f"_flint_mpn_poly_mulmid_{self.n1}x{self.n2}_{self.s}_signmag"
        if kind == "sqr_signmag":
            return f"_flint_mpn_poly_sqrmid_{self.n1}x{self.n2}_{self.s}_signmag"
        raise ValueError

    def prototype(self, kind, direction="rev"):
        if kind == "kernel":
            if direction == "strided":
                return f"void {self.name(kind, direction)}(nn_ptr res, nn_srcptr a, slong astride, nn_srcptr b, slong bstride, slong len);"
            return f"void {self.name(kind, direction)}(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len);"
        if kind in ("mul", "mul_short", "mul_signmag"):
            return f"void {self.name(kind)}(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlo, slong nhi);"
        if kind in ("sqr", "sqr_short", "sqr_signmag"):
            return f"void {self.name(kind)}(nn_ptr res, nn_srcptr f, slong flen, slong nlo, slong nhi);"
        raise ValueError

    # -- code generators -------------------------------------------------

    def gen_kernel(self, direction="rev"):
        """direction: "fwd" (b[k]), "rev" (b[len-1-k]) or "strided"
        (explicit limb strides for a and b, which may be negative).
        For n1 >= ASMKN, the products are computed with the assembly
        multiplication routines and accumulated in registers, which was
        measured to be at least as fast as the carry chains for 4 limbs."""
        n1, n2, s = self.n1, self.n2, self.s
        asm = (n1 >= ASMKN)
        out = []
        out.append("static " + self.prototype("kernel", direction)[:-1])
        out.append("{")
        out += ["    " + d for d in self.declarations(asm)]
        out.append("    ulong " + ", ".join(f"r{j}" for j in range(s)) + ";")
        if asm:
            out.append(f"    ulong t[{n1 + n2}];")
        if direction == "strided":
            out.append("    slong k;")
        else:
            out.append("    nn_srcptr aend;")
        out.append("")
        if asm:
            out.append("    " + " = ".join(f"r{j}" for j in range(s)) + " = 0;")
        else:
            out += self.init("    ")
        out.append("")
        if direction == "strided":
            out.append("    for (k = 0; k < len; k++)")
        else:
            out.append(f"    aend = a + len * {n1};")
            if direction == "rev":
                out.append(f"    b += (len - 1) * {n2};")
            out.append("")
            out.append("    for ( ; a != aend; )")
        out.append("    {")
        if not asm or self.signed:
            out += self.loads("        ", asm=asm)
            out.append("")
        if asm:
            out += self.product_accumulate_asm("        ")
        else:
            out += self.body("        ")
        out.append("")
        if direction == "strided":
            out.append("        a += astride;")
            out.append("        b += bstride;")
        else:
            out.append(f"        a += {n1};")
            out.append(f"        b {'-' if direction == 'rev' else '+'}= {n2};")
        out.append("    }")
        out.append("")
        if not asm:
            out += self.combine("    ")
            out.append("")
        for j in range(s):
            out.append(f"    res[{j}] = r{j};")
        out.append("}")
        out.append("")
        return prune_unused("\n".join(out))

    def product_accumulate(self, indent, square=False, signmag=False):
        """Compute the product of the loaded A and B (n1 <= 2) with a
        minimal number of additions and add it to r0..r_{s-1}
        (mod 2^(64 s), with sign extension for two's complement inputs).
        With square=True (B == A), symmetric terms are computed once."""
        n1, n2, s, signed = self.n1, self.n2, self.s, self.signed
        assert n1 <= 2
        out = []
        m = min(n1 + n2, s)      # product limbs P0..P_{m-1}
        P = [f"P{j}" for j in range(m)]

        def mul(hi, lo, x, y, toptop):
            out.append(f"{indent}{'smul_ppmm' if (signed and toptop) else 'umul_ppmm'}({hi}, {lo}, {x}, {y});")

        if n1 == 1:
            if m == 1:
                out.append(f"{indent}P0 = A0 * B0;")
            else:
                mul("P1", "P0", "A0", "B0", True)
        elif n2 == 1:
            # A0 B0 at offset 0, A1 B0 at offset 1 (A1 and B0 are top limbs)
            mul("P1", "P0", "A0", "B0", False)
            if m == 2:
                out.append(f"{indent}P1 += A1 * B0;")
            else:
                mul("P2", "t0", "A1", "B0", True)
                out.append(f"{indent}add_ssaaaa(P2, P1, P2, P1, UWORD(0), t0);")
            if signed:
                # mixed term A0 x B0 where B0 is the signed top limb
                if m == 2:
                    out.append(f"{indent}P1 -= (mb & A0);")
                else:
                    out.append(f"{indent}sub_ddmmss(P2, P1, P2, P1, UWORD(0), (mb & A0));")
        else:
            # 2 x 2: A0 B0 (offset 0), A0 B1 + A1 B0 (offset 1), A1 B1 (offset 2)
            mul("P1", "P0", "A0", "B0", False)
            if m == 3:
                out.append(f"{indent}P2 = A1 * B1;")
                mul("t1", "t0", "A0", "B1", False)
                out.append(f"{indent}add_ssaaaa(P2, P1, P2, P1, t1, t0);")
                if square:
                    out.append(f"{indent}add_ssaaaa(P2, P1, P2, P1, t1, t0);")
                else:
                    mul("t1", "t0", "A1", "B0", False)
                    out.append(f"{indent}add_ssaaaa(P2, P1, P2, P1, t1, t0);")
                if signed:
                    out.append(f"{indent}P2 -= (ma & B0);")
                    out.append(f"{indent}P2 -= (mb & A0);")
            else:
                mul("P3", "P2", "A1", "B1", True)
                mul("t1", "t0", "A0", "B1", False)
                out.append(f"{indent}add_sssaaaaaa(P3, P2, P1, P3, P2, P1, UWORD(0), t1, t0);")
                if square:
                    out.append(f"{indent}add_sssaaaaaa(P3, P2, P1, P3, P2, P1, UWORD(0), t1, t0);")
                else:
                    mul("t1", "t0", "A1", "B0", False)
                    out.append(f"{indent}add_sssaaaaaa(P3, P2, P1, P3, P2, P1, UWORD(0), t1, t0);")
                if signed:
                    out.append(f"{indent}sub_ddmmss(P3, P2, P3, P2, UWORD(0), (ma & B0));")
                    out.append(f"{indent}sub_ddmmss(P3, P2, P3, P2, UWORD(0), (mb & A0));")

        # accumulate
        if signmag:
            # r += (P ^ ms), i.e. r += ~P = -P - 1 when ms = -1; the
            # missing +1 completions are counted in cnt (one carry-free
            # op instead of an s-limb subtraction per term) and added
            # once after the loop
            for j in range(m):
                out.append(f"{indent}P{j} ^= ms;")
            src = P + ["ms"] * (s - m)
            emit_add(out, [f"r{j}" for j in range(s)], src, indent)
            if not square:
                out.append(f"{indent}cnt -= ms;")
            return out
        if s > m:
            if signed:
                out.append(f"{indent}sx = FLINT_SIGN_EXT(P{m - 1});")
                pad = "sx"
            else:
                pad = "UWORD(0)"
            src = P + [pad] * (s - m)
        else:
            src = P
        emit_add(out, [f"r{j}" for j in range(s)], src, indent)
        return out

    def product_accumulate_asm(self, indent, square=False, signmag=False):
        """Compute the product of the terms at a and b with the assembly
        multiplication routines into t[] and add it to r0..r_{s-1}
        (mod 2^(64 s)), correcting for two's complement inputs. A and B
        must have been loaded when signed. Requires s <= 8."""
        n1, n2, s, signed = self.n1, self.n2, self.s, self.signed
        out = []
        m = min(n1 + n2, s)
        off = 1 if signmag else 0
        ap = f"a + {off}" if signmag else "a"
        bp = f"b + {off}" if signmag else "b"
        if square:
            out.append(f"{indent}flint_mpn_sqr(t, {ap}, {n1});")
        elif n1 == n2:
            out.append(f"{indent}flint_mpn_mul_n(t, {ap}, {bp}, {n1});")
        else:
            out.append(f"{indent}flint_mpn_mul(t, {ap}, {n1}, {bp}, {n2});")
        if signmag:
            # xor in the operands (rather than in memory) to avoid
            # store-forwarding stalls; the +1 completions of the negated
            # products are counted in cnt and added once after the loop
            src = [f"(t[{j}] ^ ms)" for j in range(m)] + ["ms"] * (s - m)
            emit_add(out, [f"r{j}" for j in range(s)], src, indent)
            if not square:
                out.append(f"{indent}cnt -= ms;")
            return out
        src = [f"t[{j}]" for j in range(m)] + ["UWORD(0)"] * (s - m)
        emit_add(out, [f"r{j}" for j in range(s)], src, indent)
        if signed:
            # a b = A B - 2^(64 n1) (ma & B) - 2^(64 n2) (mb & A) + 2^(64 (n1+n2)) (ma & mb)
            reps = 2 if square else 1
            if n1 < s:
                dst = [f"r{j}" for j in range(n1, s)]
                src = [f"(ma & B{q})" if q < n2 else "UWORD(0)" for q in range(s - n1)]
                for _ in range(reps):
                    emit_sub(out, dst, src, indent)
            if not square and n2 < s:
                dst = [f"r{j}" for j in range(n2, s)]
                src = [f"(mb & A{p})" if p < n1 else "UWORD(0)" for p in range(s - n2)]
                emit_sub(out, dst, src, indent)
            if s > n1 + n2:
                out.append(f"{indent}r{n1 + n2} -= (ma & {'ma' if square else 'mb'});")
        return out

    def use_acc(self):
        """Whether to accumulate full products (computed with a minimal
        number of additions) into a single register accumulator, instead
        of keeping one chain per offset across the dot product. The former
        has less per-output overhead and fewer live registers and was
        measured to be faster for the 1- and 2-limb shapes at all lengths;
        for 3 limbs the full product needs too many additions per term."""
        return self.n1 <= ACCN

    def gen_mul_loop(self, asm=False, signmag=False):
        n1, n2, s = self.n1, self.n2, self.s
        assert not (signmag and self.signed)
        if signmag:
            asm = (n1 >= 3)
        st1, st2 = (n1 + 1, n2 + 1) if signmag else (n1, n2)
        out = []
        out.append("static " + self.prototype("mul_signmag" if signmag else ("mul_short" if asm else "mul"))[:-1])
        out.append("{")
        out.append("    slong i, top1, top2;")
        out.append("    nn_srcptr a, b, aend;")
        out += ["    " + d for d in self.declarations(asm, loop=True, signmag=signmag)]
        out.append("    ulong " + ", ".join(f"r{j}" for j in range(s)) + ";")
        if asm:
            out.append(f"    ulong t[{n1 + n2}];")
        out.append("")
        out.append("    for (i = nlo; i < nhi; i++)")
        out.append("    {")
        out.append("        top1 = FLINT_MIN(flen - 1, i);")
        out.append("        top2 = FLINT_MIN(glen - 1, i);")
        out.append(f"        a = f + (i - top2) * {st1};")
        out.append(f"        aend = a + (top1 + top2 - i + 1) * {st1};")
        out.append(f"        b = g + top2 * {st2};")
        out.append("")
        acc = asm or self.use_acc() or signmag
        if acc:
            out.append("        " + " = ".join(f"r{j}" for j in range(s)) + " = 0;")
        else:
            out += self.init("        ")
        if signmag:
            out.append("        cnt = 0;")
        out.append("")
        out.append("        for ( ; a != aend; )")
        out.append("        {")
        self.memops = (not asm and not self.signed and not signmag
                       and not self.use_acc() and self.n1 == 3)
        if (not asm or self.signed or signmag) and not self.memops:
            out += self.loads("            ", asm=asm, signmag=signmag)
            out.append("")
        if asm:
            out += self.product_accumulate_asm("            ", signmag=signmag)
        elif self.use_acc() or signmag:
            out += self.product_accumulate("            ", signmag=signmag)
        else:
            out += self.body("            ")
        self.memops = False
        out.append("")
        out.append(f"            a += {st1};")
        out.append(f"            b -= {st2};")
        out.append("        }")
        out.append("")
        if not acc:
            out += self.combine("        ")
            out.append("")
        if signmag:
            rr = [f"r{j}" for j in range(s)]
            emit_add(out, rr, ["cnt"] + ["UWORD(0)"] * (s - 1), "        ")
            out.append("")
        for j in range(s):
            out.append(f"        res[{j}] = r{j};")
        out.append(f"        res += {s};")
        out.append("    }")
        out.append("}")
        out.append("")
        return prune_unused("\n".join(out))

    def gen_sqr_loop(self, asm=False, signmag=False):
        n, s = self.n1, self.s
        assert self.n1 == self.n2
        assert not (signmag and self.signed)
        if signmag:
            asm = (n >= 3)
        st = n + 1 if signmag else n
        out = []
        out.append("static " + self.prototype("sqr_signmag" if signmag else ("sqr_short" if asm else "sqr"))[:-1])
        out.append("{")
        out.append("    slong i, start, stop;")
        out.append("    nn_srcptr a, b, aend;")
        out += ["    " + d for d in self.declarations(asm, loop=True, signmag=signmag)]
        out.append("    ulong " + ", ".join(f"r{j}" for j in range(s)) + ";")
        if asm or (not self.use_acc() and not self.signed):
            out.append(f"    ulong t[{2 * n}];")
        out.append("")
        out.append("    for (i = nlo; i < nhi; i++)")
        out.append("    {")
        out.append("        /* cross terms f_j f_{i-j} with j < i - j, counted twice */")
        out.append("        start = FLINT_MAX(0, i - flen + 1);")
        out.append("        stop = FLINT_MIN(flen - 1, (i + 1) / 2 - 1);")
        out.append("")
        acc = asm or self.use_acc() or signmag
        if acc:
            out.append("        " + " = ".join(f"r{j}" for j in range(s)) + " = 0;")
        else:
            out += self.init("        ")
        if signmag:
            out.append("        cnt = 0;")
        out.append("")
        out.append(f"        a = f + start * {st};")
        out.append(f"        aend = f + (stop + 1) * {st};")
        out.append(f"        b = f + (i - start) * {st};")
        out.append("")
        out.append("        for ( ; a < aend; )")
        out.append("        {")
        self.memops = False
        if (not asm or self.signed or signmag) and not self.memops:
            out += self.loads("            ", asm=asm, signmag=signmag)
            out.append("")
        if asm:
            out += self.product_accumulate_asm("            ", signmag=signmag)
        elif self.use_acc() or signmag:
            out += self.product_accumulate("            ", signmag=signmag)
        else:
            out += self.body("            ")
        self.memops = False
        out.append("")
        out.append(f"            a += {st};")
        out.append(f"            b -= {st};")
        out.append("        }")
        out.append("")
        if not acc:
            out += self.combine("        ")
            out.append("")
        out.append("        /* double */")
        rr = [f"r{j}" for j in range(s)]
        if signmag:
            emit_add(out, rr, ["cnt"] + ["UWORD(0)"] * (s - 1), "        ")
        emit_add(out, rr, rr, "        ")
        out.append("")
        out.append("        /* square term */")
        out.append("        if (i % 2 == 0)")
        out.append("        {")
        out.append(f"            a = f + (i / 2) * {st};")
        if signmag:
            # the square is nonnegative: ms = 0
            out += self.loads("            ", "a", "a", asm=asm, signmag=True)
            if asm:
                out += self.product_accumulate_asm("            ", square=True, signmag=True)
            else:
                out += self.product_accumulate("            ", square=True, signmag=True)
        elif asm:
            if self.signed:
                out += self.loads("            ", "a", "a", asm=True)
            out += self.product_accumulate_asm("            ", square=True)
        elif self.use_acc():
            out += self.loads("            ", "a", "a")
            out += self.product_accumulate("            ", square=True)
        elif not self.signed:
            # the assembly squaring routine is faster than inline code
            # for 3 or more limbs
            m = min(2 * n, s)
            out.append(f"            flint_mpn_sqr(t, a, {n});")
            src = [f"t[{j}]" for j in range(m)] + ["UWORD(0)"] * (s - m)
            emit_add(out, [f"r{j}" for j in range(s)], src, "            ")
        else:
            out += self.loads("            ", "a", "a")
            out += self.init("            ")
            out += self.body("            ", square=True)
            out += self.combine("            ", accumulate=True)
        out.append("        }")
        out.append("")
        for j in range(s):
            out.append(f"        res[{j}] = r{j};")
        out.append(f"        res += {s};")
        out.append("    }")
        out.append("}")
        out.append("")
        return prune_unused("\n".join(out))


HEADER = '''/*
    Copyright (C) 2024, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* This file is generated by dev/gen_mpn_dot_rev.py. Do not edit by hand. */

#include "longlong.h"
#include "mpn_extras.h"

'''


def kernel_shapes(maxn):
    for n1 in range(1, maxn + 1):
        for n2 in range(1, n1 + 1):
            for s in (n1 + n2 - 1, n1 + n2, n1 + n2 + 1):
                if s <= MAXS:
                    yield n1, n2, s


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "kernels"
    out = []

    if mode == "prototypes":
        # extern declarations of the dispatch tables, for mpn_extras.h
        out.append(f"FLINT_DLL extern const flint_mpn_dot_strided_func_t flint_mpn_dot_strided_tab[2][{MAXN + 1}][{MAXN + 1}][3];")
        out.append(f"FLINT_DLL extern const flint_mpn_dot_func_t flint_mpn_dot_tab[2][{DEDN + 1}][{DEDN + 1}][3];")
        out.append(f"FLINT_DLL extern const flint_mpn_dot_func_t flint_mpn_dot_rev_tab[2][{DEDN + 1}][{DEDN + 1}][3];")
        out.append(f"FLINT_DLL extern const flint_mpn_poly_mulmid_func_t flint_mpn_poly_mulmid_tab[2][{LOOPN + 1}][{LOOPN + 1}][3];")
        out.append(f"FLINT_DLL extern const flint_mpn_poly_sqrmid_func_t flint_mpn_poly_sqrmid_tab[2][{LOOPN + 1}][{LOOPN + 1}][3];")
        out.append(f"FLINT_DLL extern const flint_mpn_poly_mulmid_func_t flint_mpn_poly_mulmid_short_tab[2][{max(ASMN) + 1}][{max(ASMN) + 1}][3];")
        out.append(f"FLINT_DLL extern const flint_mpn_poly_sqrmid_func_t flint_mpn_poly_sqrmid_short_tab[2][{max(ASMN) + 1}][{max(ASMN) + 1}][3];")
        out.append(f"FLINT_DLL extern const flint_mpn_poly_mulmid_func_t flint_mpn_poly_mulmid_signmag_tab[{SIGNMAGN + 1}][{SIGNMAGN + 1}][3];")
        out.append(f"FLINT_DLL extern const flint_mpn_poly_sqrmid_func_t flint_mpn_poly_sqrmid_signmag_tab[{SIGNMAGN + 1}][{SIGNMAGN + 1}][3];")
        print("\n".join(out))
        return

    out.append(HEADER)

    if mode == "kernels":
        for signed in (0, 1):
            for (n1, n2, s) in kernel_shapes(MAXN):
                K = Kernel(n1, n2, s, signed)
                out.append(K.gen_kernel("strided"))
                if n1 <= DEDN:
                    out.append(K.gen_kernel("fwd"))
                    out.append(K.gen_kernel("rev"))
        # tables indexed [sgn][n1][n2][s - (n1 + n2 - 1)]
        def table(tabname, direction, maxn, dim):
            out.append(f"const {'flint_mpn_dot_strided_func_t' if direction == 'strided' else 'flint_mpn_dot_func_t'} {tabname}[2][{dim}][{dim}][3] = {{")
            for signed in (0, 1):
                out.append("    {")
                for n1 in range(dim):
                    out.append("        {")
                    for n2 in range(dim):
                        row = []
                        for i in range(3):
                            s = n1 + n2 - 1 + i
                            if 1 <= n2 <= n1 <= maxn and s <= MAXS:
                                row.append(Kernel(n1, n2, s, signed).name("kernel", direction))
                            else:
                                row.append("NULL")
                        out.append("            { " + ", ".join(row) + " },")
                    out.append("        },")
                out.append("    },")
            out.append("};")
            out.append("")
        table("flint_mpn_dot_strided_tab", "strided", MAXN, MAXN + 1)
        table("flint_mpn_dot_tab", "fwd", DEDN, DEDN + 1)
        table("flint_mpn_dot_rev_tab", "rev", DEDN, DEDN + 1)
    elif mode == "loops":
        for signed in (0, 1):
            for (n1, n2, s) in kernel_shapes(LOOPN):
                K = Kernel(n1, n2, s, signed)
                out.append(K.gen_mul_loop())
                if n1 == n2:
                    out.append(K.gen_sqr_loop())
        for signed in (0, 1):
            for (n1, n2, s) in kernel_shapes(max(ASMN)):
                if n1 not in ASMN:
                    continue
                K = Kernel(n1, n2, s, signed)
                out.append(K.gen_mul_loop(asm=True))
                if n1 == n2:
                    out.append(K.gen_sqr_loop(asm=True))
        for (n1, n2, s) in kernel_shapes(SIGNMAGN):
            K = Kernel(n1, n2, s, 0)
            out.append(K.gen_mul_loop(signmag=True))
            if n1 == n2:
                out.append(K.gen_sqr_loop(signmag=True))

        # dispatch tables indexed [sgn][n1][n2][s - (n1 + n2 - 1)]
        # (signmag tables have no sign dimension), with NULL where no
        # kernel exists
        def poly_table(tabname, kind, shapes_n, restrict, dim, sgndim, sq):
            functype = "flint_mpn_poly_sqrmid_func_t" if sq else "flint_mpn_poly_mulmid_func_t"
            out.append(f"const {functype} {tabname}{'[2]' if sgndim else ''}[{dim}][{dim}][3] = {{")
            for signed in ((0, 1) if sgndim else (0,)):
                if sgndim:
                    out.append("    {")
                pad = "        " if sgndim else "    "
                for n1 in range(dim):
                    out.append(pad + "{")
                    for n2 in range(dim):
                        row = []
                        for q in range(3):
                            s = n1 + n2 - 1 + q
                            ok = 1 <= n2 <= n1 and (n1, n2, s) in kernel_shapes(shapes_n) and restrict(n1) and (not sq or n1 == n2)
                            row.append(Kernel(n1, n2, s, signed).name(kind) if ok else "NULL")
                        out.append(pad + "    { " + ", ".join(row) + " },")
                    out.append(pad + "},")
                if sgndim:
                    out.append("    },")
            out.append("};")
            out.append("")
        poly_table("flint_mpn_poly_mulmid_tab", "mul", LOOPN, lambda n: True, LOOPN + 1, 1, 0)
        poly_table("flint_mpn_poly_sqrmid_tab", "sqr", LOOPN, lambda n: True, LOOPN + 1, 1, 1)
        poly_table("flint_mpn_poly_mulmid_short_tab", "mul_short", max(ASMN), lambda n: n in ASMN, max(ASMN) + 1, 1, 0)
        poly_table("flint_mpn_poly_sqrmid_short_tab", "sqr_short", max(ASMN), lambda n: n in ASMN, max(ASMN) + 1, 1, 1)
        poly_table("flint_mpn_poly_mulmid_signmag_tab", "mul_signmag", SIGNMAGN, lambda n: True, SIGNMAGN + 1, 0, 0)
        poly_table("flint_mpn_poly_sqrmid_signmag_tab", "sqr_signmag", SIGNMAGN, lambda n: True, SIGNMAGN + 1, 0, 1)
    else:
        raise SystemExit("unknown mode " + mode)

    sys.stdout.write("\n".join(out))


if __name__ == "__main__":
    main()
