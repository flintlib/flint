dnl
dnl Copyright (C) 2024 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

define(`rp',  `%rdi')
define(`ap',  `%rsi')
define(`bp',  `%rdx')
define(`cnt', `%rcx')

define(`tnc', `%r8')
define(`sx', `%rax')

define(`s0', `%r9')
define(`s1', `%r10')
define(`s2', `%r11')

dnl r <- a +/- 2^n b
dnl
dnl For 0 <= i < n - 1, we have
dnl
dnl     r_{i} = a_{i} +/- (b_{i} >> n + b_{i + 1} << (64 - n)),
dnl
dnl and
dnl
dnl     r_{n - 1} = a_{n - 1} +/- (b_{n - 1} >> n).

dnl The idea is the following:
dnl
dnl Assume that bp[i] is loaded in a register b0.
dnl
dnl t = b0 >> n		C shrx
dnl b1 = bp[i + 1]	C mov, and fullfills assumption for next iteration
dnl s = b1 << (64 - n)	C shlx
dnl s = s + t		C lea, carry-less
dnl if OP = add, then
dnl   s += ap[i]	C adc
dnl   rp[i] = s		C mov
dnl else
dnl   u = ap[i]		C mov
dnl   u -= s		C sbb
dnl   rp[i] = u		C mov
dnl fi

define(ALL_AORS,`
	ALIGN(16)
PROLOGUE(flint_mpn_addrsh_1)
	shrx	cnt, 0*8(bp), s0
	xor	R32(sx), R32(sx)
	add	0*8(ap), s0
	mov	s0, 0*8(rp)
	setc	R8(sx)
	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_addrsh_2)
	xor	R32(tnc), R32(tnc)
	sub	cnt, tnc

	xor	R32(sx), R32(sx)

	mov	1*8(bp), s1
C
	shrx	cnt, 0*8(bp), s0
	shlx	tnc, s1, s2
	shrx	cnt, s1, s1
	C (0, 2), 1

	adox	s2, s0
C
	adcx	0*8(ap), s0
	mov	s0, 0*8(rp)
	adox	sx, s1		C cannot overflow
	adcx	1*8(ap), s1
C
	mov	s1, 1*8(rp)

	setc	R8(sx)
	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_addrsh_3)
	xor	R32(tnc), R32(tnc)
	sub	cnt, tnc
	xor	R32(sx), R32(sx)

	shrx	cnt, 0*8(bp), s0
	mov	1*8(bp), s1
	shlx	tnc, s1, s2
	lea	(s0, s2), s2
ifelse(OP,`add',`
	add	0*8(ap), s2
	mov	s2, 0*8(rp)
',`
	mov	0*8(ap), s0
	sub	s2, s0
	mov	s0, 0*8(rp)
')
	C Used: s1

	shrx	cnt, s1, s1
	mov	2*8(bp), s2
	shlx	tnc, s2, s0
	lea	(s1, s0), s0
ifelse(OP,`add',`
	adc	1*8(ap), s0
	mov	s0, 1*8(rp)
',`
	mov	1*8(ap), s1
	sbb	s0, s1
	mov	s1, 1*8(rp)
')
	C Used: s2

	shrx	cnt, s2, s2
ifelse(OP,`add',`
	adc	2*8(ap), s2
	mov	s2, 2*8(rp)
',`
	mov	2*8(ap), s0
	sbb	s2, s0
	mov	s0, 1*8(rp)
')

	setc	R8(sx)
	ret
EPILOGUE()
')

	ALIGN(16)
PROLOGUE(flint_mpn_addrsh_4)
	xor	R32(tnc), R32(tnc)
	sub	cnt, tnc
	xor	R32(sx), R32(sx)

	shrx	cnt, 0*8(bp), s0
	mov	1*8(bp), s1
	shlx	tnc, s1, s2
	lea	(s0, s2), s2
ifelse(OP,`add',`
	add	0*8(ap), s2
	mov	s2, 0*8(rp)
',`
	mov	0*8(ap), s0
	sub	s2, s0
	mov	s0, 0*8(rp)
')
	C Used: s1

	shrx	cnt, s1, s1
	mov	2*8(bp), s2
	shlx	tnc, s2, s0
	lea	(s1, s0), s0
ifelse(OP,`add',`
	adc	1*8(ap), s0
	mov	s0, 1*8(rp)
',`
	mov	1*8(ap), s1
	sbb	s0, s1
	mov	s1, 1*8(rp)
')
	C Used: s2

C
	shrx	cnt, s1, s1
	mov	3*8(bp), s2
	shlx	tnc, s2, s0
	lea	(s1, s0), s0
ifelse(OP,`add',`
	adc	2*8(ap), s0
	mov	s0, 2*8(rp)
',`
	mov	2*8(ap), s1
	sbb	s0, s1
	mov	s1, 2*8(rp)
')
	C Used: s2
C

	shrx	cnt, s2, s2
ifelse(OP,`add',`
	adc	2*8(ap), s2
	mov	s2, 2*8(rp)
',`
	mov	2*8(ap), s0
	sbb	s2, s0
	mov	s0, 1*8(rp)
')

	setc	R8(sx)
	ret
EPILOGUE()
')

	TEXT
define(`flint_mpn_aorsrsh',`flint_mpn_addrsh_$1')
define(`OP',`add')
define(`OPC',`adc')
ALL_AORSRSH
undefine(`flint_mpn_aorsrsh')
undefine(`OP')
undefine(`OPC')

define(`flint_mpn_aorsrsh',`flint_mpn_subrsh_$1')
define(`OP',`sub')
define(`OPC',`sbb')
ALL_AORSRSH
undefine(`flint_mpn_aorsrsh')
undefine(`OP')
undefine(`OPC')
