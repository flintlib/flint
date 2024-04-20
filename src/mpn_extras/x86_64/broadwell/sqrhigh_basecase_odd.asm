dnl  X64-64 mpn_mullo_basecase optimised for Intel Broadwell.

dnl  Contributed to the GNU project by Torbjorn Granlund.

dnl  Copyright 2017 Free Software Foundation, Inc.

dnl  This file is part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can redistribute it and/or modify
dnl  it under the terms of either:
dnl
dnl    * the GNU Lesser General Public License as published by the Free
dnl      Software Foundation; either version 3 of the License, or (at your
dnl      option) any later version.
dnl
dnl  or
dnl
dnl    * the GNU General Public License as published by the Free Software
dnl      Foundation; either version 2 of the License, or (at your option) any
dnl      later version.
dnl
dnl  or both in parallel, as here.
dnl
dnl  The GNU MP Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl  for more details.
dnl
dnl  You should have received copies of the GNU General Public License and the
dnl  GNU Lesser General Public License along with the GNU MP Library.  If not,
dnl  see https://www.gnu.org/licenses/.
dnl
dnl   Copyright (C) 2024 Albin Ahlbäck
dnl
dnl   This file is part of FLINT.
dnl
dnl   FLINT is free software: you can redistribute it and/or modify it under
dnl   the terms of the GNU Lesser General Public License (LGPL) as published
dnl   by the Free Software Foundation; either version 3 of the License, or
dnl   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

define(`rp',	   `%rdi')
define(`ap_param', `%rsi')
define(`no2_param',`%rdx')

define(`ap',	   `%r8')
define(`no2',	   `%rcx')

define(`ix',	   `no2')
define(`ix_save',  `ap_param')
define(`mx',	   `%rbx')
define(`mx_final', `%rbp')

define(`jmpreg',   `%r14')

define(`rx',	   `%rax')

define(`r0',	   `%r9')
define(`r1',	   `%r10')
define(`r2',	   `%r11')
define(`r3a',	   `%r12')
define(`r3b',	   `%r13')

dnl - Use mx_final as a comparison against m to see when to do final iteration
dnl   and when to jump to diagonal.
dnl - Use mx to reset ap. Increases by 2*8 in between each loop.
dnl - Use ix for loop variable (should be rcx).
dnl - Use ix_save for saving initial n in between each loop.
dnl - jmpreg for jumping to the next label

dnl Scheme for odd n
dnl
dnl n = 9:
dnl   0 1 2 3 4 5 6 7 8
dnl 0
dnl 1
dnl 2
dnl 3
dnl 4       h d
dnl 5     h x x d
dnl 6   h x x x x d
dnl 7 h x x x x x x d
dnl 8 x x x x x x x x d
dnl
dnl n = 11:
dnl    0 1 2 3 4 5 6 7 8 9 10
dnl 0
dnl 1
dnl 2
dnl 3
dnl 4
dnl 5          h d
dnl 6        h x x d
dnl 7      h x x x x d
dnl 8    h x x x x x x d
dnl 9  h x x x x x x x x d
dnl 10 x x x x x x x x x x d

	TEXT
	ALIGN(32)
PROLOGUE(_flint_mpn_sqrhigh_basecase_odd)
	mov	R32(no2_param), R32(no2)
	lea	(ap_param,no2_param,8), ap	C ap += no2 where no2 = (n - 1) / 2

	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	C Initial triangle
	C
	C     v--< Top is (no2, no2 - 1) in matrix
	C     h
	C   h x x
	C h x x x x
define(`sx', `rx')
define(`s0', `r0')
define(`s1', `r1')
define(`s2', `r2')
define(`s3', `mx')
define(`s4', `r3a')
define(`s5', `r3b')
define(`s6', `jmpreg')
define(`s7', `mx_final')
define(`zr', `ix_save')
	mov	0*8(ap), %rdx
	mulx	-1*8(ap), sx, sx
	C x

	mov	1*8(ap), %rdx
	xor	R32(zr), R32(zr)
	mulx	-2*8(ap), s2, s2
	mulx	-1*8(ap), s3, s0
	mulx	0*8(ap), s4, s1
	add	s2, sx
	adc	zr, s0
	add	s3, sx
	adc	s4, s0
	adc	zr, s1
	C x, 0, 1

	mov	2*8(ap), %rdx
	mulx	-3*8(ap), s4, s4
	mulx	-2*8(ap), s2, s5
	mulx	-1*8(ap), s3, s6
	add	s2, s4
	adc	s3, s5
	mulx	0*8(ap), s7, s2
	mulx	1*8(ap), %rdx, s3
	adc	s7, s6
	adc	%rdx, s2
	adc	zr, s3
	add	s4, sx
	adc	s5, s0
	adc	s6, s1
	mov	s0, 0*8(rp)
	mov	s1, 1*8(rp)
	adc	zr, s2
	adc	zr, s3
	mov	s2, 2*8(rp)
	mov	s3, 3*8(rp)
	C x, 0, 1, 2, 3
undefine(`sx')
undefine(`s0')
undefine(`s1')
undefine(`s2')
undefine(`s3')
undefine(`s4')
undefine(`s5')
undefine(`s6')
undefine(`s7')
undefine(`zr')

	C Addmul chain for lower mid triangle into rp
	C
	C In the initial loop, we want to do 7 multiplications, of which one is
	C a half-multiplication (only upper result). Throughout all loops, the
	C first two multiplications contribute to rx, so these has to be
	C gathered. Moreover, the last multiplication is "disjoint" from rp, so
	C this cannot be in the loop as we have to add two zeroes here.
	C
	C mx       = -7 * 8
	C mx_final = -(no2 * 2 - 1) * 8
	C ix       = 0
	C ix_save  = 0			NOTE: It is already zero
	lea	-4*8(ap), ap
	lea	(,no2,2), R32(mx_final)
	mov	7*8(ap), %rdx
	lea	-1*8(,mx_final,8), R32(mx_final)
	mulx	0*8(ap), r1, r1
	mulx	1*8(ap), r2, r3a
	neg	mx_final
	xor	R32(ix), R32(ix)	C ix <- 0 and reset flags
	adcx	r1, rx
	movq	$-7*8, mx
	adox	r2, rx
	C jmp L(f4)

L(f4):	lea	-2*8(ap), ap
	lea	-4*8(rp), rp
	lea	L(f6)(%rip), jmpreg
	jmp	L(a4)

L(f6):	lea	-2*8(rp), rp
	mov	r3a, r3b
	lea	L(f0)(%rip), jmpreg
	lea	1(ix_save), ix_save
	jmp	L(a6)

L(end):	adox	r3a, r0
	adox	ix, r1
	adc	ix, r0
	adc	ix, r1
	cmp	R32(mx), R32(mx_final)
	mov	r0, 0*8(rp)
	mov	r1, 1*8(rp)
	lea	3*8(rp,mx), rp		C Reset rp
	jg	L(lsh)
	mov	2*8(ap), %rdx		C Load second factor
	lea	(ap,mx), ap		C Reset ap
	je	L(jmp)
	mulx	0*8(ap), ix, ix
L(jmp):	mulx	1*8(ap), r2, r3a
	test	%al, %al
	adcx	ix, rx
	lea	-2*8(mx), mx
	mov	R32(ix_save), R32(ix)
	adox	r2, rx
	jmp	*jmpreg

L(f0):	lea	2*8(ap), ap
	lea	L(f2)(%rip), jmpreg
	C jmp	L(a0)

	ALIGN(32)
L(a0):	mulx	0*8(ap), r0, r1
	jrcxz	L(end)
	mulx	1*8(ap), r2, r3b
	adox	0*8(rp), r0
	lea	-1(ix), ix
	adox	1*8(rp), r2
	adcx	r3a, r0
	adcx	r1, r2
	mov	r0, 0*8(rp)
	mov	r2, 1*8(rp)
L(a6):	mulx	2*8(ap), r0, r1
	mulx	3*8(ap), r2, r3a
	adox	2*8(rp), r0
	adox	3*8(rp), r2
	adcx	r3b, r0
	adcx	r1, r2
	mov	r0, 2*8(rp)
	mov	r2, 3*8(rp)
L(a4):	mulx	4*8(ap), r0, r1
	mulx	5*8(ap), r2, r3b
	adox	4*8(rp), r0
	adox	5*8(rp), r2
	adcx	r3a, r0
	adcx	r1, r2
	mov	r0, 4*8(rp)
	mov	r2, 5*8(rp)
	lea	8*8(ap), ap
L(a2):	mulx	-2*8(ap), r0, r1
	mulx	-1*8(ap), r2, r3a
	adox	6*8(rp), r0
	adox	7*8(rp), r2
	adcx	r3b, r0
	lea	8*8(rp), rp
	adcx	r1, r2
	mov	r0, -2*8(rp)
	mov	r2, -1*8(rp)
	jmp	L(a0)

L(f2):	lea	4*8(ap), ap
	lea	-6*8(rp), rp
	mov	r3a, r3b
	lea	L(f4)(%rip), jmpreg
	jmp	L(a2)

define(`r3', `jmpreg')
define(`ld0', `r3a')
define(`ld1', `r3b')
define(`ld2', `mx')
define(`ld3', `ix_save')

	C Left shift rp by one and add diagonal
L(lsh):	sar	$4, mx_final
	lea	L(dtab)(%rip), ld2
	lea	1*8(ap,mx_final,8), ap
	neg	mx_final
	mov	0*8(rp), ld1
	or	R32(mx_final), R32(ix)
	shr	$2, R32(ix)
	and	$3, R32(mx_final)	C Also reset flags

	mov	0*8(ap), %rdx
	mulx	%rdx, r0, r1
ifdef(`PIC',
`	movslq	(ld2,mx_final,4), mx_final
',`')
	adox	rx, rx
	adox	ld1, ld1
ifdef(`PIC',
`	lea	(mx_final,ld2), ld2
',`')
	adcx	r0, rx
	adcx	r1, ld1
ifdef(`PIC',
`	jmp	*ld2
',`
	jmp	*(ld2,mx_final,8)
')

L(dp3):	lea	-1*8(ap), ap
	lea	-2*8(rp), rp
	mov	ld1, ld3
	jmp	L(ds3)

L(dp2):	lea	2*8(ap), ap
	lea	-4*8(rp), rp
	jmp	L(ds2)

L(dp1):	lea	1*8(ap), ap
	lea	2*8(rp), rp
	mov	ld1, ld3
	jmp	L(ds1)

	ALIGN(32)
L(dp0):	C Do nothing
L(ds0):	mov	1*8(ap), %rdx
	mov	1*8(rp), ld2
	mov	2*8(rp), ld3
	mulx	%rdx, r2, r3
	lea	-1(ix), R32(ix)
	adox	ld2, ld2
	adox	ld3, ld3
	mov	ld1, 0*8(rp)
	adcx	r2, ld2
	adcx	r3, ld3
	mov	ld2, 1*8(rp)
L(ds3):	mov	2*8(ap), %rdx
	mov	3*8(rp), ld0
	mov	4*8(rp), ld1
	mov	ld3, 2*8(rp)
	mulx	%rdx, r0, r1
	adox	ld0, ld0
	adox	ld1, ld1
	lea	4*8(ap), ap
	adcx	r0, ld0
	adcx	r1, ld1
	mov	ld0, 3*8(rp)
L(ds2):	mov	-1*8(ap), %rdx
	mov	5*8(rp), ld2
	mov	6*8(rp), ld3
	mov	ld1, 4*8(rp)
	mulx	%rdx, r2, r3
	adox	ld2, ld2
	adox	ld3, ld3
	lea	8*8(rp), rp
	adcx	r2, ld2
	adcx	r3, ld3
	mov	ld2, -3*8(rp)
L(ds1):	mov	0*8(ap), %rdx
	mov	-1*8(rp), ld0
	mov	0*8(rp), ld1
	mov	ld3, -2*8(rp)
	mulx	%rdx, r0, r1
	adox	ld0, ld0
	jrcxz	L(fin)
	adox	ld1, ld1
	adcx	r0, ld0
	adcx	r1, ld1
	mov	ld0, -1*8(rp)
	jmp	L(ds0)

L(fin):	adox	ix, r1		C ix = 0
	adc	r0, ld0
	adc	ix, r1
	mov	ld0, -1*8(rp)
	mov	r1, 0*8(rp)

undefine(`r3')
undefine(`ld0')
undefine(`ld1')
undefine(`ld2')
undefine(`ld3')

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()
	JUMPTABSECT
	ALIGN(8)
L(dtab):JMPENT(	L(dp0), L(dtab))
	JMPENT(	L(dp1), L(dtab))
	JMPENT(	L(dp2), L(dtab))
	JMPENT(	L(dp3), L(dtab))
