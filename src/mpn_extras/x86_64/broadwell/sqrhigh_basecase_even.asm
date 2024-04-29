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

dnl Scheme for even n
dnl
dnl n = 10:
dnl   0 1 2 3 4 5 6 7 8 9
dnl 0
dnl 1
dnl 2
dnl 3
dnl 4         e
dnl 5       h x d
dnl 6     h x x x d
dnl 7   h x x x x x d
dnl 8 h x x x x x x x d
dnl 9 x x x x x x x x x d
dnl
dnl n = 12:
dnl    0 1 2 3 4 5 6 7 8 9 0 1
dnl  0
dnl  1
dnl  2
dnl  3
dnl  4
dnl  5           e
dnl  6         h x d
dnl  7       h x x x d
dnl  8     h x x x x x d
dnl  9   h x x x x x x x d
dnl 10 h x x x x x x x x x d
dnl 11 x x x x x x x x x x x d

	TEXT
	ALIGN(32)
PROLOGUE(_flint_mpn_sqrhigh_basecase_even)
	mov	R32(no2_param), R32(no2)
	lea	(ap_param,no2_param,8), ap	C ap += no2 where no2 = n / 2

	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	C Initial triangle
	C
	C       v--< (no2, no2 - 1) in matrix
	C     h x
	C   h x x x
	C h x x x x x
define(`sx', `rx')
define(`s0', `r0')
define(`s1', `r1')
define(`s2', `r2')
define(`s3', `mx')
define(`s4', `r3a')
define(`s5', `r3b')
define(`s6', `jmpreg')
define(`s7', `mx_final')
define(`s8', `ix_save')
	mov	0*8(ap), %rdx
	mulx	-2*8(ap), sx, sx
	mulx	-1*8(ap), s1, s0
	add	s1, sx
	adc	$0, s0
	C x, 0

	mov	1*8(ap), %rdx
	mulx	-3*8(ap), s3, s3
	mulx	-2*8(ap), s4, s5
	mulx	-1*8(ap), s6, s1
	mulx	0*8(ap), s7, s2
	add	s3, sx
	adc	s5, s0
	adc	$0, s1
	add	s4, sx
	adc	s6, s0
	adc	s7, s1
	adc	$0, s2
	C x, 0, 1, 2

	mov	2*8(ap), %rdx
	mulx	-4*8(ap), s3, s3
	mulx	-3*8(ap), s7, s4
	mulx	-2*8(ap), s8, s5
	add	s7, s3
	adc	s8, s4
	mulx	-1*8(ap), s7, s6
	adc	s7, s5
	adc	$0, s6
	xor	R32(s8), R32(s8)
	add	s3, sx
	adc	s4, s0
	adc	s5, s1
	adc	s6, s2
	mulx	0*8(ap), s5, s3
	mulx	1*8(ap), s6, s4
	adc	s6, s3
	adc	s8, s4
	adc	s5, s2
	adc	s8, s3
	adc	s8, s4
	C x, 0, 1, 2, 3, 4

	mov	s0, 0*8(rp)
	mov	s1, 1*8(rp)
	mov	s2, 2*8(rp)
	mov	s3, 3*8(rp)
	mov	s4, 4*8(rp)
undefine(`sx')
undefine(`s0')
undefine(`s1')
undefine(`s2')
undefine(`s3')
undefine(`s4')
undefine(`s5')
undefine(`s6')
undefine(`s7')
undefine(`s8')

	C Addmul chain for lower mid triangle into rp
	C
	C mx       = -7 * 8
	C mx_final = -(2 * no2 - 3) * 8
	C ix       = 0
	C ix_save  = 0			NOTE: It is already zero
	lea	(,no2,2), R32(mx_final)
	lea	-5*8(ap), ap
	lea	-3*8(,mx_final,8), R32(mx_final)
	mov	8*8(ap), %rdx
	mulx	0*8(ap), r1, r1
	mulx	1*8(ap), r2, r3a
	neg	mx_final
	xor	R32(ix), R32(ix)	C Also resets flags
	adcx	r1, rx
	movq	$-7*8, mx
	adox	r2, rx
	C jmp L(e4)

L(e4):	lea	-2*8(ap), ap
	lea	-4*8(rp), rp
	lea	L(e6)(%rip), jmpreg
	jmp	L(b4)

L(e6):	lea	-2*8(rp), rp
	mov	r3a, r3b
	lea	L(e0)(%rip), jmpreg
	lea	1(ix_save), ix_save
	jmp	L(b6)

	nop; nop; nop; nop		C Must for relative 8-bit jump to L(fnd)
L(fnd):	adox	ix, r2
	adox	ix, r3b
	adc	ix, r3b
	cmp	R32(mx), R32(mx_final)
	mov	r2, 1*8(rp)
	mov	r3b, 2*8(rp)
	lea	3*8(rp,mx), rp		C Reset rp
	jg	L(msh)
	mov	3*8(ap), %rdx		C Load second factor
	lea	(ap,mx), ap		C Reset ap
	je	L(kmp)
	mulx	0*8(ap), ix, ix
L(kmp):	mulx	1*8(ap), r2, r3a
	test	%al, %al
	adcx	ix, rx
	lea	-2*8(mx), mx
	mov	R32(ix_save), R32(ix)
	adox	r2, rx
	jmp	*jmpreg

L(e0):	lea	2*8(ap), ap
	lea	L(e2)(%rip), jmpreg
	C jmp	L(b0)

	ALIGN(32)
L(b0):	mulx	0*8(ap), r0, r1
	mulx	1*8(ap), r2, r3b
	adcx	r3a, r0
	adox	0*8(rp), r0
	adcx	r1, r2
	mov	r0, 0*8(rp)
	jrcxz	L(fnd)
	lea	-1(ix), ix
	adox	1*8(rp), r2
	mov	r2, 1*8(rp)
L(b6):	mulx	2*8(ap), r0, r1
	mulx	3*8(ap), r2, r3a
	adox	2*8(rp), r0
	adox	3*8(rp), r2
	adcx	r3b, r0
	adcx	r1, r2
	mov	r0, 2*8(rp)
	mov	r2, 3*8(rp)
L(b4):	mulx	4*8(ap), r0, r1
	mulx	5*8(ap), r2, r3b
	adox	4*8(rp), r0
	adox	5*8(rp), r2
	adcx	r3a, r0
	adcx	r1, r2
	mov	r0, 4*8(rp)
	mov	r2, 5*8(rp)
	lea	8*8(ap), ap
L(b2):	mulx	-2*8(ap), r0, r1
	mulx	-1*8(ap), r2, r3a
	adox	6*8(rp), r0
	adox	7*8(rp), r2
	adcx	r3b, r0
	lea	8*8(rp), rp
	adcx	r1, r2
	mov	r0, -2*8(rp)
	mov	r2, -1*8(rp)
	jmp	L(b0)

L(e2):	lea	4*8(ap), ap
	lea	-6*8(rp), rp
	mov	r3a, r3b
	lea	L(e4)(%rip), jmpreg
	jmp	L(b2)

define(`r3', `jmpreg')
define(`ld0', `r3a')
define(`ld1', `r3b')
define(`ld2', `mx')
define(`ld3', `ix_save')


	C Left shift rp by one and add diagonal
L(msh):	sar	$4, mx_final
	lea	L(etab)(%rip), ld2
	lea	1*8(ap,mx_final,8), ap	C Reset ap
	neg	R32(mx_final)
	mov	0*8(rp), ld0
	mov	1*8(rp), ld1
	or	R32(mx_final), R32(ix)
	shr	$2, R32(ix)
	and	$3, R32(mx_final)	C Also reset flags

	mov	0*8(ap), %rdx
	mulx	%rdx, r1, r1
ifdef(`PIC',
`	movslq	(ld2,mx_final,4), mx_final
	lea	(mx_final,ld2), ld2
',`')
	mov	1*8(ap), %rdx
	mulx	%rdx, r2, r3
	adox	rx, rx
	adox	ld0, ld0
	adcx	r1, rx
	adox	ld1, ld1
	adcx	r2, ld0
	adcx	r3, ld1
	mov	ld0, 0*8(rp)
ifdef(`PIC',
`	jmp	*ld2
',`
	jmp	*(ld2,mx_final,8)
')

L(ep3):	lea	-1*8(ap), ap
	lea	-2*8(rp), rp
	mov	ld1, ld3
	jmp	L(es3)

L(ep2):	lea	2*8(ap), ap
	lea	-4*8(rp), rp
	jmp	L(es2)

L(ep1):	lea	1*8(ap), ap
	lea	2*8(rp), rp
	mov	ld1, ld3
	jmp	L(es1)

	ALIGN(32)
L(ep0):	# Do nothing
L(es0):	mov	2*8(ap), %rdx
	mov	2*8(rp), ld2
	mov	3*8(rp), ld3
	mulx	%rdx, r2, r3
	lea	-1(ix), R32(ix)
	adox	ld2, ld2
	adox	ld3, ld3
	mov	ld1, 1*8(rp)
	adcx	r2, ld2
	adcx	r3, ld3
	mov	ld2, 2*8(rp)
L(es3):	mov	3*8(ap), %rdx
	mov	4*8(rp), ld0
	mov	5*8(rp), ld1
	mov	ld3, 3*8(rp)
	mulx	%rdx, r0, r1
	adox	ld0, ld0
	adox	ld1, ld1
	lea	4*8(ap), ap
	adcx	r0, ld0
	adcx	r1, ld1
	mov	ld0, 4*8(rp)
L(es2):	mov	0*8(ap), %rdx
	mov	6*8(rp), ld2
	mov	7*8(rp), ld3
	mov	ld1, 5*8(rp)
	mulx	%rdx, r2, r3
	adox	ld2, ld2
	adox	ld3, ld3
	lea	8*8(rp), rp
	adcx	r2, ld2
	adcx	r3, ld3
	mov	ld2, -2*8(rp)
L(es1):	mov	1*8(ap), %rdx
	mov	0*8(rp), ld0
	mov	1*8(rp), ld1
	mov	ld3, -1*8(rp)
	mulx	%rdx, r0, r1
	adox	ld0, ld0
	jrcxz	L(gin)
	adox	ld1, ld1
	adcx	r0, ld0
	adcx	r1, ld1
	mov	ld0, 0*8(rp)
	jmp	L(es0)

L(gin):	adox	ix, r1		C ix = 0
	adc	r0, ld0
	adc	ix, r1
	mov	ld0, 0*8(rp)
	mov	r1, 1*8(rp)

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
L(etab):JMPENT(	L(ep0), L(etab))
	JMPENT(	L(ep1), L(etab))
	JMPENT(	L(ep2), L(etab))
	JMPENT(	L(ep3), L(etab))
