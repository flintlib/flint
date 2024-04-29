dnl
dnl Copyright 2017 Free Software Foundation, Inc.
dnl Contributed to the GNU project by Torbjorn Granlund.
dnl Copyright (C) 2024 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

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

include(`config.m4')

define(`rp',	   `%rdi')
define(`ap',	   `%rsi')
define(`bp_param', `%rdx')
define(`n',	   `%rcx')

define(`bp',	`%r11')

define(`jmpreg',`%r15')
define(`nn',    `%rbp')
define(`mm',    `%rbx')

define(`s0', `%r8')
define(`s1', `%r9')
define(`s2', `%r10')
define(`s3', `%r12')
define(`s4', `%r13')
define(`s5', `%r14')

define(`sx', `%rax')

dnl Scheme:
dnl   0 1 2 3 4 5 6 7 8
dnl 0 x x x x x x x x x
dnl 1 x x x x x x x x l
dnl 2 x x x x x x x l
dnl 3 x x x x x x l
dnl 4 x x x x x l
dnl 5 x x x x l
dnl 6 x x x l
dnl 7 x x l
dnl 8 x l

dnl NOTE: Requires n > 8.

	TEXT
	ALIGN(32)
PROLOGUE(flint_mpn_mullow_basecase)
	push	mm
	push	nn
	push	s3
	push	s4
	push	s5
	push	jmpreg

	lea	-2(n), R32(nn)
	lea	1*8(bp_param), bp	C Prepare bp for addmul_1
	mov	0*8(bp_param), %rdx	C Load rdx for mul_1

	mov	R32(n), R32(sx)
	shr	$3, R32(n)
	and	$7, R32(sx)		C clear OF, CF as side-effect
	lea	L(mtab)(%rip), s0
ifdef(`PIC',
`	movslq	(s0,sx,4), sx
	lea	(sx,s0), s0
	jmp	*s0
',`
	jmp	*(s0,sx,8)
')

L(mf0):	mulx	0*8(ap), s0, s2
	lea	7*8(ap), ap
	lea	-1*8(rp), rp
	lea	L(f0)(%rip), jmpreg
	jmp	L(mb0)

L(mf3):	mulx	0*8(ap), s1, sx
	lea	2*8(ap), ap
	lea	2*8(rp), rp
	inc	R32(n)
	lea	L(f3)(%rip), jmpreg
	jmp	L(mb3)

L(mf4):	mulx	0*8(ap), s0, s2
	lea	3*8(ap), ap
	lea	3*8(rp), rp
	inc	R32(n)
	lea	L(f4)(%rip), jmpreg
	jmp	L(mb4)

L(mf5):	mulx	0*8(ap), s1, sx
	lea	4*8(ap), ap
	lea	4*8(rp), rp
	inc	R32(n)
	lea	L(f5)(%rip), jmpreg
	jmp	L(mb5)

L(mf6):	mulx	0*8(ap), s0, s2
	lea	5*8(ap), ap
	lea	5*8(rp), rp
	inc	R32(n)
	lea	L(f6)(%rip), jmpreg
	jmp	L(mb6)

L(mf7):	mulx	0*8(ap), s1, sx
	lea	6*8(ap), ap
	lea	6*8(rp), rp
	inc	R32(n)
	lea	L(f7)(%rip), jmpreg
	jmp	L(mb7)

L(mf1):	mulx	0*8(ap), s1, sx
	lea	L(f1)(%rip), jmpreg
	jmp	L(mb1)

L(mf2):	mulx	0*8(ap), s0, s2
	lea	1*8(ap), ap
	lea	1*8(rp), rp
	lea	L(f2)(%rip), jmpreg
	mulx	0*8(ap), s1, sx

	ALIGN(32)
L(mtop):mov	s0, -1*8(rp)
	adc	s2, s1
L(mb1):	mulx	1*8(ap), s0, s2
	adc	sx, s0
	lea	8*8(ap), ap
	mov	s1, 0*8(rp)
L(mb0):	mov	s0, 1*8(rp)
	mulx	-6*8(ap), s1, sx
	lea	8*8(rp), rp
	adc	s2, s1
L(mb7):	mulx	-5*8(ap), s0, s2
	mov	s1, -6*8(rp)
	adc	sx, s0
L(mb6):	mov	s0, -5*8(rp)
	mulx	-4*8(ap), s1, sx
	adc	s2, s1
L(mb5):	mulx	-3*8(ap), s0, s2
	mov	s1, -4*8(rp)
	adc	sx, s0
L(mb4):	mulx	-2*8(ap), s1, sx
	mov	s0, -3*8(rp)
	adc	s2, s1
L(mb3):	mulx	-1*8(ap), s0, s2
	adc	sx, s0
	mov	s1, -2*8(rp)
	dec	R32(n)
	mulx	0*8(ap), s1, sx
	jnz	L(mtop)

	lea	1*8(,nn,8), R32(mm)
	mov	s0, -1*8(rp)
	adc	s2, s1
	mov	0*8(bp), %rdx
	mov	s1, 0*8(rp)
	adc	n, sx			C n = 0
	shr	$3, R32(nn)
	neg	mm
	or	R32(nn), R32(n)		C Reset n, clear OF and CF
	lea	1*8(rp,mm), rp		C Reset rp
	lea	(ap,mm), ap		C Reset ap
	jmp	*jmpreg

	nop;nop;nop;nop;nop;nop;nop;nop;nop;nop;nop;nop;nop;nop;nop;nop;nop;
L(f7):	mulx	0*8(ap), s0, s2
	lea	5*8(ap), ap
	lea	-3*8(rp), rp
	lea	L(f6)(%rip), jmpreg
	jmp	L(b7)

L(f6):	mulx	0*8(ap), s1, s3
	lea	4*8(ap), ap
	lea	-4*8(rp), rp
	lea	L(f5)(%rip), jmpreg
	jmp	L(b6)

L(end):	adox	0*8(rp), s1
	mulx	1*8(ap), s0, s2		C Only s0 is used
	lea	1*8(mm), mm
	adox	s3, sx
	mov	s1, 0*8(rp)
	lea	1*8(bp), bp
	adc	s0, sx
	lea	(ap,mm), ap		C Reset ap
	mov	0*8(bp), %rdx
	or	R32(nn), R32(n)		C Reset count, clear CF and OF
	lea	1*8(rp,mm), rp		C Reset rp
	jmp	*jmpreg

L(f0):	mulx	0*8(ap), s1, s3
	lea	-2*8(ap), ap
	lea	-2*8(rp), rp
	lea	L(f7)(%rip), jmpreg
	jmp	L(b0)

L(f3):	mulx	0*8(ap), s0, s2
	lea	1*8(ap), ap
	lea	1*8(rp), rp
	mulx	0*8(ap), s1, s3
	lea	L(f2)(%rip), jmpreg

	ALIGN(32)
L(top):	adox	-1*8(rp), s0
	adcx	s2, s1
	mov	s0, -1*8(rp)
	jrcxz	L(end)
L(b2):	mulx	1*8(ap), s0, s2
	adox	0*8(rp), s1
	lea	-1(n), R32(n)
	mov	s1, 0*8(rp)
	adcx	s3, s0
L(b1):	mulx	2*8(ap), s1, s3
	adcx	s2, s1
	adox	1*8(rp), s0
	mov	s0, 1*8(rp)
L(b0):	mulx	3*8(ap), s0, s2
	lea	8*8(ap), ap
	adcx	s3, s0
	adox	2*8(rp), s1
	mov	s1, 2*8(rp)
L(b7):	mulx	-4*8(ap), s1, s3
	adox	3*8(rp), s0
	adcx	s2, s1
	mov	s0, 3*8(rp)
L(b6):	mulx	-3*8(ap), s0, s2
	adcx	s3, s0
	adox	4*8(rp), s1
	mov	s1, 4*8(rp)
L(b5):	mulx	-2*8(ap), s1, s3
	adox	5*8(rp), s0
	adcx	s2, s1
	mov	s0, 5*8(rp)
L(b4):	adox	6*8(rp), s1
	mulx	-1*8(ap), s0, s2
	mov	s1, 6*8(rp)
	lea	8*8(rp), rp
	adcx	s3, s0
	mulx	0*8(ap), s1, s3
	jmp	L(top)

L(f5):	mulx	0*8(ap), s0, s2
	lea	3*8(ap), ap
	lea	-5*8(rp), rp
	lea	L(f4)(%rip), jmpreg
	jmp	L(b5)

L(f4):	mulx	0*8(ap), s1, s3
	lea	2*8(ap), ap
	lea	-6*8(rp), rp
	lea	L(f3)(%rip), jmpreg
	jmp	L(b4)

L(f2):	mulx	0*8(ap), s1, s3
	lea	-1(nn), R32(nn)
	lea	L(f1)(%rip), jmpreg
	jmp	L(b2)

L(f1):	mulx	0*8(ap), s0, s2
	jrcxz	L(cor)
	lea	-1*8(ap), ap
	lea	-1*8(rp), rp
	lea	L(f0)(%rip), jmpreg
	jmp	L(b1)

define(`t0', `s0')
define(`t1', `s2')
define(`t2', `s1')
define(`t3', `s3')
define(`t4', `s4')
define(`t5', `s5')
define(`t6', `jmpreg')
define(`t7', `mm')
define(`t8', `nn')
define(`t9', `n')
define(`tx', `sx')
L(cor):	mulx	1*8(ap), t2, t3
	mulx	2*8(ap), t4, t5
	adcx	0*8(rp), t0
	adox	t1, t2
	mulx	3*8(ap), t6, t7
	adcx	1*8(rp), t2
	adox	t3, t4
	mulx	4*8(ap), t8, t9
	adcx	2*8(rp), t4
	adox	t5, t6
	mulx	5*8(ap), t1, t3
	adcx	3*8(rp), t6
	adox	t7, t8
	mulx	6*8(ap), t5, t7
	mov	t0, 0*8(rp)
	adcx	4*8(rp), t8
	adox	t9, t1
	mulx	7*8(ap), t0, t9
	adcx	5*8(rp), t1
	adox	t3, t5
	mulx	8*8(ap), t3, %rdx	C %rdx unused
	adcx	6*8(rp), t5
	adox	t7, t0
	adcx	7*8(rp), t0
	adox	t9, tx
	C 2, 4, 6, 8, 1, 5, 0, x

	mov	1*8(bp), %rdx
	mulx	0*8(ap), t7, t9
	adc	t3, tx
	test	%al, %al		C Reset OF and CF
	adcx	t7, t2
	adox	t9, t4
	mov	t2, 1*8(rp)
	mulx	1*8(ap), t7, t9
	mulx	2*8(ap), t2, t3
	adcx	t7, t4
	adox	t9, t6
	mulx	3*8(ap), t7, t9
	adcx	t2, t6
	adox	t3, t8
	mulx	4*8(ap), t2, t3
	adcx	t7, t8
	adox	t9, t1
	mulx	5*8(ap), t7, t9
	adcx	t2, t1
	adox	t3, t5
	mulx	6*8(ap), t2, t3
	adcx	t7, t5
	adox	t9, t0
	mulx	7*8(ap), t7, t9		C t9 unused
	adcx	t2, t0
	adox	t3, tx
	C 4, 6, 8, 1, 5, 0, x

	mov	2*8(bp), %rdx
	mulx	0*8(ap), t2, t3
	adc	t7, tx
	test	%al, %al
	mulx	1*8(ap), t7, t9
	adcx	t2, t4
	adox	t3, t6
	mov	t4, 2*8(rp)
	mulx	2*8(ap), t2, t3
	adcx	t7, t6
	adox	t9, t8
	mulx	3*8(ap), t7, t9
	adcx	t2, t8
	adox	t3, t1
	mulx	4*8(ap), t2, t3
	adcx	t7, t1
	adox	t9, t5
	mulx	5*8(ap), t7, t9
	adcx	t2, t5
	adox	t3, t0
	mulx	6*8(ap), t2, t3		C t3 unused
	adcx	t7, t0
	adox	t9, tx
	C 6, 8, 1, 5, 0, x

	mov	3*8(bp), %rdx
	mulx	0*8(ap), t4, t7
	adc	t2, tx
	test	%al, %al
	mulx	1*8(ap), t2, t3
	adcx	t4, t6
	adox	t7, t8
	mov	t6, 3*8(rp)
	mulx	2*8(ap), t4, t7
	mulx	3*8(ap), t6, t9
	adcx	t2, t8
	adox	t3, t1
	mulx	4*8(ap), t2, t3
	adcx	t4, t1
	adox	t7, t5
	mulx	5*8(ap), t4, t7		C t7 unused
	adcx	t6, t5
	adox	t9, t0
	adcx	t2, t0
	adox	t3, tx
	C 8, 1, 5, 0, x

	mov	4*8(bp), %rdx
	mulx	0*8(ap), t2, t3
	adc	t4, tx
	test	%al, %al
	mulx	1*8(ap), t4, t6
	mulx	2*8(ap), t7, t9
	adcx	t2, t8
	adox	t3, t1
	mulx	3*8(ap), t2, t3
	adcx	t4, t1
	adox	t6, t5
	mulx	4*8(ap), t4, t6		C t6 unused
	adcx	t7, t5
	adox	t9, t0
	mov	t8, 4*8(rp)
	adox	t3, tx
	adc	t2, t0
	C 1, 5, 0, x

	mov	5*8(bp), %rdx
	mulx	0*8(ap), t2, t3
	adc	t4, tx
	mulx	1*8(ap), t4, t6
	mulx	2*8(ap), t7, t8
	imul	3*8(ap), %rdx
	add	t2, t1
	adc	t3, t5
	mov	t1, 5*8(rp)
	adc	t6, t0
	adc	t8, tx
	add	t4, t5
	adc	t7, t0
	adc	%rdx, tx
	C 5, 0, x

	mov	6*8(bp), %rdx
	mulx	0*8(ap), t1, t2
	mulx	1*8(ap), t3, t4
	imul	2*8(ap), %rdx
	add	t1, t5
	adc	t2, t0
	mov	t5, 6*8(rp)
	adc	t4, tx
	add	t3, t0
	adc	%rdx, tx
	C 0, x

	mov	7*8(bp), %rdx
	mulx	0*8(ap), t1, t2
	imul	1*8(ap), %rdx

	pop	jmpreg
	pop	s5
	pop	s4
	pop	s3
	pop	nn
	pop	mm

	add	t1, t0
	adc	t2, tx
	mov	t0, 7*8(rp)
	add	%rdx, tx

	ret
EPILOGUE()
	JUMPTABSECT
	ALIGN(8)
L(mtab):JMPENT(	L(mf0), L(mtab))
	JMPENT(	L(mf1), L(mtab))
	JMPENT(	L(mf2), L(mtab))
	JMPENT(	L(mf3), L(mtab))
	JMPENT(	L(mf4), L(mtab))
	JMPENT(	L(mf5), L(mtab))
	JMPENT(	L(mf6), L(mtab))
	JMPENT(	L(mf7), L(mtab))
