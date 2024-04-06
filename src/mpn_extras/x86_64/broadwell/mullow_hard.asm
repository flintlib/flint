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

define(`rp',       `%rdi')
define(`ap',       `%rsi')
define(`bp_param', `%rdx')

define(`bp',	   `%r8')

define(`s0', `%rcx')
define(`s1', `%r9')
define(`s2', `%r10')
define(`s3', `%r11')
define(`s4', `%rbx')
define(`s5', `%rbp')
define(`s6', `%r12')
define(`s7', `%r13')
define(`s8', `%r14')
define(`s9', `%r15')

define(`sx', `%rax')

	TEXT

	ALIGN(16)
PROLOGUE(flint_mpn_mullow_1)
	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), s0, sx		C a0 b0
	mov	s0, 0*8(rp)
	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mullow_2)
	mov	1*8(bp_param), sx
	mov	0*8(bp_param), %rdx

	mulx	0*8(ap), s0, s1		C a0 b0
	mulx	1*8(ap), s2, s3		C a1 b0
	mov	sx, %rdx
	mulx	0*8(ap), bp_param, bp	C a0 b1
	imul	1*8(ap), sx		C L(a1 b1)
	C s0, (s1, s2, bp_param), (s3, bp, sx)

	mov	s0, 0*8(rp)

	add	s2, s1
	adc	s3, bp
	C (s1, bp_param), (bp, sx)

	add	bp_param, s1
	adc	bp, sx
	C s1, sx

	mov	s1, 1*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mullow_3)
	push	s4
	push	s5

	mov	bp_param, bp

	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b0
	mulx	1*8(ap), s2, s3		C   a1 b0
	mulx	2*8(ap), s4, sx		C   a2 b0
	C 0, (1, 2), (3, 4), x

	mov	s0, 0*8(rp)
	add	s2, s1

	mov	1*8(bp), %rdx
	mulx	0*8(ap), s0, s2		C   a0 b1
	adc	s4, s3
	adc	$0, sx
	mulx	1*8(ap), s4, s5		C   a1 b1
	imul	2*8(ap), %rdx		C L(a2 b1)
	C 0, (2, 4), (5, rdx)

	add	s0, s1
	adc	s2, s3
	mov	s1, 1*8(rp)
	adc	%rdx, sx

	mov	2*8(bp), %rdx
	mulx	0*8(ap), s0, s2		C   a0 b2
	imul	1*8(ap), %rdx		C L(a1 b2)
	C 0, (2, rdx)

	add	s4, s3
	adc	s5, sx

	pop	s5
	pop	s4

	add	s2, sx
	add	s0, s3
	adc	%rdx, sx
	mov	s3, 2*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mullow_4)
	mov	bp_param, bp

	push	s4
	push	s5
	push	s6
	push	s7
	push	s8

	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b0
	mulx	1*8(ap), s2, s3		C   a1 b0
	mulx	2*8(ap), s4, s5		C   a2 b0
	mulx	3*8(ap), s6, sx		C   a3 b0
	C 0, (1, 2), (3, 4), (5, 6), x

	mov	s0, 0*8(rp)
	add	s2, s1
	adc	s4, s3
	adc	s6, s5
	adc	$0, sx
	C 1, 3, 5, x

	mov	1*8(bp), %rdx
	mulx	0*8(ap), s0, s2		C   a0 b1
	mulx	1*8(ap), s4, s6		C   a1 b1
	mulx	2*8(ap), s7, s8		C   a2 b1
	imul	3*8(ap), %rdx		C L(a3 b1)
	C 0, (2, 4), (6, 7), (8, rdx)

	add	s0, s1
	adc	s2, s3
	adc	s6, s5
	adc	s8, sx
	mov	s1, 1*8(rp)
	add	s4, s3
	adc	s7, s5
	adc	%rdx, sx
	C 3, 5, x

	mov	2*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b2
	mulx	1*8(ap), s2, s4		C   a1 b2
	imul	2*8(ap), %rdx		C L(a2 b2)
	C 0, (1, 2), (4, rdx)

	add	s0, s3
	adc	s1, s5
	mov	s3, 2*8(rp)
	adc	s4, sx
	add	s5, s2
	adc	%rdx, sx
	C 2, x

	mov	3*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b3
	imul	1*8(ap), %rdx		C L(a1 b3)
	C 0, (1, rdx)

	pop	s8
	pop	s7
	pop	s6
	pop	s5
	pop	s4

	add	s0, s2
	adc	s1, sx
	mov	s2, 3*8(rp)
	add	%rdx, sx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mullow_5)
	mov	bp_param, bp

	push	s4
	push	s5
	push	s6
	push	s7
	push	s8
	push	s9

	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b0
	mulx	1*8(ap), s2, s3		C   a1 b0
	mulx	2*8(ap), s4, s5		C   a2 b0
	mulx	3*8(ap), s6, s7		C   a3 b0
	mulx	4*8(ap), s8, sx		C   a4 b0
	C 0, (1, 2), (3, 4), (5, 6), (7, 8), x

	mov	s0, 0*8(rp)
	add	s2, s1
	adc	s4, s3
	adc	s6, s5
	adc	s8, s7
	adc	$0, sx
	C 1, 3, 5, 7, x

	mov	1*8(bp), %rdx
	mulx	0*8(ap), s0, s2		C   a0 b1
	mulx	1*8(ap), s4, s6		C   a1 b1
	mulx	2*8(ap), s8, s9		C   a2 b1
	add	s0, s1
	adc	s2, s3
	mov	s1, 1*8(rp)
	mulx	3*8(ap), s0, s2		C   a3 b1
	mulx	4*8(ap), %rdx, s1	C   a4 b1, but we only use %rdx
	C -, (-, 4), (6, 8), (9, 0), (2, %rdx)

	adc	s6, s5
	adc	s9, s7
	adc	s2, sx
	add	s4, s3
	adc	s8, s5
	adc	s0, s7
	adc	%rdx, sx
	C 3, 5, 7, x

	mov	2*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b2
	mulx	1*8(ap), s2, s4		C   a1 b2
	mulx	2*8(ap), s6, s8		C   a2 b2
	imul	3*8(ap), %rdx		C L(a3 b2)
	C 0, (1, 2), (4, 6), (8, rdx)

	add	s0, s3
	adc	s1, s5
	adc	s4, s7
	adc	s8, sx
	add	s2, s5
	adc	s6, s7
	adc	%rdx, sx
	mov	s3, 2*8(rp)
	C 5, 7, x

	mov	3*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b3
	mulx	1*8(ap), s2, s3		C   a1 b3
	imul	2*8(ap), %rdx		C L(a2 b3)
	C 0, (1, 2), (3, rdx)

	add	s0, s5
	adc	s1, s7
	mov	s5, 3*8(rp)
	adc	s3, sx
	add	s7, s2
	adc	%rdx, sx
	C 2, x

	mov	4*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b4
	imul	1*8(ap), %rdx		C L(a1 b4)
	C 0, (1, rdx)

	pop	s9
	pop	s8
	pop	s7
	pop	s6
	pop	s5
	pop	s4

	add	s0, s2
	adc	s1, sx
	mov	s2, 4*8(rp)
	add	%rdx, sx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mullow_6)
	mov	bp_param, bp

	push	s4
	push	s5
	push	s6
	push	s7
	push	s8
	push	s9

	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b0
	mulx	1*8(ap), s2, s3		C   a1 b0
	mulx	2*8(ap), s4, s5		C   a2 b0
	mulx	3*8(ap), s6, s7		C   a3 b0
	mulx	4*8(ap), s8, s9		C   a4 b0
	mulx	5*8(ap), %rdx, sx	C   a5 b0
	C 0, (1, 2), (3, 4), (5, 6), (7, 8), (9, rdx), x

	add	s2, s1
	adc	s4, s3
	adc	s6, s5
	mov	s0, 0*8(rp)
	adc	s8, s7
	adc	%rdx, s9
	adc	$0, sx
	C 1, 3, 5, 7, 9, x

	test	%al, %al
	mov	1*8(bp), %rdx
	mulx	0*8(ap), s2, s4		C   a0 b1
	mulx	1*8(ap), s6, s8		C   a1 b1
	adcx	s2, s1
	adox	s4, s3
	mov	s1, 1*8(rp)
	mulx	2*8(ap), s2, s4		C   a2 b1
	adcx	s6, s3
	adox	s8, s5
	mulx	3*8(ap), s6, s8		C   a3 b1
	adcx	s2, s5
	adox	s4, s7
	mulx	4*8(ap), s2, s4		C   a4 b1
	adcx	s6, s7
	adox	s8, s9
	mulx	5*8(ap), %rdx, s6	C   a5 b1, but we only use %rdx
	C -, (-, -), (-, -), (-, -), (-, 2), (4, %rdx)

	adcx	s2, s9
	adox	s4, sx
	adc	%rdx, sx
	C 3, 5, 7, 9, x

	mov	2*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b2
	mulx	1*8(ap), s2, s4		C   a1 b2
	mulx	2*8(ap), s6, s8		C   a2 b2
	add	s0, s3
	adc	s1, s5
	mov	s3, 2*8(rp)
	mulx	3*8(ap), s0, s1		C   a3 b2
	mulx	4*8(ap), %rdx, s3	C   a4 b2, but we only use %rdx
	C -, (-, 2), (4, 6), (8, 0), (1, rdx)

	adc	s4, s7
	adc	s8, s9
	adc	s1, sx
	add	s2, s5
	adc	s6, s7
	adc	s0, s9
	adc	%rdx, sx
	C 5, 7, 9, x

	mov	3*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b3
	mulx	1*8(ap), s2, s3		C   a1 b3
	mulx	2*8(ap), s4, s6		C   a2 b3
	imul	3*8(ap), %rdx		C L(a3 b3)
	C 0, (1, 2), (3, 4), (6, rdx)

	add	s0, s5
	adc	s1, s7
	adc	s3, s9
	mov	s5, 3*8(rp)
	adc	s6, sx
	add	s2, s7
	adc	s4, s9
	adc	%rdx, sx
	C 7, 9, x

	mov	4*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b4
	mulx	1*8(ap), s2, s3		C   a1 b4
	imul	2*8(ap), %rdx		C L(a2 b4)
	C 0, (1, 2), (3, rdx)

	add	s0, s7
	adc	s1, s9
	adc	s3, sx
	mov	s7, 4*8(rp)
	add	s9, s2
	adc	%rdx, sx
	C 2, x

	mov	5*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b5
	imul	1*8(ap), %rdx		C L(a1 b5)
	C 0, (1, rdx)

	pop	s9
	pop	s8
	pop	s7
	pop	s6
	pop	s5
	pop	s4

	add	s0, s2
	adc	s1, sx
	mov	s2, 5*8(rp)
	add	%rdx, sx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mullow_7)
	mov	bp_param, bp

	push	s4
	push	s5
	push	s6
	push	s7
	push	s8
	push	s9

	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b0
	mulx	1*8(ap), s2, s3		C   a1 b0
	mulx	2*8(ap), s4, s5		C   a2 b0
	mulx	3*8(ap), s6, s7		C   a3 b0
	mov	s0, 0*8(rp)
	mulx	4*8(ap), s8, s9		C   a4 b0
	add	s2, s1
	adc	s4, s3
	mulx	5*8(ap), s0, s2		C   a5 b0
	adc	s6, s5
	adc	s8, s7
	mulx	6*8(ap), s4, sx		C   a6 b0
	C 0, (1, -), (3, -), (5, -), (7, -), (9, 0), (2, 4), x

	adc	s0, s9
	adc	s4, s2
	mov	1*8(bp), %rdx
	adc	$0, sx
	C 1, 3, 5, 7, 9, 2, x

	test	%al, %al
	mulx	0*8(ap), s0, s4		C   a0 b1
	mulx	1*8(ap), s6, s8		C   a1 b1
	adcx	s0, s1
	adox	s4, s3
	mulx	2*8(ap), s0, s4		C   a2 b1
	adcx	s6, s3
	adox	s8, s5
	mulx	3*8(ap), s6, s8		C   a3 b1
	adcx	s0, s5
	adox	s4, s7
	mulx	4*8(ap), s0, s4		C   a4 b1
	adcx	s6, s7
	adox	s8, s9
	mulx	5*8(ap), s6, s8		C   a5 b1
	adcx	s0, s9
	adox	s4, s2
	mulx	6*8(ap), %rdx, s0	C   a6 b1, but we only use %rdx
	adcx	s6, s2
	adox	s8, sx
	mov	s1, 1*8(rp)
	adc	%rdx, sx
	C 3, 5, 7, 9, 2, x

	mov	2*8(bp), %rdx
	test	%al, %al
	mulx	0*8(ap), s0, s1		C   a0 b2
	mulx	1*8(ap), s4, s6		C   a1 b2
	adcx	s0, s3
	adox	s1, s5
	mov	s3, 2*8(rp)
	mulx	2*8(ap), s0, s1		C   a2 b2
	adcx	s4, s5
	adox	s6, s7
	mulx	3*8(ap), s4, s6		C   a3 b2
	adcx	s0, s7
	adox	s1, s9
	mulx	4*8(ap), s0, s1		C   a4 b2
	adcx	s4, s9
	adox	s6, s2
	mulx	5*8(ap), %rdx, s3	C   a5 b2, but we only use %rdx
	adcx	s0, s2
	adox	s1, sx
	adc	%rdx, sx
	C 5, 7, 9, 2, x

	mov	3*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b3
	mulx	1*8(ap), s3, s4		C   a1 b3
	mulx	2*8(ap), s6, s8		C   a2 b3
	add	s0, s5
	adc	s3, s1
	adc	$0, s4
	mulx	3*8(ap), s0, s3		C   a3 b3
	imul	4*8(ap), %rdx		C L(a4 b3)
	C -, (1, -), (4, 6), (8, 0), (3, rdx)

	mov	s5, 3*8(rp)
	add	s1, s7
	adc	s4, s9
	adc	s8, s2
	adc	s3, sx
	add	s6, s9
	adc	s0, s2
	adc	%rdx, sx
	C 7, 9, 2, x

	mov	4*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b4
	mulx	1*8(ap), s3, s4		C   a1 b4
	mulx	2*8(ap), s5, s6		C   a2 b4
	imul	3*8(ap), %rdx		C L(a3 b4)
	C 0, (1, 3), (4, 5), (6, rdx)

	add	s0, s7
	adc	s1, s9
	adc	s4, s2
	mov	s7, 4*8(rp)
	adc	s6, sx
	add	s9, s3
	adc	s5, s2
	adc	%rdx, sx
	C 3, 2, x

	mov	5*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b5
	mulx	1*8(ap), s4, s5		C   a1 b5
	imul	2*8(ap), %rdx		C L(a2 b5)
	C 0, (1, 4), (5, rdx)

	pop	s9
	pop	s8

	add	s0, s3
	adc	s1, s2
	adc	s5, sx
	mov	s3, 5*8(rp)
	add	s4, s2
	adc	%rdx, sx
	C 2, x

	mov	6*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b6
	imul	1*8(ap), %rdx		C L(a1 b6)
	C 0, (1, rdx)

	pop	s7
	pop	s6
	pop	s5
	pop	s4

	add	s0, s2
	adc	s1, sx
	mov	s2, 6*8(rp)
	add	%rdx, sx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mullow_8)
	mov	bp_param, bp

	push	s4
	push	s5
	push	s6
	push	s7
	push	s8
	push	s9

	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b0
	mulx	1*8(ap), s2, s3		C   a1 b0
	mulx	2*8(ap), s4, s5		C   a2 b0
	mulx	3*8(ap), s6, s7		C   a3 b0
	mov	s0, 0*8(rp)
	mulx	4*8(ap), s8, s9		C   a4 b0
	add	s2, s1
	adc	s4, s3
	mulx	5*8(ap), s0, s2		C   a5 b0
	adc	s6, s5
	adc	s8, s7
	mulx	6*8(ap), s4, s6		C   a6 b0
	mulx	7*8(ap), s8, sx		C   a7 b0
	C -, (1, -), (3, -), (5, -), (7, -), (9, 0), (2, 4), (6, 8), x

	adc	s0, s9
	adc	s4, s2
	adc	s8, s6
	mov	1*8(bp), %rdx
	adc	$0, sx
	C 1, 3, 5, 7, 9, 2, 6, x

	test	%al, %al
	mulx	0*8(ap), s0, s4		C   a0 b1
	adcx	s0, s1
	adox	s4, s3
	mov	s1, 1*8(rp)
	mulx	1*8(ap), s8, s0		C   a1 b1
	mulx	2*8(ap), s4, s1		C   a2 b1
	adcx	s8, s3
	adox	s0, s5
	mulx	3*8(ap), s8, s0		C   a3 b1
	adcx	s4, s5
	adox	s1, s7
	mulx	4*8(ap), s4, s1		C   a4 b1
	adcx	s8, s7
	adox	s0, s9
	mulx	5*8(ap), s8, s0		C   a5 b1
	adcx	s4, s9
	adox	s1, s2
	mulx	6*8(ap), s4, s1		C   a6 b1
	adcx	s8, s2
	adox	s0, s6
	mulx	7*8(ap), s8, s0		C   a7 b1, but we only use s8
	adcx	s4, s6
	adox	s1, sx
	adc	s8, sx
	C 3, 5, 7, 9, 2, 6, x

	mov	2*8(bp), %rdx
	test	%al, %al
	mulx	0*8(ap), s0, s1		C   a0 b2
	mulx	1*8(ap), s4, s8		C   a1 b2
	adcx	s0, s3
	adox	s1, s5
	mov	s3, 2*8(rp)
	mulx	2*8(ap), s0, s1		C   a2 b2
	adcx	s4, s5
	adox	s8, s7
	mulx	3*8(ap), s4, s8		C   a3 b2
	adcx	s0, s7
	adox	s1, s9
	mulx	4*8(ap), s0, s1		C   a4 b2
	adcx	s4, s9
	adox	s8, s2
	mulx	5*8(ap), s4, s8		C   a5 b2
	adcx	s0, s2
	adox	s1, s6
	mulx	6*8(ap), s0, s1		C   a6 b2, but we only use s0
	adcx	s4, s6
	adox	s8, sx
	adc	s0, sx
	C 5, 7, 9, 2, 6, x

	mov	3*8(bp), %rdx
	test	%al, %al
	mulx	0*8(ap), s0, s1		C   a0 b3
	mulx	1*8(ap), s3, s4		C   a1 b3
	adcx	s0, s5
	adox	s1, s7
	mulx	2*8(ap), s0, s1		C   a2 b3
	adcx	s3, s7
	adox	s4, s9
	mulx	3*8(ap), s3, s4		C   a3 b3
	adcx	s0, s9
	adox	s1, s2
	mulx	4*8(ap), s0, s1		C   a4 b3
	adcx	s3, s2
	adox	s4, s6
	mulx	5*8(ap), s3, s4		C   a5 b3, but we only use s3
	mov	s5, 3*8(rp)
	adcx	s0, s6
	adox	s1, sx
	adc	s3, sx
	C 7, 9, 2, 6, x

	mov	4*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b4
	mulx	1*8(ap), s3, s4		C   a1 b4
	mulx	2*8(ap), s5, s8		C   a2 b4
	add	s0, s7
	adc	$0, s1
	mov	s7, 4*8(rp)
	mulx	3*8(ap), s0, s7		C   a3 b4
	imul	4*8(ap), %rdx		C L(a4 b4)
	C -, (1, 3), (4, 5), (8, 0), (7, rdx)

	add	s1, s9
	adc	s4, s2
	adc	s8, s6
	adc	s7, sx
	add	s3, s9
	adc	s5, s2
	adc	s0, s6
	adc	%rdx, sx
	C 9, 2, 6, x

	mov	5*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b5
	mulx	1*8(ap), s5, s4		C   a1 b5
	mulx	2*8(ap), s3, s7		C   a2 b5
	imul	3*8(ap), %rdx		C L(a3 b5)
	C 0, (1, 5), (4, 3), (7, rdx)

	add	s0, s9
	adc	s1, s2
	adc	s4, s6
	adc	s7, sx
	mov	s9, 5*8(rp)
	add	s5, s2
	adc	s6, s3
	adc	%rdx, sx
	C 2, 3, x

	mov	6*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b6
	mulx	1*8(ap), s4, s5		C   a1 b6
	imul	2*8(ap), %rdx		C L(a2 b6)
	C 0, (1, 4), (5, rdx)

	pop	s9
	pop	s8

	add	s0, s2
	adc	s1, s3
	adc	s5, sx
	mov	s2, 6*8(rp)
	add	s4, s3
	adc	%rdx, sx
	C 3, x

	mov	7*8(bp), %rdx
	mulx	0*8(ap), s0, s1		C   a0 b7
	imul	1*8(ap), %rdx		C L(a1 b7)
	C 0, (1, rdx)

	pop	s7
	pop	s6
	pop	s5
	pop	s4

	add	s0, s3
	adc	s1, sx
	mov	s3, 7*8(rp)
	add	%rdx, sx

	ret
EPILOGUE()
