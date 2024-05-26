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

dnl TODO:
dnl  * Redo n = 6, 7, 8, just like n < 6.

define(`rp',	   `%rdi')
define(`ap',	   `%rsi')
define(`bp_param', `%rdx')

define(`bp',	   `%r8')

define(`s0',  `%rax')
define(`s1',  `%rcx')
define(`s2',  `%r9')
define(`s3',  `%r10')
define(`s4',  `%r11')
define(`s5',  `%rbx')
define(`s6',  `%rbp')
define(`s7',  `%r12')
define(`s8',  `%r13')
define(`s9',  `%r14')
define(`s10', `%r15')

	TEXT
	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_1)
	C   0
	C 0 x
	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), s0, s1
	xor	%edx, %edx
	test	s1, s1
	js	L(1)
	add	s0, s0
	adc	s1, s1
	inc	%edx
L(1):	mov	s1, 0*8(rp)
	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_2)
	C   0 1
	C 0 h x
	C 1 x x
	mov	1*8(bp_param), s2
	mov	0*8(bp_param), %rdx

	mulx	0*8(ap), s0, s0		C a0 b0
	xor	R32(s1), R32(s1)
	mulx	1*8(ap), s3, s4		C a1 b0
	adox	s3, s0
	C 0, 4 + o

	mov	s2, %rdx
	mulx	0*8(ap), bp, s3		C a0 b1
	adcx	bp, s0
	adox	s3, s4
	mulx	1*8(ap), s2, ap		C a1 b1
	adcx	s2, s4
	adox	s1, ap
	C 0, 4, ap + c

	mov	$0, %edx
	adc	s1, ap
	js	L(2)
	add	s0, s0
	adc	s4, s4
	inc	%edx
	adc	ap, ap
L(2):	mov	s4, 0*8(rp)
	mov	ap, 1*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_3)
	C   0 1 2
	C 0   h x
	C 1 h x x
	C 2 x x x
	mov	bp_param, bp
	mov	0*8(bp_param), %rdx

	mulx	1*8(ap), s0, s0		C a1 b0
	push	s5
	xor	R32(s1), R32(s1)
	mulx	2*8(ap), s2, s3		C a2 b0
	push	s6
	adcx	s2, s0
	C 0, 3 + c

	mov	1*8(bp), %rdx
	mulx	1*8(ap), s4, s5		C a1 b1
	adox	s4, s0
	adcx	s5, s3
	mulx	2*8(ap), s6, s2		C a2 b1
	adox	s6, s3
	adcx	s1, s2
	mulx	0*8(ap), s4, s4		C a0 b1
	adox	s1, s2
	C (0, 4), 3, 2

	mov	2*8(bp), %rdx
	mulx	0*8(ap), bp, s5		C a0 b2
	adox	s4, s0
	mulx	1*8(ap), s6, s4		C a1 b2
	adcx	bp, s0
	adox	s5, s3
	mulx	2*8(ap), bp, ap		C a2 b2
	adcx	s6, s3
	adox	s4, s2
	pop	s6
	adcx	bp, s2
	adox	s1, ap
	C 0, 3, 2, ap + c

	mov	$0, %edx
	pop	s5
	adc	s1, ap
	C 0, 3, 2, ap

	js	L(3)
	add	s0, s0
	adc	s3, s3
	inc	%edx
	adc	s2, s2
	adc	ap, ap
L(3):	mov	s3, 0*8(rp)
	mov	s2, 1*8(rp)
	mov	ap, 2*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_4)
	C   0 1 2 3
	C 0     h x
	C 1   h x x
	C 2 h x x x
	C 3 x x x x
	mov	bp_param, bp
	mov	0*8(bp_param), %rdx

	mulx	2*8(ap), s0, s0		C a2 b0
	push	s5
	xor	R32(s1), R32(s1)
	mulx	3*8(ap), s2, s3		C a3 b0
	push	s6
	adcx	s2, s0
	C 0, 3 + c

	mov	1*8(bp), %rdx
	mulx	2*8(ap), s4, s5		C a2 b1
	push	s7
	adox	s4, s0
	adcx	s5, s3
	mulx	3*8(ap), s6, s2		C a3 b1
	adox	s6, s3
	adcx	s1, s2
	mulx	1*8(ap), s4, s4		C a1 b1
	adox	s1, s2
	C (0, 4), 3, 2

	mov	2*8(bp), %rdx
	mulx	1*8(ap), s7, s5		C a1 b2
	adox	s4, s0
	mulx	2*8(ap), s6, s4		C a2 b2
	adcx	s7, s0
	adox	s5, s3
	mulx	3*8(ap), s7, s5		C a3 b2
	adcx	s6, s3
	adox	s4, s2
	mulx	0*8(ap), s6, s6		C a0 b2
	adcx	s7, s2
	adox	s1, s5
	C (0, 6), 3, 2, 5 + c

	mov	3*8(bp), %rdx
	mulx	0*8(ap), bp, s4		C a0 b3
	adcx	s1, s5
	adox	s6, s0
	mulx	1*8(ap), s7, s6		C a1 b3
	adcx	bp, s0
	adox	s4, s3
	mulx	2*8(ap), bp, s4		C a2 b3
	adcx	s7, s3
	adox	s6, s2
	mulx	3*8(ap), s7, ap		C a3 b3
	adcx	bp, s2
	adox	s5, s4
	adcx	s7, s4
	mov	$0, %edx
	adox	s1, ap
	pop	s7
	adc	s1, ap
	C 0, 3, 2, 4, ap

	js	L(4)
	add	s0, s0
	adc	s3, s3
	adc	s2, s2
	inc	%edx
	adc	s4, s4
	adc	ap, ap
L(4):	mov	s3, 0*8(rp)
	pop	s6
	mov	s2, 1*8(rp)
	pop	s5
	mov	s4, 2*8(rp)
	mov	ap, 3*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_5)
	C   0 1 2 3 4
	C 0       h x
	C 1     h x x
	C 2   h x x x
	C 3 h x x x x
	C 4 x x x x x

	mov	bp_param, bp
	mov	0*8(bp_param), %rdx

	mulx	3*8(ap), s0, s0		C a3 b0
	push	s5
	xor	R32(s1), R32(s1)
	mulx	4*8(ap), s2, s3		C a4 b0
	push	s6
	adox	s2, s0
	C 0, 3 + o

	mov	1*8(bp), %rdx
	mulx	3*8(ap), s4, s5		C a3 b1
	push	s7
	adcx	s4, s0
	adox	s5, s3
	mulx	4*8(ap), s6, s2		C a4 b1
	push	s8
	adcx	s6, s3
	adox	s1, s2
	mulx	2*8(ap), s4, s4		C a2 b1
	adcx	s1, s2
	C (0, 4), 3, 2

	mov	2*8(bp), %rdx
	mulx	2*8(ap), s7, s8		C a2 b2
	adox	s4, s0
	mulx	3*8(ap), s5, s6		C a3 b2
	adcx	s7, s0
	adox	s8, s3
	mulx	4*8(ap), s4, s7		C a4 b2
	adcx	s5, s3
	adox	s6, s2
	mulx	1*8(ap), s8, s8		C a1 b2
	adcx	s4, s2
	adox	s1, s7
	C (0, 8), 3, 2, 7 + c

	mov	3*8(bp), %rdx
	mulx	1*8(ap), s5, s6		C a1 b3
	adcx	s1, s7
	adox	s8, s0
	mulx	2*8(ap), s4, s8		C a2 b3
	adcx	s5, s0
	adox	s6, s3
	mulx	3*8(ap), s5, s6		C a3 b3
	adcx	s4, s3
	adox	s8, s2
	mulx	4*8(ap), s4, s8		C a4 b3
	adcx	s5, s2
	adox	s6, s7
	mulx	0*8(ap), s5, s5		C a0 b3
	adcx	s4, s7
	adox	s1, s8
	C (0, 5), 3, 2, 7, 8 + c

	mov	4*8(bp), %rdx
	mulx	0*8(ap), s6, s4		C a0 b4
	adcx	s1, s8
	adox	s5, s0
	mulx	1*8(ap), bp, s5		C a1 b4
	adcx	s6, s0
	adox	s4, s3
	mulx	2*8(ap), s6, s4		C a2 b4
	adcx	bp, s3
	adox	s5, s2
	mulx	3*8(ap), bp, s5		C a3 b4
	adcx	s6, s2
	adox	s4, s7
	mulx	4*8(ap), s6, s4		C a4 b4
	adcx	s7, bp
	adox	s5, s8
	C 0, 3, 2, bp, (8 + c, 6), 4 + o

	adcx	s8, s6
	adox	s1, s4
	mov	$0, %edx
	pop	s8
	adc	s1, s4
	C 0, 3, 2, bp, 6, 4

	js	L(5)
	add	s0, s0
	adc	s3, s3
	adc	s2, s2
	inc	%edx
	adc	bp, bp
	adc	s6, s6
	adc	s4, s4
L(5):	mov	s3, 0*8(rp)
	pop	s7
	mov	s2, 1*8(rp)
	mov	bp, 2*8(rp)
	mov	s6, 3*8(rp)
	pop	s6
	mov	s4, 4*8(rp)
	pop	s5

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_6)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	4*8(%rsi), %r8, %rax
	mulx	5*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %rbx
	mulx	4*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	5*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %rbp
	mulx	3*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	4*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	5*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r12
	mulx	2*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	3*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r13
	mulx	1*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	2*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r13
	adcx	%r8, %rax
	adcx	%r13, %r10
	mulx	1*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adox	%r13, %r11
	mulx	2*8(%rsi), %r8, %r13
	adcx	%r8, %r11
	adcx	%r13, %rbx
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adox	%r13, %rbp
	mulx	4*8(%rsi), %r8, %r13
	adcx	%r8, %rbp
	adcx	%r13, %r12
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	$0, %rdx
	test	%r13, %r13
	setns	%dl
	js	L(6)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
L(6):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_7)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	5*8(%rsi), %r8, %rax
	mulx	6*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %rbx
	mulx	5*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	6*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %rbp
	mulx	4*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	5*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	6*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r12
	mulx	3*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r13
	mulx	2*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r14
	mulx	1*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	2*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r14
	adcx	%r8, %rax
	adcx	%r14, %r10
	mulx	1*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adox	%r14, %r11
	mulx	2*8(%rsi), %r8, %r14
	adcx	%r8, %r11
	adcx	%r14, %rbx
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adox	%r14, %rbp
	mulx	4*8(%rsi), %r8, %r14
	adcx	%r8, %rbp
	adcx	%r14, %r12
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adox	%r14, %r13
	mulx	6*8(%rsi), %r8, %r14
	adcx	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	$0, %rdx
	test	%r14, %r14
	setns	%dl
	js	L(7)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	adc	%r14, %r14
L(7):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)
	mov	%r14, 6*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_8)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	6*8(%rsi), %r8, %rax
	mulx	7*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	5*8(%rsi), %r8, %rbx
	mulx	6*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	7*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %rbp
	mulx	5*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	6*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	7*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %r12
	mulx	4*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	7*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r13
	mulx	3*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	7*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r14
	mulx	2*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	7*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r15
	mulx	1*8(%rsi), %r8, %r14
	adcx	%r15, %rax
	adox	%r8, %rax
	adcx	%r14, %r10
	mulx	2*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adcx	%r14, %r11
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %r11
	adcx	%r14, %rbx
	mulx	4*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adcx	%r14, %rbp
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %rbp
	adcx	%r14, %r12
	mulx	6*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adcx	%r14, %r13
	mulx	7*8(%rsi), %r8, %r14
	adox	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	7*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r15
	adcx	%r8, %rax
	adcx	%r15, %r10
	mulx	1*8(%rsi), %r8, %r15
	adox	%r8, %r10
	adox	%r15, %r11
	mulx	2*8(%rsi), %r8, %r15
	adcx	%r8, %r11
	adcx	%r15, %rbx
	mulx	3*8(%rsi), %r8, %r15
	adox	%r8, %rbx
	adox	%r15, %rbp
	mulx	4*8(%rsi), %r8, %r15
	adcx	%r8, %rbp
	adcx	%r15, %r12
	mulx	5*8(%rsi), %r8, %r15
	adox	%r8, %r12
	adox	%r15, %r13
	mulx	6*8(%rsi), %r8, %r15
	adcx	%r8, %r13
	adcx	%r15, %r14
	mulx	7*8(%rsi), %r8, %r15
	adox	%r8, %r14
	adcx	%r9, %r15
	adox	%r9, %r15

	mov	$0, %rdx
	test	%r15, %r15
	setns	%dl
	js	L(8)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	adc	%r14, %r14
	adc	%r15, %r15
L(8):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)
	mov	%r14, 6*8(%rdi)
	mov	%r15, 7*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_9)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	push	%rdi

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	7*8(%rsi), %r8, %rax
	mulx	8*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	6*8(%rsi), %r8, %rbx
	mulx	7*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	8*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	5*8(%rsi), %r8, %rbp
	mulx	6*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	7*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	8*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %r12
	mulx	5*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	7*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	8*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %r13
	mulx	4*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	7*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	8*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r14
	mulx	3*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	7*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	8*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r15
	mulx	2*8(%rsi), %r8, %r14
	adcx	%r15, %rax
	adox	%r8, %rax
	adcx	%r14, %r10
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adcx	%r14, %r11
	mulx	4*8(%rsi), %r8, %r14
	adox	%r8, %r11
	adcx	%r14, %rbx
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adcx	%r14, %rbp
	mulx	6*8(%rsi), %r8, %r14
	adox	%r8, %rbp
	adcx	%r14, %r12
	mulx	7*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adcx	%r14, %r13
	mulx	8*8(%rsi), %r8, %r14
	adox	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	7*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	mulx	1*8(%rsi), %r8, %r15
	adcx	%rdi, %rax
	adox	%r8, %rax
	adcx	%r15, %r10
	mulx	2*8(%rsi), %r8, %r15
	adox	%r8, %r10
	adcx	%r15, %r11
	mulx	3*8(%rsi), %r8, %r15
	adox	%r8, %r11
	adcx	%r15, %rbx
	mulx	4*8(%rsi), %r8, %r15
	adox	%r8, %rbx
	adcx	%r15, %rbp
	mulx	5*8(%rsi), %r8, %r15
	adox	%r8, %rbp
	adcx	%r15, %r12
	mulx	6*8(%rsi), %r8, %r15
	adox	%r8, %r12
	adcx	%r15, %r13
	mulx	7*8(%rsi), %r8, %r15
	adox	%r8, %r13
	adcx	%r15, %r14
	mulx	8*8(%rsi), %r8, %r15
	adox	%r8, %r14
	adcx	%r9, %r15
	adox	%r9, %r15

	mov	8*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	adcx	%r8, %rax
	adcx	%rdi, %r10
	mulx	1*8(%rsi), %r8, %rdi
	adox	%r8, %r10
	adox	%rdi, %r11
	mulx	2*8(%rsi), %r8, %rdi
	adcx	%r8, %r11
	adcx	%rdi, %rbx
	mulx	3*8(%rsi), %r8, %rdi
	adox	%r8, %rbx
	adox	%rdi, %rbp
	mulx	4*8(%rsi), %r8, %rdi
	adcx	%r8, %rbp
	adcx	%rdi, %r12
	mulx	5*8(%rsi), %r8, %rdi
	adox	%r8, %r12
	adox	%rdi, %r13
	mulx	6*8(%rsi), %r8, %rdi
	adcx	%r8, %r13
	adcx	%rdi, %r14
	mulx	7*8(%rsi), %r8, %rdi
	adox	%r8, %r14
	adox	%rdi, %r15
	mulx	8*8(%rsi), %r8, %rdi
	adcx	%r8, %r15
	adcx	%r9, %rdi
	pop	%r8
	adox	%r9, %rdi

	mov	$0, %rdx
	test	%rdi, %rdi
	setns	%dl
	js	L(9)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	adc	%r14, %r14
	adc	%r15, %r15
	adc	%rdi, %rdi
L(9):
	mov	%r10, 0*8(%r8)
	mov	%r11, 1*8(%r8)
	mov	%rbx, 2*8(%r8)
	mov	%rbp, 3*8(%r8)
	mov	%r12, 4*8(%r8)
	mov	%r13, 5*8(%r8)
	mov	%r14, 6*8(%r8)
	mov	%r15, 7*8(%r8)
	mov	%rdi, 8*8(%r8)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()
