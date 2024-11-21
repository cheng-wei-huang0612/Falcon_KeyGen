	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 0	sdk_version 15, 1
	.globl	_Signed_ModMul                  ; -- Begin function Signed_ModMul
	.p2align	2
_Signed_ModMul:                         ; @Signed_ModMul
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #48
	.cfi_def_cfa_offset 48
	str	x0, [sp, #40]
	str	x1, [sp, #32]
	str	x2, [sp, #24]
	str	x3, [sp, #16]
	str	wzr, [sp, #12]
	b	LBB0_1
LBB0_1:                                 ; =>This Inner Loop Header: Depth=1
	ldr	w8, [sp, #12]
	subs	w8, w8, #4
	cset	w8, ge
	tbnz	w8, #0, LBB0_6
	b	LBB0_2
LBB0_2:                                 ;   in Loop: Header=BB0_1 Depth=1
	ldr	x8, [sp, #40]
	ldrsw	x9, [sp, #12]
	ldrsw	x8, [x8, x9, lsl #2]
	ldr	x9, [sp, #32]
	ldrsw	x10, [sp, #12]
	ldrsw	x9, [x9, x10, lsl #2]
	mul	x8, x8, x9
	str	x8, [sp]
	ldr	x8, [sp]
	ldr	x9, [sp, #24]
	ldrsw	x10, [sp, #12]
	ldrsw	x10, [x9, x10, lsl #2]
	sdiv	x9, x8, x10
	mul	x9, x9, x10
	subs	x8, x8, x9
	str	x8, [sp]
	ldr	x8, [sp]
	sxtw	x8, w8
	str	x8, [sp]
	ldr	x8, [sp]
	ldr	x9, [sp, #24]
	ldrsw	x10, [sp, #12]
	ldr	w9, [x9, x10, lsl #2]
	asr	w9, w9, #2
	subs	x8, x8, w9, sxtw
	cset	w8, le
	tbnz	w8, #0, LBB0_4
	b	LBB0_3
LBB0_3:                                 ;   in Loop: Header=BB0_1 Depth=1
	ldr	x8, [sp]
	ldr	x9, [sp, #24]
	ldrsw	x10, [sp, #12]
	ldrsw	x9, [x9, x10, lsl #2]
	subs	x8, x8, x9
	str	x8, [sp]
	b	LBB0_4
LBB0_4:                                 ;   in Loop: Header=BB0_1 Depth=1
	ldr	x8, [sp]
	ldr	x9, [sp, #16]
	ldrsw	x10, [sp, #12]
                                        ; kill: def $w8 killed $w8 killed $x8
	str	w8, [x9, x10, lsl #2]
	b	LBB0_5
LBB0_5:                                 ;   in Loop: Header=BB0_1 Depth=1
	ldr	w8, [sp, #12]
	add	w8, w8, #1
	str	w8, [sp, #12]
	b	LBB0_1
LBB0_6:
	add	sp, sp, #48
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_NEON_Signed_Barrett_ModMul     ; -- Begin function NEON_Signed_Barrett_ModMul
	.p2align	2
_NEON_Signed_Barrett_ModMul:            ; @NEON_Signed_Barrett_ModMul
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #480
	stp	x28, x27, [sp, #464]            ; 16-byte Folded Spill
	.cfi_def_cfa_offset 480
	.cfi_offset w27, -8
	.cfi_offset w28, -16
	str	x0, [sp, #296]
	str	x1, [sp, #288]
	str	x2, [sp, #280]
	str	x3, [sp, #272]
	str	x4, [sp, #264]
	ldr	x8, [sp, #296]
	ldr	q0, [x8]
	str	q0, [sp, #224]
	ldr	q0, [sp, #224]
	str	q0, [sp, #208]
	ldr	q0, [sp, #208]
	str	q0, [sp, #240]
	ldr	x8, [sp, #288]
	ldr	q0, [x8]
	str	q0, [sp, #176]
	ldr	q0, [sp, #176]
	str	q0, [sp, #160]
	ldr	q0, [sp, #160]
	str	q0, [sp, #192]
	ldr	x8, [sp, #280]
	ldr	q0, [x8]
	str	q0, [sp, #128]
	ldr	q0, [sp, #128]
	str	q0, [sp, #112]
	ldr	q0, [sp, #112]
	str	q0, [sp, #144]
	ldr	x8, [sp, #272]
	ldr	q0, [x8]
	str	q0, [sp, #80]
	ldr	q0, [sp, #80]
	str	q0, [sp, #64]
	ldr	q0, [sp, #64]
	str	q0, [sp, #96]
	ldr	q1, [sp, #240]
	ldr	q0, [sp, #192]
	str	q1, [sp, #336]
	str	q0, [sp, #320]
	ldr	q0, [sp, #336]
	ldr	q1, [sp, #320]
	mul.4s	v0, v0, v1
	str	q0, [sp, #304]
	ldr	q0, [sp, #304]
	str	q0, [sp, #48]
	ldr	q1, [sp, #240]
	ldr	q0, [sp, #144]
	str	q1, [sp, #384]
	str	q0, [sp, #368]
	ldr	q0, [sp, #384]
	ldr	q1, [sp, #368]
	sqrdmulh.4s	v0, v0, v1
	str	q0, [sp, #352]
	ldr	q0, [sp, #352]
	str	q0, [sp, #32]
	ldr	q2, [sp, #48]
	ldr	q1, [sp, #32]
	ldr	q0, [sp, #96]
	str	q2, [sp, #448]
	str	q1, [sp, #432]
	str	q0, [sp, #416]
	ldr	q0, [sp, #448]
	ldr	q1, [sp, #432]
	ldr	q2, [sp, #416]
	mul.4s	v1, v1, v2
	sub.4s	v0, v0, v1
	str	q0, [sp, #400]
	ldr	q0, [sp, #400]
	str	q0, [sp, #16]
	ldr	q0, [sp, #16]
	str	q0, [sp]
	ldr	x8, [sp, #264]
	ldr	q0, [sp]
	str	q0, [x8]
	ldp	x28, x27, [sp, #464]            ; 16-byte Folded Reload
	add	sp, sp, #480
	ret
	.cfi_endproc
                                        ; -- End function
	.section	__TEXT,__literal8,8byte_literals
	.p2align	3, 0x0                          ; -- Begin function Mu
lCPI2_0:
	.quad	0x41e0000000000000              ; double 2147483648
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_Mu
	.p2align	2
_Mu:                                    ; @Mu
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #64
	.cfi_def_cfa_offset 64
	str	x0, [sp, #56]
	str	x1, [sp, #48]
	str	x2, [sp, #40]
	str	wzr, [sp, #36]
	b	LBB2_1
LBB2_1:                                 ; =>This Inner Loop Header: Depth=1
	ldr	w8, [sp, #36]
	subs	w8, w8, #4
	cset	w8, ge
	tbnz	w8, #0, LBB2_4
	b	LBB2_2
LBB2_2:                                 ;   in Loop: Header=BB2_1 Depth=1
	ldr	x8, [sp, #56]
	ldrsw	x9, [sp, #36]
	ldr	s1, [x8, x9, lsl #2]
                                        ; implicit-def: $d0
	fmov	s0, s1
	sshll.2d	v0, v0, #0
                                        ; kill: def $d0 killed $d0 killed $q0
	scvtf	d0, d0
	adrp	x8, lCPI2_0@PAGE
	ldr	d1, [x8, lCPI2_0@PAGEOFF]
	fmul	d0, d0, d1
	str	d0, [sp, #24]
	ldr	d1, [sp, #24]
	ldr	x8, [sp, #48]
	ldrsw	x9, [sp, #36]
	ldr	s2, [x8, x9, lsl #2]
                                        ; implicit-def: $d0
	fmov	s0, s2
	sshll.2d	v0, v0, #0
                                        ; kill: def $d0 killed $d0 killed $q0
	scvtf	d2, d0
	fmov	d0, #2.00000000
	fmul	d2, d0, d2
	fdiv	d1, d1, d2
	str	d1, [sp, #16]
	ldr	d1, [sp, #16]
	fmov	d2, #0.50000000
	fadd	d1, d1, d2
	frintm	d1, d1
	str	d1, [sp, #8]
	ldr	d1, [sp, #8]
	fmul	d0, d0, d1
	fcvtzs	w8, d0
	str	w8, [sp, #4]
	ldr	w8, [sp, #4]
	ldr	x9, [sp, #40]
	ldrsw	x10, [sp, #36]
	str	w8, [x9, x10, lsl #2]
	b	LBB2_3
LBB2_3:                                 ;   in Loop: Header=BB2_1 Depth=1
	ldr	w8, [sp, #36]
	add	w8, w8, #1
	str	w8, [sp, #36]
	b	LBB2_1
LBB2_4:
	add	sp, sp, #64
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_main                           ; -- Begin function main
	.p2align	2
_main:                                  ; @main
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #192
	stp	x29, x30, [sp, #176]            ; 16-byte Folded Spill
	add	x29, sp, #176
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	adrp	x8, ___stack_chk_guard@GOTPAGE
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
	ldr	x8, [x8]
	stur	x8, [x29, #-8]
	str	wzr, [sp, #76]
	adrp	x8, l___const.main.a@PAGE
	add	x8, x8, l___const.main.a@PAGEOFF
	ldr	q0, [x8]
	sub	x8, x29, #32
	str	x8, [sp, #64]                   ; 8-byte Folded Spill
	stur	q0, [x29, #-32]
	adrp	x8, l___const.main.b@PAGE
	add	x8, x8, l___const.main.b@PAGEOFF
	ldr	q0, [x8]
	sub	x0, x29, #48
	str	x0, [sp, #32]                   ; 8-byte Folded Spill
	stur	q0, [x29, #-48]
	adrp	x8, l___const.main.q@PAGE
	add	x8, x8, l___const.main.q@PAGEOFF
	ldr	q0, [x8]
	sub	x1, x29, #64
	str	x1, [sp, #48]                   ; 8-byte Folded Spill
	stur	q0, [x29, #-64]
	adrp	x8, l___const.main.mu@PAGE
	add	x8, x8, l___const.main.mu@PAGEOFF
	ldr	q0, [x8]
	sub	x2, x29, #80
	str	x2, [sp, #40]                   ; 8-byte Folded Spill
	stur	q0, [x29, #-80]
	add	x8, sp, #80
	str	x8, [sp, #56]                   ; 8-byte Folded Spill
	str	xzr, [sp, #80]
	str	xzr, [sp, #88]
	bl	_Mu
	adrp	x0, l_.str@PAGE
	add	x0, x0, l_.str@PAGEOFF
	bl	_printf
	ldur	w8, [x29, #-80]
                                        ; implicit-def: $x12
	mov	x12, x8
	ldur	w8, [x29, #-76]
                                        ; implicit-def: $x11
	mov	x11, x8
	ldur	w8, [x29, #-72]
                                        ; implicit-def: $x10
	mov	x10, x8
	ldur	w9, [x29, #-68]
                                        ; implicit-def: $x8
	mov	x8, x9
	mov	x9, sp
	str	x12, [x9]
	str	x11, [x9, #8]
	str	x10, [x9, #16]
	str	x8, [x9, #24]
	adrp	x0, l_.str.1@PAGE
	add	x0, x0, l_.str.1@PAGEOFF
	bl	_printf
	ldr	x1, [sp, #32]                   ; 8-byte Folded Reload
	ldr	x2, [sp, #48]                   ; 8-byte Folded Reload
	ldr	x3, [sp, #56]                   ; 8-byte Folded Reload
	ldr	x0, [sp, #64]                   ; 8-byte Folded Reload
	bl	_Signed_ModMul
	adrp	x0, l_.str.2@PAGE
	add	x0, x0, l_.str.2@PAGEOFF
	bl	_printf
	ldr	w8, [sp, #80]
                                        ; implicit-def: $x12
	mov	x12, x8
	ldr	w8, [sp, #84]
                                        ; implicit-def: $x11
	mov	x11, x8
	ldr	w8, [sp, #88]
                                        ; implicit-def: $x10
	mov	x10, x8
	ldr	w9, [sp, #92]
                                        ; implicit-def: $x8
	mov	x8, x9
	mov	x9, sp
	str	x12, [x9]
	str	x11, [x9, #8]
	str	x10, [x9, #16]
	str	x8, [x9, #24]
	adrp	x0, l_.str.3@PAGE
	add	x0, x0, l_.str.3@PAGEOFF
	bl	_printf
	adrp	x0, l_.str.4@PAGE
	add	x0, x0, l_.str.4@PAGEOFF
	bl	_printf
	ldr	x1, [sp, #32]                   ; 8-byte Folded Reload
	ldr	x2, [sp, #40]                   ; 8-byte Folded Reload
	ldr	x3, [sp, #48]                   ; 8-byte Folded Reload
	ldr	x4, [sp, #56]                   ; 8-byte Folded Reload
	ldr	x0, [sp, #64]                   ; 8-byte Folded Reload
	bl	_NEON_Signed_Barrett_ModMul
	ldr	w8, [sp, #80]
                                        ; implicit-def: $x12
	mov	x12, x8
	ldr	w8, [sp, #84]
                                        ; implicit-def: $x11
	mov	x11, x8
	ldr	w8, [sp, #88]
                                        ; implicit-def: $x10
	mov	x10, x8
	ldr	w9, [sp, #92]
                                        ; implicit-def: $x8
	mov	x8, x9
	mov	x9, sp
	str	x12, [x9]
	str	x11, [x9, #8]
	str	x10, [x9, #16]
	str	x8, [x9, #24]
	adrp	x0, l_.str.5@PAGE
	add	x0, x0, l_.str.5@PAGEOFF
	bl	_printf
	ldur	x9, [x29, #-8]
	adrp	x8, ___stack_chk_guard@GOTPAGE
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
	ldr	x8, [x8]
	subs	x8, x8, x9
	cset	w8, eq
	tbnz	w8, #0, LBB3_2
	b	LBB3_1
LBB3_1:
	bl	___stack_chk_fail
LBB3_2:
	mov	w0, #0                          ; =0x0
	ldp	x29, x30, [sp, #176]            ; 16-byte Folded Reload
	add	sp, sp, #192
	ret
	.cfi_endproc
                                        ; -- End function
	.section	__TEXT,__literal16,16byte_literals
	.p2align	2, 0x0                          ; @__const.main.a
l___const.main.a:
	.long	123                             ; 0x7b
	.long	7234                            ; 0x1c42
	.long	7345                            ; 0x1cb1
	.long	745622                          ; 0xb6096

	.p2align	2, 0x0                          ; @__const.main.b
l___const.main.b:
	.long	7234                            ; 0x1c42
	.long	7345                            ; 0x1cb1
	.long	7456                            ; 0x1d20
	.long	756722                          ; 0xb8bf2

	.p2align	2, 0x0                          ; @__const.main.q
l___const.main.q:
	.long	33550337                        ; 0x1fff001
	.long	33540097                        ; 0x1ffc801
	.long	33538049                        ; 0x1ffc001
	.long	33533953                        ; 0x1ffb001

	.p2align	2, 0x0                          ; @__const.main.mu
l___const.main.mu:
	.long	463032                          ; 0x710b8
	.long	470280                          ; 0x72d08
	.long	477418                          ; 0x748ea
	.long	48459784                        ; 0x2e37008

	.section	__TEXT,__cstring,cstring_literals
l_.str:                                 ; @.str
	.asciz	"The mu(multiplier) values are\n"

l_.str.1:                               ; @.str.1
	.asciz	"Mu = %d %d %d %d\n\n"

l_.str.2:                               ; @.str.2
	.asciz	"The standard result is\n"

l_.str.3:                               ; @.str.3
	.asciz	"r = %d %d %d %d\n\n"

l_.str.4:                               ; @.str.4
	.asciz	"Barrett with NEON intrinsics\n"

l_.str.5:                               ; @.str.5
	.asciz	"r = %d %d %d %d\n"

.subsections_via_symbols
