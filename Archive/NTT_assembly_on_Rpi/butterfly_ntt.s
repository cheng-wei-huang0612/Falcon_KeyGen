
    .global butterfly_two_layers_ntt
    .text

butterfly_two_layers_ntt:
    // Function parameters:
    // x0: int32_t *f      (pointer to array f)
    // w1: uint32_t i      (index i)
    // w2: uint32_t d      (stride d)
    // w3: int32_t w1      (twiddle factor w1)
    // w4: int32_t w2      (twiddle factor w2)
    // w5: int32_t w3      (twiddle factor w3)
    // w6: int32_t mod     (modulus P)
    // x7: int64_t R       (precomputed inverse R)

    // Compute indices
    // idx_d = i + d
    ADD     w9, w1, w2               // w9 = idx_d = i + d

    // idx_2d = i + 2*d
    ADD     w10, w1, w2, LSL #1      // w10 = idx_2d = i + 2*d

    // idx_3d = i + 3*d
    ADD     w11, w10, w2             // w11 = idx_3d = i + 3*d

    // Load values from f[]
    LDR     w12, [x0, w1, UXTW #2]   // w12  = f[i]
    LDR     w13, [x0, w9, UXTW #2]   // w13 = f[i + d]
    LDR     w14, [x0, w10, UXTW #2]  // w14 = f[i + 2*d]
    LDR     w15, [x0, w11, UXTW #2]  // w15 = f[i + 3*d]

    // First layer twiddle multiplication for f[i + 2*d]
    // A = w14, B = w3 (twiddle factor w1)
    SMULL   x14, w14, w3             // x14 = A * B mod 2^64
    MUL     x14, x14, x7             // x14 = (A * B) * R mod 2^64
    LSR     x14, x14, #32            // x14 = ((A * B * R) >> 32)
    ADD     x14, x14, #16            // x14 = e
    SMULL   x14, w14, w6             // x14 = e * P mod 2^64
    LSR     x14, x14, #32            // x14 = (e * P) >> 32

    // First layer twiddle multiplication for f[i + 3*d]
    // A = w15, B = w3
    SMULL   x15, w15, w3             // x15 = A * B mod 2^64
    MUL     x15, x15, x7             // x15 = (A * B) * R mod 2^64
    LSR     x15, x15, #32            // x15 = ((A * B * R) >> 32)
    ADD     x15, x15, #16            // x15 = e
    SMULL   x15, w15, w6             // x15 = e * P mod 2^64
    LSR     x15, x15, #32            // w15 = (e * P) >> 32

    // First layer add and subtract (in-place updates)

    ADD     w12, w12, w14

    SUB     w14, w12, w14, LSL #1

    ADD     w13, w13, w15

    SUB     w15, w13, w15, LSL #1


    // Second layer twiddle multiplication for f[i + d]
    // A = w13, B = w4 (twiddle factor w2)
    SMULL   x13, w13, w4             // x13 = A * B mod 2^64
    MUL     x13, x13, x7             // x13 = (A * B) * R mod 2^64
    LSR     x13, x13, #32            // x13 = ((A * B * R) >> 32)
    ADD     x13, x13, #16            // x13 = e
    SMULL   x13, w13, w6             // x13 = e * P mod 2^64
    LSR     x13, x13, #32            // w13 = (e * P) >> 32

    // Second layer twiddle multiplication for f[i + 3*d]
    // A = w15, B = w5 (twiddle factor w3)
    SMULL   x15, w15, w5             // x15 = A * B mod 2^64
    MUL     x15, x15, x7             // x15 = (A * B) * R mod 2^64
    LSR     x15, x15, #32            // x15 = ((A * B * R) >> 32)
    ADD     x15, x15, #16            // x15 = e
    SMULL   x15, w15, w6             // x15 = e * P mod 2^64
    LSR     x15, x15, #32            // w15 = (e * P) >> 32

    // Second layer add and subtract (in-place updates)

    ADD     w12, w12, w13

    SUB     w13, w12, w13, LSL #1

    ADD     w14, w14, w15

    SUB     w15, w14, w15, LSL #1

    // Store final results back to f[]
    STR     w12, [x0, w1, UXTW #2]    // f[i] = w12
    STR     w13, [x0, w9, UXTW #2]   // f[i + d] = w13
    STR     w14, [x0, w10, UXTW #2]  // f[i + 2*d] = w14
    STR     w15, [x0, w11, UXTW #2]  // f[i + 3*d] = w15

    RET
