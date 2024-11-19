// barrett_mul_neon.S

    .text
    .global barrett_mul_neon
    .type barrett_mul_neon, %function

// Function prototype:
// void barrett_mul_neon(uint32_t *a, uint32_t *b, uint32_t *N, uint32_t *mu, uint32_t *z);

barrett_mul_neon:
    // Input registers:
    // x0 - pointer to array 'a' (uint32_t a[4])
    // x1 - pointer to array 'b' (uint32_t b[4])
    // x2 - pointer to array 'N' (uint32_t N[4])
    // x3 - pointer to array 'mu' (uint32_t mu[4])  // Precomputed floor(bR/N)_2
    // x4 - pointer to array 'z' (uint32_t z[4])    // Result

    // Load operands into vector registers
    ld1    {v0.4s}, [x0]           // v0 = a[0..3]
    ld1    {v1.4s}, [x1]           // v1 = b[0..3]
    ld1    {v2.4s}, [x2]           // v2 = N[0..3]
    ld1    {v3.4s}, [x3]           // v3 = mu[0..3]

    // Step 1: Multiply z = a * b
    mul    v4.4s, v0.4s, v1.4s     // v4 = v0 * v1

    // Step 2: Compute t = sqrdmulh(a, mu)
    sqrdmulh v5.4s, v0.4s, v3.4s   // v5 = sqrdmulh(v0, v3)

    // Step 3: Compute z = z - t * N using Multiply-Subtract
    mls    v6.4s, v5.4s, v2.4s, v4.4s   // v6 = v4 - (v5 * v2)

    // Store the result
    st1    {v6.4s}, [x4]           // Store z[0..3]

    ret
