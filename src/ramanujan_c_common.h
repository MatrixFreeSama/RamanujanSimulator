#ifndef RAMANUJAN_C_COMMON_H
#define RAMANUJAN_C_COMMON_H

#include <gmp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    int ix;
    int iy;
    const char *coeff_dec;
} RjConstTerm;

typedef struct {
    int ix;
    int iy;
    mpz_t coeff;
} RjTerm;

typedef struct {
    int degree_x;
    int degree_y;
    size_t nterms;
    RjTerm *terms;
} RjPoly;

typedef struct {
    mpf_t z;
    mpf_t u;
    mpf_t v;
} RjState;

extern const int RJ_SEQUENCE_10K4[32];
extern const size_t RJ_SEQUENCE_10K4_LEN;

int rj_poly_init_from_const(RjPoly *out, int dx, int dy, const RjConstTerm *src, size_t n);
void rj_poly_clear(RjPoly *p);
int rj_get_phi_poly(int p, RjPoly *out);
int rj_get_z_poly(int p, RjPoly *out);
int rj_generate_phi_prime(int p, int target_hi, RjPoly *out);
int rj_z_polynomial_from_phi(const RjPoly *phi, RjPoly *out);
int rj_poly_equal(const RjPoly *a, const RjPoly *b);

void rj_state_init(RjState *s, mp_bitcnt_t bits);
void rj_state_clear(RjState *s);
void rj_initial_state(RjState *s);

void rj_eval_poly_all(const RjPoly *poly, const mpf_t x, const mpf_t y,
                      mpf_t P, mpf_t Px, mpf_t Py,
                      mpf_t Pxx, mpf_t Pxy, mpf_t Pyy);

int rj_transform_step(const RjPoly *poly, int p, int newton_iters,
                      const RjState *in, RjState *out);

int rj_hyper_F_theta(const mpf_t z, unsigned long working_digits,
                     mpf_t F, mpf_t thetaF, unsigned long max_terms);

void rj_inverse_pi_from_state(const RjState *s, unsigned long working_digits,
                              mpf_t inv_pi, mpf_t F, mpf_t thetaF);

/* Exact q-series residual check of Phi_p(j(q),j(q^p)) through q^hi. */
int rj_verify_phi_qseries_exact(const RjPoly *phi, int p, int hi,
                                int *out_lo, size_t *out_nonzero);

/* Decimal output helper. digits_after_decimal excludes the leading '3'. */
int rj_write_decimal_mpf(FILE *fp, const mpf_t x, size_t digits_after_decimal);

/* Utility: decimal digits -> binary precision with guard digits. */
mp_bitcnt_t rj_bits_for_decimal(uint64_t decimal_digits, uint64_t guard_digits);

#endif
