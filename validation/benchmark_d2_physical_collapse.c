#define main original_d2_basic_main
#include "benchmark_original_d2_basic_vs_chud.c"
#undef main

/*
 * Hybrid physical-embedding collapse for the D=-427, h(Delta)=2 formula.
 *
 * The lower product-tree levels remain exact in Q(omega),
 *   omega=(1+sqrt(61))/2, omega^2=omega+15.
 * Each exact chunk is then evaluated once at the physical embedding and the
 * upper tree is collapsed to one scalar mpf lane:
 *
 *   P_total <- P_total * P_chunk
 *   W_total <- W_total + P_total_before * W_chunk
 *
 * This deliberately removes the slow conjugate coordinate from the largest
 * product-tree levels while preserving exact algebraic arithmetic inside each
 * chunk.  It is a numerical physical-branch evaluator, so byte equality to
 * Chudnovsky-BS is the validation gate.
 */

static const char *HC_Z0 = "-59818419102592333";
static const char *HC_Z1 = "13579278976889262";
static const char *HC_ZD = "609541351191872000";
static const char *HC_B0 = "320004750671";
static const char *HC_B1 = "-56552805760";
static const char *HC_BD = "766923569391";

typedef struct {
    mpz_t a;
    mpz_t b;
} hc_ring_t;

typedef struct {
    mpz_t q;
    hc_ring_t p;
    hc_ring_t w;
} hc_seg_t;

typedef struct {
    mpz_t num, den, g, weight0;
    mpz_t m0, m1, m2, t0, t1;
    hc_ring_t rw;
} hc_ws_t;

static void hc_ring_init(hc_ring_t *r) {
    mpz_init(r->a);
    mpz_init(r->b);
}

static void hc_ring_clear(hc_ring_t *r) {
    mpz_clear(r->a);
    mpz_clear(r->b);
}

static void hc_seg_init(hc_seg_t *s) {
    mpz_init(s->q);
    hc_ring_init(&s->p);
    hc_ring_init(&s->w);
}

static void hc_seg_clear(hc_seg_t *s) {
    mpz_clear(s->q);
    hc_ring_clear(&s->p);
    hc_ring_clear(&s->w);
}

static void hc_ws_init(hc_ws_t *ws) {
    mpz_inits(ws->num, ws->den, ws->g, ws->weight0,
              ws->m0, ws->m1, ws->m2, ws->t0, ws->t1, NULL);
    hc_ring_init(&ws->rw);
}

static void hc_ws_clear(hc_ws_t *ws) {
    hc_ring_clear(&ws->rw);
    mpz_clears(ws->num, ws->den, ws->g, ws->weight0,
               ws->m0, ws->m1, ws->m2, ws->t0, ws->t1, NULL);
}

static void hc_ring_mul(hc_ring_t *out,
                        const hc_ring_t *x,
                        const hc_ring_t *y,
                        hc_ws_t *ws) {
    mpz_mul(ws->m0, x->a, y->a);
    mpz_mul(ws->m1, x->b, y->b);
    mpz_add(ws->t0, x->a, x->b);
    mpz_add(ws->t1, y->a, y->b);
    mpz_mul(ws->m2, ws->t0, ws->t1);

    mpz_set(out->a, ws->m0);
    mpz_addmul_ui(out->a, ws->m1, 15UL);
    mpz_sub(out->b, ws->m2, ws->m0);
}

static void hc_leaf(hc_seg_t *out,
                    unsigned long n,
                    const mpz_t z0,
                    const mpz_t z1,
                    const mpz_t zd,
                    const mpz_t b0,
                    const mpz_t b1,
                    const mpz_t bd,
                    hc_ws_t *ws) {
    mpz_set_ui(ws->num, 6UL * n + 1UL);
    mpz_mul_ui(ws->num, ws->num, 2UL * n + 1UL);
    mpz_mul_ui(ws->num, ws->num, 6UL * n + 5UL);

    mpz_set_ui(ws->den, n + 1UL);
    mpz_mul_ui(ws->den, ws->den, n + 1UL);
    mpz_mul_ui(ws->den, ws->den, n + 1UL);
    mpz_mul_ui(ws->den, ws->den, 72UL);
    mpz_mul(ws->den, ws->den, zd);

    mpz_gcd(ws->g, ws->num, ws->den);
    mpz_divexact(ws->num, ws->num, ws->g);
    mpz_divexact(out->q, ws->den, ws->g);

    mpz_mul(out->p.a, z0, ws->num);
    mpz_mul(out->p.b, z1, ws->num);

    mpz_mul_ui(ws->weight0, bd, 6UL * n);
    mpz_add(ws->weight0, ws->weight0, b0);
    mpz_mul(out->w.a, ws->weight0, out->q);
    mpz_mul(out->w.b, b1, out->q);
}

static void hc_build(hc_seg_t *out,
                     unsigned long a,
                     unsigned long b,
                     const mpz_t z0,
                     const mpz_t z1,
                     const mpz_t zd,
                     const mpz_t b0,
                     const mpz_t b1,
                     const mpz_t bd,
                     hc_ws_t *ws) {
    if (b - a == 1UL) {
        hc_leaf(out, a, z0, z1, zd, b0, b1, bd, ws);
        return;
    }

    unsigned long m = a + (b - a) / 2UL;
    hc_seg_t left, right;
    hc_seg_init(&left);
    hc_seg_init(&right);
    hc_build(&left, a, m, z0, z1, zd, b0, b1, bd, ws);
    hc_build(&right, m, b, z0, z1, zd, b0, b1, bd, ws);

    mpz_mul(out->q, left.q, right.q);
    hc_ring_mul(&out->p, &left.p, &right.p, ws);
    hc_ring_mul(&ws->rw, &left.p, &right.w, ws);

    mpz_mul(out->w.a, left.w.a, right.q);
    mpz_add(out->w.a, out->w.a, ws->rw.a);
    mpz_mul(out->w.b, left.w.b, right.q);
    mpz_add(out->w.b, out->w.b, ws->rw.b);

    hc_seg_clear(&left);
    hc_seg_clear(&right);
}

static void hc_project(const hc_seg_t *seg,
                       const mpf_t omega,
                       const mpf_t fbd,
                       mpf_t p,
                       mpf_t w,
                       mpf_t t0,
                       mpf_t t1,
                       mpf_t fq) {
    mpf_set_z(t0, seg->p.a);
    mpf_set_z(t1, seg->p.b);
    mpf_mul(t1, t1, omega);
    mpf_add(p, t0, t1);

    mpf_set_z(fq, seg->q);
    mpf_div(p, p, fq);

    mpf_set_z(t0, seg->w.a);
    mpf_set_z(t1, seg->w.b);
    mpf_mul(t1, t1, omega);
    mpf_add(w, t0, t1);
    mpf_div(w, w, fq);
    mpf_div(w, w, fbd);
}

static int hc_pi(uint64_t digits,
                 uint64_t guard,
                 unsigned long chunk_terms,
                 const char *path,
                 double *seconds,
                 unsigned long *terms_out) {
    double t0s = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);

    const long double digits_per_term =
        24.95589965765426673091470594098813130578L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard + 30.0L) /
                            digits_per_term) + 2UL;
    if (chunk_terms == 0UL) chunk_terms = 1024UL;

    mpz_t z0, z1, zd, b0, b1, bd;
    mpz_inits(z0, z1, zd, b0, b1, bd, NULL);
    if (mpz_set_str(z0, HC_Z0, 10) != 0 ||
        mpz_set_str(z1, HC_Z1, 10) != 0 ||
        mpz_set_str(zd, HC_ZD, 10) != 0 ||
        mpz_set_str(b0, HC_B0, 10) != 0 ||
        mpz_set_str(b1, HC_B1, 10) != 0 ||
        mpz_set_str(bd, HC_BD, 10) != 0) {
        mpz_clears(z0, z1, zd, b0, b1, bd, NULL);
        return 0;
    }

    mpf_t omega, root61, fbd, p_total, w_total, p_chunk, w_chunk,
          tmp0, tmp1, fq, tmp2, z, K, inv_pi, pi;
    mpf_init2(omega, bits); mpf_init2(root61, bits);
    mpf_init2(fbd, bits); mpf_init2(p_total, bits);
    mpf_init2(w_total, bits); mpf_init2(p_chunk, bits);
    mpf_init2(w_chunk, bits); mpf_init2(tmp0, bits);
    mpf_init2(tmp1, bits); mpf_init2(fq, bits);
    mpf_init2(tmp2, bits); mpf_init2(z, bits);
    mpf_init2(K, bits); mpf_init2(inv_pi, bits);
    mpf_init2(pi, bits);

    mpf_set_ui(root61, 61UL);
    mpf_sqrt(root61, root61);
    mpf_add_ui(omega, root61, 1UL);
    mpf_div_ui(omega, omega, 2UL);
    mpf_set_z(fbd, bd);

    mpf_set_ui(p_total, 1UL);
    mpf_set_ui(w_total, 0UL);

    hc_ws_t ws;
    hc_ws_init(&ws);
    for (unsigned long a = 0UL; a < terms; a += chunk_terms) {
        unsigned long b = a + chunk_terms;
        if (b > terms) b = terms;

        hc_seg_t seg;
        hc_seg_init(&seg);
        hc_build(&seg, a, b, z0, z1, zd, b0, b1, bd, &ws);
        hc_project(&seg, omega, fbd, p_chunk, w_chunk, tmp0, tmp1, fq);

        /* W_total += P_total * W_chunk; P_total *= P_chunk. */
        mpf_mul(tmp2, p_total, w_chunk);
        mpf_add(w_total, w_total, tmp2);
        mpf_mul(p_total, p_total, p_chunk);

        hc_seg_clear(&seg);
    }
    hc_ws_clear(&ws);

    /* z = (Z0 + Z1*omega)/ZD for the final prefactor. */
    mpf_set_z(tmp0, z0);
    mpf_set_z(tmp1, z1);
    mpf_mul(tmp1, tmp1, omega);
    mpf_add(z, tmp0, tmp1);
    mpf_set_z(fq, zd);
    mpf_div(z, z, fq);

    mpf_set_ui(K, 427UL);
    mpf_sqrt(K, K);
    mpf_div_ui(K, K, 6UL);
    mpf_ui_sub(tmp0, 1UL, z);
    mpf_sqrt(tmp0, tmp0);
    mpf_mul(K, K, tmp0);

    mpf_mul(inv_pi, K, w_total);
    mpf_ui_div(pi, 1UL, inv_pi);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0s;
    if (terms_out) *terms_out = terms;

    mpf_clears(omega, root61, fbd, p_total, w_total, p_chunk, w_chunk,
               tmp0, tmp1, fq, tmp2, z, K, inv_pi, pi, NULL);
    mpz_clears(z0, z1, zd, b0, b1, bd, NULL);
    return ok;
}

static int run_hc_case(uint64_t digits, unsigned long chunk_terms) {
    const uint64_t guard = 256ULL;
    char ref_path[256], hc_path[256];
    snprintf(ref_path, sizeof(ref_path), "build/hc_ref_%llu_%lu.txt",
             (unsigned long long)digits, chunk_terms);
    snprintf(hc_path, sizeof(hc_path), "build/hc_%llu_%lu.txt",
             (unsigned long long)digits, chunk_terms);

    double tc = 0.0, th = 0.0;
    unsigned long nc = 0UL, nh = 0UL;
    if (!chud_bs_pi(digits, guard, ref_path, &tc, &nc) ||
        !hc_pi(digits, guard, chunk_terms, hc_path, &th, &nh)) return 0;
    int eq = files_equal(ref_path, hc_path);

    printf("%llu,%lu,%.9f,%lu,%.9f,%lu,%.6f,%d\n",
           (unsigned long long)digits,
           chunk_terms,
           tc, nc,
           th, nh,
           th / tc,
           eq);
    fflush(stdout);
    return eq;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s DIGITS CHUNK_TERMS\n", argv[0]);
        return 2;
    }
    uint64_t digits = strtoull(argv[1], NULL, 10);
    unsigned long chunk_terms = strtoul(argv[2], NULL, 10);
    printf("digits,chunk_terms,chud_seconds,chud_terms,collapsed_seconds,d2_terms,collapsed_over_chud,equal\n");
    if (digits == 0ULL || chunk_terms == 0UL) return 2;
    return run_hc_case(digits, chunk_terms) ? 0 : 1;
}
