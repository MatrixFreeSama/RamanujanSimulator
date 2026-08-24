#include "ramanujan_c_common.h"
#include <json-c/json.h>
#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Constructive certificate materializer with a geometric precision ladder.
 *
 * The mathematical path is unchanged: the exact seed and the complete
 * certificate transform sequence are replayed.  What changes is the arithmetic
 * schedule.  Instead of executing every modular-root Newton correction at the
 * final precision from the beginning, the path is rebuilt at successively
 * doubled precisions.  The previous stage supplies the small-branch root seed
 * for the next stage.  In the quadratic regime one Newton correction therefore
 * approximately doubles the number of correct bits.
 */

static void usage(const char *argv0) {
    fprintf(stderr,
            "Usage: %s CERTIFICATE [--digits N] [--guard N] [--output FILE] [--trace-newton]\n",
            argv0);
}

static const char *jstr(struct json_object *o, const char *key) {
    struct json_object *v = NULL;
    if (!json_object_object_get_ex(o, key, &v)) return NULL;
    return json_object_get_string(v);
}

static int poly_from_json(struct json_object *node, RjPoly *out) {
    struct json_object *terms = NULL, *dx = NULL, *dy = NULL;
    if (!json_object_object_get_ex(node, "terms", &terms) ||
        !json_object_is_type(terms, json_type_array)) return 0;
    if (!json_object_object_get_ex(node, "degree_x", &dx) ||
        !json_object_object_get_ex(node, "degree_y", &dy)) return 0;

    size_t n = json_object_array_length(terms);
    memset(out, 0, sizeof(*out));
    out->degree_x = json_object_get_int(dx);
    out->degree_y = json_object_get_int(dy);
    out->nterms = n;
    out->terms = (RjTerm *)calloc(n, sizeof(RjTerm));
    if (!out->terms) return 0;

    for (size_t k = 0; k < n; ++k) {
        struct json_object *t = json_object_array_get_idx(terms, k);
        struct json_object *ix = NULL, *iy = NULL, *c = NULL;
        if (!json_object_object_get_ex(t, "x_pow", &ix) ||
            !json_object_object_get_ex(t, "y_pow", &iy) ||
            !json_object_object_get_ex(t, "coeff", &c)) {
            rj_poly_clear(out);
            return 0;
        }
        out->terms[k].ix = json_object_get_int(ix);
        out->terms[k].iy = json_object_get_int(iy);
        mpz_init(out->terms[k].coeff);
        if (mpz_set_str(out->terms[k].coeff,
                        json_object_get_string(c), 10) != 0) {
            rj_poly_clear(out);
            return 0;
        }
    }
    return 1;
}

static int seed_from_json(struct json_object *rec, RjState *s) {
    struct json_object *seed = NULL, *z0 = NULL, *alpha0 = NULL;
    if (!json_object_object_get_ex(rec, "exact_seed", &seed)) return 0;
    if (!json_object_object_get_ex(seed, "z0", &z0) ||
        !json_object_object_get_ex(seed, "alpha0", &alpha0)) return 0;

    const char *zn = jstr(z0, "numerator");
    const char *zd = jstr(z0, "denominator");
    const char *an = jstr(alpha0, "numerator");
    const char *ad = jstr(alpha0, "denominator");
    if (!zn || !zd || !an || !ad) return 0;

    mp_bitcnt_t bits = mpf_get_prec(s->z);
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    mpf_t a, C, tmp;
    mpf_init2(a, bits);
    mpf_init2(C, bits);
    mpf_init2(tmp, bits);

    mpz_set_str(n, zn, 10);
    mpz_set_str(d, zd, 10);
    mpf_set_z(s->z, n);
    mpf_set_z(tmp, d);
    mpf_div(s->z, s->z, tmp);

    mpz_set_str(n, an, 10);
    mpz_set_str(d, ad, 10);
    mpf_set_z(a, n);
    mpf_set_z(tmp, d);
    mpf_div(a, a, tmp);

    mpf_ui_sub(tmp, 1, s->z);
    mpf_sqrt(tmp, tmp);
    mpf_set_ui(C, 163);
    mpf_sqrt(C, C);
    mpf_div_ui(C, C, 6);
    mpf_mul(C, C, tmp);

    mpf_ui_sub(tmp, 1, a);
    mpf_mul(s->u, C, tmp);
    mpf_mul_ui(s->v, C, 6);

    mpz_clear(n);
    mpz_clear(d);
    mpf_clear(a);
    mpf_clear(C);
    mpf_clear(tmp);
    return 1;
}

typedef struct {
    RjState *state;
    size_t count;
    mp_bitcnt_t bits;
} RjStatePath;

static int path_init(RjStatePath *path, size_t count, mp_bitcnt_t bits) {
    memset(path, 0, sizeof(*path));
    path->state = (RjState *)calloc(count, sizeof(RjState));
    if (!path->state) return 0;
    path->count = count;
    path->bits = bits;
    for (size_t i = 0; i < count; ++i) rj_state_init(&path->state[i], bits);
    return 1;
}

static void path_clear(RjStatePath *path) {
    if (!path || !path->state) return;
    for (size_t i = 0; i < path->count; ++i) rj_state_clear(&path->state[i]);
    free(path->state);
    memset(path, 0, sizeof(*path));
}

static long double decimal_capacity(mp_bitcnt_t bits) {
    return (long double)bits * 0.30102999566398119521373889472449L;
}

static long double minus_log10_abs(const mpf_t x) {
    if (mpf_sgn(x) == 0) return INFINITY;
    mpf_t a;
    mpf_init2(a, mpf_get_prec(x));
    mpf_abs(a, x);
    mp_exp_t e2 = 0;
    double m = mpf_get_d_2exp(&e2, a);
    mpf_clear(a);
    return -(log10l(fabsl((long double)m)) +
             (long double)e2 * log10l(2.0L));
}

/*
 * Same modular state transport as rj_transform_step, but the modular-root
 * Newton solve may start from the previous precision stage.  The cap is only
 * a guard.  Normally the warm-started stages stabilize in one or two steps.
 */
static int transform_step_refined(const RjPoly *poly, int p,
                                  const RjState *in,
                                  const mpf_t previous_y,
                                  int have_previous_y,
                                  int newton_cap,
                                  RjState *out,
                                  int *used_iterations) {
    mp_bitcnt_t bits = mpf_get_prec(in->z);
    mpf_t y, oldy, P, Px, Py, Pxx, Pxy, Pyy, delta;
    mpf_t r, s, M, L, tmp, tmp2, den;
    mpf_init2(y, bits);       mpf_init2(oldy, bits);
    mpf_init2(P, bits);       mpf_init2(Px, bits);
    mpf_init2(Py, bits);      mpf_init2(Pxx, bits);
    mpf_init2(Pxy, bits);     mpf_init2(Pyy, bits);
    mpf_init2(delta, bits);   mpf_init2(r, bits);
    mpf_init2(s, bits);       mpf_init2(M, bits);
    mpf_init2(L, bits);       mpf_init2(tmp, bits);
    mpf_init2(tmp2, bits);    mpf_init2(den, bits);

    if (have_previous_y) {
        mpf_set(y, previous_y);
    } else {
        /* Small-branch predictor y0 = x^p / 1728^(p-1). */
        mpf_pow_ui(y, in->z, (unsigned long)p);
        mpf_set_ui(den, 1728);
        mpf_pow_ui(den, den, (unsigned long)(p - 1));
        mpf_div(y, y, den);
    }

    int used = 0;
    for (int it = 0; it < newton_cap; ++it) {
        rj_eval_poly_all(poly, in->z, y, P, Px, Py, Pxx, Pxy, Pyy);
        if (mpf_sgn(Py) == 0) goto fail;
        mpf_div(delta, P, Py);
        mpf_set(oldy, y);
        mpf_sub(y, y, delta);
        used = it + 1;
        if (mpf_cmp(y, oldy) == 0) break;
    }

    rj_eval_poly_all(poly, in->z, y, P, Px, Py, Pxx, Pxy, Pyy);
    if (mpf_sgn(Py) == 0 || mpf_sgn(y) == 0) goto fail;

    /* r = -Px/Py */
    mpf_div(r, Px, Py);
    mpf_neg(r, r);

    /* s = -(Pxx + 2 Pxy r + Pyy r^2)/Py */
    mpf_mul(tmp, Pxy, r);
    mpf_mul_ui(tmp, tmp, 2);
    mpf_add(tmp, tmp, Pxx);
    mpf_mul(tmp2, r, r);
    mpf_mul(tmp2, tmp2, Pyy);
    mpf_add(tmp, tmp, tmp2);
    mpf_div(s, tmp, Py);
    mpf_neg(s, s);

    /* M = (1/p)(x/y)r sqrt((1-x)/(1-y)) */
    mpf_div(M, in->z, y);
    mpf_mul(M, M, r);
    mpf_div_ui(M, M, (unsigned long)p);
    mpf_ui_sub(tmp, 1, in->z);
    mpf_ui_sub(tmp2, 1, y);
    mpf_div(tmp, tmp, tmp2);
    mpf_sqrt(tmp, tmp);
    mpf_mul(M, M, tmp);
    if (mpf_sgn(M) == 0 || mpf_sgn(r) == 0) goto fail;

    /* L = 1-xr/y+xs/r-x/(2(1-x))+xr/(2(1-y)) */
    mpf_set_ui(L, 1);
    mpf_mul(tmp, in->z, r); mpf_div(tmp, tmp, y); mpf_sub(L, L, tmp);
    mpf_mul(tmp, in->z, s); mpf_div(tmp, tmp, r); mpf_add(L, L, tmp);
    mpf_ui_sub(tmp2, 1, in->z); mpf_mul_ui(tmp2, tmp2, 2);
    mpf_div(tmp, in->z, tmp2); mpf_sub(L, L, tmp);
    mpf_ui_sub(tmp2, 1, y); mpf_mul_ui(tmp2, tmp2, 2);
    mpf_mul(tmp, in->z, r); mpf_div(tmp, tmp, tmp2); mpf_add(L, L, tmp);

    /* u'=(u-vL)/M ; v'=v*x*r/(y*M) */
    mpf_mul(tmp, in->v, L);
    mpf_sub(tmp, in->u, tmp);
    mpf_div(out->u, tmp, M);
    mpf_mul(tmp, in->v, in->z);
    mpf_mul(tmp, tmp, r);
    mpf_mul(tmp2, y, M);
    mpf_div(out->v, tmp, tmp2);
    mpf_set(out->z, y);

    if (used_iterations) *used_iterations = used;
    mpf_clear(y); mpf_clear(oldy); mpf_clear(P); mpf_clear(Px); mpf_clear(Py);
    mpf_clear(Pxx); mpf_clear(Pxy); mpf_clear(Pyy); mpf_clear(delta);
    mpf_clear(r); mpf_clear(s); mpf_clear(M); mpf_clear(L);
    mpf_clear(tmp); mpf_clear(tmp2); mpf_clear(den);
    return 1;

fail:
    mpf_clear(y); mpf_clear(oldy); mpf_clear(P); mpf_clear(Px); mpf_clear(Py);
    mpf_clear(Pxx); mpf_clear(Pxy); mpf_clear(Pyy); mpf_clear(delta);
    mpf_clear(r); mpf_clear(s); mpf_clear(M); mpf_clear(L);
    mpf_clear(tmp); mpf_clear(tmp2); mpf_clear(den);
    return 0;
}

/*
 * Build the whole certificate path at one arithmetic precision.  When a lower
 * precision path is supplied, each modular root is warm-started from the
 * corresponding previous-layer z value.
 */
static int build_path_stage(struct json_object *rec,
                            struct json_object *seq,
                            const RjPoly *p2,
                            const RjPoly *p5,
                            mp_bitcnt_t bits,
                            const RjStatePath *previous,
                            RjStatePath *out,
                            int *total_newton) {
    size_t L = json_object_array_length(seq);
    if (!path_init(out, L + 1, bits)) return 0;
    if (!seed_from_json(rec, &out->state[0])) {
        path_clear(out);
        return 0;
    }

    int total = 0;
    for (size_t k = 0; k < L; ++k) {
        int p = json_object_get_int(json_object_array_get_idx(seq, k));
        const RjPoly *poly = (p == 2) ? p2 : ((p == 5) ? p5 : NULL);
        int used = 0;
        int have_previous = previous && previous->state &&
                            previous->count == L + 1;
        if (!poly ||
            !transform_step_refined(poly, p, &out->state[k],
                                    have_previous ? previous->state[k + 1].z
                                                  : out->state[k].z,
                                    have_previous,
                                    32, &out->state[k + 1], &used)) {
            path_clear(out);
            return 0;
        }
        total += used;
    }

    if (total_newton) *total_newton = total;
    return 1;
}

/*
 * Reciprocal by Newton-Schulz:
 *
 *     x_{n+1} = x_n (2 - a x_n)
 *
 * Once the seed is in the quadratic basin, the relative error squares.  The
 * working precision is doubled after each correction so arithmetic effort and
 * available correct bits grow together rather than paying final precision at
 * every early iteration.
 */
static int reciprocal_precision_doubling(const mpf_t a,
                                         mp_bitcnt_t target_bits,
                                         mpf_t out,
                                         int trace) {
    double ad = mpf_get_d(a);
    if (ad == 0.0 || !isfinite(ad)) return 0;

    mp_bitcnt_t bits = target_bits < 128 ? target_bits : 128;
    mpf_t x;
    mpf_init2(x, bits);
    mpf_set_d(x, 1.0 / ad);

    /* Fill the initial 128-bit stage from the ~53-bit hardware seed. */
    {
        mpf_t aa, t;
        mpf_init2(aa, bits);
        mpf_init2(t, bits);
        mpf_set(aa, a);
        for (int i = 0; i < 2; ++i) {
            mpf_mul(t, aa, x);
            mpf_ui_sub(t, 2, t);
            mpf_mul(x, x, t);
        }
        if (trace) {
            mpf_mul(t, aa, x);
            mpf_ui_sub(t, 1, t);
            long double d = minus_log10_abs(t);
            if (isinf(d)) d = decimal_capacity(bits);
            fprintf(stderr,
                    "RECIP_NEWTON bits=%lu capacity=%.1Lf residual_digits>=%.1Lf\n",
                    (unsigned long)bits, decimal_capacity(bits), d);
        }
        mpf_clear(aa);
        mpf_clear(t);
    }

    while (bits < target_bits) {
        mp_bitcnt_t next_bits = bits > target_bits / 2 ? target_bits : bits * 2;
        mpf_t next, aa, t;
        mpf_init2(next, next_bits);
        mpf_init2(aa, next_bits);
        mpf_init2(t, next_bits);
        mpf_set(next, x);
        mpf_set(aa, a);

        /* One correction at doubled precision approximately doubles accuracy. */
        mpf_mul(t, aa, next);
        mpf_ui_sub(t, 2, t);
        mpf_mul(next, next, t);

        if (trace) {
            mpf_mul(t, aa, next);
            mpf_ui_sub(t, 1, t);
            long double d = minus_log10_abs(t);
            if (isinf(d)) d = decimal_capacity(next_bits);
            fprintf(stderr,
                    "RECIP_NEWTON bits=%lu capacity=%.1Lf residual_digits>=%.1Lf\n",
                    (unsigned long)next_bits,
                    decimal_capacity(next_bits), d);
        }

        mpf_clear(x);
        mpf_init2(x, next_bits);
        mpf_set(x, next);
        mpf_clear(next);
        mpf_clear(aa);
        mpf_clear(t);
        bits = next_bits;
    }

    mpf_set(out, x);
    mpf_clear(x);
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 2;
    }

    const char *cert_path = argv[1];
    const char *out_path = NULL;
    uint64_t digits = 1000;
    uint64_t guard = 100;
    int trace_newton = 0;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--digits") == 0 && i + 1 < argc) {
            digits = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--guard") == 0 && i + 1 < argc) {
            guard = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--trace-newton") == 0) {
            trace_newton = 1;
        } else {
            usage(argv[0]);
            return 2;
        }
    }

    if (digits == 0 || digits > (uint64_t)SIZE_MAX - 2) {
        fprintf(stderr, "digits too large for this build\n");
        return 2;
    }

    struct json_object *root = json_object_from_file(cert_path);
    if (!root) {
        fprintf(stderr, "cannot parse certificate\n");
        return 1;
    }

    struct json_object *rec = NULL, *polys = NULL;
    struct json_object *p2n = NULL, *p5n = NULL, *seq = NULL;
    if (!json_object_object_get_ex(root, "constructive_recovery", &rec) ||
        !json_object_object_get_ex(rec, "modular_z_polynomials", &polys) ||
        !json_object_object_get_ex(polys, "2", &p2n) ||
        !json_object_object_get_ex(polys, "5", &p5n) ||
        !json_object_object_get_ex(rec, "transform_sequence", &seq)) {
        fprintf(stderr, "certificate missing constructive fields\n");
        json_object_put(root);
        return 1;
    }

    RjPoly p2 = {0}, p5 = {0};
    if (!poly_from_json(p2n, &p2) || !poly_from_json(p5n, &p5)) {
        fprintf(stderr, "bad polynomial table\n");
        json_object_put(root);
        return 1;
    }

    mp_bitcnt_t target_bits = rj_bits_for_decimal(digits, guard);
    mp_bitcnt_t stage_bits = target_bits < 512 ? target_bits : 512;
    RjStatePath previous = {0};
    int stage_index = 0;

    for (;;) {
        RjStatePath current = {0};
        int stage_newton = 0;
        if (!build_path_stage(rec, seq, &p2, &p5, stage_bits,
                              previous.state ? &previous : NULL,
                              &current, &stage_newton)) {
            fprintf(stderr,
                    "precision-doubling transform failed at stage %d (%lu bits)\n",
                    stage_index, (unsigned long)stage_bits);
            path_clear(&previous);
            rj_poly_clear(&p2);
            rj_poly_clear(&p5);
            json_object_put(root);
            return 1;
        }

        if (trace_newton) {
            fprintf(stderr,
                    "PATH_NEWTON stage=%d bits=%lu decimal_capacity=%.1Lf corrections=%d\n",
                    stage_index,
                    (unsigned long)stage_bits,
                    decimal_capacity(stage_bits),
                    stage_newton);
        }

        path_clear(&previous);
        previous = current;

        if (stage_bits == target_bits) break;
        stage_bits = stage_bits > target_bits / 2 ? target_bits : stage_bits * 2;
        ++stage_index;
    }

    size_t L = json_object_array_length(seq);
    RjState *final_state = &previous.state[L];
    mpf_t F, T, R, pi;
    mpf_init2(F, target_bits);
    mpf_init2(T, target_bits);
    mpf_init2(R, target_bits);
    mpf_init2(pi, target_bits);

    if (!rj_hyper_F_theta(final_state->z,
                          (unsigned long)(digits + guard),
                          F, T, 10000000UL)) {
        fprintf(stderr, "hypergeometric observer failed\n");
        return 1;
    }

    mpf_mul(R, final_state->u, F);
    mpf_mul(pi, final_state->v, T);
    mpf_add(R, R, pi);

    if (!reciprocal_precision_doubling(R, target_bits, pi, trace_newton)) {
        fprintf(stderr, "Newton reciprocal failed\n");
        return 1;
    }

    FILE *fp = stdout;
    if (out_path) {
        fp = fopen(out_path, "wb");
        if (!fp) {
            perror("fopen");
            return 1;
        }
    }

    if (!rj_write_decimal_mpf(fp, pi, (size_t)digits)) {
        fprintf(stderr, "decimal serialization failed\n");
        return 1;
    }
    fputc('\n', fp);

    if (out_path) {
        fclose(fp);
        printf("WROTE %s (%llu digits after decimal)\n",
               out_path, (unsigned long long)digits);
    }

    mpf_clear(F);
    mpf_clear(T);
    mpf_clear(R);
    mpf_clear(pi);
    path_clear(&previous);
    rj_poly_clear(&p2);
    rj_poly_clear(&p5);
    json_object_put(root);
    return 0;
}
