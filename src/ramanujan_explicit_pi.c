#include "ramanujan_c_common.h"

#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(const char *argv0) {
    fprintf(stderr,
            "Usage: %s [--digits N] [--guard N] [--output FILE] [--quiet]\n",
            argv0);
}

int main(int argc, char **argv) {
    uint64_t digits = 1000;
    uint64_t guard = 120;
    const char *out_path = NULL;
    int quiet = 0;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--digits") == 0 && i + 1 < argc) {
            digits = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--guard") == 0 && i + 1 < argc) {
            guard = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--quiet") == 0) {
            quiet = 1;
        } else {
            usage(argv[0]);
            return 2;
        }
    }

    if (digits == 0 || digits > (uint64_t)SIZE_MAX - 2) {
        fprintf(stderr, "invalid digit count\n");
        return 2;
    }

    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    mpf_t z, C, tmp, term, sum, inv_pi, pi, den;
    mpf_init2(z, bits);
    mpf_init2(C, bits);
    mpf_init2(tmp, bits);
    mpf_init2(term, bits);
    mpf_init2(sum, bits);
    mpf_init2(inv_pi, bits);
    mpf_init2(pi, bits);
    mpf_init2(den, bits);

    /* Exact explicit seed:
     * z0 = -1/53360^3
     * alpha0 = 77265280/90856689
     * C0 = sqrt(163)/6 * sqrt(1-z0)
     * 1/pi = C0 * sum c_n z0^n (1-alpha0+6n)
     */
    mpf_set_ui(den, 53360);
    mpf_pow_ui(den, den, 3);
    mpf_ui_div(z, 1, den);
    mpf_neg(z, z);

    mpf_ui_sub(tmp, 1, z);
    mpf_sqrt(tmp, tmp);
    mpf_set_ui(C, 163);
    mpf_sqrt(C, C);
    mpf_mul(C, C, tmp);
    mpf_div_ui(C, C, 6);

    /*
     * 1-alpha0+6n = (13591409 + 545140134 n)/90856689.
     * Accumulating the integer numerator avoids one full-precision multiply
     * by alpha on every term without changing the explicit series.
     */
    mpf_set_ui(term, 1);
    mpf_set_ui(sum, 13591409UL);

    const long double digits_per_term = 14.1816474627254776555L;
    unsigned long terms =
        (unsigned long)(((long double)(digits + guard)) / digits_per_term) + 6UL;

    for (unsigned long n = 0; n + 1 < terms; ++n) {
        const unsigned long np1 = n + 1;

        /* c_{n+1} z^{n+1} from c_n z^n. */
        mpf_mul_ui(term, term, 6 * n + 1);
        mpf_mul_ui(term, term, 2 * n + 1);
        mpf_mul_ui(term, term, 6 * n + 5);
        mpf_div_ui(term, term, 72UL);
        mpf_div_ui(term, term, np1);
        mpf_div_ui(term, term, np1);
        mpf_div_ui(term, term, np1);
        mpf_mul(term, term, z);

        const unsigned long weight = 13591409UL + 545140134UL * np1;
        mpf_mul_ui(tmp, term, weight);
        mpf_add(sum, sum, tmp);
    }

    mpf_mul(inv_pi, C, sum);
    mpf_div_ui(inv_pi, inv_pi, 90856689UL);
    mpf_ui_div(pi, 1, inv_pi);

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
    if (out_path) fclose(fp);

    if (!quiet) {
        fprintf(stderr,
                "EXPLICIT_SEED digits=%llu terms=%lu bits=%lu\n",
                (unsigned long long)digits,
                terms,
                (unsigned long)bits);
    }

    mpf_clear(z);
    mpf_clear(C);
    mpf_clear(tmp);
    mpf_clear(term);
    mpf_clear(sum);
    mpf_clear(inv_pi);
    mpf_clear(pi);
    mpf_clear(den);
    return 0;
}
