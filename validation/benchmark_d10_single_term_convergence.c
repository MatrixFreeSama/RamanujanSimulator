#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define DEG_H 10
#define DEG_C 10
#define DEG_F 30

static const char *H_ASC[11] = {
"-370379678802534904801278783784712516669511657545292686541317920733750417605390874497797114616248919946333551664716964637548220313449661671973306049203399881120521650176000000000000000000000000000000",
"-16384491475280569857609606791335695552666625484765086387174938658663813020359754659184361553122188956373339449603305081745136824483583488172583077501846833379783380828160000000000000000000000000000",
"-217350279991328190761985550176942873049694056194412714001325339136414424489817793868040748634553400591278443383092632547675780889787677168496425167768288712702245816238080000000000000000000000000",
"418369010729471937388837995958372278400439494046430042571202929122465271027869228122325714946016657143974277552123652771875407817385041791004748472603827319008838266388480000000000000000000000",
"-284990521388995246316162827656061115574939004181220209497083465422269713552869415073979491221639144085206093288268428426324235462148160575299769337129153682332826491420672000000000000000000",
"70847587572824648936525224372543186878973589986261303386779554786605638312391735378070425335097090463401517210305516918216139834018597205809483527019979311359260098560000000000000000000",
"-257203233837035525461655402777787102338440698835254775521328981264810710382429441484236329698930431430416431736284956301679573229169891425668900992511329873883243216896000000000000",
"3095576677826745001791647813730045774163789379693210329707194097504345347048291101401257306468618558671368494783528572132060746016311647242448231212195585331759677440000000000",
"-373347885499824232804010412882330497087294334996397135698951537102944942755876372620309978804222318847935128447068247937919760067465510041830754488642742834429952000000",
"33685656633958669107888967803661931714647776290935715365741200718405786535119228031508290986188379856550182868818065771167002034269118770355222717680986890240000",
"1"
};

static const char *C_ASC[11] = {
"-718150912949642523148693959180027408505991659454736629760000000000",
"148269606635177168317423665367562949221038833187015360512000000000",
"-63536734918606054543228298377769047547857440465492954316800000000",
"421486683565974884866530269709770250139307349495722803200000000",
"-211565635591557038546781996408514792373254059354608369664000000",
"33483032339905697583284254695534390058767103248033329971200000",
"1234312959090914569645822713213259492862969867888586915840000",
"17531495856827170677348882327567868131596697367321055232000",
"53975830800870124552573137895158808579232875400442163200",
"322959704590412332725686274638226069771593878764442080",
"1"
};

static int is_squarefree(unsigned n) {
    if (n == 0U) return 0;
    for (unsigned p = 2U; (uint64_t)p * (uint64_t)p <= n; ++p) {
        unsigned pp = p * p;
        if (n % pp == 0U) return 0;
    }
    return 1;
}

static int is_fundamental_negative_abs(unsigned d) {
    if ((d & 3U) == 3U) return is_squarefree(d);
    if ((d & 3U) == 0U) {
        unsigned n = d / 4U;
        unsigned r = n & 3U;
        return (r == 1U || r == 2U) && is_squarefree(n);
    }
    return 0;
}

static unsigned gcd_u(unsigned a, unsigned b) {
    while (b) { unsigned t = a % b; a = b; b = t; }
    return a;
}

static unsigned class_number_negative_discriminant(int D) {
    unsigned d = (unsigned)(-D);
    unsigned max_a = (unsigned)(sqrt((double)d / 3.0) + 2.0);
    unsigned count = 0U;
    for (unsigned A = 1U; A <= max_a; ++A) {
        for (int B = -(int)A; B <= (int)A; ++B) {
            long long num = (long long)B * (long long)B - (long long)D;
            long long den = 4LL * (long long)A;
            if (num % den != 0LL) continue;
            long long C = num / den;
            if ((long long)A > C) continue;
            unsigned g = gcd_u(A, (unsigned)llabs((long long)B));
            g = gcd_u(g, (unsigned)C);
            if (g != 1U) continue;
            if ((((unsigned)abs(B) == A) || ((long long)A == C)) && B < 0) continue;
            ++count;
        }
    }
    return count;
}

static int verify_exact_factorization(void) {
    mpz_t H[11], P[11], F[31], q, tmp;
    for (int i = 0; i <= DEG_H; ++i) { mpz_init(H[i]); mpz_set_str(H[i], H_ASC[i], 10); }
    for (int i = 0; i <= DEG_C; ++i) { mpz_init(P[i]); mpz_set_str(P[i], C_ASC[i], 10); }
    for (int i = 0; i <= DEG_F; ++i) mpz_init(F[i]);
    mpz_init(q); mpz_init(tmp);

    for (int i = 0; i <= DEG_H; ++i) mpz_set(F[3 * i], H[i]);

    int ok = mpz_cmp_ui(P[DEG_C], 1U) == 0;
    for (int d = DEG_F; ok && d >= DEG_C; --d) {
        if (mpz_sgn(F[d]) == 0) continue;
        mpz_set(q, F[d]);
        int shift = d - DEG_C;
        for (int i = 0; i <= DEG_C; ++i) {
            mpz_mul(tmp, q, P[i]);
            mpz_sub(F[shift + i], F[shift + i], tmp);
        }
    }
    for (int i = 0; ok && i <= DEG_F; ++i) if (mpz_sgn(F[i]) != 0) ok = 0;

    for (int i = 0; i <= DEG_H; ++i) mpz_clear(H[i]);
    for (int i = 0; i <= DEG_C; ++i) mpz_clear(P[i]);
    for (int i = 0; i <= DEG_F; ++i) mpz_clear(F[i]);
    mpz_clear(q); mpz_clear(tmp);
    return ok;
}

static void eval_poly_and_derivative(mpf_t p, mpf_t dp, const mpf_t x, mpz_t P[11], mpf_t ztmp) {
    mpf_set_z(p, P[DEG_C]);
    mpf_set_ui(dp, 0U);
    for (int i = DEG_C - 1; i >= 0; --i) {
        mpf_mul(dp, dp, x);
        mpf_add(dp, dp, p);
        mpf_mul(p, p, x);
        mpf_set_z(ztmp, P[i]);
        mpf_add(p, p, ztmp);
    }
}

static double principal_c_abs(void) {
    mpz_t P[11];
    for (int i = 0; i <= DEG_C; ++i) { mpz_init(P[i]); mpz_set_str(P[i], C_ASC[i], 10); }
    mpf_set_default_prec(1024U);
    mpf_t x, p, dp, step, ztmp;
    mpf_inits(x, p, dp, step, ztmp, NULL);
    mpf_set_str(x, "-322959704590412332725686274638226069771593878764440000", 10);
    for (int it = 0; it < 24; ++it) {
        eval_poly_and_derivative(p, dp, x, P, ztmp);
        mpf_div(step, p, dp);
        mpf_sub(x, x, step);
    }
    double v = fabs(mpf_get_d(x));
    mpf_clears(x, p, dp, step, ztmp, NULL);
    for (int i = 0; i <= DEG_C; ++i) mpz_clear(P[i]);
    return v;
}

int main(void) {
    const unsigned target_d = 13843U;
    unsigned h = class_number_negative_discriminant(-(int)target_d);
    unsigned maxd[11] = {0U};
    unsigned counts[11] = {0U};
    for (unsigned d = 3U; d <= target_d; ++d) {
        if (!is_fundamental_negative_abs(d)) continue;
        unsigned ch = class_number_negative_discriminant(-(int)d);
        if (ch >= 1U && ch <= 10U) { ++counts[ch]; if (d > maxd[ch]) maxd[ch] = d; }
    }

    int factor_ok = verify_exact_factorization();
    double cabs = principal_c_abs();
    double digits_per_term = 3.0 * log10(cabs) - log10(1728.0);

    printf("D10 convergence certificate\n");
    printf("discriminant=-%u class_number=%u reduced_forms=%u\n", target_d, h, h);
    printf("exact_factor_H_of_x3_by_C_polynomial=%s\n", factor_ok ? "PASS" : "FAIL");
    printf("principal_abs_C=%.17e\n", cabs);
    printf("asymptotic_digits_per_term=%.15f\n", digits_per_term);
    printf("terms_for_100k=%llu\n", (unsigned long long)ceil(100000.0 / digits_per_term));
    printf("terms_for_200k=%llu\n", (unsigned long long)ceil(200000.0 / digits_per_term));
    printf("terms_for_1M=%llu\n", (unsigned long long)ceil(1000000.0 / digits_per_term));
    printf("\nmax fundamental |D| by class number through h=10:\n");
    for (unsigned k = 1U; k <= 10U; ++k)
        printf("h=%u maxD=%u count=%u\n", k, maxd[k], counts[k]);

    if (h != 10U || maxd[10] != 13843U || counts[10] != 87U || !factor_ok) return 1;
    if (!(digits_per_term > 157.28 && digits_per_term < 157.30)) return 1;
    return 0;
}
