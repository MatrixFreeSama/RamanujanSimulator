#define main original_d2_benchmark_main
#include "benchmark_original_d2_basic_vs_chud.c"
#undef main

/* Classical iterative references only.  The bare direct D2/Chud series are
 * compiled only to reuse common helpers and are never executed here. */

static unsigned long ceil_log2_u64_fast(uint64_t n) {
    unsigned long k = 0UL;
    uint64_t x = 1ULL;
    while (x < n && k < 63UL) {
        x <<= 1;
        ++k;
    }
    return k;
}

static int agm_pi_only(uint64_t digits, uint64_t guard, const char *path,
                       double *seconds, unsigned long *iterations_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    mpf_t a,b,t,p,an,bn,delta,tmp,sum,pi;
    mpf_init2(a,bits); mpf_init2(b,bits); mpf_init2(t,bits);
    mpf_init2(p,bits); mpf_init2(an,bits); mpf_init2(bn,bits);
    mpf_init2(delta,bits); mpf_init2(tmp,bits); mpf_init2(sum,bits);
    mpf_init2(pi,bits);

    mpf_set_ui(a,1UL);
    mpf_set_ui(tmp,2UL); mpf_sqrt(tmp,tmp); mpf_ui_div(b,1UL,tmp);
    mpf_set_ui(t,1UL); mpf_div_ui(t,t,4UL);
    mpf_set_ui(p,1UL);

    unsigned long iterations = ceil_log2_u64_fast(digits + guard + 1ULL) + 4UL;
    for (unsigned long i=0UL;i<iterations;++i) {
        mpf_add(an,a,b); mpf_div_ui(an,an,2UL);
        mpf_mul(tmp,a,b); mpf_sqrt(bn,tmp);
        mpf_sub(delta,a,an); mpf_mul(delta,delta,delta);
        mpf_mul(tmp,p,delta); mpf_sub(t,t,tmp);
        mpf_mul_ui(p,p,2UL);
        mpf_set(a,an); mpf_set(b,bn);
    }
    mpf_add(sum,a,b); mpf_mul(sum,sum,sum);
    mpf_mul_ui(tmp,t,4UL); mpf_div(pi,sum,tmp);

    int ok = write_pi(path,pi,digits);
    if (seconds) *seconds = now_seconds()-t0;
    if (iterations_out) *iterations_out=iterations;
    mpf_clears(a,b,t,p,an,bn,delta,tmp,sum,pi,NULL);
    return ok;
}

static int borwein4_pi_only(uint64_t digits, uint64_t guard, const char *path,
                            double *seconds, unsigned long *iterations_out) {
    double t0=now_seconds();
    mp_bitcnt_t bits=rj_bits_for_decimal(digits,guard);
    mpf_t sqrt2,y,yn,a,y2,y4,root,num,den,one,oneplus,oneplus2,oneplus4,
          poly,tmp,factor,pi;
    mpf_init2(sqrt2,bits); mpf_init2(y,bits); mpf_init2(yn,bits);
    mpf_init2(a,bits); mpf_init2(y2,bits); mpf_init2(y4,bits);
    mpf_init2(root,bits); mpf_init2(num,bits); mpf_init2(den,bits);
    mpf_init2(one,bits); mpf_init2(oneplus,bits); mpf_init2(oneplus2,bits);
    mpf_init2(oneplus4,bits); mpf_init2(poly,bits); mpf_init2(tmp,bits);
    mpf_init2(factor,bits); mpf_init2(pi,bits);

    mpf_set_ui(one,1UL);
    mpf_set_ui(sqrt2,2UL); mpf_sqrt(sqrt2,sqrt2);
    mpf_sub_ui(y,sqrt2,1UL);
    mpf_mul_ui(a,sqrt2,4UL); mpf_ui_sub(a,6UL,a);

    unsigned long log2n=ceil_log2_u64_fast(digits+guard+1ULL);
    unsigned long iterations=(log2n+1UL)/2UL+3UL;
    unsigned long factor_ui=8UL;
    for (unsigned long i=0UL;i<iterations;++i) {
        mpf_mul(y2,y,y); mpf_mul(y4,y2,y2);
        mpf_ui_sub(root,1UL,y4); mpf_sqrt(root,root); mpf_sqrt(root,root);
        mpf_ui_sub(num,1UL,root); mpf_add_ui(den,root,1UL); mpf_div(yn,num,den);
        mpf_add_ui(oneplus,yn,1UL);
        mpf_mul(oneplus2,oneplus,oneplus); mpf_mul(oneplus4,oneplus2,oneplus2);
        mpf_mul(y2,yn,yn); mpf_add(poly,one,yn); mpf_add(poly,poly,y2);
        mpf_mul(poly,poly,yn);
        mpf_set_ui(factor,factor_ui); mpf_mul(poly,poly,factor);
        mpf_mul(tmp,a,oneplus4); mpf_sub(a,tmp,poly);
        mpf_set(y,yn);
        if (factor_ui<=ULONG_MAX/4UL) factor_ui*=4UL;
    }
    mpf_ui_div(pi,1UL,a);

    int ok=write_pi(path,pi,digits);
    if (seconds) *seconds=now_seconds()-t0;
    if (iterations_out) *iterations_out=iterations;
    mpf_clears(sqrt2,y,yn,a,y2,y4,root,num,den,one,oneplus,oneplus2,
               oneplus4,poly,tmp,factor,pi,NULL);
    return ok;
}

static int run_classic_only(uint64_t digits) {
    const uint64_t guard=192ULL;
    char ref[256],agmp[256],bqp[256];
    snprintf(ref,sizeof(ref),"build/d2bs_ref_%llu.txt",(unsigned long long)digits);
    snprintf(agmp,sizeof(agmp),"build/d2bs_agm_%llu.txt",(unsigned long long)digits);
    snprintf(bqp,sizeof(bqp),"build/d2bs_bq_%llu.txt",(unsigned long long)digits);
    double ta=0.0,tb=0.0; unsigned long ia=0UL,ib=0UL;
    if (!agm_pi_only(digits,guard,agmp,&ta,&ia) ||
        !borwein4_pi_only(digits,guard,bqp,&tb,&ib)) return 0;
    int eqa=files_equal(ref,agmp), eqb=files_equal(ref,bqp);
    printf("%llu,%.9f,%lu,%.9f,%lu,%d,%d\n",
           (unsigned long long)digits,ta,ia,tb,ib,eqa,eqb);
    fflush(stdout);
    return eqa && eqb;
}

int main(int argc,char **argv) {
    printf("digits,agm_seconds,agm_iterations,borwein4_seconds,borwein4_iterations,equal_agm,equal_borwein4\n");
    if (argc<=1) return 1;
    for (int i=1;i<argc;++i) {
        uint64_t digits=strtoull(argv[i],NULL,10);
        if (!digits || !run_classic_only(digits)) return 1;
    }
    return 0;
}
