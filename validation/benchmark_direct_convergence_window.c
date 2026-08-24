#define main original_d2_direct_benchmark_main
#include "benchmark_original_d2_basic_vs_chud.c"
#undef main

/* Convergence-Windowed Direct Recurrence (CWDR).
 * Keeps the original scalar/direct evaluator, adds tapered precision and
 * short magnitude-window accumulation. No binary splitting is used here.
 */
static mp_bitcnt_t cw_active_bits(mp_bitcnt_t max_bits,
                                  long double digits_per_term,
                                  unsigned long n) {
    long double dropped = digits_per_term * (long double)n * 3.32192809488736234787L;
    long double keep = (long double)max_bits - dropped + 192.0L;
    if (keep < 320.0L) keep = 320.0L;
    if (keep > (long double)max_bits) keep = (long double)max_bits;
    return (mp_bitcnt_t)ceill(keep);
}

static mp_bitcnt_t cw_window_bits(long double digits_per_term,
                                  unsigned long block_len,
                                  mp_bitcnt_t max_bits) {
    long double span = digits_per_term * (long double)(block_len + 2UL) *
                       3.32192809488736234787L + 320.0L;
    if (span > (long double)max_bits) span = (long double)max_bits;
    if (span < 512.0L) span = 512.0L;
    return (mp_bitcnt_t)ceill(span);
}

static void cw_series(const mpf_t z,
                      long double digits_per_term,
                      unsigned long terms,
                      mp_bitcnt_t bits,
                      unsigned long block_terms,
                      mpf_t F,
                      mpf_t T) {
    mpf_t term, zwork, bF, bT, nt;
    mpf_init2(term, bits); mpf_init2(zwork, bits);
    mpf_init2(bF, bits);   mpf_init2(bT, bits); mpf_init2(nt, bits);
    mpf_set(zwork, z);
    mpf_set_ui(term, 1UL); mpf_set_ui(F, 1UL); mpf_set_ui(T, 0UL);

    for (unsigned long a = 1UL; a < terms; a += block_terms) {
        unsigned long b = a + block_terms;
        if (b > terms) b = terms;
        mp_bitcnt_t wbits = cw_window_bits(digits_per_term, b - a, bits);
        mpf_set_prec_raw(bF, wbits); mpf_set_prec_raw(bT, wbits);
        mpf_set_prec_raw(nt, wbits);
        mpf_set_ui(bF, 0UL); mpf_set_ui(bT, 0UL);

        for (unsigned long n = a; n < b; ++n) {
            mp_bitcnt_t abits = cw_active_bits(bits, digits_per_term, n);
            mpf_set_prec_raw(term, abits);
            mpf_set_prec_raw(zwork, abits);
            hyper_term_advance(term, zwork, n);
            mpf_add(bF, bF, term);
            mpf_mul_ui(nt, term, n);
            mpf_add(bT, bT, nt);
        }
        mpf_add(F, F, bF);
        mpf_add(T, T, bT);
    }

    mpf_set_prec_raw(term, bits); mpf_set_prec_raw(zwork, bits);
    mpf_set_prec_raw(bF, bits); mpf_set_prec_raw(bT, bits); mpf_set_prec_raw(nt, bits);
    mpf_clears(term, zwork, bF, bT, nt, NULL);
}

static int d2_cw_pi(uint64_t digits, uint64_t guard, unsigned long block_terms,
                    const char *path, double *seconds, unsigned long *terms_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    mpf_t y,v,z,alpha,beta,K,F,T,tmp,inv_pi,pi;
    mpf_init2(y,bits); mpf_init2(v,bits); mpf_init2(z,bits); mpf_init2(alpha,bits);
    mpf_init2(beta,bits); mpf_init2(K,bits); mpf_init2(F,bits); mpf_init2(T,bits);
    mpf_init2(tmp,bits); mpf_init2(inv_pi,bits); mpf_init2(pi,bits);
    if (!tiny_quadratic_root(D2_H1,D2_H2,y,bits) ||
        !tiny_quadratic_root(D2_L1,D2_L2,v,bits)) return 0;
    mpf_mul_ui(z,y,1728UL);
    mpf_ui_sub(tmp,1UL,z); mpf_mul(tmp,tmp,v); mpf_mul_ui(tmp,tmp,427UL);
    mpf_div(alpha,y,tmp); mpf_ui_sub(beta,1UL,alpha);
    mpf_set_ui(K,427UL); mpf_sqrt(K,K); mpf_div_ui(K,K,6UL);
    mpf_ui_sub(tmp,1UL,z); mpf_sqrt(tmp,tmp); mpf_mul(K,K,tmp);
    const long double q = 24.95589965765426673091470594098813130578L;
    unsigned long terms=(unsigned long)ceill(((long double)digits+(long double)guard+20.0L)/q)+2UL;
    cw_series(z,q,terms,bits,block_terms,F,T);
    mpf_mul(inv_pi,beta,F); mpf_mul_ui(tmp,T,6UL); mpf_add(inv_pi,inv_pi,tmp);
    mpf_mul(inv_pi,inv_pi,K); mpf_ui_div(pi,1UL,inv_pi);
    int ok=write_pi(path,pi,digits);
    if(seconds)*seconds=now_seconds()-t0; if(terms_out)*terms_out=terms;
    mpf_clears(y,v,z,alpha,beta,K,F,T,tmp,inv_pi,pi,NULL);
    return ok;
}

static int chud_cw_pi(uint64_t digits, uint64_t guard, unsigned long block_terms,
                      const char *path, double *seconds, unsigned long *terms_out) {
    double t0=now_seconds(); mp_bitcnt_t bits=rj_bits_for_decimal(digits,guard);
    mpf_t z,alpha,beta,K,F,T,tmp,den,inv_pi,pi;
    mpf_init2(z,bits); mpf_init2(alpha,bits); mpf_init2(beta,bits); mpf_init2(K,bits);
    mpf_init2(F,bits); mpf_init2(T,bits); mpf_init2(tmp,bits); mpf_init2(den,bits);
    mpf_init2(inv_pi,bits); mpf_init2(pi,bits);
    mpf_set_si(z,-1); mpf_set_ui(den,53360UL); mpf_pow_ui(den,den,3UL); mpf_div(z,z,den);
    mpf_set_ui(alpha,77265280UL); mpf_div_ui(alpha,alpha,90856689UL); mpf_ui_sub(beta,1UL,alpha);
    mpf_set_ui(K,163UL); mpf_sqrt(K,K); mpf_div_ui(K,K,6UL);
    mpf_ui_sub(tmp,1UL,z); mpf_sqrt(tmp,tmp); mpf_mul(K,K,tmp);
    const long double q=14.1816474627254776555255216782L;
    unsigned long terms=(unsigned long)ceill(((long double)digits+(long double)guard+20.0L)/q)+2UL;
    cw_series(z,q,terms,bits,block_terms,F,T);
    mpf_mul(inv_pi,beta,F); mpf_mul_ui(tmp,T,6UL); mpf_add(inv_pi,inv_pi,tmp);
    mpf_mul(inv_pi,inv_pi,K); mpf_ui_div(pi,1UL,inv_pi);
    int ok=write_pi(path,pi,digits);
    if(seconds)*seconds=now_seconds()-t0; if(terms_out)*terms_out=terms;
    mpf_clears(z,alpha,beta,K,F,T,tmp,den,inv_pi,pi,NULL);
    return ok;
}

static int run_cw(uint64_t digits,unsigned long block_terms){
    const uint64_t guard=192ULL; char ref[256],db[256],dcw[256],cb[256],ccw[256];
    snprintf(ref,sizeof(ref),"build/cw_ref_%llu_%lu.txt",(unsigned long long)digits,block_terms);
    snprintf(db,sizeof(db),"build/cw_d2_base_%llu_%lu.txt",(unsigned long long)digits,block_terms);
    snprintf(dcw,sizeof(dcw),"build/cw_d2_%llu_%lu.txt",(unsigned long long)digits,block_terms);
    snprintf(cb,sizeof(cb),"build/cw_chud_base_%llu_%lu.txt",(unsigned long long)digits,block_terms);
    snprintf(ccw,sizeof(ccw),"build/cw_chud_%llu_%lu.txt",(unsigned long long)digits,block_terms);
    double tr=0,tdb=0,tdcw=0,tcb=0,tccw=0; unsigned long nr=0,ndb=0,ndcw=0,ncb=0,nccw=0;
    if(!chud_bs_pi(digits,guard,ref,&tr,&nr)||!d2_basic_pi(digits,guard,db,&tdb,&ndb)||
       !d2_cw_pi(digits,guard,block_terms,dcw,&tdcw,&ndcw)||
       !chud_direct_pi(digits,guard,cb,&tcb,&ncb)||!chud_cw_pi(digits,guard,block_terms,ccw,&tccw,&nccw))return 0;
    int e1=files_equal(ref,db),e2=files_equal(ref,dcw),e3=files_equal(ref,cb),e4=files_equal(ref,ccw);
    printf("%llu,%lu,%.9f,%.9f,%.9f,%.9f,%.9f,%.6f,%.6f,%d,%d,%d,%d\n",
      (unsigned long long)digits,block_terms,tr,tdb,tdcw,tcb,tccw,tdcw/tdb,tdcw/tccw,e1,e2,e3,e4);
    fflush(stdout); return e1&&e2&&e3&&e4;
}

int main(int argc,char**argv){
    if(argc!=3){fprintf(stderr,"usage: %s DIGITS BLOCK_TERMS\n",argv[0]);return 2;}
    uint64_t digits=strtoull(argv[1],NULL,10); unsigned long b=strtoul(argv[2],NULL,10);
    printf("digits,block_terms,chud_bs_s,d2_direct_s,d2_cw_s,chud_direct_s,chud_cw_s,d2_cw_over_d2_direct,d2_cw_over_chud_cw,eq_d2_base,eq_d2_cw,eq_chud_base,eq_chud_cw\n");
    return (digits&&b&&run_cw(digits,b))?0:1;
}
