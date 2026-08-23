int rj_z_polynomial_from_phi(const RjPoly *phi,RjPoly *out){
    int d=phi->degree_x;if(phi->degree_y!=d)return 0;memset(out,0,sizeof(*out));out->degree_x=d;out->degree_y=d;out->nterms=phi->nterms;out->terms=(RjTerm*)calloc(out->nterms,sizeof(RjTerm));if(!out->terms)return 0;
    mpz_t pow1728,g,ab;mpz_init(pow1728);mpz_init_set_ui(g,0);mpz_init(ab);
    for(size_t k=0;k<phi->nterms;k++){
        out->terms[k].ix=d-phi->terms[k].ix;out->terms[k].iy=d-phi->terms[k].iy;mpz_init(out->terms[k].coeff);mpz_ui_pow_ui(pow1728,1728,(unsigned long)(phi->terms[k].ix+phi->terms[k].iy));mpz_mul(out->terms[k].coeff,phi->terms[k].coeff,pow1728);mpz_abs(ab,out->terms[k].coeff);if(mpz_sgn(ab)!=0){if(mpz_sgn(g)==0)mpz_set(g,ab);else mpz_gcd(g,g,ab);}
    }
    if(mpz_sgn(g)==0){mpz_clear(pow1728);mpz_clear(g);mpz_clear(ab);rj_poly_clear(out);return 0;}for(size_t k=0;k<out->nterms;k++)mpz_divexact(out->terms[k].coeff,out->terms[k].coeff,g);qsort(out->terms,out->nterms,sizeof(RjTerm),poly_term_cmp_desc);mpz_clear(pow1728);mpz_clear(g);mpz_clear(ab);return 1;
}

int rj_poly_equal(const RjPoly *a,const RjPoly *b){
    if(a->degree_x!=b->degree_x||a->degree_y!=b->degree_y||a->nterms!=b->nterms)return 0;
    for (size_t k = 0; k < a->nterms; k++) {
        if (a->terms[k].ix != b->terms[k].ix ||
            a->terms[k].iy != b->terms[k].iy ||
            mpz_cmp(a->terms[k].coeff, b->terms[k].coeff) != 0)
            return 0;
    }
    return 1;
}

int rj_verify_phi_qseries_exact(const RjPoly *phi, int p, int hi,
                                int *out_lo, size_t *out_nonzero) {
    int dx=phi->degree_x,dy=phi->degree_y;
    int lo=-(dx+p*dy)-5;
    int internal_hi=hi+dx+p*dy+30;
    QSeries J={0},Jsmall={0},K={0},total={0};
    if(!qs_init(&J,-1,internal_hi+10)||!qs_init(&Jsmall,-1,internal_hi/p+20)||!qs_init(&K,-p,internal_hi)||!qs_init(&total,lo,hi))return 0;
    if(!j_series_exact(internal_hi+10,&J)||!j_series_exact(internal_hi/p+20,&Jsmall)){qs_clear(&J);qs_clear(&Jsmall);qs_clear(&K);qs_clear(&total);return 0;}
    qs_qsub(&K,&Jsmall,p);
    QSeries *Jp=(QSeries*)calloc((size_t)dx+1,sizeof(QSeries));QSeries *Kp=(QSeries*)calloc((size_t)dy+1,sizeof(QSeries));
    if(!Jp||!Kp)return 0;
    for(int i=0;i<=dx;i++)qs_init(&Jp[i],lo,internal_hi);
    for(int j=0;j<=dy;j++)qs_init(&Kp[j],lo,internal_hi);
    mpq_set_ui(qs_at(&Jp[0],0),1,1);mpq_set_ui(qs_at(&Kp[0],0),1,1);
    for(int i=1;i<=dx;i++)qs_conv(&Jp[i],&Jp[i-1],&J);
    for(int j=1;j<=dy;j++)qs_conv(&Kp[j],&Kp[j-1],&K);
    QSeries term={0};qs_init(&term,lo,hi);mpq_t c,t;mpq_init(c);mpq_init(t);
    for(size_t k=0;k<phi->nterms;k++){
        int i=phi->terms[k].ix,j=phi->terms[k].iy;qs_conv(&term,&Jp[i],&Kp[j]);mpq_set_z(c,phi->terms[k].coeff);
        for(int e=lo;e<=hi;e++){mpq_mul(t,qs_at(&term,e),c);mpq_add(qs_at(&total,e),qs_at(&total,e),t);}
    }
    size_t nz=0;for(int e=lo;e<=hi;e++)if(mpq_sgn(qs_at(&total,e))!=0)nz++;
    if (out_lo) *out_lo = lo;
    if (out_nonzero) *out_nonzero = nz;
    mpq_clear(c);mpq_clear(t);qs_clear(&term);
    for (int i = 0; i <= dx; i++) qs_clear(&Jp[i]);
    for (int j = 0; j <= dy; j++) qs_clear(&Kp[j]);
    free(Jp);
    free(Kp);
    qs_clear(&J);qs_clear(&Jsmall);qs_clear(&K);qs_clear(&total);return nz==0;
}

mp_bitcnt_t rj_bits_for_decimal(uint64_t decimal_digits, uint64_t guard_digits) {
    long double d=(long double)(decimal_digits+guard_digits+16ULL);
    long double bits=ceill(d*3.32192809488736234787L)+64.0L;
    if(bits>(long double)(~(mp_bitcnt_t)0)) return (mp_bitcnt_t)(~(mp_bitcnt_t)0);
    return (mp_bitcnt_t)bits;
}

int rj_write_decimal_mpf(FILE *fp, const mpf_t x, size_t digits_after_decimal) {
    if(!fp)return 0;
    size_t sig=digits_after_decimal+1;
    mp_exp_t exp10=0;
    char *mant=mpf_get_str(NULL,&exp10,10,sig,x);
    if(!mant)return 0;
    int neg=(mant[0]=='-');char *m=mant+neg;size_t len=strlen(m);
    if(neg)fputc('-',fp);
    if(exp10<=0){fputs("0.",fp);for(mp_exp_t i=0;i<-exp10;i++)fputc('0',fp);fwrite(m,1,len,fp);}
    else {
        size_t e=(size_t)exp10;
        if(e>=len){fwrite(m,1,len,fp);for(size_t i=len;i<e;i++)fputc('0',fp);if(digits_after_decimal)fputc('.',fp);}
        else {fwrite(m,1,e,fp);if(digits_after_decimal)fputc('.',fp);fwrite(m+e,1,len-e,fp);}
    }
    /* mpf_get_str may omit trailing zeros; pad fractional tail deterministically. */
    if(exp10>0){size_t produced_frac = len>(size_t)exp10 ? len-(size_t)exp10 : 0;for(size_t i=produced_frac;i<digits_after_decimal;i++)fputc('0',fp);} 
    free(mant);return 1;
}
