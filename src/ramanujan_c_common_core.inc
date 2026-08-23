const int RJ_SEQUENCE_10K4[32] = {
    2,2,2,2,5,5,5,5,
    2,2,2,2,5,5,5,5,
    2,2,2,2,5,5,5,5,
    2,2,2,2,5,5,5,5
};
const size_t RJ_SEQUENCE_10K4_LEN = 32;

/* Exact modular polynomials Phi_p(X,Y), frozen from the exact q-series generator. */
static const RjConstTerm PHI2[] = {
    {3,0,"1"},{2,2,"-1"},{2,1,"1488"},{2,0,"-162000"},
    {1,2,"1488"},{1,1,"40773375"},{1,0,"8748000000"},
    {0,3,"1"},{0,2,"-162000"},{0,1,"8748000000"},{0,0,"-157464000000000"}
};
static const RjConstTerm PHI3[] = {
    {4,0,"1"},{3,3,"-1"},{3,2,"2232"},{3,1,"-1069956"},{3,0,"36864000"},
    {2,3,"2232"},{2,2,"2587918086"},{2,1,"8900222976000"},{2,0,"452984832000000"},
    {1,3,"-1069956"},{1,2,"8900222976000"},{1,1,"-770845966336000000"},
    {1,0,"1855425871872000000000"},{0,4,"1"},{0,3,"36864000"},
    {0,2,"452984832000000"},{0,1,"1855425871872000000000"}
};
static const RjConstTerm PHI5[] = {
    {6,0,"1"},{5,5,"-1"},{5,4,"3720"},{5,3,"-4550940"},{5,2,"2028551200"},
    {5,1,"-246683410950"},{5,0,"1963211489280"},{4,5,"3720"},
    {4,4,"1665999364600"},{4,3,"107878928185336800"},{4,2,"383083609779811215375"},
    {4,1,"128541798906828816384000"},{4,0,"1284733132841424456253440"},
    {3,5,"-4550940"},{3,4,"107878928185336800"},{3,3,"-441206965512914835246100"},
    {3,2,"26898488858380731577417728000"},{3,1,"-192457934618928299655108231168000"},
    {3,0,"280244777828439527804321565297868800"},{2,5,"2028551200"},
    {2,4,"383083609779811215375"},{2,3,"26898488858380731577417728000"},
    {2,2,"5110941777552418083110765199360000"},
    {2,1,"36554736583949629295706472332656640000"},
    {2,0,"6692500042627997708487149415015068467200"},{1,5,"-246683410950"},
    {1,4,"128541798906828816384000"},{1,3,"-192457934618928299655108231168000"},
    {1,2,"36554736583949629295706472332656640000"},
    {1,1,"-264073457076620596259715790247978782949376"},
    {1,0,"53274330803424425450420160273356509151232000"},{0,6,"1"},
    {0,5,"1963211489280"},{0,4,"1284733132841424456253440"},
    {0,3,"280244777828439527804321565297868800"},
    {0,2,"6692500042627997708487149415015068467200"},
    {0,1,"53274330803424425450420160273356509151232000"},
    {0,0,"141359947154721358697753474691071362751004672000"}
};

/* z-polynomials P_p(x,y) after X=1728/x, Y=1728/y and primitive normalization. */
static const RjConstTerm P2[] = {
    {3,3,"-1953125"},{3,2,"187500"},{3,1,"-6000"},{3,0,"64"},
    {2,3,"187500"},{2,2,"1510125"},{2,1,"95232"},
    {1,3,"-6000"},{1,2,"95232"},{1,1,"-110592"},{0,3,"64"}
};
static const RjConstTerm P3[] = {
    {4,3,"262144000000000"},{4,2,"110592000000"},{4,1,"15552000"},{4,0,"729"},
    {3,4,"262144000000000"},{3,3,"-188194816000000"},{3,2,"3754781568000"},
    {3,1,"-779997924"},{2,4,"110592000000"},{2,3,"3754781568000"},
    {2,2,"1886592284694"},{2,1,"2811677184"},{1,4,"15552000"},
    {1,3,"-779997924"},{1,2,"2811677184"},{1,1,"-2176782336"},{0,4,"729"}
};
static const RjConstTerm P5[] = {
    {6,6,"5309626171273360722362368000"},{6,5,"3457795560648760910413824000"},
    {6,4,"750608416927050074633011200"},{6,3,"54313315434020926285414400"},
    {6,2,"430254526762844160"},{6,1,"1136117760"},{6,0,"1"},
    {5,6,"3457795560648760910413824000"},{5,5,"-29617595563122405481849552896"},
    {5,4,"7084552847250663218872320000"},{5,3,"-64453772899964735127552000"},
    {5,2,"74387615108118528000"},{5,1,"-246683410950"},
    {4,6,"750608416927050074633011200"},{4,5,"7084552847250663218872320000"},
    {4,4,"1711644060233550509015040000"},{4,3,"15566255126377738181376000"},
    {4,2,"383083609779811215375"},{4,1,"3505336473600"},
    {3,6,"54313315434020926285414400"},{3,5,"-64453772899964735127552000"},
    {3,4,"15566255126377738181376000"},{3,3,"-441206965512914835246100"},
    {3,2,"186414787904261990400"},{3,1,"-13589034024960"},
    {2,6,"430254526762844160"},{2,5,"74387615108118528000"},
    {2,4,"383083609779811215375"},{2,3,"186414787904261990400"},
    {2,2,"4974647446705766400"},{2,1,"19194382909440"},
    {1,6,"1136117760"},{1,5,"-246683410950"},{1,4,"3505336473600"},
    {1,3,"-13589034024960"},{1,2,"19194382909440"},{1,1,"-8916100448256"},
    {0,6,"1"}
};

int rj_poly_init_from_const(RjPoly *out, int dx, int dy, const RjConstTerm *src, size_t n) {
    memset(out, 0, sizeof(*out));
    out->degree_x = dx;
    out->degree_y = dy;
    out->nterms = n;
    out->terms = (RjTerm*)calloc(n, sizeof(RjTerm));
    if (!out->terms) return 0;
    for (size_t k=0; k<n; ++k) {
        out->terms[k].ix = src[k].ix;
        out->terms[k].iy = src[k].iy;
        mpz_init(out->terms[k].coeff);
        if (mpz_set_str(out->terms[k].coeff, src[k].coeff_dec, 10) != 0) {
            rj_poly_clear(out);
            return 0;
        }
    }
    return 1;
}

void rj_poly_clear(RjPoly *p) {
    if (!p || !p->terms) return;
    for (size_t k=0; k<p->nterms; ++k) mpz_clear(p->terms[k].coeff);
    free(p->terms);
    memset(p, 0, sizeof(*p));
}

int rj_get_phi_poly(int p, RjPoly *out) {
    switch (p) {
        case 2: return rj_poly_init_from_const(out,3,3,PHI2,sizeof(PHI2)/sizeof(PHI2[0]));
        case 3: return rj_poly_init_from_const(out,4,4,PHI3,sizeof(PHI3)/sizeof(PHI3[0]));
        case 5: return rj_poly_init_from_const(out,6,6,PHI5,sizeof(PHI5)/sizeof(PHI5[0]));
        default: return 0;
    }
}

int rj_get_z_poly(int p, RjPoly *out) {
    switch (p) {
        case 2: return rj_poly_init_from_const(out,3,3,P2,sizeof(P2)/sizeof(P2[0]));
        case 3: return rj_poly_init_from_const(out,4,4,P3,sizeof(P3)/sizeof(P3[0]));
        case 5: return rj_poly_init_from_const(out,6,6,P5,sizeof(P5)/sizeof(P5[0]));
        default: return 0;
    }
}

void rj_state_init(RjState *s, mp_bitcnt_t bits) {
    mpf_init2(s->z,bits); mpf_init2(s->u,bits); mpf_init2(s->v,bits);
}
void rj_state_clear(RjState *s) { mpf_clear(s->z); mpf_clear(s->u); mpf_clear(s->v); }

void rj_initial_state(RjState *s) {
    mpf_t den, alpha, C, tmp, one;
    mp_bitcnt_t bits = mpf_get_prec(s->z);
    mpf_init2(den,bits); mpf_init2(alpha,bits); mpf_init2(C,bits); mpf_init2(tmp,bits); mpf_init2(one,bits);
    mpf_set_ui(one,1);

    mpf_set_ui(den,53360); mpf_pow_ui(den,den,3);
    mpf_ui_div(s->z,1,den); mpf_neg(s->z,s->z);

    mpf_set_ui(alpha,77265280); mpf_div_ui(alpha,alpha,90856689);

    mpf_ui_sub(tmp,1,s->z);      /* 1-z */
    mpf_sqrt(tmp,tmp);
    mpf_set_ui(C,163); mpf_sqrt(C,C); mpf_div_ui(C,C,6); mpf_mul(C,C,tmp);

    mpf_ui_sub(tmp,1,alpha); mpf_mul(s->u,C,tmp);
    mpf_mul_ui(s->v,C,6);

    mpf_clear(den); mpf_clear(alpha); mpf_clear(C); mpf_clear(tmp); mpf_clear(one);
}

void rj_eval_poly_all(const RjPoly *poly, const mpf_t x, const mpf_t y,
                      mpf_t P, mpf_t Px, mpf_t Py,
                      mpf_t Pxx, mpf_t Pxy, mpf_t Pyy) {
    mp_bitcnt_t bits = mpf_get_prec(P);
    int dx=poly->degree_x, dy=poly->degree_y;
    mpf_t *xp=(mpf_t*)malloc((dx+1)*sizeof(mpf_t));
    mpf_t *yp=(mpf_t*)malloc((dy+1)*sizeof(mpf_t));
    for(int i=0;i<=dx;i++) mpf_init2(xp[i],bits);
    for(int j=0;j<=dy;j++) mpf_init2(yp[j],bits);
    mpf_set_ui(xp[0],1); for(int i=1;i<=dx;i++) mpf_mul(xp[i],xp[i-1],x);
    mpf_set_ui(yp[0],1); for(int j=1;j<=dy;j++) mpf_mul(yp[j],yp[j-1],y);
    mpf_set_ui(P,0); mpf_set_ui(Px,0); mpf_set_ui(Py,0); mpf_set_ui(Pxx,0); mpf_set_ui(Pxy,0); mpf_set_ui(Pyy,0);
    mpf_t t,c;
    mpf_init2(t,bits); mpf_init2(c,bits);
    for(size_t k=0;k<poly->nterms;k++) {
        int i=poly->terms[k].ix, j=poly->terms[k].iy;
        mpf_set_z(c,poly->terms[k].coeff);
        mpf_mul(t,c,xp[i]); mpf_mul(t,t,yp[j]); mpf_add(P,P,t);
        if(i) { mpf_mul_ui(t,c,(unsigned long)i); mpf_mul(t,t,xp[i-1]); mpf_mul(t,t,yp[j]); mpf_add(Px,Px,t); }
        if(j) { mpf_mul_ui(t,c,(unsigned long)j); mpf_mul(t,t,xp[i]); mpf_mul(t,t,yp[j-1]); mpf_add(Py,Py,t); }
        if(i>=2) { mpf_mul_ui(t,c,(unsigned long)(i*(i-1))); mpf_mul(t,t,xp[i-2]); mpf_mul(t,t,yp[j]); mpf_add(Pxx,Pxx,t); }
        if(i && j) { mpf_mul_ui(t,c,(unsigned long)(i*j)); mpf_mul(t,t,xp[i-1]); mpf_mul(t,t,yp[j-1]); mpf_add(Pxy,Pxy,t); }
        if(j>=2) { mpf_mul_ui(t,c,(unsigned long)(j*(j-1))); mpf_mul(t,t,xp[i]); mpf_mul(t,t,yp[j-2]); mpf_add(Pyy,Pyy,t); }
    }
    mpf_clear(t); mpf_clear(c);
    for(int i=0;i<=dx;i++) mpf_clear(xp[i]);
    for(int j=0;j<=dy;j++) mpf_clear(yp[j]);
    free(xp); free(yp);
}

int rj_transform_step(const RjPoly *poly, int p, int newton_iters,
                      const RjState *in, RjState *out) {
    mp_bitcnt_t bits=mpf_get_prec(in->z);
    mpf_t y,P,Px,Py,Pxx,Pxy,Pyy,delta,r,s,M,L,tmp,tmp2,den;
    mpf_init2(y,bits); mpf_init2(P,bits); mpf_init2(Px,bits); mpf_init2(Py,bits);
    mpf_init2(Pxx,bits); mpf_init2(Pxy,bits); mpf_init2(Pyy,bits); mpf_init2(delta,bits);
    mpf_init2(r,bits); mpf_init2(s,bits); mpf_init2(M,bits); mpf_init2(L,bits);
    mpf_init2(tmp,bits); mpf_init2(tmp2,bits); mpf_init2(den,bits);

    /* y0=x^p/1728^(p-1) */
    mpf_pow_ui(y,in->z,(unsigned long)p);
    mpf_set_ui(den,1728); mpf_pow_ui(den,den,(unsigned long)(p-1)); mpf_div(y,y,den);

    for(int it=0;it<newton_iters;it++) {
        rj_eval_poly_all(poly,in->z,y,P,Px,Py,Pxx,Pxy,Pyy);
        if (mpf_sgn(Py)==0) goto fail;
        mpf_div(delta,P,Py); mpf_sub(y,y,delta);
    }
    rj_eval_poly_all(poly,in->z,y,P,Px,Py,Pxx,Pxy,Pyy);
    if (mpf_sgn(Py)==0 || mpf_sgn(y)==0) goto fail;

    /* r=-Px/Py */
    mpf_div(r,Px,Py); mpf_neg(r,r);

    /* s=-(Pxx+2*Pxy*r+Pyy*r^2)/Py */
    mpf_mul(tmp,Pxy,r); mpf_mul_ui(tmp,tmp,2); mpf_add(tmp,tmp,Pxx);
    mpf_mul(tmp2,r,r); mpf_mul(tmp2,tmp2,Pyy); mpf_add(tmp,tmp,tmp2);
    mpf_div(s,tmp,Py); mpf_neg(s,s);

    /* M=(1/p)*(x/y)*r*sqrt((1-x)/(1-y)) */
    mpf_div(M,in->z,y); mpf_mul(M,M,r); mpf_div_ui(M,M,(unsigned long)p);
    mpf_ui_sub(tmp,1,in->z); mpf_ui_sub(tmp2,1,y); mpf_div(tmp,tmp,tmp2); mpf_sqrt(tmp,tmp); mpf_mul(M,M,tmp);
    if (mpf_sgn(M)==0 || mpf_sgn(r)==0) goto fail;

    /* L=1-x*r/y+x*s/r-x/(2(1-x))+x*r/(2(1-y)) */
    mpf_set_ui(L,1);
    mpf_mul(tmp,in->z,r); mpf_div(tmp,tmp,y); mpf_sub(L,L,tmp);
    mpf_mul(tmp,in->z,s); mpf_div(tmp,tmp,r); mpf_add(L,L,tmp);
    mpf_ui_sub(tmp2,1,in->z); mpf_mul_ui(tmp2,tmp2,2); mpf_div(tmp,in->z,tmp2); mpf_sub(L,L,tmp);
    mpf_ui_sub(tmp2,1,y); mpf_mul_ui(tmp2,tmp2,2); mpf_mul(tmp,in->z,r); mpf_div(tmp,tmp,tmp2); mpf_add(L,L,tmp);

    /* u'=(u-vL)/M ; v'=v*x*r/(y*M) */
    mpf_mul(tmp,in->v,L); mpf_sub(tmp,in->u,tmp); mpf_div(out->u,tmp,M);
    mpf_mul(tmp,in->v,in->z); mpf_mul(tmp,tmp,r); mpf_mul(tmp2,y,M); mpf_div(out->v,tmp,tmp2);
    mpf_set(out->z,y);

    mpf_clear(y);mpf_clear(P);mpf_clear(Px);mpf_clear(Py);mpf_clear(Pxx);mpf_clear(Pxy);mpf_clear(Pyy);
    mpf_clear(delta);mpf_clear(r);mpf_clear(s);mpf_clear(M);mpf_clear(L);mpf_clear(tmp);mpf_clear(tmp2);mpf_clear(den);
    return 1;
fail:
    mpf_clear(y);mpf_clear(P);mpf_clear(Px);mpf_clear(Py);mpf_clear(Pxx);mpf_clear(Pxy);mpf_clear(Pyy);
    mpf_clear(delta);mpf_clear(r);mpf_clear(s);mpf_clear(M);mpf_clear(L);mpf_clear(tmp);mpf_clear(tmp2);mpf_clear(den);
    return 0;
}

int rj_hyper_F_theta(const mpf_t z, unsigned long working_digits,
                     mpf_t F, mpf_t thetaF, unsigned long max_terms) {
    mp_bitcnt_t bits=mpf_get_prec(z);
    mpf_t term,ratio,tmp,num,den,cutoff,ten;
    mpf_init2(term,bits);mpf_init2(ratio,bits);mpf_init2(tmp,bits);mpf_init2(num,bits);mpf_init2(den,bits);
    mpf_init2(cutoff,bits);mpf_init2(ten,bits);
    mpf_set_ui(term,1); mpf_set(F,term); mpf_set_ui(thetaF,0);
    mpf_set_ui(ten,10);
    unsigned long cd = working_digits > 20 ? working_digits-20 : 1;
    mpf_pow_ui(cutoff,ten,cd); mpf_ui_div(cutoff,1,cutoff);
    unsigned long n=0;
    for(;;) {
        /* ((6n+1)(2n+1)(6n+5))/(72(n+1)^3) */
        mpf_set_ui(num,6*n+1); mpf_mul_ui(num,num,2*n+1); mpf_mul_ui(num,num,6*n+5);
        unsigned long np1=n+1;
        mpf_set_ui(den,np1); mpf_pow_ui(den,den,3); mpf_mul_ui(den,den,72);
        mpf_div(ratio,num,den); mpf_mul(term,term,ratio); mpf_mul(term,term,z);
        n++;
        mpf_add(F,F,term); mpf_mul_ui(tmp,term,n); mpf_add(thetaF,thetaF,tmp);
        mpf_abs(tmp,term);
        if (mpf_cmp(tmp,cutoff)<0) break;
        if (n>=max_terms) { mpf_clear(term);mpf_clear(ratio);mpf_clear(tmp);mpf_clear(num);mpf_clear(den);mpf_clear(cutoff);mpf_clear(ten); return 0; }
    }
    mpf_clear(term);mpf_clear(ratio);mpf_clear(tmp);mpf_clear(num);mpf_clear(den);mpf_clear(cutoff);mpf_clear(ten);
    return 1;
}

void rj_inverse_pi_from_state(const RjState *s, unsigned long working_digits,
                              mpf_t inv_pi, mpf_t F, mpf_t thetaF) {
    mp_bitcnt_t bits=mpf_get_prec(s->z);
    mpf_t t;
    mpf_init2(t,bits);
    rj_hyper_F_theta(s->z,working_digits,F,thetaF,10000000UL);
    mpf_mul(inv_pi,s->u,F); mpf_mul(t,s->v,thetaF); mpf_add(inv_pi,inv_pi,t);
    mpf_clear(t);
}
