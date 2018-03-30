data {
int<lower=0> N_obs;
int<lower=0> N_SN;
int<lower=0> N_filt;
vector[N_obs] t;
vector[N_obs] fL;
vector[N_obs] dfL;
vector[N_SN] z;
vector[N_SN] t0_mean; int<lower=1,upper=N_filt> J[N_obs]; int<lower=1,upper=N_SN> SNid[N_obs]; int<lower=0> Kcor_N;
real Kcor[N_SN, N_filt,Kcor_N]; real<lower=0> fluxscale; vector<lower=0,upper=1>[N_SN] duringseason;
}
transformed data {
vector[N_filt] prior_t_hF[4]; vector[N_filt] prior_t_hF_s[4]; vector[N_filt] prior_r_hF[5]; vector[N_filt] prior_r_hF_s[5]; for (i in 1:N_filt) {
prior_t_hF[1,i] <- 0;
prior_t_hF_s[1,i] <- 0.1; }
prior_t_hF[2,1] <- -1;
prior_t_hF[2,2] <- -0.5;
prior_t_hF[2,3] <- 0;
prior_t_hF[2,4] <- 0.5;
prior_t_hF[2,5] <- 1;
for (i in 1:N_filt) {prior_t_hF_s[2,i] <- 0.1;} for (i in 1:N_filt) {
prior_t_hF[3,i] <- 0;
prior_t_hF_s[3,i] <- 0.1; }
for (i in 1:N_filt) { prior_t_hF[4,i] <- 0; prior_t_hF_s[4,i] <- 0.1;
}
for (i in 1:N_filt) {
prior_r_hF[1,i] <- 0;
prior_r_hF_s[1,i] <- 0.1; }
prior_r_hF[2,1] <- 2;
prior_r_hF[2,2] <- 1;
prior_r_hF[2,3] <- 0;
prior_r_hF[2,4] <- -0.5;
prior_r_hF[2,5] <- -1;
for (i in 1:N_filt) {prior_r_hF_s[2,i] <- 0.1;} prior_r_hF[3,1] <- 1;
prior_r_hF[3,2] <- 0.3;
prior_r_hF[3,3] <- 0;
prior_r_hF[3,4] <- -1;
prior_r_hF[3,5] <- -1;
for (i in 1:N_filt) {prior_r_hF_s[3,i] <- 0.1;} for (i in 1:N_filt) {
prior_r_hF[4,i] <- 0;
prior_r_hF_s[4,i] <- 0.1; }
for (i in 1:N_filt) { prior_r_hF[5,i] <- 0; prior_r_hF_s[5,i] <- 0.1;
}
}
parameters {
vector[4] t_hP;
vector<lower=0>[4] sig_t_hP;
vector[N_filt] t_hF[4]; vector<lower=0>[N_filt] sig_t_hF[4]; vector[N_SN * N_filt] t_hSNF[4]; vector<lower=0>[N_SN * N_filt] sig_t_hSNF[4]; vector[5] r_hP;
vector<lower=0>[5] sig_r_hP;
vector[N_filt] r_hF[5];
vector<lower=0>[5] sig_r_hF[5];
vector[N_SN * N_filt] r_hSNF[5]; vector<lower=0>[N_SN * N_filt] sig_r_hSNF[5]; real M_h;
real<lower=0> sig_M_h;
vector[N_filt] M_hF;
vector<lower=0>[N_filt] sig_M_hF;
vector[N_SN * N_filt] M_hSNF; vector<lower=0>[N_SN * N_filt] sig_M_hSNF; real Y_h;
real<lower=0> sig_Y_h;
vector[N_SN * N_filt] Y_hSNF; vector<lower=0>[N_SN * N_filt] sig_Y_hSNF; real t0s_h;
real<lower=0> sig_t0s_h;
vector[N_SN] t0s_hSN;
vector<lower=0>[N_SN] sig_t0s_hSN;
real t0l_h;
real<lower=0> sig_t0l_h;
vector[N_SN] t0l_hSN;
vector<lower=0>[N_SN] sig_t0l_hSN; real<lower=0> V_h;
vector<lower=0>[N_filt] V_hF; vector<lower=0>[N_SN * N_filt] V_hSNF;
}
transformed parameters {
vector[N_obs] mm;
vector[N_obs] dm; vector<upper=0>[N_SN] pt0; matrix<lower=0>[N_SN, N_filt] t1; matrix<lower=0>[N_SN, N_filt] t2; matrix<lower=0>[N_SN, N_filt] td; matrix<lower=0>[N_SN, N_filt] tp; matrix[N_SN, N_filt] lalpha; matrix[N_SN, N_filt] lbeta1; matrix[N_SN, N_filt] lbeta2; matrix[N_SN, N_filt] lbetadN; matrix[N_SN, N_filt] lbetadC; matrix[N_SN, N_filt] Mp; matrix[N_SN, N_filt] Yb; matrix<lower=0>[N_SN, N_filt] V; matrix<lower=0>[N_SN, N_filt] M1; matrix<lower=0>[N_SN, N_filt] M2; matrix<lower=0>[N_SN, N_filt] Md; for (l in 1:N_SN) {
if (duringseason[l] == 1) {
pt0[l] <- -exp( t0s_h + sig_t0s_h * ( t0s_hSN[l] .* sig_t0s_hSN[l] ));
} else {
pt0[l] <- -exp( t0l_h + sig_t0l_h * ( t0l_hSN[l] .* sig_t0l_hSN[l] ));
} }
for (i in 1:N_filt) {
for (j in 1:N_SN) {
t1[j,i] <- exp( log(1) + t_hP[1] + sig_t_hP[1] * (
                   t_hF[1,i] * sig_t_hF[1,i]
                 + sig_t_hSNF[1,(i-1)*N_SN+j] * t_hSNF[1,(i-1)*N_SN+j]
));
tp[j,i] <- exp( log(10) + t_hP[2] + sig_t_hP[2] * ( t_hF[2,i] * sig_t_hF[2,i]
                 + sig_t_hSNF[2,(i-1)*N_SN+j] * t_hSNF[2,(i-1)*N_SN+j]
                   ));
t2[j,i] <- exp( log(100) + t_hP[3] + sig_t_hP[3] * ( t_hF[3,i] * sig_t_hF[3,i]
                 + sig_t_hSNF[3,(i-1)*N_SN+j] * t_hSNF[3,(i-1)*N_SN+j]
                   ));
td[j,i] <- exp( log(10) + t_hP[4] + sig_t_hP[4] * ( t_hF[4,i] * sig_t_hF[4,i]
                 + sig_t_hSNF[4,(i-1)*N_SN+j] * t_hSNF[4,(i-1)*N_SN+j]
                   ));
lalpha[j,i] <- -1 + ( r_hP[1] + sig_r_hP[1] * ( r_hF[1,i] * sig_r_hF[1,i]
                     + sig_r_hSNF[1,(i-1)*N_SN+j] * r_hSNF[1,(i-1)*N_SN+j]
                       ));
lbeta1[j,i] <- -4 + ( r_hP[2] + sig_r_hP[2] * ( r_hF[2,i] * sig_r_hF[2,i]
                     + sig_r_hSNF[2,(i-1)*N_SN+j] * r_hSNF[2,(i-1)*N_SN+j]
                       ));
lbeta2[j,i] <- -4 + ( r_hP[3] + sig_r_hP[3] * ( r_hF[3,i] * sig_r_hF[3,i]
                     + sig_r_hSNF[3,(i-1)*N_SN+j] * r_hSNF[3,(i-1)*N_SN+j]
                       ));
lbetadN[j,i] <- -3 + ( r_hP[4] + sig_r_hP[4] * ( r_hF[4,i] * sig_r_hF[4,i]
                      + sig_r_hSNF[4,(i-1)*N_SN+j] * r_hSNF[4,(i-1)*N_SN+j]
                        ));
lbetadC[j,i] <- -5 + ( r_hP[5] + sig_r_hP[5] * ( r_hF[5,i] * sig_r_hF[5,i]
                      + sig_r_hSNF[5,(i-1)*N_SN+j] * r_hSNF[5,(i-1)*N_SN+j]
                        ));
Mp[j,i] <- exp(M_h + sig_M_h * ( M_hF[i] * sig_M_hF[i]
                     + sig_M_hSNF[(i-1)*N_SN+j] * M_hSNF[(i-1)*N_SN+j]
                       ));
Yb[j,i] <- Y_h + sig_Y_h * (Y_hSNF[(i-1)*N_SN+j] .* sig_Y_hSNF[(i-1)*N_SN+j]);
V[j,i] <- V_h * V_hF[i] * V_hSNF[(i-1)*N_SN+j]; }
}
M1 <- Mp ./ exp( exp(lbeta1) .* tp ); M2 <- Mp .* exp( -exp(lbeta2) .* t2 ); Md <- M2 .* exp( -exp(lbetadN) .* td ); for (n in 1:N_obs) {
    real N_SNc;
    int Kc_up;
    int Kc_down;
    real t_exp;
    int j;
int k; real mm_1; real mm_2; real mm_3; real mm_4; real mm_5; real mm_6; j <- J[n];
k <- SNid[n];
t_exp <- ( t[n] - (t0_mean[k] + pt0[k]) ) / (1 + z[k]); if (t_exp<0) {
mm_1 } else { mm_1
}
if ((t_exp>=0) mm_2 } else { mm_2
<- Yb[k,j]; <- 0;
&& (t_exp < t1[k,j])) {
<- Yb[k,j] + M1[k,j] * pow(t_exp / t1[k,j] , exp(lalpha[k,j]));
<- 0;
}
if ((t_exp >= t1[k,j]) && (t_exp < t1[k,j] + tp[k,j])) {
mm_3 <- Yb[k,j] + M1[k,j] * exp(exp(lbeta1[k,j]) * (t_exp - t1[k,j])); } else {
mm_3 <- 0; }
if ((t_exp >= t1[k,j] + tp[k,j]) && (t_exp < t1[k,j] + tp[k,j] + t2[k,j])) {
mm_4 <- Yb[k,j] + Mp[k,j] * exp(-exp(lbeta2[k,j]) * (t_exp - t1[k,j] - tp[k,j]));
} else {
mm_4 <- 0;
}
if ((t_exp >= t1[k,j] + tp[k,j] + t2[k,j]) && (t_exp < t1[k,j] + tp[k,j] + t2[k,j] + td[k,j])) {
mm_5 <- Yb[k,j] + M2[k,j] * exp(-exp(lbetadN[k,j]) * (t_exp - t1[k,j] - tp[k,j] - t2[k,j]));
} else {
mm_5 <- 0;
}
if (t_exp >= t1[k,j] + tp[k,j] + t2[k,j] + td[k,j]) {
mm_6 <- Yb[k,j] + Md[k,j] * exp(-exp(lbetadC[k,j]) * (t_exp - t1[k,j] - tp[k,j] - t2[k,j] - td[k,j]));
} else {
mm_6 <- 0;
}
dm[n] <- sqrt(pow(dfL[n],2) + pow(V[k,j],2)); if (t_exp<0) {
N_SNc <- 0;
} else if (t_exp<Kcor_N-2){
Kc_down <- 0;
while ((Kc_down+1) < t_exp) { Kc_down <- Kc_down + 1;
}
Kc_up <- Kc_down+1;
N_SNc <- Kcor[k,j,Kc_down+1] + (t_exp - floor(t_exp)) * (Kcor[k,j,Kc_up+1]
-Kcor[k,j,Kc_down+1]); } else {
N_SNc <- Kcor[k,j,Kcor_N]; }
mm[n] <- (mm_1+mm_2+mm_3+mm_4+mm_5+mm_6) / (pow(10, N_SNc/(-2.5))); }
}
model {
    t0s_h   normal(0, 0.5);
    sig_t0s_h   cauchy(0, 0.1);
    t0l_h   normal(log(100), 1);
    sig_t0l_h   cauchy(0, 0.1);
    V_h   cauchy(0, 0.001);
    Y_h   normal(0, 0.1);
    sig_Y_h   cauchy(0, 0.01);
    M_h   normal(0, 1);
    sig_M_h   cauchy(0, 0.1);
    t_hP   normal(0,0.1);
    sig_t_hP cauchy(0, 0.1); for (i in 1:4) {
  t_hF[i]   normal(prior_t_hF[i], prior_t_hF_s[i]);
  sig_t_hF[i]   cauchy(0, 0.1);
  t_hSNF[i]   normal(0,1);
  sig_t_hSNF[i]   cauchy(0, 0.1);
}
r_hP normal(0,1); sig_r_hP cauchy(0, 0.1); for (i in 1:5) {
  r_hF[i]   normal(prior_r_hF[i], prior_r_hF_s[i]);
  sig_r_hF[i]   cauchy(0, 0.1);
  r_hSNF[i]   normal(0,1);
  sig_r_hSNF[i]   cauchy(0, 0.1);
}
M_hF   normal(0,1);
sig_M_hF   cauchy(0, 0.1);
M_hSNF   normal(0,1);
sig_M_hSNF   cauchy(0, 0.1);
Y_hSNF   normal(0,1);
sig_Y_hSNF   cauchy(0, 0.1);
      V_hF
      V_hSNF
      t0s_hSN
      sig_t0s_hSN   cauchy(0, 0.1);
      t0l_hSN   normal(0,1);
      sig_t0l_hSN   cauchy(0, 0.1);
      fL   normal(mm,dm);
}
