sigma = 2;
beta = 0.99;
phi_pi = 1.5;
phi_y = 0.125;

sigmaepsilon = 0.2;
sigmas = 10;
sigmaeta = 0.4;
sigmau = 2.5;

alpha = 0;
epsilon = 2;
phi = 0.75;
eta = 2;
psi_c = 0.48;
psi_d = 1 - psi_c;

delta = 0.05;

lambda = 0.1;
zl = 1.5;
zh = 2;
zbar = 1.75;
p0 = 0.65;

y_ss = delta*psi_d;


pi_c_ss = 0;

Idss = delta;

D_ss = 1;

%y_ss = zetayminus1(5);
%pi_c_ss = zetapiminus1(5);
%Idss = zetaIdminus1(5);
%D_ss = zetaDtminus1(5);

M = 15;

zetayminus1 = zeros(M+1);

zetapiminus1 = zeros(M+1);

zetaIdminus1 = zeros(M+1);

zetaDtminus1 = zeros(M+1);



zetayminus1(1) = y_ss;

zetapiminus1(1) = pi_c_ss;

zetaIdminus1(1) = Idss;

zetaDtminus1(1) = D_ss;



Theta = (1 - alpha)/(1-alpha + alpha*epsilon);
q = ((1 - phi*beta)*(1-phi)*Theta)/(1 - phi);
kappa = q*(sigma + (eta + alpha)/(1 - alpha));


mu = 0.1;

I_D = delta;

markup = (epsilon - 1)/epsilon;

zeta_d = 0.6;

c1_d = (markup/1-alpha)*((zeta_d + alpha)/(1-alpha));
c1tild_d = c1_d*(1-mu);
c2_y = (markup/1-alpha)*((1-zeta_d + sigma*(1-alpha))/(1-alpha));

chi = (c1_d - (c2_y*(1-psi_c))/(psi_c));

%a2 = -pdf(shat)*(partialdshat/partialDt-1);







psi = [1 0; 0 -1; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0];

pi = [0 0 0; 1 0 -psi_d; (-kappa/psi_c) 1 (kappa*psi_d/psi_c); c2_y 0 chi; 0 0 I_D; 0 0 0; 1 0 0;
    0 1 0;
    0 0 1]; 

%M = 5;
N = 10;

ztild = zeros(M+1);
ztild(1) = p0*zh + (1-p0)*zl;
p = zeros(M+1);
m = zeros(M+1);
pr = zeros(M+1);
p(1) = p0;

c = zeros(M+1);

signalcutoff = zeros(M);

c(1) = ((-psi_c)/sigma)*(lambda*ztild(1) + (1-lambda)*zbar);


e = [0; c(1); 0; 0; 0; 0; 0; 0; 0];




mminus = zeros(M+1);

output = zeros(M);
interestrate = zeros(M);
nondurableinflation = zeros(M);
nondurableconsumption = zeros(M);
durablepricelevel = zeros(M);
durablegoodlevelD = zeros(M);
durablegoodinvestmentI = zeros(M);

shat = 10;

a2_old = 0.3;
a_3_2_old = 0.2;
a_3_3_old = 0.3;
a_3_4_old = 0.3;
a_3_6_old = 0.7;

a2_new = 0.3;
a_3_2_new = 0.2;
a_3_3_new = 0.3;
a_3_6_new = 0.7;

g1_c_y = zeros(M+1);
g1_c_pi = zeros(M+1);
g1_c_D = zeros(M+1);

a_c_y = zeros(M+1);
a_c_Id = zeros(M+1);

Dtminus1 = D_ss;


diff = 0.5;
tol = 10^(-3);
k = 1;

eta0 = -2*sigmaeta;



ratio = zeros(M+1);

for j = 1:1:M

 while (diff > tol) && (k <= 1000)
g0 = [1 0 0 0 0 0 0 0 0; 
    (-psi_c/sigma) 1 (psi_c/sigma) -psi_d 0 0 0 0 0;
    0 0 beta 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 a2_old*(1-mu) 0 0 0;
    0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 1];

g1 = [0 phi_y phi_pi 0 0 0 0 0 0; 
    0 1 0 -psi_d 0 0 0 0 0;
    0 (-kappa/psi_c) 1 (kappa*psi_d/psi_c) 0 0 0 0 0;
    0 c2_y 0 chi 0 0 0 0 0;
    0 0 0 I_D 0 (1-delta) 0 0 0;
    0 a_3_2_old a_3_3_old a_3_4_old 0 a_3_6_old 0 0 0;
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0];
[G1,C,impact,fmat,fwt,ywt,gev,eu]= gensys(g0,g1,e,psi,pi);
Q = G1;
Ay = [(G1(7,2)/c(1)) (G1(7,3)/c(1)) (G1(7,4)/c(1)) (G1(7,6)/c(1))]; % output
AId = [(G1(9,2)/c(1)) (G1(9,3)/c(1)) (G1(7,4)/c(1)) (G1(9,6)/c(1))]; % Id


D = Ay - psi_d*AId;
G1_c = (D)/psi_c;

a_c_y(j) = C(7)/c(1);
a_c_Id(j) = C(9)/c(1);

c2_c = c2_y;

%Dtminus1 = D_ss;

h = 0.1;
g1_c_y(j) = G1_c(1);
g1_c_pi(j) = G1_c(2);
g1_c_Id(j) = G1_c(3);
g1_c_D(j) = G1_c(4);


y_ss = zetayminus1(j);
pi_c_ss = zetapiminus1(j);
Idss = zetaIdminus1(j);
D_ss = zetaDtminus1(j);
g_y = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_yplush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss + h, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_yminush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss - h, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

%g_pi = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_piplush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss+h, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_piminush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss-h, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

%g_Id = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_Idplush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss + h, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_Idminush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss - h, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

%g_Dtminus = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_Dtminusplush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss+h, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_Dtminusminush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss-h, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

%g_Dtminusdelta = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_Dtminusdeltaplush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1+h)

g_Dtminusdeltaminush = G_dynamicinvestormass_FP(j, mminus, N, p0, a_c_y(j), a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y(j), g1_c_pi(j), g1_c_Id(j), g1_c_D(j), psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1-h)

dshatdy = (g_yplush(1) - g_yminush(1))/2*h;

dshatdpi = (g_piplush(1) - g_piminush(1))/2*h;

dshatdId = (g_Idplush(1) - g_Idminush(1))/2*h;

dshatddtminus = (g_Dtminusplush(1) - g_Dtminusminush(1))/2*h;

dshatddtminusdelta = (g_Dtminusdeltaplush(1) - g_Dtminusdeltaminush(1))/2*h;

a_3_2_new = -dshatdy*pdf('Normal', g_y(1), zl + eta0, sigmas);

a_3_3_new = -dshatdpi*pdf('Normal', g_y(1), zl + eta0, sigmas);

a_3_4_new = -dshatdId*pdf('Normal', g_y(1), zl + eta0, sigmas);

a_3_6_new = -dshatddtminus*pdf('Normal', g_y(1), zl + eta0, sigmas);

a2_new = -dshatddtminusdelta*pdf('Normal', g_y(1), zl + eta0, sigmas);

diff_3_2 = abs(a_3_2_new - a_3_2_old);
diff_3_3 = abs(a_3_3_new - a_3_3_old);
diff_3_4 = abs(a_3_4_new - a_3_4_old);
diff_3_6 = abs(a_3_6_new - a_3_6_old);
diff_2 = abs(a2_new - a2_old);
diff_array = [diff_3_2 diff_3_3 diff_3_4 diff_3_6 diff_2];
diff = max(diff_array);


if diff > tol
    a_3_2_old = a_3_2_new;
    a_3_3_old = a_3_3_new;
    a_3_4_old = a_3_4_new;
    a_3_6_old = a_3_6_new;
    a2_old = a2_new;
    
else
    if g_y(2) == 0
        sh = g_y(1);
    else
        break
    end
    m(j) = (1-mu)*normcdf((sh - (zl + eta0))/(sigmas));
    dm = m(j);
   if j == 1
       mminus(j) = 0;
   else
       mminus(j) = m(j-1);
   end
  templ = simpsonpublicbelief(j, N, zl, zh, sh,zl,dm,p0, mminus, sigmaeta, sigmas, sigmaepsilon);
  temph = simpsonpublicbelief(j, N, zl, zh, sh,zh,dm,p0, mminus, sigmaeta, sigmas, sigmaepsilon);
 %if templ == 0 & temph == 0
 %    templ = 1;
     %temph = 1;
 %elseif templ == 0 & temph~=0
   %  templ = 1;
    % temph = 1;
 %end
 ratio(j) = templ/temph;
 pr(j) = (p(j)*pdf('Normal', zl - zh, 0, sigmau))/(p(j)*pdf('Normal', zl - zh, 0, sigmau) + (1-p(j))*pdf('Normal', zl-zl, 0, sigmau));
  p(j+1) = 1/( 1 + ((1-pr(j))*templ)/(pr(j)*temph));
  p0 = p(j+1);
  ztild(j+1) = p(j+1)*zh + (1-p(j+1))*zl;
  c(j+1) = ((-psi_c)/sigma)*(lambda*ztild(j+1) + (1-lambda)*zbar);
dd = c(j+1);
e = [0; dd; 0; 0; 0; 0; 0; 0; 0];
    g0 = [1 0 0 0 0 0 0 0 0; 
    (-psi_c/sigma) 1 (psi_c/sigma) -psi_d 0 0 0 0 0;
    0 0 beta 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 a2_new*(1-mu) 0 0 0;
    0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 1];
g1 = [0 phi_y phi_pi 0 0 0 0 0 0; 
    0 1 0 -psi_d 0 0 0 0 0;
    0 (-kappa/psi_c) 1 (kappa*psi_d/psi_c) 0 0 0 0 0;
    0 c2_y 0 chi 0 0 0 0 0;
    0 0 0 I_D 0 (1-delta) 0 0 0;
    0 a_3_2_new a_3_3_new a_3_4_new 0 a_3_6_new 0 0 0;
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0];
[G1,C,impact,fmat,fwt,ywt,gev,eu]= gensys(g0,g1,e,psi,pi);
Azetay = [G1(2,2) G1(2,3) G1(2,4) G1(2,6)]; % zetaoutput
 czetay = C(2);
Azetapi = [G1(3,2) G1(3,3) G1(3,4) G1(3,6)]; % zetainflation
czetapi = C(3);
AzetaId = [G1(4,2) G1(4,3) G1(4,4) G1(4,6)/c(1)]; % zetaId
czetaId = C(4);
AzetaDt = [G1(6,2) G1(6,3) G1(6,4) G1(6,6)]; % zetaDt
czetaDt = C(6);
zetayminus1(j+1) = G1(2,2)*zetayminus1(j) + G1(2,3)*zetapiminus1(j) + G1(2,4)*zetaIdminus1(j) + G1(2,6)*zetaDtminus1(j) + C(2);
zetapiminus1(j+1) = G1(3,2)*zetayminus1(j) + G1(3,3)*zetapiminus1(j) + G1(3,4)*zetaIdminus1(j) + G1(3,6)*zetaDtminus1(j) + C(3);
zetaIdminus1(j+1) = G1(4,2)*zetayminus1(j) + G1(4,3)*zetapiminus1(j) + G1(4,4)*zetaIdminus1(j) + G1(4,6)*zetaDtminus1(j) + C(4);
zetaDtminus1(j+1) = G1(6,2)*zetayminus1(j) + G1(6,3)*zetapiminus1(j) + G1(6,4)*zetaIdminus1(j) + G1(6,6)*zetaDtminus1(j) + C(6);
%Impulse Responses Here
output(j) = G1(7,2)*zetayminus1(j) + G1(7,3)*zetapiminus1(j) + G1(7,4)*zetaIdminus1(j) + G1(7,6)*zetaDtminus1(j) + C(7);
nondurableinflation(j) = G1(8,2)*zetayminus1(j) + G1(8,3)*zetapiminus1(j) + G1(8,4)*zetaIdminus1(j) + G1(8,6)*zetaDtminus1(j) + C(8);
interestrate(j) = G1(1,2)*zetayminus1(j) + G1(1,3)*zetapiminus1(j) + C(1);
durablegoodinvestmentI(j) = G1(9,2)*zetayminus1(j) + G1(9,3)*zetapiminus1(j) + G1(9,4)*zetaIdminus1(j) + G1(9,6)*zetaDtminus1(j) + C(9);
durablepricelevel(j) = G1(5,2)*zetayminus1(j) + G1(5,3)*zetapiminus1(j) + G1(5,4)*zetaIdminus1(j) + G1(5,6)*zetaDtminus1(j) + C(5);
durablegoodlevelD(j) = G1(6,2)*zetayminus1(j) + G1(6,3)*zetapiminus1(j) + G1(6,4)*zetaIdminus1(j) + G1(6,6)*zetaDtminus1(j) + C(6);
nondurableconsumption(j) = (output(j) - psi_d*durablegoodinvestmentI(j))/psi_c;
%Simulations Here 
end
k = k+1;
 end
 k = 1;
 diff = 0.5;
 if g_y(2) == 0
        sh = g_y(1);
        signalcutoff(j) = sh;
  else
        break
  end
 
end


plot(p);

