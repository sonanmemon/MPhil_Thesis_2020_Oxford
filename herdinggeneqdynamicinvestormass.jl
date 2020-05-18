using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.4.0")
using LinearAlgebra, Statistics
using DataFrames, Parameters, Plots, Printf, Random
gr(fmt = :png);
using Distributions
using LaTeXStrings
using .LRESolve, Test






let
sigma = 2
beta = 0.99
phi_pi = 1.5
phi_y = 0.125
pcurrent = 0

g_y = 0

sigmaepsilon = 5
sigmas = 15
sigmaeta = 0.3
sigmau = 2.5

alpha = 0
epsilon = 2
phi = 0.75
eta = 2
psi_c = 0.48
psi_d = 1 - psi_c

delta = 0.07
gamma = 1.02

lambda = 0.1
zl = 1.5
zh = 2
zbar = 1.75
p0 = 0.65



rho = 0.9



pi_c_ss = 0

#Idss = delta


D_ss = 2
Idss = delta
y_ss = delta*psi_d;

psi = 1

mu_y = ((psi - alpha)*(1-alpha) - sigma)/psi_c

mu_Id = (sigma*psi_d - psi_d*(psi - alpha)*(1-alpha))/psi_c




dshatdy = 0







#y_ss = zetayminus1(5)
#pi_c_ss = zetapiminus1(5)
#Idss = zetaIdminus1(5)
#D_ss = zetaDtminus1(5)

M = 25

zetayminus1 = zeros(M+1)

zetapiminus1 = zeros(M+1)

zetaIdminus1 = zeros(M+1)

zetaDtminus1 = zeros(M+1)

irfoutputgap = zeros(M+1)

irfnaturaloutput = zeros(M+1)
irfoutput = zeros(M+1)
irfnondurableconsumption = zeros(M+1)
irfaggregatehours = zeros(M+1)

ztild = zeros(M+1)
ztild[1] = p0*zh + (1-p0)*zl
p = zeros(M+1)
m = zeros(M+1)
pr = zeros(M+1)
p[1] = p0

c = zeros(M+1)

signalcutoff = zeros(M)

c[1] = ((-psi_c)/sigma)*(lambda*ztild[1] + (1-lambda)*zbar)

zetayminus1[1] = 0

zetapiminus1[1] = 0

zetaIdminus1[1] = 0

zetaDtminus1[1] = 0






Theta = (1 - alpha)/(1-alpha + alpha*epsilon)
LAMBDA = ((1 - phi*beta)*(1-phi)*Theta)/(1 - phi)
#kappa = q*(sigma + (eta + alpha)/(1 - alpha))

kappa1 = -LAMBDA*mu_y


mu = 0.1

I_D = delta/1

markup = (epsilon - 1)/epsilon

zeta_d = 0.6

c1_d = (markup/1-alpha)*((zeta_d + alpha)/(1-alpha))
c1tild_d = c1_d*(1-mu)
c2_y = (markup/1-alpha)*((1-zeta_d + sigma*(1-alpha))/(1-alpha))

chi = (c1_d - (c2_y*(1-psi_c))/(psi_c))

#a2 = -pdf(shat)*(partialdshat/partialDt-1);







psi = [1 0; 0 -1; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]




N = 10

ztild = zeros(M+1)
ztild[1] = p0*zh + (1-p0)*zl
p = zeros(M+1)
m = zeros(M+1)
pr = zeros(M+1)
p[1] = p0

c = zeros(M+1)

signalcutoff = zeros(M)



c[1] = ((-psi_c)/sigma)*(lambda*ztild[1] + (1-lambda)*zbar)


e = [0; c[1]; 0; 0; 0; 0; 0; 0; 0]






mminus = zeros(M+1)

outputgap = zeros(M+1)
outputgaprelated = zeros(M+1)
output = zeros(M+1)
durableinvestment = zeros(M+1)
naturaloutput = zeros(M+1)
aggregatehours = zeros(M+1)
interestrate = zeros(M+1)
nondurableinflation = zeros(M+1)
nondurableconsumption = zeros(M+1)
durablepricelevel = zeros(M+1)
durablegoodlevelD = zeros(M+1)
durablegoodinvestmentI = zeros(M+1)

shat = 10



a_6_2_new = -0.0
a_6_3_new = -0.0
a_6_4_new = -0.0
a_new = -3.4739270016290274e-10
a_6_6_new = 5.174623512762339e-9
g0 = [1 0 0 0 0 0 0 0 0;
    (-psi_c/sigma) 1 (psi_c/sigma) -psi_d 0 0 0 0 0;
    0 0 beta 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 a_new*(1-mu) 0 0 0;
    0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 1]

g1 = [0 phi_y phi_pi -(phi_y*mu_Id)/mu_y 0 0 0 0 0;
    0 1 0 (-psi_d + (mu_Id/mu_y)) 0 0 0 0 0;
    0 -kappa1 1 0 0 0 0 0 0;
    0 c2_y 0 (chi - ((c2_y*mu_Id)/mu_y)) 0 0 0 0 0;
    0 0 0 I_D 0 (1-delta) 0 0 0;
    0 a_6_2_new a_6_3_new a_6_4_new 0 a_6_6_new 0 0 0;
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0]

    pi = [phi_y phi_pi -(phi_y*mu_Id)/mu_y;
    1 0 (-psi_d + (mu_Id/mu_y));
    -kappa1 1 0;
    c2_y 0 (chi - ((c2_y*mu_Id)/mu_y));
    0 0 I_D;
    a_6_2_new a_6_3_new a_6_4_new;
    1 0 0;
    0 1 0;
    0 0 1]


M0 = ModelSims(g0, g1, e, psi, pi)


C, G0, G1 = solve_sims(M0)
Cprior = C
Gprior = G1

C11 = C
G11 = G1

X1 = zeros(9, M+1)

Xp = zeros(9, M+1)

irf = zeros(9, M+1)

a_old = 0.3
a_6_2_old = 0.2
a_6_3_old = 0.3
a_6_4_old = 0.3
a_6_6_old = 0.7

a_new = 0.3
a_6_2_new = 0.2
a_6_3_new = 0.3
a_6_4_new = 0.3
a_6_6_new = 0.7

g1_c_y = zeros(M+1)
g1_c_pi = zeros(M+1)
g1_c_Id = zeros(M+1)
g1_c_D = zeros(M+1)

a_c_y = zeros(M+1)
a_c_Id = zeros(M+1)

Dtminus1 = D_ss


D_SS = zeros(M)

D_SS[1] = 2


diff = 0.5
tol = 10^(-3)
k = 1

eta0 = 2*sigmaeta



ratio = zeros(M+1)

for j in 1:M

 while (diff > tol) && (k <= 1000)
g0 = [1 0 0 0 0 0 0 0 0;
    (-psi_c/sigma) 1 (psi_c/sigma) -psi_d 0 0 0 0 0;
    0 0 beta 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 a_old*(1-mu) 0 0 0;
    0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 1]



g1 = [0 phi_y phi_pi -(phi_y*mu_Id)/mu_y 0 0 0 0 0;
    0 1 0 (-psi_d + (mu_Id/mu_y)) 0 0 0 0 0;
    0 -kappa1 1 0 0 0 0 0 0;
    0 c2_y 0 (chi - ((c2_y*mu_Id)/mu_y)) 0 0 0 0 0;
    0 0 0 I_D 0 (1-delta) 0 0 0;
    0 a_6_2_old a_6_3_old a_6_4_old 0 a_6_6_old 0 0 0;
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0]

    pi = [phi_y phi_pi -(phi_y*mu_Id)/mu_y;
    1 0 (-psi_d + (mu_Id/mu_y));
    -kappa1 1 0;
    c2_y 0 (chi - ((c2_y*mu_Id)/mu_y));
    0 0 I_D;
    a_6_2_old a_6_3_old a_6_4_old;
    1 0 0;
    0 1 0;
    0 0 1]

    e = [0; c[1]; 0; 0; 0; 0; 0; 0; 0]
    M0 = ModelSims(g0, g1, e, psi, pi)


    C, G0, G1 = solve_sims(M0)
Q = G1
Ay = [(G1[7,2]/c[1]) (G1[7,3]/c[1]) (G1[7,4]/c[1]) (G1[7,6]/c[1])]; # output
AId = [(G1[9,2]/c[1]) (G1[9,3]/c[1]) (G1[7,4]/c[1]) (G1[9,6]/c[1])]; # Id


D = Ay - psi_d*AId


G1_c = (D)/psi_c

a_c_y[j] = C[7]/c[1]
a_c_Id[j] = C[9]/c[1]

c2_c = c2_y

#Dtminus1 = D_ss

h = 0.001
g1_c_y[j] = G1_c[1]
g1_c_pi[j] = G1_c[2]
g1_c_Id[j] = G1_c[3]
g1_c_D[j] = G1_c[4]

if j == 1
    pi_c_ss = 0
    D_ss = 2
    Idss = delta
    y_ss = delta*psi_d
else
y_ss = delta*psi_d*(1 + zetayminus1[j])
pi_c_ss = zetapiminus1[j]
Idss = m[j]
D_SS[j] = (1 - delta)*D_SS[j-1] + m[j]
D_ss = D_SS[j]
end

g_y = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, gamma, Dtminus1)

g_yplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss + h, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, gamma, Dtminus1)

g_yminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss - h, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

#g_pi = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id(j), c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_piplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss+h, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

g_piminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss-h, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

#g_Id = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_Idplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss + h, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

g_Idminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss - h, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

#g_Dtminus = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_Dtminusplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss, D_ss+h, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

g_Dtminusminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss, D_ss-h, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

#g_Dtminusdelta = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

g_Dtminusdeltaplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1+h)

g_Dtminusdeltaminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], c2_c, y_ss, pi_c_ss, Idss, D_ss, c1tild_d, shat, g1_c_y[j], g1_c_pi[j], g1_c_Id[j], g1_c_D[j], psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1-h)

dshatdy = (g_yplush[1] - g_yminush[1])/2*h





dshatdpi = (g_piplush[1] - g_piminush[1])/2*h

dshatdId = (g_Idplush[1] - g_Idminush[1])/2*h

dshatddtminus = (g_Dtminusplush[1] - g_Dtminusminush[1])/2*h

dshatddtminusdelta = (g_Dtminusdeltaplush[1] - g_Dtminusdeltaminush[1])/2*h


w = Normal(zl + eta0, sigmas)
a_6_2_new = -dshatdy*pdf(w, g_y[1])

a_6_3_new = -dshatdpi*pdf(w, g_y[1])

a_6_4_new = -dshatdId*pdf(w, g_y[1])

a_6_6_new = -dshatddtminus*pdf(w, g_y[1])

a_new = -dshatddtminusdelta*pdf(w, g_y[1])

diff_6_2 = abs(a_6_2_new - a_6_2_old)
diff_6_3 = abs(a_6_3_new - a_6_3_old)
diff_6_4 = abs(a_6_4_new - a_6_4_old)
diff_6_6 = abs(a_6_6_new - a_6_6_old)
diff_a = abs(a_new - a_old)
diff_array = [diff_6_2 diff_6_3 diff_6_4 diff_6_6 diff_a]
diff = findmax(diff_array)[1]


if diff > tol
    a_6_2_old = a_6_2_new
    a_6_3_old = a_6_3_new
    a_6_4_old = a_6_4_new
    a_6_6_old = a_6_6_new
    a_old = a_new

else
    if g_y[2] == 0
        sh = g_y[1]
    else
        break
    end
    x = Normal(0,1)
    #if j == 1
    #    eta0[j] = 2*sigmaeta
    #else
    #    eta0[j] = rho*eta0[j-1]
    #end
    m[j] = 1 - ((1-mu)*cdf(x, (sh - (zl + eta0))/(sigmas)))
    dm = m[j]
   if j == 1
       mminus[j] = 0
   else
       mminus[j] = m[j-1]
   end
  templ = simpsonpublicbelief(j, N, zl, zh, sh,zl,dm,p0, mminus, sigmaeta, sigmas, sigmaepsilon)
  temph = simpsonpublicbelief(j, N, zl, zh, sh,zh,dm,p0, mminus, sigmaeta, sigmas, sigmaepsilon)
 #if templ == 0 & temph == 0
 #    templ = 1;
     #temph = 1;
 # elseif templ == 0 & temph~=0
   #  templ = 1;
    # temph = 1;
 #end
 ratio[j] = templ/temph
 u = Normal(0, sigmau)
 #pr[j] = (p[j]*pdf(u, zl - zh))/(p[j]*pdf(u, zl - zh) + (1-p[j])*pdf(u, zl-zl))
  p[j+1] = 1/( 1 + (((1-p[j])*templ)/(p[j]*temph)))
  d = p[j+1]
  if isnan(d)
      p[j+1] = 0
  end
  #p0 = p[j+1]
  ztild[j+1] = p[j+1]*zh + (1-p[j+1])*zl
  c[j+1] = ((-psi_c)/sigma)*(lambda*ztild[j+1] + (1-lambda)*zbar)
dd = c[j+1]
e = [0; dd; 0; 0; 0; 0; 0; 0; 0];
g0 = [1 0 0 0 0 0 0 0 0;
    (-psi_c/sigma) 1 (psi_c/sigma) -psi_d 0 0 0 0 0;
    0 0 beta 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 a_new*(1-mu) 0 0 0;
    0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 1]

g1 = [0 phi_y phi_pi -(phi_y*mu_Id)/mu_y 0 0 0 0 0;
    0 1 0 (-psi_d + (mu_Id/mu_y)) 0 0 0 0 0;
    0 -kappa1 1 0 0 0 0 0 0;
    0 c2_y 0 (chi - ((c2_y*mu_Id)/mu_y)) 0 0 0 0 0;
    0 0 0 I_D 0 (1-delta) 0 0 0;
    0 a_6_2_new a_6_3_new a_6_4_new 0 a_6_6_new 0 0 0;
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0]

    pi = [phi_y phi_pi -(phi_y*mu_Id)/mu_y;
    1 0 (-psi_d + (mu_Id/mu_y));
    -kappa1 1 0;
    c2_y 0 (chi - ((c2_y*mu_Id)/mu_y));
    0 0 I_D;
    a_6_2_new a_6_3_new a_6_4_new;
    1 0 0;
    0 1 0;
    0 0 1]
M0 = ModelSims(g0, g1, e, psi, pi)
C, G0, G1 = solve_sims(M0)
if j == 16
    #print(a_6_2_new)
    #print(a_6_3_new)
    #print(a_6_4_new)
    #print(a_new)
end
Azetay = [G1[2,2] G1[2,3] G1[2,4] G1[2,6]] # zetaoutput
 czetay = C[2]
Azetapi = [G1[3,2] G1[3,3] G1[3,4] G1[3,6]] # zetainflation
czetapi = C[3]
AzetaId = [G1[4,2] G1[4,3] G1[4,4] G1[4,6]/c[1]] # zetaId
czetaId = C[4]
AzetaDt = [G1[6,2] G1[6,3] G1[6,4] G1[6,6]] # zetaDt
czetaDt = C[6]
zetayminus1[j+1] = G1[2,2]*zetayminus1[j] + G1[2,3]*zetapiminus1[j] + G1[2,4]*zetaIdminus1[j] + G1[2,6]*zetaDtminus1[j] + C[2]
zetapiminus1[j+1] = G1[3,2]*zetayminus1[j] + G1[3,3]*zetapiminus1[j] + G1[3,4]*zetaIdminus1[j] + G1[3,6]*zetaDtminus1[j] + C[3]
zetaIdminus1[j+1] = G1[4,2]*zetayminus1[j] + G1[4,3]*zetapiminus1[j] + G1[4,4]*zetaIdminus1[j] + G1[4,6]*zetaDtminus1[j] + C[4]
zetaDtminus1[j+1] = G1[6,2]*zetayminus1[j] + G1[6,3]*zetapiminus1[j] + G1[6,4]*zetaIdminus1[j] + G1[6,6]*zetaDtminus1[j] + C[6]
#Impulse Responses Here
outputgap[j+1] = G1[7,2]*zetayminus1[j] + G1[7,3]*zetapiminus1[j] + G1[7,4]*zetaIdminus1[j] + G1[7,6]*zetaDtminus1[j] + C[7]

if j == 1
    C11 = C
    G11 = G1
end
#if j == 1
#irfoutputgap[j+1] = G1[7,2]*zetayminus1[j] + G1[7,3]*zetapiminus1[j] + G1[7,4]*zetaIdminus1[j] + G1[7,6]*zetaDtminus1[j] + C[7] - (
#G1[7,2]*zetayminus1[j] + G1[7,3]*zetapiminus1[j] + G1[7,4]*zetaIdminus1[j] + G1[7,6]*zetaDtminus1[j] + Cprior[7])
#print(irfoutputgap[j+1])
#end
X1[:,j+1] = G11*X1[:,j] + C11
Xp[:, j+1] = Gprior*Xp[:,j] + Cprior
#if j >=6
#    irf[:,j+1] = rho.*irf[:,j]
#else
    irf[:,j+1] = (X1[:,j+1] - Xp[:,j+1])
#end
irfnaturaloutput[j+1] = -(mu_Id*irf[9,j+1])/mu_y
irfoutput[j+1] = irf[7,j+1] + irfnaturaloutput[j+1]
irfnondurableconsumption[j+1] = (irfoutput[j] - psi_d*irf[9, j])/psi_c
irfaggregatehours[j+1] = nondurableconsumption[j+1] + m[j]
#outputgap[j+1] = outputgaprelated[j+1] + c[j+1]
#if j == 16
#    print(outputgap[j])
#end
nondurableinflation[j+1] = G1[8,2]*zetayminus1[j] + G1[8,3]*zetapiminus1[j] + G1[8,4]*zetaIdminus1[j] + G1[8,6]*zetaDtminus1[j] + C[8]
interestrate[j+1] = G1[1,2]*zetayminus1[j] + G1[1,3]*zetapiminus1[j] + C[1];
#durablegoodinvestmentI[j+1] = G1[9,2]*zetayminus1[j] + G1[9,3]*zetapiminus1[j] + G1[9,4]*zetaIdminus1[j] + G1[9,6]*zetaDtminus1[j] + C[9]
durablepricelevel[j+1] = G1[5,2]*zetayminus1[j] + G1[5,3]*zetapiminus1[j] + G1[5,4]*zetaIdminus1[j] + G1[5,6]*zetaDtminus1[j] + C[5]
durablegoodlevelD[j+1] = G1[6,2]*zetayminus1[j] + G1[6,3]*zetapiminus1[j] + G1[6,4]*zetaIdminus1[j] + G1[6,6]*zetaDtminus1[j] + C[6]
#nondurableconsumption[j+1] = (output[j] - psi_d*durablegoodinvestmentI[j])/psi_c
durableinvestment[j+1] = G1[9,2]*zetayminus1[j] + G1[9,3]*zetapiminus1[j] + G1[9,4]*zetaIdminus1[j] + G1[9,6]*zetaDtminus1[j] + C[9]
naturaloutput[j+1] = -(mu_Id*durableinvestment[j+1])/mu_y
output[j+1] = outputgap[j+1] + naturaloutput[j+1]
nondurableconsumption[j+1] = (output[j+1] - psi_d*durableinvestment[j+1])/psi_c
aggregatehours[j+1] = nondurableconsumption[j+1] + m[j]
#Simulations Here
end
k = k+1
 end
 #if j == 1
#     print(dshatdy)
 #end
 k = 1
 diff = 0.5
 if g_y[2] == 0
        sh = g_y[1]
        signalcutoff[j] = sh
  else
        break
  end

end



p[isnan.(p)] .= 0
plot(p, color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Beliefs}", title = L"\textbf{Evolution of Public Beliefs}", label = L"\textbf{Beliefs}", linewidth = 3, leg = false)
savefig("GEpublicbeliefsFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
output[isnan.(output)] .= 0
plot(irfoutput, color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviation From SS}", title = L"\textbf{Response of Output}", linewidth = 3, label = L"\textbf{Output}")
savefig("GEoutputFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
outputgap[isnan.(outputgap)] .= 0
plot(irf[7, :], color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviation From SS}", title = L"\textbf{Response of Output Gap}", linewidth = 3, label = L"\textbf{Output Gap}")
savefig("GEoutputgapFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
aggregatehours[isnan.(aggregatehours)] .= 0
plot(irfaggregatehours, color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviation From SS}", title = L"\textbf{Response of Aggregate Hours}", linewidth = 3, label = L"\textbf{Hours}")
savefig("GEaggregatehoursFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
durableinvestment[isnan.(durableinvestment)] .= 0
plot(irf[9, :], color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviation From SS}", title = L"\textbf{Response of Duarable Investment}", linewidth = 3, label = L"\textbf{Durable Investment}")
savefig("GEdurableinvestmentFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
signalcutoff[isnan.(signalcutoff)] .= 0
plot(signalcutoff, color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Signal Cutoff}", title = L"\textbf{Response of Signal Cutoff}", linewidth = 3, label = L"\textbf{Cutoff}")
savefig("GEsignalcutoffFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
interestrate[isnan.(interestrate)] .= 0
plot(irf[1,:], color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviations From SS}", title = L"\textbf{Response of Interest Rate}", linewidth = 3, label = L"\textbf{Interest Rate}")
savefig("GEinterestrateFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
durablepricelevel[isnan.(durablepricelevel)] .= 0
plot(irf[5,:], color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviations From SS}", title = L"\textbf{Response of Durable Good Pirce Level}", label = L"\textbf{Durable Price}", linewidth = 3)
savefig("GEdurablepricelevelFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
nondurableinflation[isnan.(nondurableinflation)] .= 0
plot(irf[8,:], color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviations From SS}", title = L"\textbf{Response of Non-Durable Good Inflation}", linewidth = 3, label = L"\textbf{Non-duable Inflation}")
savefig("GEnondurableinflationFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
nondurableconsumption[isnan.(nondurableconsumption)] .= 0
plot(irfnondurableconsumption, color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviations From SS}", title = L"\textbf{Response of Non-durable Consumption}", linewidth = 3, label = L"\textbf{Non-durable Consumption}")
savefig("GEnondurableconsumptionFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
m[isnan.(m)] .= 0
plot(m, color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Level of Mass}", title = L"\textbf{Evolution of Mass of Investors}", linewidth = 3, label = L"\textbf{Mass}")
savefig("GEmassofinvestorsFP2sdsigmas15zl1pt5sigmaepsilon518may.pdf")
plot(irfoutput, color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Log Deviation From SS}", title = L"\textbf{Response of Output}", linewidth = 3, label = L"\textbf{Output}")
end
