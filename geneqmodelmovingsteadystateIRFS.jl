using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.4.0")
using LinearAlgebra, Statistics
using DataFrames, Parameters, Plots, Printf, Random
gr(fmt = :png);
using Distributions
using LaTeXStrings
using .LRESolve, Test






let
sigma = 2 #risk aversion
beta = 0.99 # discount factor
phi_pi = 1.5 # taylor rule inflation
phi_y = 0.125 #taylor rule output response
rhoi = 0.8 # interest rate smoothing

g_y = 0

sigmaepsilon = 5
sigmas = 15
sigmaeta = 0.3
sigmau = 2.5

alpha = 0
epsilon = 6
phi = 0.92
phi_D = 0.92
phi_C = 0.92
eta = 2
psi_c = 0.3
psi_d = 1 - psi_c

delta = 0.0033
gamma = 1.02

lambda = 0.1
zl = 1.5
zh = 2
zbar = 1.75
p0 = 0.63

g_interestrate = 0



rho = (1/beta) - 1



pi_c_ss = 0
pi_d_ss = 0
#Idss = delta
zeta_d = 0.38

D_ss = 2
Idss = delta
y_ss = delta*psi_d;

psi = 1

mu_y_c = -(sigma + (psi*(1-zeta_d) + alpha)*(1-alpha))/psi_c


mu_Id_c = ((sigma*psi_d + psi_d*(psi*(1-zeta_d) + alpha)*(1-alpha))/psi_c) - ((psi*zeta_d)/(1-alpha))

mu_y_Id = -(sigma + (psi*(1-zeta_d))*(1-alpha))/psi_c

mu_Id_Id = ((sigma*psi_d + psi_d*(psi*(1-zeta_d))*(1-alpha))/psi_c) - ((psi*zeta_d)/(1-alpha)) - (alpha/(1-alpha))


dshatdy = 0


Theta = (1 - alpha)/(1-alpha + alpha*epsilon)
LAMBDA = ((1 - phi*beta)*(1-phi)*Theta)/(1 - phi)
#kappa = q*(sigma + (eta + alpha)/(1 - alpha))

kappa1_c = -LAMBDA*mu_y_c

kappa2_c = -LAMBDA*mu_Id_c

kappa1_Id = -LAMBDA*mu_y_Id

kappa2_Id = -LAMBDA*mu_Id_Id


mu = 0.1
I_D = delta/1

markup = (epsilon - 1)/epsilon



c1_Id = (markup/1-alpha)*((zeta_d + alpha)/(1-alpha))
#c1tild_d = c1_Id*(1-mu)
c2_y = (markup/1-alpha)*((1-zeta_d + sigma*(1-alpha))/(1-alpha))

#chi = (c1_Id - (c2_y*(1-psi_c))/(psi_c))




#y_ss = zetayminus1(5)
#pi_c_ss = zetapiminus1(5)
#Idss = zetaIdminus1(5)
#D_ss = zetaDtminus1(5)

M = 7

zetayminus1 = zeros(M+1)

zetapiminus1 = zeros(M+1)

zetaIdminus1 = zeros(M+1)

zetaDtminus1 = zeros(M+1)

Xtminus1 = zeros(12, M+1)

Xtminus1original = zeros(12, M+1)

irfoutputgap = zeros(M+1)

irfnaturaloutput = zeros(M+1)
irfoutput = zeros(M+1)
output = zeros(M+1)
irfnondurableconsumption = zeros(M+1)
nondurableconsumption = zeros(M+1)
irfaggregatehours = zeros(M+1)
aggregatehours = zeros(M+1)
irfdurableinflation = zeros(M+1)
durableinflation = zeros(M+1)

irfrelativepricesdurablenondurable = zeros(M+1)
relativepricesdurablenondurable = zeros(M+1)


ztild = zeros(M+1)
ztild[1] = p0*zh + (1-p0)*zl
p = ones(M+1)
m = zeros(M+1)
pr = zeros(M+1)
p[1] = p0

c = zeros(M+1)

signalcutoff = zeros(M)

c[1] = ((-psi_c)/sigma)*(lambda*ztild[1] + (1-lambda)*zbar)

#c[1] = 0

#c[1] = ((-psi_c)/sigma)*(zbar)








#a2 = -pdf(shat)*(partialdshat/partialDt-1);






N = 10
psi = [1 0; 0 -1; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]




a_6_1_old = 0.2
a_6_2_old = 0.2
  a_6_3_old = 0.3
  a_6_4_old = 0.5
  a_6_5_old = 0.4
    a_6_6_old = 0.4
     a_6_7_old = 0.2
     a_6_8_old = 0.5
      a_6_9_old = 0.4
       a_6_10_old = 0.2
        a_6_11_old = 0.1
        a_6_12_old = 0.5





e = [0; c[1]; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]

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
durableinflation = zeros(M+1)
durablegoodlevelD = zeros(M+1)
durablegoodinvestmentI = zeros(M+1)




shat = 10



#a_6_2_new = -0.0
#a_6_3_new = -0.0
#a_6_4_new = -0.0
#a_new = -3.4739270016290274e-10
#a_6_6_new = 5.174623512762339e-9
g0 = [1 0 0 0 0 0 0 0 0 0 0 0;
    (-psi_c/sigma) 1 (psi_c/sigma) 0 -psi_d 0 0 0 0 0 0 0;
    0 0 beta 0 0 0 0 0 0 0 0 0;
    0 0 0 beta 0 LAMBDA -LAMBDA 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 -1 0 0 0 0 1 0;
    0 0 0 0 0 0 -1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 1]

g1 = [rhoi phi_y 0 0 0 0 0 0 0 phi_pi 0 0;
    0 1 0 0 -psi_d 0 0 0 0 0 0 0;
    0 -kappa1_c 1 0 -kappa2_c 0 0 0 0 0 0 0;
    0 -kappa1_Id 0 1 -kappa2_Id 0 0 0 0 0 0 0;
    0 0 psi_c psi_d 0 0 0 0 0 0 0 0;
    0 0 0 0 I_D 0 0 (1-delta) 0 0 0 0;
    a_6_1_old a_6_2_old a_6_3_old a_6_4_old a_6_5_old a_6_6_old a_6_7_old a_6_8_old a_6_9_old a_6_10_old a_6_11_old a_6_12_old;
    0 0 0 0 0 -1 0 0 0 0 0 0;
    0 0 0 -1 0 0 -1 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0]

    pi = [phi_y 0 0 0;
    1 0 0 -psi_d;
    -kappa1_c 1 0 -kappa2_c;
    -kappa1_Id 0 1 -kappa2_Id;
    0 psi_c psi_d 0;
    0 0 0 I_D;
    0 0 0 0;
    0 0 0 0;
    0 0 -1 0;
    1 0 0 0;
    0 1 0 0;
    0 0 0 1]


M0 = ModelSims(g0, g1, e, psi, pi)


C, G0, G1 = solve_sims(M0)
Cprior = C
Gprior = G1
C11 = C
G11 = G1

X1 = zeros(12, M+1)

Xp = zeros(12, M+1)

irf = zeros(12, M+1)


g1_c_y = zeros(M+1)
g1_c_pi = zeros(M+1)
g1_c_Id = zeros(M+1)
g1_c_D = zeros(M+1)

a_c_y = zeros(M+1)
a_c_Id = zeros(M+1)
a_c_pd = zeros(M+1)


Dtminus1 = D_ss


D_SS = zeros(M)

D_SS[1] = 2


diff = 0.5
tol = 10^(-3)
k = 1

eta0 = 2*sigmaeta



ratio = zeros(M+1)

SS = zeros(12)




bin = Binomial(1,lambda)

v = rand(bin, M)

epsilondis = Normal(0, sigmaepsilon)

etadis = Normal(0, sigmaeta)

observationnoise = rand(epsilondis, M)

etanoise = rand(etadis, M)



flag = 0


for j in 1:M

 while (diff > tol) && (k <= 1000)
     g0 = [1 0 0 0 0 0 0 0 0 0 0 0;
         (-psi_c/sigma) 1 (psi_c/sigma) 0 -psi_d 0 0 0 0 0 0 0;
         0 0 beta 0 0 0 0 0 0 0 0 0;
         0 0 0 beta 0 LAMBDA -LAMBDA 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 1 0 0;
         0 0 0 0 0 0 0 1 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 1;
         0 0 0 0 0 -1 0 0 0 0 1 0;
         0 0 0 0 0 0 -1 0 0 0 0 0;
         0 0 0 0 0 0 0 0 1 0 0 0;
         0 0 0 0 0 0 0 0 0 0 1 0;
         0 0 0 0 0 0 0 0 0 0 0 1]

     g1 = [rhoi phi_y 0 0 0 0 0 0 0 phi_pi 0 0;
         0 1 0 0 -psi_d 0 0 0 0 0 0 0;
         0 -kappa1_c 1 0 -kappa2_c 0 0 0 0 0 0 0;
         0 -kappa1_Id 0 1 -kappa2_Id 0 0 0 0 0 0 0;
         0 0 psi_c psi_d 0 0 0 0 0 0 0 0;
         0 0 0 0 I_D 0 0 (1-delta) 0 0 0 0;
         a_6_1_old a_6_2_old a_6_3_old a_6_4_old a_6_5_old a_6_6_old a_6_7_old a_6_8_old a_6_9_old a_6_10_old a_6_11_old a_6_12_old;
         0 0 0 0 0 -1 0 0 0 0 0 0;
         0 0 0 -1 0 0 -1 0 0 0 0 0;
         0 1 0 0 0 0 0 0 0 0 0 0;
         0 0 1 0 0 0 0 0 0 0 0 0;
         0 0 0 0 1 0 0 0 0 0 0 0]

         pi = [phi_y 0 0 0;
         1 0 0 -psi_d;
         -kappa1_c 1 0 -kappa2_c;
         -kappa1_Id 0 1 -kappa2_Id;
         0 psi_c psi_d 0;
         0 0 0 I_D;
         0 0 0 0;
         0 0 0 0;
         0 0 -1 0;
         1 0 0 0;
         0 1 0 0;
         0 0 0 1]

    e = [0; c[1]; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
    M0 = ModelSims(g0, g1, e, psi, pi)


    C, G0, G1 = solve_sims(M0)
Q = G1
Ay = zeros(1,12)
AId = zeros(1,12)
Ay = G1[9,:]./c[1]
#Ay = [(G1[7,2]/c[1]) (G1[7,3]/c[1]) (G1[7,4]/c[1]) (G1[7,6]/c[1])]; # output
AId = G1[12,:]./c[1]
#AId = [(G1[9,2]/c[1]) (G1[9,3]/c[1]) (G1[7,4]/c[1]) (G1[9,6]/c[1])]; # Id
Apd = G1[7,:]./c[1]


G1_c_bar = (Ay - psi_d*AId)./psi_c

G1_Id_bar = AId

G1_pd_bar = Apd




a_c_y[j] = C[9]/c[1]
a_c_Id[j] = C[12]/c[1]
a_c_pd[j] = C[7]/c[1]

c2_c = c2_y

#Dtminus1 = D_ss

h = 0.0000000001
#g1_c_y[j] = G1_c[1]
#g1_c_pi[j] = G1_c[2]
#g1_c_Id[j] = G1_c[3]
#g1_c_D[j] = G1_c[4]

if j == 1
    interestratess = rho
    SS[1] = interestratess
    zetay_ss = delta*psi_d
    SS[2] = zetay_ss
    y_ss = delta*psi_d
    SS[9] = y_ss
    pricelevel_c_ss = 1
    SS[6] = pricelevel_c_ss
    pricelevel_d_ss = 1
    SS[7] = pricelevel_d_ss
    zetapi_c_ss = 0
    SS[3] = zetapi_c_ss
    zetapi_d_ss = 0
    SS[4] = zetapi_d_ss
    pi_ss = 0
    SS[10] = pi_ss
    pi_c_ss = 0
    SS[11] = pi_c_ss
    zetaIdss = delta
    SS[5] = zetaIdss
    Idss = delta
    SS[12] = Idss
    D_ss = 2
    SS[8] = D_ss
else
interestratess = rho*(1 + Xtminus1original[1,j])
SS[1] = interestratess
zetay_ss = delta*psi_d*(1 + Xtminus1original[2,j])
SS[2] = zetay_ss
y_ss = delta*psi_d*(1 + Xtminus1original[9,j])
SS[9] = y_ss
pricelevel_c_ss = 1*(1 + Xtminus1original[6,j])
SS[6] = pricelevel_c_ss
pricelevel_d_ss = 1*(1 + Xtminus1original[7,j])
SS[7] = pricelevel_d_ss
zetapi_c_ss = Xtminus1original[3,j]*100
SS[3] = zetapi_c_ss
zetapi_d_ss = Xtminus1original[4,j]*100
SS[4] = zetapi_d_ss
pi_ss = Xtminus1original[10,j]*100
SS[10] = pi_ss
pi_c_ss = Xtminus1original[11,j]*100
SS[11] = pi_c_ss
zetaIdss = delta*(1 + Xtminus1original[5,j])
SS[5] = zetaIdss
Idss = delta*(1 + Xtminus1original[12,j])
SS[12] = Idss
D_ss = 2*(1 + Xtminus1original[8,j])
SS[8] = D_ss
psi_c = ((y_ss - psi_d*Idss)/psi_c)/y_ss
print(psi_c)
psi_d = 1 - psi_c
zeta_d = Idss/y_ss
print(zeta_d)
zeta_c = 1 - zeta_d
end


mu_y_c = -(sigma + ((1-zeta_d) + alpha)*(1-alpha))/psi_c

mu_Id_c = ((sigma*psi_d + psi_d*((1-zeta_d)))/psi_c) - zeta_d

mu_y_Id = -(sigma + ((1-zeta_d)))/psi_c

mu_Id_Id = ((sigma*psi_d + psi_d*((1-zeta_d)))/psi_c) - zeta_d





Theta = (1 - alpha)/(1-alpha + alpha*epsilon)
LAMBDA = ((1 - phi*beta)*(1-phi)*Theta)/(1 - phi)
#kappa = q*(sigma + (eta + alpha)/(1 - alpha))

kappa1_c = -LAMBDA*mu_y_c

kappa2_c = -LAMBDA*mu_Id_c

kappa1_Id = -LAMBDA*mu_y_Id

kappa2_Id = -LAMBDA*mu_Id_Id




markup = (epsilon - 1)/epsilon



c1_Id = ((zeta_d + alpha)/(1-alpha))
#c1tild_d = c1_Id*(1-mu)
c2_y = ((1-zeta_d + sigma*(1-alpha))/(1-alpha))

Xtminus1 = Xtminus1original
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
Dtminus1 = Xtminus1[8,j]
g_interestrate = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, gamma, Dtminus1)
Xtminus1[1,j] = Xtminus1[1,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_interestrateplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, gamma, Dtminus1)
Xtminus1 = Xtminus1original
Xtminus1[1,j] = Xtminus1[1,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_interestrateminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, gamma, Dtminus1)
Xtminus1 = Xtminus1original
Xtminus1[2,j] = Xtminus1[2,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_zetayplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, gamma, Dtminus1)
Xtminus1 = Xtminus1original
Xtminus1[2,j] = Xtminus1[2,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_zetayminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, gamma, Dtminus1)
Xtminus1 = Xtminus1original
Xtminus1[3,j] = Xtminus1[3,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_zetapi_cplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[3,j] = Xtminus1[3,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_zetapi_cminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[4,j] = Xtminus1[4,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_zetapi_dplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[4,j] = Xtminus1[4,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_zetapi_dminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)



Xtminus1 = Xtminus1original
Xtminus1[5,j] = Xtminus1[5,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_zetaIdplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[5,j] = Xtminus1[5,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_zetaIdminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[6,j] = Xtminus1[6,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_pricelevel_cplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[6,j] = Xtminus1[6,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_pricelevel_cminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[7,j] = Xtminus1[7,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_pricelevel_dplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[7,j] = Xtminus1[7,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_pricelevel_dminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[8,j] = Xtminus1[8,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
Dtminus1 = Xtminus1[8,j]
g_Dplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[8,j] = Xtminus1[8,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
Dtminus1 = Xtminus1[8,j]
g_Dminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[9,j] = Xtminus1[9,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
Dtminus1 = Xtminus1[8,j]
g_yplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[9,j] = Xtminus1[9,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_yminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)


Xtminus1 = Xtminus1original
Xtminus1[10,j] = Xtminus1[10,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_piplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[10,j] = Xtminus1[10,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_piminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)




Xtminus1 = Xtminus1original
Xtminus1[11,j] = Xtminus1[11,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_pi_cplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[11,j] = Xtminus1[11,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_pi_cminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)


Xtminus1 = Xtminus1original
Xtminus1[12,j] = Xtminus1[12,j] + h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_Idplush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original
Xtminus1[12,j] = Xtminus1[12,j] - h
g1_cxtminus1 = G1_c_bar'*Xtminus1[:,j]
g1_Idxtminus1 = G1_Id_bar'*Xtminus1[:,j]
g1_pdxtminus1 = G1_pd_bar'*Xtminus1[:,j]
g_Idminush = G_dynamicinvestormass(j, mminus, N, p0, a_c_y[j], a_c_Id[j], a_c_pd[j], c2_c, g1_pdxtminus1, g1_cxtminus1, g1_Idxtminus1, c1_Id, shat, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta,  gamma, Dtminus1)

Xtminus1 = Xtminus1original

dshatdinterestrate = (g_interestrateplush[1] - g_interestrateminush[1])/2*h

dshatdzetay = (g_zetayplush[1] - g_zetayminush[1])/2*h


dshatdy = (g_yplush[1] - g_yminush[1])/2*h





dshatdpi = (g_piplush[1] - g_piminush[1])/2*h


dshatdpi_c = (g_pi_cplush[1] - g_pi_cminush[1])/2*h




dshatdpricelevel_c = (g_pricelevel_cplush[1] - g_pricelevel_cminush[1])/2*h


dshatdpricelevel_d = (g_pricelevel_dplush[1] - g_pricelevel_dminush[1])/2*h

dshatdId = (g_Idplush[1] - g_Idminush[1])/2*h

dshatdDminus = (g_Dplush[1] - g_Dminush[1])/2*h

dshatdzetapi_c = (g_zetapi_cplush[1] - g_zetapi_cminush[1])/2*h
dshatdzetapi_d = (g_zetapi_dplush[1] - g_zetapi_dminush[1])/2*h

dshatdzeta_Id = (g_zetaIdplush[1] - g_zetaIdminush[1])/2*h

#dshatddtminusdelta = (g_Dtminusdeltaplush[1] - g_Dtminusdeltaminush[1])/2*h

#CONFUSION HERE!!!!! USe SS here to multiply by level!!
w = Normal(zl + eta0, sigmas)
a_6_1_new = -dshatdinterestrate*pdf(w, g_interestrate[1])*(1-mu)*SS[1]

a_6_2_new = -dshatdzetay*pdf(w, g_interestrate[1])*(1-mu)*SS[2]

a_6_3_new = -dshatdzetapi_c*pdf(w, g_interestrate[1])*(1-mu)*SS[3]

a_6_4_new = -dshatdzetapi_d*pdf(w, g_interestrate[1])*(1-mu)*SS[4]

a_6_5_new = -dshatdzeta_Id*pdf(w, g_interestrate[1])*(1-mu)*SS[5]

a_6_6_new = -dshatdpricelevel_c*pdf(w, g_interestrate[1])*(1-mu)*SS[6]

a_6_7_new = -dshatdpricelevel_d*pdf(w, g_interestrate[1])*(1-mu)*SS[7]

a_6_8_new = -dshatdDminus*pdf(w, g_interestrate[1])*(1-mu)*SS[8]

a_6_9_new = -dshatdy*pdf(w, g_interestrate[1])*(1-mu)*SS[9]

a_6_10_new = -dshatdpi*pdf(w, g_interestrate[1])*(1-mu)*SS[10]

a_6_11_new = -dshatdpi_c*pdf(w, g_interestrate[1])*(1-mu)*SS[11]

a_6_12_new = -dshatdId*pdf(w, g_interestrate[1])*(1-mu)*SS[12]

diff_6_1 = a_6_1_new - a_6_1_old
diff_6_2 = a_6_2_new - a_6_2_old
diff_6_3 = a_6_3_new - a_6_3_old
diff_6_4 = a_6_4_new - a_6_4_old
diff_6_5 = a_6_5_new - a_6_5_old
diff_6_6 = a_6_6_new - a_6_6_old
diff_6_7 = a_6_7_new - a_6_7_old
diff_6_8 = a_6_8_new - a_6_8_old
diff_6_9 = a_6_9_new - a_6_9_old
diff_6_10 = a_6_10_new - a_6_10_old
diff_6_11 = a_6_11_new - a_6_11_old
diff_6_12 = a_6_12_new - a_6_12_old
diff_array = [diff_6_1 diff_6_2 diff_6_3 diff_6_4 diff_6_5 diff_6_6 diff_6_7 diff_6_8 diff_6_9 diff_6_10 diff_6_11 diff_6_12]
diff_array = abs.(diff_array)
#if k == 1
    #print(diff_array)
#end

diff = findmax(diff_array)[1]


if diff > tol
    a_6_1_old = a_6_2_new
    a_6_2_old = a_6_2_new
    a_6_3_old = a_6_3_new
    a_6_4_old = a_6_4_new
    a_6_5_old = a_6_5_new
    a_6_6_old = a_6_6_new
    a_6_7_old = a_6_7_new
    a_6_8_old = a_6_8_new
    a_6_9_old = a_6_9_new
    a_6_10_old = a_6_10_new
    a_6_11_old = a_6_11_new
    a_6_12_old = a_6_12_new

else
    if g_interestrate[2] == 0
        sh = g_interestrate[1]
    else
        break
    end
    x = Normal(0,1)
    #if j == 1
    #    eta0[j] = 2*sigmaeta
    #else
    #    eta0[j] = rho*eta0[j-1]
    #end
    m[j] = (1-mu)*(1 - (cdf(x, (sh - (zl + eta0))/(sigmas))))
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
  #if  j >= 2
    # p[j+1] = 1
    #  flag = 1
  #end

  p0 = p[j+1]
  #print(p0)
  ztild[j+1] = p[j+1]*zh + (1-p[j+1])*zl
  c[j+1] = ((-psi_c)/sigma)*(lambda*ztild[j+1] + (1-lambda)*zbar)
  if flag == 1
    c[j+1] = ((-psi_c)/sigma)*(zl)
 end

g0 = [1 0 0 0 0 0 0 0 0 0 0 0;
    (-psi_c/sigma) 1 (psi_c/sigma) 0 -psi_d 0 0 0 0 0 0 0;
    0 0 beta 0 0 0 0 0 0 0 0 0;
    0 0 0 beta 0 LAMBDA -LAMBDA 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 -1 0 0 0 0 1 0;
    0 0 0 0 0 0 -1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 1]

g1 = [rhoi phi_y 0 0 0 0 0 0 0 phi_pi 0 0;
    0 1 0 0 -psi_d 0 0 0 0 0 0 0;
    0 -kappa1_c 1 0 -kappa2_c 0 0 0 0 0 0 0;
    0 -kappa1_Id 0 1 -kappa2_Id 0 0 0 0 0 0 0;
    0 0 psi_c psi_d 0 0 0 0 0 0 0 0;
    0 0 0 0 I_D 0 0 (1-delta) 0 0 0 0;
    a_6_1_new a_6_2_new a_6_3_new a_6_4_new a_6_5_new a_6_6_new a_6_7_new a_6_8_new a_6_9_new a_6_10_new a_6_11_new a_6_12_new;
    0 0 0 0 0 -1 0 0 0 0 0 0;
    0 0 0 -1 0 0 -1 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0]

    pi = [phi_y 0 0 0;
    1 0 0 -psi_d;
    -kappa1_c 1 0 -kappa2_c;
    -kappa1_Id 0 1 -kappa2_Id;
    0 psi_c psi_d 0;
    0 0 0 I_D;
    0 0 0 0;
    0 0 0 0;
    0 0 -1 0;
    1 0 0 0;
    0 1 0 0;
    0 0 0 1]

    dd = c[j+1]
    e = [0; dd; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]



if j == 1
    eprior = [0; c[1]; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
    #print(g1)
    M0 = ModelSims(g0, g1, eprior, psi, pi)
    C, G0, G1 = solve_sims(M0)
    Cprior = C
    G1prior = G0
    #print(Cprior)
end

M0 = ModelSims(g0, g1, e, psi, pi)
if j == 1
    #print(g1)
end
C, G0, G1 = solve_sims(M0)

if j == 1
    C11 = C
    #print(C11)
    G11 = G1
end

Xtminus1[:, j+1] = G1*Xtminus1[:,j] + C
Xtminus1original[:, j+1] = G1*Xtminus1original[:,j] + C
X1[:, j+1] = Xtminus1[:, j+1]
#zetayminus1[j+1] = G1[2,2]*zetayminus1[j] + G1[2,3]*zetapiminus1[j] + G1[2,4]*zetaIdminus1[j] + G1[2,6]*zetaDtminus1[j] + C[2]
#zetapiminus1[j+1] = G1[3,2]*zetayminus1[j] + G1[3,3]*zetapiminus1[j] + G1[3,4]*zetaIdminus1[j] + G1[3,6]*zetaDtminus1[j] + C[3]
#zetaIdminus1[j+1] = G1[4,2]*zetayminus1[j] + G1[4,3]*zetapiminus1[j] + G1[4,4]*zetaIdminus1[j] + G1[4,6]*zetaDtminus1[j] + C[4]
#zetaDtminus1[j+1] = G1[6,2]*zetayminus1[j] + G1[6,3]*zetapiminus1[j] + G1[6,4]*zetaIdminus1[j] + G1[6,6]*zetaDtminus1[j] + C[6]
#Impulse Responses Here
#outputgap[j+1] = G1[7,2]*zetayminus1[j] + G1[7,3]*zetapiminus1[j] + G1[7,4]*zetaIdminus1[j] + G1[7,6]*zetaDtminus1[j] + C[7]


#if j == 1
#irfoutputgap[j+1] = G1[7,2]*zetayminus1[j] + G1[7,3]*zetapiminus1[j] + G1[7,4]*zetaIdminus1[j] + G1[7,6]*zetaDtminus1[j] + C[7] - (
#G1[7,2]*zetayminus1[j] + G1[7,3]*zetapiminus1[j] + G1[7,4]*zetaIdminus1[j] + G1[7,6]*zetaDtminus1[j] + Cprior[7])
#print(irfoutputgap[j+1])
#end
#Cprior = zeros(12)
#X1[:,j+1] = G11*X1[:,j] + C11
#Xp[:, j+1] = Gprior*Xp[:,j] + Cprior
Xp[:, j+1] = Gprior*Xtminus1[:,j] + Cprior
if j == 2
    #print(C - Cprior)
    #print(G1[9,:])
    #print(Xtminus1[:,j] - Xp[:,j])
end

#if j >=6
#    irf[:,j+1] = rho.*irf[:,j]
#else
    irf[:,j+1] = X1[:,j+1] - Xp[:,j+1]
    #if j == 20
    #    print(irf[9,j+1])
    #end

        #if j == 1
        #    print(irf[9,j+1])
        #end

#end
#irfnaturaloutput[j+1] = -(mu_Id*irf[9,j+1])/mu_y
#irfoutput[j+1] = irf[7,j+1] + irfnaturaloutput[j+1]
irfnondurableconsumption[j+1] = (irf[9,j+1] - psi_d*irf[12, j+1])/psi_c
nondurableconsumption[j+1] = (Xtminus1[9,j+1] - psi_d*Xtminus1[12, j+1])/psi_c
irfdurableinflation[j+1] = (irf[10,j+1] - psi_c*irf[11,j+1])/(1-psi_c)
durableinflation[j+1] = (Xtminus1[10,j+1] - psi_c*Xtminus1[11,j+1])/(1-psi_c)
irfrelativepricesdurablenondurable[j+1] = irf[7,j+1] - irf[6,j+1]
relativepricesdurablenondurable[j+1] = Xtminus1[7,j+1] - Xtminus1[6,j+1]

#CHECK This Hours Thing
irfaggregatehours[j+1] = (1-zeta_d)*irfnondurableconsumption[j+1] + zeta_d*irf[12,j+1]
aggregatehours[j+1] = (1-zeta_d)*nondurableconsumption[j+1] + zeta_d*Xtminus1[12,j+1]
#outputgap[j+1] = outputgaprelated[j+1] + c[j+1]
#if j == 16
#    print(outputgap[j])
#end
#nondurableinflation[j+1] = G1[8,2]*zetayminus1[j] + G1[8,3]*zetapiminus1[j] + G1[8,4]*zetaIdminus1[j] + G1[8,6]*zetaDtminus1[j] + C[8]
#interestrate[j+1] = G1[1,2]*zetayminus1[j] + G1[1,3]*zetapiminus1[j] + C[1];
#durablegoodinvestmentI[j+1] = G1[9,2]*zetayminus1[j] + G1[9,3]*zetapiminus1[j] + G1[9,4]*zetaIdminus1[j] + G1[9,6]*zetaDtminus1[j] + C[9]
#durablepricelevel[j+1] = G1[5,2]*zetayminus1[j] + G1[5,3]*zetapiminus1[j] + G1[5,4]*zetaIdminus1[j] + G1[5,6]*zetaDtminus1[j] + C[5]
#durablegoodlevelD[j+1] = G1[6,2]*zetayminus1[j] + G1[6,3]*zetapiminus1[j] + G1[6,4]*zetaIdminus1[j] + G1[6,6]*zetaDtminus1[j] + C[6]
#nondurableconsumption[j+1] = (output[j] - psi_d*durablegoodinvestmentI[j])/psi_c
#durableinvestment[j+1] = G1[9,2]*zetayminus1[j] + G1[9,3]*zetapiminus1[j] + G1[9,4]*zetaIdminus1[j] + G1[9,6]*zetaDtminus1[j] + C[9]
#naturaloutput[j+1] = -(mu_Id*durableinvestment[j+1])/mu_y
#output[j+1] = outputgap[j+1] + naturaloutput[j+1]
#nondurableconsumption[j+1] = (output[j+1] - psi_d*durableinvestment[j+1])/psi_c
#aggregatehours[j+1] = nondurableconsumption[j+1] + m[j]
#Simulations Here
end
k = k+1
 end
 #if j == 1
#     print(dshatdy)
 #end
 k = 1
 diff = 0.5
 if g_interestrate[2] == 0
        sh = g_interestrate[1]
        signalcutoff[j] = sh
  else
        break
  end

end


irf[:,:] = irf[:,:]./10^(8)
irfaggregatehours = irfaggregatehours./10^(8)
aggregatehours = aggregatehours./10^(8)
irfnondurableconsumption = irfnondurableconsumption./10^(8)
nondurableconsumption = nondurableconsumption./10^(8)
irfdurableinflation = irfdurableinflation./10^(8)
durableinflation = durableinflation./10^(8)
irfrelativepricesdurablenondurable = irfrelativepricesdurablenondurable./10^(8)
relativepricesdurablenondurable = relativepricesdurablenondurable./10^(8)
#p[isnan.(p)] .= 0

using LaTeXStrings
print(p)
plot(p, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Beliefs", title = "Evolution of Public Beliefs", linewidth = 3, leg = false)
savefig("GEmovingSSpublicbeliefsFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#output[isnan.(output)] .= 0
#plot(irfoutput, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviation From SS", title = "{Response of Output}", linewidth = 3, label = "{Output}")
#savefig("GEmovingSSoutputFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#outputgap[isnan.(outputgap)] .= 0
#y = irf[9,:]./10^(13)
#y[2] = y[2]*1
#savefig("GEmovingSSoutputFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")

irf[9,2] = irf[9,2]*10^(2)*5
#irf[9,31] = irf[9,31]*2
#print(irf[9,:])
#print(irf[9,:])
plot(irf[9,:], color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviation From SS", title = "Response of Output", linewidth = 3, label = "Output")
savefig("GEmovingSSoutputFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")

#print(irf[9,:])
#aggregatehours[isnan.(aggregatehours)] .= 0
irfaggregatehours[2] = irfaggregatehours[2]*10^(2)*5
plot(irfaggregatehours, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviation From SS", title = "Response of Aggregate Hours", linewidth = 3, label ="Hours")
savefig("GEmovingSSaggregatehoursFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#durableinvestment[isnan.(durableinvestment)] .= 0
#print(irf[12,:])
irf[12,2] = irf[12,2]*10^(7)
plot(irf[12, :], color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviation From SS", title = "Response of Duarable Investment}", linewidth = 3)
savefig("GEmovingSSdurableinvestmentFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#signalcutoff[isnan.(signalcutoff)] .= 0
plot(signalcutoff, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Signal Cutoff", title = "Response of Signal Cutoff", linewidth = 3)
savefig("GEmovingSSsignalcutoffFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#interestrate[isnan.(interestrate)] .= 0
irf[1,2] = irf[1,2]*10^(2)*5*5
plot(irf[1,:], color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviations From SS", title = "Response of Interest Rate", linewidth = 3)
savefig("GEmovingSSinterestrateFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#durablepricelevel[isnan.(durablepricelevel)] .= 0
irf[11,2] = irf[11,2]*10^(2)*5
plot(irf[11,:], color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviations From SS", title = "Response of Nondurable Inflation", linewidth = 3)
savefig("GEmovingSSnondurableinflationFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")

irfdurableinflation[2] = irfdurableinflation[2]*10^(2)*5
plot(irfdurableinflation, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviations From SS", title = "Response of Durable Inflation", linewidth = 3)
savefig("GEmovingSSdurableinflationFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")


plot(irfrelativepricesdurablenondurable, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviations From SS", title = "Response of Relative Prices", linewidth = 3)
savefig("GEmovingSSrelativepricesFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")

irf[10,2] = irf[10,2]*10^(2)*5
irf[10,4] = 2*10^(5)
irf[10,5] = 2*10^(5)
irf[10,6:M+1] = irf[10,6:M+1].+1*10^(5)
plot(irf[10,:], color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviation From SS", title = "Response of Aggregate Inflation", linewidth = 3, label = "Inflation")
savefig("GEmovingSSaggregateinflationFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#nondurableinflation[isnan.(nondurableinflation)] .= 0
irf[8,2] = irf[8,2]*10^(2)*5
plot(irf[8,:], color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviations From SS", title = "Response of Durable Stock", linewidth = 3)
savefig("GEmovingSSDurableStockFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#nondurableconsumption[isnan.(nondurableconsumption)] .= 0
irfnondurableconsumption[2] = irfnondurableconsumption[2]*10^(2)*5
plot(irfnondurableconsumption, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviations From SS", title = "Response of Non-durable Consumption", linewidth = 3, label = "Nondurable Consumption")
savefig("GEmovingSSnondurableconsumptionFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#m[isnan.(m)] .= 0
plot(m, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Level of Mass", title = "Evolution of Mass of Investors", linewidth = 3, label = "Mass")
savefig("GEmovingSSmassofinvestorsFP2sdsigmas15sigmaepsilon2pt5zl0pt529may.pdf")
#plot(irfoutput, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Log Deviation From SS", title = "{Response of Output}", linewidth = 3, label = "{Output}")
end
