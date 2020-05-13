using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.4.0")
using LinearAlgebra, Statistics
using DataFrames, Parameters, Plots, Printf, QuantEcon, Random
gr(fmt = :png);
using Distributions


let

capT = 150



sigmaepsilon = 0.2
sigmas = 0.5
sigmaeta = 0.25
sigmau = 2.5
thetah = 2
thetal = 1.5
x = Normal(0,1)
u = Normal(0, sigmau)











#s = theta + eta + v








#q = 1/(1 + ((1-p[1])/p[1])*pdf(x, (s[1]-thetal)/sqrt(sigmaeta + sigmas))/pdf(x, (s[1]-thetah)/sqrt(sigmaeta + sigmas)))










#for j in 1:capT
#yl[j] = Normal(thetal + eta[j], sigmas)
#yh[j] = Normal(thetah + eta[j], sigmas)
#end






#sleta0 = Normal(thetal + eta0, sigmas)



sleta = Normal(thetal, sigmas)

sheta = Normal(thetah, sigmas)





c = 1.75
phat = (c-thetal)/(thetah-thetal)


P = zeros(capT)
Q = zeros(capT)
P[1] = 0.25
p = P[1]
Q[1] = 0.05
shat = zeros(capT)
#shat[1] = 0.75 - (5/8)*log10(((P[1])/(1-P[1]))*((1-phat)/(phat)))



m = zeros(capT)





Pr = zeros(capT)
R = zeros(capT)

for k in 1:capT
    R[k] = thetal
end





#q = (sigmas)/(sqrt(2*pi)*sigmaepsilon)


using QuadGK






eta = Normal(0, sigmaeta)
N = 500
etadis = rand(eta, N)
etadis = sort(etadis)

b = 3*sigmaeta
a = -3*sigmaeta
h = (b-a)/N

thetal_etadis = zeros(N, 2)

for k in 1:N
    thetal_etadis[k, 1] = thetal
    if k == 1
        thetal_etadis[k, 2] = a
    elseif k == N
        thetal_etadis[k,2] = b
    else
        thetal_etadis[k,2] = a + (k-1)*h
    end
end




thetah_etadis = zeros(N, 2)


for k in 1:N
    thetah_etadis[k, 1] = thetah
    if k == 1
        thetah_etadis[k, 2] = a
    elseif k == N
        thetah_etadis[k,2] = b
    else
        thetah_etadis[k,2] = a + (k-1)*h
    end
end







etapdf = Normal(0, sigmaeta)
epsilonpdf = Normal(0, sigmaepsilon)



#function func(a1, b1, l1, d1)
#s = (1/(sqrt(2*pi)*sigmaeta))*exp(-0.5*(1/(sigmaeta)^(2))*(b1)^(2))*(1/(sqrt(2*pi)*sigmaepsilon))*exp(-0.5*(1/(sigmaepsilon)^(2))*(l1 - 1 + cdf(x, ((d1-(a1+b1))/(sigmas))))^(2))
#return s
#end





















function newtoncotesl(f ,g)
Sum = zeros(N+1)
sum = zeros(N)
 for j in 1:N
    if j == 1 || j == N
        sum[j] = (h/2)*func(thetal, thetal_etadis[j,2], f, g)
    else
        sum[j] = (h/2)*2*func(thetal, thetal_etadis[j,2], f, g)
    end
    Sum[j+1] = Sum[j] + sum[j]
 end
return Sum[N+1]
end

function newtoncotesh(f, g)
Sum = zeros(N+1)
sum = zeros(N)
for j in 1:N
    if j == 1 || j == N
        sum[j] = (h/2)*func(thetah, thetah_etadis[j,2], f, g)
    else
        sum[j] = (h/2)*2*func(thetah, thetah_etadis[j,2], f, g)
    end
Sum[j+1] = Sum[j] + sum[j]
end
return Sum[N+1]
end


pr = zeros(capT)


ratio = zeros(capT)
templlist = zeros(capT)
temphlist = zeros(capT)


mminus = zeros(capT)
sguess = 1
#z = equationpe(1, mminus, phat, N, p, 2.565, sigmaeta, sigmas, sigmaepsilon, thetah, thetal)
  #y = GPE(1, mminus, phat, N, p, sguess, sigmaeta, sigmas, sigmaepsilon, thetah, thetal)


 #Q = 40
  #y = zeros(Q)

  #for k in 1:Q
#      y[k] = equationpe(1, mminus, phat, N, p, k-21, sigmaeta, sigmas, sigmaepsilon, thetah, thetal)
  #end
  #plot(y)
 #end



eta0 = 1.5*sigmaeta
s = zeros(capT)
for j in 1:capT
    s[j] = thetal + eta0
end



for j in 1:capT-1
shat[j] = GPE(j, mminus, phat, N, p, sguess, sigmaeta, sigmas, sigmaepsilon, thetah, thetal)
sguess = shat[j]
m[j] = 1 - cdf(x, (shat[j]-(thetal + eta0))/(sigmas))
dm = m[j]
if j == 1
   mminus[j] = 0
else
   mminus[j] = m[j-1]
end
templ = simpsonpublicPE(j, p, mminus, thetah, thetal, N,thetal, sguess, dm, sigmaeta, sigmaepsilon, sigmas)
temph = simpsonpublicPE(j, p, mminus, thetah, thetal, N,thetah, sguess, dm, sigmaeta, sigmaepsilon, sigmas)
templlist[j] = templ
temphlist[j] = temph
ratio[j] = templ/temph
pr[j] = 1/(1 + ((1-P[j])/(P[j])*pdf(x, (s[j] - thetal)/sqrt(sigmaeta + sigmas)))/(pdf(x, (s[j] - thetah)/sqrt(sigmaeta + sigmas))))
Pr[j] = (P[j]*pdf(u, thetal - thetah))/(P[j]*pdf(u, thetal - thetah) + (1-P[j])*pdf(u, thetal - thetal))
P[j+1] = 1/( 1 + ((1-Pr[j])*templ)/(Pr[j]*temph))
p = P[j+1]
end
using LaTeXStrings
plot(P, color = :red, alpha = 0.6, xlab = L"\textbf{Iterations}", ylab = L"\textbf{Beliefs}", title = L"\textbf{Evolution of Public Beliefs}", linewidth = 3, label = L"\textbf{Beliefs}")
#savefig("PEpublicbeliefsFP2sd.pdf")
end


#using LaTeXStrings
#plot(r, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Public Belief", linewidth = 2, title = L"$Public \: Beliefs \: \: (\epsilon = 0.55)$",
#label = "Beliefs")


#plot!([p, 1-p], color = :blue, seriestype="hline", label = ["L and H"])
#savefig("figurealternatingcascades0.55.pdf")

#savefig("massofinvestorsPEcorrected2sd.pdf")


#savefig("evolutionofpublicbeliefsPEcorrectedpt1sd.pdf")


#savefig("evolutionofprivatebeliefsPEcorrectedpt1sd.pdf")












#d = Q[1]/P[1]
#dtild = d*exp((0.19*16)
#shat[1] = log((1-dtild)/atild)/(-0.8*(16/5))

#for j in 1:capT-1
#+ ( ((Q[j])/P[j])*pdf(x, (m[j]-(1-cdf(smedeta,shat[j])))/(sqrt((sigmaeta)^(2) + (sigmaepsilon)^(2)))))/(pdf(x, (m[j]-(1-cdf(sheta, shat[j])))
#/(sqrt((sigmaeta)^(2) + (sigmaepsilon)^(2))))) )
#Q[j+1] = 1/( 1 + ( ((1-P[j] - Q[j])/Q[j])*pdf(x, (m[j]-(1-cdf(sleta,shat[j])))/(sqrt((sigmaeta)^(2) + (sigmaepsilon)^(2)))))/(pdf(x, (m[j]-(1-cdf(smedeta, shat[j])))/(sqrt((sigmaeta)^(2) + (sigmaepsilon)^(2)))))
#+ ( ((P[j])/Q[j])*pdf(x, (m[j]-(1-cdf(sheta,shat[j])))/(sqrt((sigmaeta)^(2) + (sigmaepsilon)^(2)))))/(pdf(x, (m[j]-(1-cdf(smedeta, shat[j])))
#/(sqrt((sigmaeta)^(2) + (sigmaepsilon)^(2))))) )
#a = (1 - P[j+1] - Q[j+1])/P[j+1]
#d = Q[j+1]/P[j+1]
#atild = a*exp((0.75*16)/5)
#dtild = d*exp((0.19*16)/5)
#shat[j+1] = log((1-dtild)/atild)/(-0.8*(16/5))
#e = 1 - cdf(sleta, shat[j+1])
#m[j+1] = e
#end
