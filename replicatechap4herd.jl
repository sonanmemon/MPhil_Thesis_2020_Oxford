using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.4.0")
using LinearAlgebra, Statistics
using DataFrames, Parameters, Plots, Printf, QuantEcon, Random
gr(fmt = :png);
using Distributions



capT = 10000
γ = zeros(capT)
γ[1] = log(3)
theta1 = 1
theta0 = 0
sigmaeta = 0.5
sigmaepsilon = 2
c = 0.5
mu_1 = 0.75
theta = theta0

#eta = zeros(capT)
y = zeros(capT)
d1 = zeros(capT)
d0 = zeros(capT)
q = Normal(0, sigmaepsilon)
function s(b)
temp = ((theta0 + theta1)/2) - (sigmaeta/(theta1 - theta0))*b
return temp
end
q = Normal(0, sigmaepsilon)
r = Normal(0, sigmaeta)
eta = rand(r, capT)
d1[1] =  s(γ[1]) - theta1
d0[1] = s(γ[1]) - theta0
y[1] = 1 - cdf(q, d0[1]) + eta[1]

for t in 1:capT-1
γ[t+1] = γ[t] + log( (pdf( r, y[t] - (1 - cdf(q, d1[t]))))/(pdf(r, y[t] - (1 - cdf(q, d0[t])))) )
d1[t+1] = s(γ[t+1]) - theta1
d0[t+1] = s(γ[t+1]) - theta0
y[t+1] = 1 - cdf(q, d0[t+1]) + eta[t+1]
end


using LaTeXStrings
#plot(γ)
#plot!(ylab = "LLR", xlab= "Iteration", label = "LLR", color = :red, title = L"$\theta = \theta_{0} = 1, \: \eta_{t} = 0$")
#savefig("figurelhighthetafixedeta.pdf")
#plot(γ, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Public LLR", label = "LLR," linewidth = 2, title = L"$\theta = \theta_{0} = 1, \: \eta_{t} = 0$")
plot(γ, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Public LLR", linewidth = 2, title = L"$Public \: LLR \: \: (\theta = 1, \: \eta_{t} = \: Random)$",
label = "LLR")

#savefig("figurehighthetarandometa.pdf")
