
using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.4.0")
using LinearAlgebra, Statistics
using DataFrames, Parameters, Plots, Printf, QuantEcon, Random
gr(fmt = :png);
using Distributions

 function simpsonpublicbelief(j::Number, N::Number, zl::Number, zh::Number, shat::Number, zk::Number, m::Number, p::Number, mminus::Array, sigmaeta::Number, sigmas::Number, sigmaepsilon::Number)


x = Normal(0,1)



function etapdfupdated1(b1::Number)
         d = Normal(0, sigmaeta)
         e = pdf(d, b1)
         return e
end

 function etapdfupdated(j::Number, b1::Number)
     y = Normal(1 - cdf(x, ((shat - (zh + b1))/sigmas)), sigmaeta)
     z = Normal(1 - cdf(x, ((shat - (zl + b1))/sigmas)), sigmaeta)
 likelihood_mcondeta = p*pdf(y, mminus[j]) +(1-p)*pdf(z, mminus[j]);
 if j == 2
     post = etapdfupdated1(b1)*likelihood_mcondeta
 else
     post = etapdfupdated(j-1, b1)*likelihood_mcondeta
 end
 return post
 end








     function funcetapdfupdated(j::Number, b1::Number)
     if j == 1
     r = etapdfupdated1(b1)
    else
        r = etapdfupdated(j,b1)
     end
     return r
    end








 function funcpublicbelief(zk::Number, eta::Number, s::Number, m::Number)
mu = 0.1
k = Normal(0, sigmaepsilon)
f = funcetapdfupdated(j,eta)*pdf(k, (1/mu)*(m - cdf(x, (s - (zk+eta))/(sigmas))*(1-mu)))
return f
 end


function  simpsonpublic(zk::Number, shat::Number, m::Number)
 a_eta = -0.75
 h_eta = 0.015
 b_eta = 0.75
nodes= zeros(N+1)
     for i in 1:N+1
     nodes[i] = a_eta + (i-1)*h_eta
     end
     sum = funcpublicbelief(zk, nodes[1], shat, m) + funcpublicbelief(zk, nodes[N], shat, m)
     for i in 3:2:N-1
     sum=sum+2*funcpublicbelief(zk, nodes[i], shat, m)
     end
     for i in 2:2:N
     sum=sum+4*funcpublicbelief(zk, nodes[i], shat, m)
     end
     s = (sum*h_eta)/3
     return s
end








s = simpsonpublic(zk, shat, m)

return s


end
