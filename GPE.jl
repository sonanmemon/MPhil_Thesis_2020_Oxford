using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.4.0")
using LinearAlgebra, Statistics
using DataFrames, Parameters, Plots, Printf, QuantEcon, Random
gr(fmt = :png);
using Distributions

function GPE(j::Number, mminus::Array, phat::Number, N::Number, p::Number, sguess::Number, sigmaeta::Number, sigmas::Number, sigmaepsilon::Number, thetah::Number, thetal::Number)



function equationpe(s)

etadis = zeros(N)

b_eta = 3*sigmaeta
a_eta = -3*sigmaeta
h_eta = (b_eta-a_eta)/N

    for k in 1:N
        if k == 1
            etadis[k] = a_eta
        elseif k == N
            etadis[k] = b_eta
        else
            etadis[k] = a_eta + (k-1)*h_eta
        end
    end

x = Normal(0,1)


    function etapdfupdated1(b1)
        etapdf = Normal(0, sigmaeta)
        res = pdf(etapdf, b1)
        return res
    end


    function etapdfupdated(j, b1)
    etapdfh = Normal(1 - cdf(x, ((s - (thetah + b1))/sigmas)), sigmaeta)
    etapdfl = Normal(1 - cdf(x, ((s - (thetal + b1))/sigmas)), sigmaeta)
    likelihood_mcondeta = p*pdf(etapdfh, mminus[j]) + (1-p)*pdf(etapdfl, mminus[j])
    if j == 2
        post = etapdfupdated1(b1)*likelihood_mcondeta
    else
        post = etapdfupdated(j-1, b1)*likelihood_mcondeta
    end
    return post
    end

    function funcetapdfupdated(b1)
        if j == 1
        r = etapdfupdated1(b1)
       else
           r = etapdfupdated(j,b1)
       end
           return r
       end

       function etapdfupdatedprivate(j::Number,b1::Number, s::Number)
           yprivate = Normal(thetah + b1, sigmas)
           zprivate = Normal(thetal + b1, sigmas)
           likelihood_scondeta = p*pdf(yprivate, s) +(1-p)*pdf(zprivate, s)
           if j == 1
               postprivate = etapdfupdated1(b1)*likelihood_scondeta
           else
            postprivate = etapdfupdated(j,b1)*likelihood_scondeta
               end
       return postprivate
       end







        function funcetapdfupdatedprivate(b1::Number, s::Number)
             r = etapdfupdatedprivate(j,b1,s)
          return r
          end





    function func(thetak, eta, s)
        q = Normal(thetak + eta, sigmas)
    s = funcetapdfupdatedprivate(eta, s)*pdf(q, s)
    return s
    end


    function simpson(thetak, s)
    sum = 0
         nodes=Array{Float64}(undef,N+1)
         for i in 1:N+1
         nodes[i] = etadis[1] + (i-1)*h_eta
         end
         sum = func(thetak, etadis[1], s) + func(thetak, etadis[N], s)
         for i in 3:2:N-1
         sum=sum+2*func(thetak, nodes[i], s)
         end
         for i in 2:2:N
         sum=sum+4*func(thetak, nodes[i], s)
         end
         return(sum*h_eta/3)
    end


c = 1/(1 + (((1-p)/p)*simpson(thetal, s))/(simpson(thetah, s)) ) - phat
return c
end

function  H(s::Number)
    h = s - equationpe(s)
    h = convert(Float64,h)
    return h
end

function answer(shat::Number)
tol = 10^(-5)
maxiter = 1000
i = 1
#c = [1.0 1.0]
c = 1
diff = 5.0
ddiff = zeros(1001)
changediff = 5.0
s = zeros(maxiter+1)
s[i] = shat
while (diff >= tol) & (i <= maxiter)
         s[i+1] = H(s[i])
        diff = abs(s[i+1] - s[i])
        ddiff[i+1] = abs(diff - ddiff[i])
        changediff = ddiff[i+1]
        if changediff < tol
            changediff = 0
        end
        #c = [s[i+1] changediff]
        c = s[i+1]
        i = i+1
end
#c = convert(Float64, c)
return c
end


r = answer(sguess)
return r

end
