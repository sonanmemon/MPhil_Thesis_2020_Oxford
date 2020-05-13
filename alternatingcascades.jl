using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.4.0")
using LinearAlgebra, Statistics
using DataFrames, Parameters, Plots, Printf, QuantEcon, Random
gr(fmt = :png);
using Distributions



capT = 200

q = zeros(capT+1)
a = zeros(capT+1)
sigma = zeros(capT+1)




mutable struct Sigma
x::Float64
end



u = capT # I want the array to be of length 10
# create an uninitialized 1D array of size `x`
arrayOfStructs = Array{Sigma, 1}(undef, u) # or `Vector{Coords}(undef, x)`

# initialize the array elements using default constructor
# `Coords(x, y, z)`
for i = 1:u
arrayOfStructs[i] = Sigma(rand())
# or you can use splatting if you already have values
# arrayOfStructs[i] = Coords(rand(3)...)
end




typeof(Sigma)
epsilon = 0.55
P = [1 - epsilon  epsilon
epsilon 1 - epsilon] # stochastic matrix
mc = MarkovChain(P)

omega_vals = Array{Vector}(undef, 1) # sample paths holder
omega = simulate_indices(mc, capT, init = 1)

#omega = cumsum(omega .== 1) ./ (1:capT) # compute state fraction. ./ required for precedence
#omega_vals = omega .-0
q[1] = 0.5


p = 0.6






for k = 1:capT
    if omega[k] == 1
        d = Bernoulli(p)
        sigma[k] = rand(d)
    else
        d = Bernoulli(1-p)
        sigma[k] = rand(d)
    end
end

# convert(Int, aux[i])

# a = dot(M[:,1], M[:, 1])



let
track = 1
for i = 1:capT
    if q[track] <= p || q[track] >= 1-p
                if sigma[track] == 0
                    if q[track] <= p
                      a[track] = 0
                    else
                      a[track] = 1
                  end
               elseif q[track] < 1 - p
                      a[track] = 0
                  else
                      a[track] = 1
                  end
     end
          if a[track] == 0
              q[track+1] = ((1 - epsilon)*(1 - p)*q[track] + epsilon*p*(1-q[track]))/((1-p)*q[track] + (p)*(1-q[track]))
          else
              q[track+1] = ((1 - epsilon)*(p)*q[track] + epsilon*(1 -p)*(1-q[track]))/(p*q[track] + (1-p)*(1-q[track]))
          end
          if q[track+1] > p || q[track+1] < 1-p
              for j in track+1:capT
              q[j+1] = (1 - epsilon)*q[j] + epsilon*(1-q[j])
              track = j+1
              if q[j+1] <= p || q[j+1] >= 1-p || track >= capT
                  break
              end
            end
         else
             track = track +1
        end

        if track >= capT
            break
        end
    end
end

r = zeros(capT)
for j in 1:capT
    r[j] = q[j]
end

k = zeros(capT)
for j in 1:capT
    k[j] = omega[j]
end




t = zeros(12,2)
for j in 1:12
        t[j,1] = r[j]
        t[j,2] = k[j]
    end

using LaTeXStrings
plot(r, color = :red, alpha = 0.6, xlab = "Iterations", ylab = "Public Belief", linewidth = 2, title = L"$Public \: Beliefs \: \: (\epsilon = 0.55)$",
label = "Beliefs")


plot!([p, 1-p], color = :blue, seriestype="hline", label = ["L and H"])
savefig("figurealternatingcascades0.55.pdf")
