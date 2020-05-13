using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.4.0")
using LinearAlgebra, Statistics
using DataFrames, Parameters, Plots, Printf, QuantEcon, Random
gr(fmt = :png);
using Distributions


function G_dynamicinvestormass(j::Number, mminus::Array, N::Number, p::Number, a_c_y::Number, a_c_Id::Number, c2_c::Number, y_ss::Number, pi_ss::Number, Idss::Number, D_ss::Number, c1tild_d::Number, shat::Number, g1_c_y::Number, g1_c_pi::Number, g1_c_Id::Number, g1_c_D::Number, psi_c::Number, psi_d::Number, lambda::Number, sigma::Number, sigmaeta::Number, sigmas::Number, sigmaepsilon::Number, mu::Number, zh::Number, zl::Number, zbar, delta, Dtminus1)



x = Normal(0,1)

function  equation(s::Number)

    b_eta = 3*sigmaeta
    a_eta = -3*sigmaeta
    h_eta = (b_eta-a_eta)/N

    etadis = zeros(N)

    etatilddis = zeros(N)






    for k in 1:N
        if k == 1
            etadis[k] = a_eta
        elseif k == N
            etadis[k] = b_eta
        else
            etadis[k] = a_eta + (k-1)*h_eta
        end
    end


    for k in 1:N
        if k == 1
            etatilddis[k] = a_eta
        elseif k == N
            etatilddis[k] = b_eta
        else
            etatilddis[k] = a_eta + (k-1)*h_eta
        end
    end



    b_epsilon = 3*sigmaepsilon
    a_epsilon = -3*sigmaepsilon
    h_epsilon = (b_epsilon-a_epsilon)/N
    epsilondis = zeros(N)

    for k in 1:N
        if k == 1
            epsilondis[k] = a_epsilon
        elseif k == N
            epsilondis[k] = b_epsilon
        else
            epsilondis[k] = a_epsilon + (k-1)*h_epsilon
        end
    end







    function dtild(zk::Number, s::Number, etatild::Number)
    mu = 0.1
    d = (1-mu)*cdf(x, (s - (zk+etatild))/(sigmas))
    return d
    end



    function etapdfupdated1(b1::Number)
        d = Normal(0, sigmaeta)
        e = pdf(d, b1)
        return e
    end


    function etapdfupdated(j::Number, b1::Number)
        y = Normal(1 - cdf(x, ((s - (zh + b1))/sigmas)), sigmaeta)
        z = Normal(1 - cdf(x, ((s - (zl + b1))/sigmas)), sigmaeta)
    likelihood_mcondeta = p*pdf(y, mminus[j]) +(1-p)*pdf(z, mminus[j])
    if j == 2
        postpublic = etapdfupdated1(b1)*likelihood_mcondeta

    else
        postpublic = etapdfupdated(j-1, b1)*likelihood_mcondeta

    end
    return postpublic
    end




    function etapdfupdatedprivate(j::Number,b1::Number, s::Number)
        yprivate = Normal(zh + b1, sigmas)
        zprivate = Normal(zl + b1, sigmas)
        likelihood_scondeta = p*pdf(yprivate, s) +(1-p)*pdf(zprivate, s)
        if j == 1
            postprivate = etapdfupdated1(b1)*likelihood_scondeta
        else
         postprivate = etapdfupdated(j,b1)*likelihood_scondeta
            end
    return postprivate
    end






        function funcetapdfupdated(b1::Number)
        if j == 1
        r = etapdfupdated1(b1)
       else
           r = etapdfupdated(j,b1)
        end
        return r
       end


       function funcetapdfupdatedprivate(b1::Number, s::Number)
          r = etapdfupdatedprivate(j,b1,s)
       return r
       end




    function func(zk::Number, eta::Number, etatild::Number, dk::Number, epsilon::Number, s::Number)
    mu = 0.1
    q = Normal(0, sigmaepsilon)
    s = funcetapdfupdated(eta)*pdf(q, (1/mu)*(dtild(dk, s, etatild) + mu*epsilon - cdf(x, (s - (zk+eta))/(sigmas))*(1-mu)))
    return s
    end









    function simpson(p::Number, s::Number, dk::Number)
    dist = Normal(0, sigmaepsilon)
    Sum = 0
    Sum2 = 0
    sum1 = 0
    sum2 = 0
    sum3 = zeros(N+1)
    denom = zeros(N+1, N+1)
    whole = zeros(N+1, N+1)
      nodes_eta = zeros(N+1)
     nodes_etatild = zeros(N+1)
    nodes_epsilon = zeros(N+1)
         for i in 1:N+1
         nodes_eta[i] = etadis[1] + (i-1)*h_eta
         end
         for i in 1:N+1
          nodes_etatild[i] = etatilddis[1] + (i-1)*h_eta
          end
         for i in 1:N+1
          nodes_epsilon[i] = epsilondis[1] + (i-1)*h_epsilon
         end
          for k in 1:N+1
             for l in 1:N+1
              sum1 = func(zl, nodes_eta[1], nodes_etatild[k], dk, nodes_epsilon[l], s) + func(zl, nodes_eta[N], nodes_etatild[k], dk, nodes_epsilon[l], s)
               for i in 3:2:N-1
               sum1= sum1 +2*func(zl, nodes_eta[i], nodes_etatild[k], dk, nodes_epsilon[l], s)
               end
              for i in 2:2:N
              sum1 = sum1 +4*func(zl, nodes_eta[i], nodes_etatild[k], dk, nodes_epsilon[l], s)
              end
              sum2 = func(zh, nodes_eta[1], nodes_etatild[k], dk, nodes_epsilon[l], s) + func(zh, nodes_eta[N], nodes_etatild[k], dk, nodes_epsilon[l], s)
              for i in 3:2:N-1
              sum2 = sum2 + 2*func(zh, nodes_eta[i], nodes_etatild[k], dk, nodes_epsilon[l], s)
              end
              for i in 2:2:N
              sum2 = sum2 + 4*func(zh, nodes_eta[i], nodes_etatild[k], dk, nodes_epsilon[l], s)
              end
              denom[l, k] = ((sum1*h_eta/3)*(1-p))/((sum2*h_eta/3)*(p)) + 1
              whole[l, k] = (1/denom[l,k])*pdf(dist, nodes_epsilon[l])
             end
             Sum = whole[1, k] + whole[N, k]
             for i in 3:2:N-1
             Sum = Sum+2*whole[i,k]
             end
             for i in 2:2:N
             Sum = Sum + 4*whole[i,k]
             end
             sum3[k] = (Sum*h_epsilon/3)*funcetapdfupdatedprivate(nodes_etatild[k], s)
          end
           Sum2 = sum3[1] + sum3[N]
           for i in 3:2:N-1
             Sum2 = Sum2 + 2*sum3[i]
            end
           for i in 2:2:N
             Sum2 = Sum2 + 4*sum3[i]
           end
             q = Sum2*h_eta/3
             return q
    end





    function  funcprivatebelief(zk::Number, eta::Number, s::Number)
    g = Normal(zk+eta, sigmas)
    f = funcetapdfupdatedprivate(eta,s)*pdf(g, s)
    return f
    end


    function  simpsonprivatebelief(s::Number,zk::Number)
     a_eta = -0.75
     h_eta = 0.015
     b_eta = 0.75
    sum = zeros(N)
         nodes= zeros(N+1)
         for i in 1:N+1
         nodes[i] = a_eta + (i-1)*h_eta
         end
         sum = funcprivatebelief(zk, nodes[1], s) + funcprivatebelief(zk, nodes[N], s)
         for i in 3:2:N-1
         sum=sum+2*funcprivatebelief(zk, nodes[i], s)
         end
         for i in 2:2:N
         sum=sum+4*funcprivatebelief(zk, nodes[i], s)
         end
         q = (sum*h_eta)/3
         return q
     end




     function  psi(p::Number,s::Number)
     v = sqrt((sigmaeta)^(2) + (sigmas)^(2))
     standl = (s - zl)/v
     standh = (s - zh)/v
     if j == 1
         #if pdf(x, standl) == 0 & pdf(x, standh) == 0
             # q = 1/( 1 + ((1-p))/(p) )
         #else
         q = 1/( 1 + ((1-p)*pdf(x, standl))/(p*pdf(x, standh )))
     else
         templ = simpsonprivatebelief(s,zl)
         temph = simpsonprivatebelief(s,zh)
         q = 1/( 1 + ((1-p)*templ)/(p*temph))
     end
     return q = convert(Float64, q)
     end


     function  Gamma(p::Number, s::Number)
      gamma = psi(p,s)*simpson(p, s, zh) + (1-psi(p,s))*simpson(p, s, zl)
      gamma = convert(Float64, gamma)
      return gamma
      end


      function  funcdf(s::Number, zk::Number, eta::Number)
      s = 1 - cdf(x, (s - (zk+eta))/(sigmas))
      return s
      end



      function simpsoncdf(s::Number,zk::Number)
      sum = zeros(N)
           nodes= zeros(N+1)
           for i in 1:N+1
           nodes[i] = etadis[1] + (i-1)*h_eta
           end
           sum = funcdf(s, zk, nodes[1]) + funcdf(s, zk, nodes[N])
           for i in 3:2:N-1
           sum=sum+2*funcdf(s, zk, nodes[i])
           end
           for i in 2:2:N
           sum=sum+4*funcdf(s, zk, nodes[i])
           end
           s = (sum*h_eta)/3
           return s
      end



      t = psi(p,shat)*simpsoncdf(shat, zh) + (1-psi(p,shat))*simpsoncdf(shat, zl)
      ztild = Gamma(p,shat)*zh + (1-Gamma(p,shat))*zl
      cons = (-psi_c/sigma)*(lambda*ztild + (1-lambda)*zbar)
      m = (a_c_y*cons -psi_d*a_c_Id*cons)/psi_c
      g1_cxminus1 = g1_c_y*cons*y_ss + g1_c_pi*cons*pi_ss + g1_c_Id*cons*Idss + g1_c_D*cons*D_ss

     c = (c1tild_d*(t) + c2_c*g1_cxminus1 + c2_c*(m) + (1-delta)*Dtminus1 - zl)/(zh - zl)


      #c = (c1tild*(t) + c2*gc_css + atild*(m) - r)/(zh - zl)


      ansl = simpsonprivatebelief(s, zl)
      ansh = simpsonprivatebelief(s, zh)
           result = ((1- c)/c)*(p/(1-p)) - (ansl/ansh)
           result = convert(Float64, result)
           return result
       end



 function  H(s::Number)
     h = s - equation(s)
     h = convert(Float64,h)
     return h
 end

 function answer(shat::Number)
 tol = 10^(-3)
 maxiter = 1000
 i = 1
 c = [1.0 1.0]
 diff = 5.0
 ddiff = zeros(1001)
 changediff = 5.0
 s = zeros(maxiter+1)
 s[i] = shat
 while (diff >= tol) & (changediff > 0) & (i <= maxiter)
          s[i+1] = H(s[i])
         diff = abs(s[i+1] - s[i])
         ddiff[i+1] = abs(diff - ddiff[i])
         changediff = ddiff[i+1]
         if changediff < tol
             changediff = 0
         end
         c = [s[i+1] changediff]
         i = i+1
 end
 #c = convert(Float64, c)
 return c
 end


r = answer(shat)
return r

 end
