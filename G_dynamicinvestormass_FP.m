
function r = G_dynamicinvestormass_FP(j, mminus, N, p, a_c_y, a_c_Id, c2_c, y_ss, pi_ss, Idss, D_ss, c1tild_d, shat, g1_c_y, g1_c_pi, g1_c_Id, g1_c_D, psi_c, psi_d, lambda, sigma, sigmaeta, sigmas, sigmaepsilon, mu, zh, zl, zbar, delta, Dtminus1)

    % Solves for the coefficients associated to the Chebychev polynomials. 

    %options = optimset('Display','Iter','TolFun',10^(-3));
    %if equation(shat) == -1000
        %shat = -1000;
    %else
        %shat = fzero(@equation,shat);
        %shat = bisection(@equation, shat + 0.04, shat+0.05, 10^(-5));
    %end
    
    r = answer(shat);
    
    
    
function a = answer(shat)
tol = 1e-3;
maxiter = 1000;
i = 1;
diff = 1000;
ddiff = zeros(maxiter+1);
changediff = 5;
s = zeros(maxiter+1);
s(i) = shat;
while (diff >= tol) && (changediff > 0) && (i <= maxiter)
        s(i+1) = H(s(i));
        diff = abs(s(i+1) - s(i));
        ddiff(i+1) = abs(diff - ddiff(i));
        changediff = ddiff(i+1);
        if changediff < 1e-3
            changediff = 0
        end
        a = [s(i+1) changediff];
        i = i+1;
end
end
    
    
    
    function g = H(shat)
        g = shat - equation(shat);
    end
    
        
    
  
    
 function result = equation(shat)



b_eta = 3*sigmaeta;
a_eta = -3*sigmaeta;
h_eta = (b_eta-a_eta)/N;

etadis = zeros(N);

etatilddis = zeros(N);






for k = 1:1:N
    if k == 1
        etadis(k) = a_eta
    elseif k == N
        etadis(k) = b_eta
    else
        etadis(k) = a_eta + (k-1)*h_eta
    end
end


for k = 1:1:N
    if k == 1
        etatilddis(k) = a_eta
    elseif k == N
        etatilddis(k) = b_eta
    else
        etatilddis(k) = a_eta + (k-1)*h_eta
    end
end



b_epsilon = 3*sigmaepsilon;
a_epsilon = -3*sigmaepsilon;
h_epsilon = (b_epsilon-a_epsilon)/N;
epsilondis = zeros(N);

for k= 1:1:N
    if k == 1
        epsilondis(k) = a_epsilon
    elseif k == N
        epsilondis(k) = b_epsilon
    else
        epsilondis(k) = a_epsilon + (k-1)*h_epsilon
    end
end







function d = dtild(zk, s, etatild)
d = (1-mu)*normcdf((s - (zk+etatild))/(sigmas))
end


function postpublic = etapdfupdated(j, b1)
likelihood_mcondeta = p*pdf('Normal', mminus(j), 1 - normcdf(((shat - (zh + b1))/sigmas)), sigmaeta) +(1-p)*pdf('Normal', mminus(j), 1 - normcdf(((shat - (zl + b1))/sigmas)), sigmaeta);
if j == 2
    postpublic = etapdfupdated1(b1)*likelihood_mcondeta;
else
    postpublic = etapdfupdated(j-1, b1)*likelihood_mcondeta;
end
end


function d = etapdfupdated1(b1)
    d = pdf('Normal', b1, 0, sigmaeta);
end




    function r = funcetapdfupdated(b1)
    if j == 1
    r = etapdfupdated1(b1);
   else
       r = etapdfupdated(j,b1);
    end
    end

   
    
    function postprivate = etapdfupdatedprivate(j, b1, s)
        likelihood_scondeta = p*pdf('Normal', s, zh + b1, sigmas) +(1-p)*pdf('Normal', s, zl + b1, sigmas);
        if j == 1
            postprivate = etapdfupdated1(b1)*likelihood_scondeta;
        else
         postprivate = etapdfupdated(j,b1)*likelihood_scondeta;
        end
    end
    
    
    function r = funcetapdfupdatedprivate(b1, s)
          r = etapdfupdatedprivate(j,b1,s);
    end
       
    
    
 function s = func(zk, eta, etatild, dk, epsilon, s)
s = funcetapdfupdated(eta)*pdf('Normal', (1/mu)*(dtild(dk, s, etatild) + mu*epsilon - normcdf((s - (zk+eta))/(sigmas))*(1-mu)), 0, sigmaepsilon);
end







function s = simpson(p, s, dk)
Sum = 0;
Sum2 = 0;
sum1 = 0;
sum2 = 0;
sum3 = zeros(N+1);
denom = zeros(N+1, N+1);
whole = zeros(N+1, N+1);
  nodes_eta = zeros(N+1);
 nodes_etatild = zeros(N+1);
nodes_epsilon = zeros(N+1);
     for i = 1:1:N+1
     nodes_eta(i) = etadis(1) + (i-1)*h_eta;
     end
     for i = 1:1:N+1
      nodes_etatild(i) = etatilddis(1) + (i-1)*h_eta;
      end
     for i = 1:1:N+1
      nodes_epsilon(i) = epsilondis(1) + (i-1)*h_epsilon;
     end
      for k = 1:1:N+1
         for l = 1:1:N+1
          sum1 = func(zl, nodes_eta(1), nodes_etatild(k), dk, nodes_epsilon(l), s) + func(zl, nodes_eta(N), nodes_etatild(k), dk, nodes_epsilon(l), s);
           for i = 3:2:N-1
           sum1= sum1 +2*func(zl, nodes_eta(i), nodes_etatild(k), dk, nodes_epsilon(l), s);
           end
          for i = 2:2:N
          sum1 = sum1 +4*func(zl, nodes_eta(i), nodes_etatild(k), dk, nodes_epsilon(l), s);
          end
          sum2 = func(zh, nodes_eta(1), nodes_etatild(k), dk, nodes_epsilon(l), s) + func(zh, nodes_eta(N), nodes_etatild(k), dk, nodes_epsilon(l), s);
          for i = 3:2:N-1
          sum2 = sum2 + 2*func(zh, nodes_eta(i), nodes_etatild(k), dk, nodes_epsilon(l), s);
          end
          for i = 2:2:N
          sum2 = sum2 + 4*func(zh, nodes_eta(i), nodes_etatild(k), dk, nodes_epsilon(l), s);
          end
          denom(l, k) = ((sum1*h_eta/3)*(1-p))/((sum2*h_eta/3)*(p)) + 1;
          whole(l, k) = (1/denom(l,k))*pdf('Normal', nodes_epsilon(l), 0, sigmaepsilon);
         end
         Sum = whole(1, k) + whole(N, k);
         for i = 3:2:N-1
         Sum = Sum+2*whole(i,k);
         end
         for i = 2:2:N
         Sum = Sum + 4*whole(i,k);
         end
         sum3(k) = (Sum*h_epsilon/3)*funcetapdfupdatedprivate(nodes_etatild(k), s);
      end
       Sum2 = sum3(1) + sum3(N);
       for i = 3:2:N-1
         Sum2 = Sum2 + 2*sum3(i);
        end
       for i = 2:2:N
         Sum2 = Sum2 + 4*sum3(i);
       end
         s = Sum2*h_eta/3;
end




v = sqrt((sigmaeta)^(2) + (sigmas)^(2));


function q = psi(p,shat)
if j == 1
    if normpdf((shat - zl)/v) == 0 & normpdf((shat - zh)/v) == 0
         q = 1/( 1 + ((1-p))/(p) );
    else
    q = 1/( 1 + ((1-p)*normpdf((shat - zl)/v))/(p*normpdf((shat - zh)/v)) );
    end
else
    q = 1/( 1 + ((1-p)*simpsonprivatebelief(shat,zl))/(p*simpsonprivatebelief(shat,zh)) );
end
end


function f = funcprivatebelief(zk, eta, s)
f = funcetapdfupdatedprivate(eta,s)*pdf('Normal', s, zk+eta, sigmas);
end
 

function s = simpsonprivatebelief(s,zk)
 a_eta = -0.75;
 h_eta = 0.015;
 b_eta = 0.75;
sum = zeros(N);
     nodes= zeros(N+1);
     for i = 1:1:N+1
     nodes(i)= a_eta + (i-1)*h_eta;
     end
     sum = funcprivatebelief(zk, nodes(1), s) + funcprivatebelief(zk, nodes(N), s);
     for i = 3:2:N-1
     sum=sum+2*funcprivatebelief(zk, nodes(i), s);
     end
     for i = 2:2:N
     sum=sum+4*funcprivatebelief(zk, nodes(i), s);
     end
     s = (sum*h_eta)/3;
 end
 





function gamma = Gamma(p, shat)
 gamma = psi(p,shat)*simpson(p, shat, zh) + (1-psi(p,shat))*simpson(p, shat, zl);
 end





function s = funcdf(s, zk, eta)
s = (1 - normcdf((s - (zk+eta))/(sigmas)))*funcetapdfupdatedprivate(eta,s);
end



function s = simpsoncdf(s,zk)
sum = zeros(N);
     nodes= zeros(N+1);
     for i = 1:1:N+1
     nodes(i)= etadis(1) + (i-1)*h_eta;
     end
     sum = funcdf(s, zk, nodes(1)) + funcdf(s, zk, nodes(N));
     for i = 3:2:N-1
     sum=sum+2*funcdf(s, zk, nodes(i));
     end
     for i = 2:2:N
     sum=sum+4*funcdf(s, zk, nodes(i));
     end
     s = (sum*h_eta)/3;
end



%r = a*zl + zl;
t = psi(p,shat)*simpsoncdf(shat, zh) + (1-psi(p,shat))*simpsoncdf(shat, zl);
%t = 0;
ztild = Gamma(p,shat)*zh + (1-Gamma(p,shat))*zl;
cons = (-psi_c/sigma)*(lambda*ztild + (1-lambda)*zbar);
m = (a_c_y*cons -psi_d*a_c_Id*cons)/psi_c;
%m = Gamma(p,shat)*zh + (1-Gamma(p, shat))*zl;
g1_cxminus1 = g1_c_y*cons*y_ss + g1_c_pi*cons*pi_ss + g1_c_Id*cons*Idss + g1_c_D*cons*D_ss;



function c = ctild(p, shat)
c = (c1tild_d*(t) + c2_c*g1_cxminus1 + c2_c*(m) + (1-delta)*Dtminus1 - zl)/(zh - zl);
end

    
   

 %if j == 1
    % if ctild(p,shat) > 0 || ctild(p,shat) < 1
    %result = (log( (1- ctild(p, shat))/(ctild(p, shat))*(p/(1-p))) - d)*e - shat;
    % else
     %result = -1000;
     %end
 %else
  g = simpsonprivatebelief(shat, zl);
 h = simpsonprivatebelief(shat, zh);
 c1 = ctild(p,shat);
 if g == 0 & h == 0
     g = 1;
     h = 1;
     result = ((1- c1)/c1)*(p/(1-p)) - (g/h);
 elseif g == 0 & h~=0
     g = 1;
     h = 1;
     result = ((1- c1)/c1)*(p/(1-p)) - (g/h);
 else
     result = ((1- c1)/c1)*(p/(1-p)) - (g/h);
 end
     
 %end
 
 
 end
end
























