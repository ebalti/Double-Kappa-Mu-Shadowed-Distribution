%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DoubleKappaMuShadowedCDF.m 
%
% Created June, 2021
% Elyes Balti
% The University of Texas at Austin
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the following paper: 
%
% N. Mensi, D. B. Rawat and E. Balti
% "PLS for V2I Communications Using Friendly Jammer and Double kappa-mu Shadowed Fading,"
% 2021 IEEE International Conference on Communications (ICC)
%
% Contact email: ebalti@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
% This script returns the analytical expression of the Cumulative
% Distribution Function (CDF) of the Double Kappa-Mu Shadowed distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = DoubleKappaMuShadowedCDF(kappa,mu,ms,md,P,L,x)
K = mu*(kappa+1);
K1 = K/(md+mu*kappa);
F = (md/(md+kappa*mu))^md*(K*x/(P*(ms-1))).^mu;

S = 0;
for ii=0:L
   F1 = (K1*mu*kappa*x./(P*(ms-1))).^ii;
   Md = gamma(md+ii)/gamma(md);
   Mimu = gamma(ii+mu+ms)/gamma(ii+mu);
   F2 = Md*Mimu/(factorial(ii)*gamma(ms)*(ii+mu));
   F3 = hypergeom([ii+mu ii+mu+ms],ii+mu+1,-K*x/(P*(ms-1)));
   S = S + F1.*F2.* F3;
end
output = F.*S;
end
