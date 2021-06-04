%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m 
%
% Created June, 2021
% Elyes Balti
% The University of Texas at Austin
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the following papers: 
%
% N. Mensi, D. B. Rawat and E. Balti
% "PLS for V2I Communications Using Friendly Jammer and Double kappa-mu Shadowed Fading,"
% 2021 IEEE International Conference on Communications (ICC)
%
% Contact email: ebalti@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
% This script produces the analytical and Monte Carlo simulations of the Cumulative Distribution Function (CDF)
% of Double Kappa-Mu Shadowed distribution. 
%% Parameters 
% T: Threshold
% p,q: vectors containing  the mean values of the in-phase and the
% quadrature phase components of the multipath clusters
% mu: number of multipath clusters
% Mc: number of Monte Carlo iterations
% L: Upper limit of the summation in the CDF expression. Originally the
% summation in the CDF expression is infinite, however, we select an upper
% limit to approximate the CDF and evaluate a finite summation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;close all
T=0:2:20;
T = db2pow(T);
Mc = 1e5;
L = 20;
mu = 4;
p = ones(mu,1);
q = ones(mu,1);
d = sqrt(sum(p.^2 + q.^2));
kappa = d^2/(2*mu);
%% Gaussian Components (X and Y are independent)
X = randn(mu,Mc);
Y = randn(mu,Mc);

%% Inverse Nakagami-m
ms=3;
A = gamrnd(ms,1/ms,Mc,1); A = 1./A; A = A/mean(A); A = sqrt(A);

%% Nakagami-m
md=2;
G = gamrnd(md,1/md,Mc,1); G = sqrt(G);

%% Fading power
R = zeros(Mc,1);
for kk=1:Mc
R(kk) = A(kk)^2*sum( (X(:,kk) + G(kk)*p  ).^2 + (Y(:,kk) + G(kk)*q ).^2);
end
rhat = mean(R);%average fading power |h|^2, assuming TX power = 1W

%% Initialize the CDF vectors
CDFm = zeros(length(T),1);
CDFa = zeros(length(T),1);

for ii=1:length(T)
%% Monte Carlo    
    for kk=1:Mc

        if R(kk) <T(ii)
            CDFm(ii) = CDFm(ii) + 1;
        end
    end
%% Analytical    
    CDFa(ii)=DoubleKappaMuShadowedCDF(kappa,mu,ms,md,rhat,L,T(ii));
end

CDFm = CDFm/Mc;% Averaging over the number of Monte Carlo iterations
T = pow2db(T);


figure; hold on
plot(T,CDFa,'linewidth',2);
plot(T,CDFm,'*','linewidth',2);
xlabel('Threshold (dB)')
ylabel('Cumulative Distribution Function (CDF)')
legend('Analytical','Monte Carlo','location','northwest')
title('Double \kappa-\mu Shadowed Distribution')

