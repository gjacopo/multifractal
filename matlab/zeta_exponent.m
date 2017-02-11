%% ZETA_EXPONENT - Zeta exponents model.
%
%% Description
% Comparison between function models for zeta exponents: log-normal, log-Poisson 
% and cascade models.
% 
%% Syntax
%       zeta = zeta_exponent()

%% Function implementation
function zeta = zeta_exponent()

p=0:0.5:20;
l=length(p);

Kolm=1; logP=2; logN=3; cast=4;
zeta=zeros(4,l);

figure, hold on;

% Modele de Kolmogorov K41
zeta(Kolm,:) = p./3.;
plot(p,zeta(Kolm,:), '--k','LineWidth',2);

% Modele de log-Poisson de She-Leveque
zeta(logP,:) = p./9. +2* (1-(2./3.).^(p/3.));
plot(p,zeta(logP,:), '-rd');


% Modele log-normal de KO62
k=0.022;
zeta(logN,:) = (1/3+3*k/2.).*p - k/2.*p.^2;
plot(p,zeta(logN,:), '-bs');

% Modele ("temperature") de cascades de Castaing
T=0.05;
zeta(cast,:) = p/3. .*(1+3*T) ./ (1+p.*T);
plot(p,zeta(cast,:), '-mo');

axis ([0 20 0 5]);
 % axis square;
 h = legend('Kolmogorov','log-Poisson','log-normal','cascade',4); 
xlabel('p'), ylabel('\zeta(p)');
