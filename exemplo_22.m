% NASBHT - Numerical and Analytical Solutions of BioHeat Transfer problems
% 
% Copyright (C) 2018. Hugo Fernando Maia Milan
%
% File:   exemplo_22.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% NASBHT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% NASBHT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with NASBHT.  If not, see <https://www.gnu.org/licenses/>.

% Algoritmo que resolve o exemplo 22.
% Encontrar a temperatura de superfície com o modelo analítico utilizando
% toda a dedução
% 1 músculo
% 2 Gordura

pb=1085;%peso específico do sangue
cb=3680;%calor específico do sangue
Tb=37;%temperatura do sangue

TC1=37;%temperatura na fronteira interna
qL2=250;%fluxo de calor para o meio externo

L1=0.02;L2=0.01;%comprimento do tecido
k1=0.5;k2=0.4;%condutividade térmica do tecido
wb1=0.0005;wb2=0.00035;%perfusão sanguínea
qmet1=700;qmet2=350;%produção de calor metabólico

l1=sqrt(wb1*pb*cb/k1);
l2=sqrt(wb2*pb*cb/k2);

n1=Tb+qmet1/(wb1*pb*cb);
n2=Tb+qmet2/(wb2*pb*cb);


g1 = k1*l1;
g2 = k2*l2;

s1 = sinh(l1*L1); s2 = sinh(l2*L2);
c1 = cosh(l1*L1); c2 = cosh(l2*L2);
t1 = tanh(l1*L1); t2 = tanh(l2*L2);

fat1 = (g1/(c2*s2) + t2*g1 + g2*t1);
fat2 = (g1 + g2*t1/t2)/fat1;
TL2 = - qL2/g2*fat2 ...
   + TC1/s2*g1/c1/fat1 ...
   + n1/s2*g1*(1 - 1/c1)/fat1 ...
   + n2/c2*(1/s2 - 1/t2)*g1/fat1 ...
   + n2*t2*fat2;

TC2 = TL2*c2 + qL2*s2/g2 - n2*(c2 - 1);
qC2 = -TL2*g2/s2 + TC2*g2/t2 + n2*g2*(1/s2 - 1/t2);
TL1 = TC2;
qC1 = -TL1*g1/s1 + TC1*g1/t1 + n1*g1*(1/s1 - 1/t1);

x1 = 0:L1/100:L1;
x2 = 0:L2/100:L2;

T1 = TL1*sinh(l1*x1)/s1 ...
    + TC1*sinh(l1*(L1 - x1))/s1 ...
    + n1*(1 - (sinh(l1*x1) + sinh(l1*(L1 - x1)))/s1);

T2 = TL2*sinh(l2*x2)/s2 ...
    + TC2*sinh(l2*(L2 - x2))/s2 ...
    + n2*(1 - (sinh(l2*x2) + sinh(l2*(L2 - x2)))/s2);

plot([x1 (x2+L1)],[T1 T2])
