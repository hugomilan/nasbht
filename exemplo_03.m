% NASBHT - Numerical and Analytical Solutions of BioHeat Transfer problems
% 
% Copyright (C) 2018. Hugo Fernando Maia Milan
%
% File:   exemplo_03.m
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

% Algoritmo que resolve o exemplo 3.

wb=0.0005;%perfusão sanguínea
pb=1085;%peso específico do sangue
cb=3680;%calor específico do sangue
Tb=37;%temperatura do sangue
k=0.3;%condutividade térmica do tecido
qmet = 700;%produção de calor metabólico
L=0.03;%comprimento do tecido

TC=37;%temperatura na fronteira interna
Tar = 20;% temperatura do ar
rc = 0.02;% resistência

l = sqrt(wb*pb*cb/k);

n = Tb + qmet/(wb*pb*cb);


g = k*l;

s = sinh(l*L);
c = cosh(l*L);
t = tanh(l*L);


x = 0:L/100:L;

TL = (Tar/rc + TC*g/s + n*g*(1/t - 1/s))/(g/t + 1/rc);

xa = x;

Ta = TL*sinh(l*xa)/s ...
    + TC*sinh(l*(L - xa))/s ...
    + n*(1 - (sinh(l*xa) + sinh(l*(L - xa)))/s);

qC = -TL*g/s + TC*g/t + n*g*(1/s - 1/t);

qconv = (TL - Tar)/rc;
plot(xa,Ta)

%% efeito do fluxo de calor
qLq = -2500:2500;
TLq = -qLq*t/g + TC/c + n*(1 - 1/c);
plot(qLq,TLq)