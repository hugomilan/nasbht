% NASBHT - Numerical and Analytical Solutions of BioHeat Transfer problems
% 
% Copyright (C) 2018. Hugo Fernando Maia Milan
%
% File:   exemplo_14.m
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

% Algoritmo que resolve o exemplo 14.
% duas camadas. Músculo (1 e 2) e gordura (3 e 4).

% dados iniciais:
pb = 1085; % peso especifico do sangue
cb = 3680; % calor especifico do sangue
Tb = 37; % temperatura do sangue

TC = 37; % temperatura na condicao de fronteira interna
qL = 250; % fluxo de calor na condicao de fronteira externa

deltat = 0.1;

Ltotal = [20 10]*1e-3; % comprimento da camada
k = [0.5 0.5 0.4 0.4]; % condutividade termica
wb = [0.0005 0.0005 0.00035 0.00035]; % perfusão sanguinea
qmet = [700 700 350 350]; % geracao de calor metabolico
p = [1000 1000 1000 1000]; % peso especifico
c = [3600 3600 3000 3000]; % calor especifico

deltax = [2 1 0.75 0.5]*1e-3; % discretizacao espacial dos meios
Lmeios = [9 11 5.25 4.75]*1e-3; % comprimento dos meios

% número de passos-de-tempo e tempo total da simulação
ktotal = 100000;
ttotal = k*deltat;

% calculos inicias
NNOS(1) = Lmeios(1)/deltax(1) - 0.5; % Eq. 70
NNOS(2) = Lmeios(2)/deltax(2); % Eq. 74
NNOS(3) = Lmeios(3)/deltax(3); % Eq. 74
NNOS(4) = Lmeios(4)/deltax(4) + 0.5; % Eq. 72
NNOSTOTAL = sum(NNOS); %numero total de nos
meios = [ones(1,NNOS(1)) 2*ones(1,NNOS(2))...
    3*ones(1,NNOS(3)) 4*ones(1,NNOS(4))];
% vetor que contem o meio de cada no

% vetor que contém a posição dos nós
xpos = [deltax(1):deltax(1):(deltax(1)*NNOS(1))];

xpos = [xpos (xpos(end) + deltax(1)/2 + deltax(2)/2):deltax(2):...
    (xpos(end) + deltax(1)/2 + deltax(2)/2 + deltax(2)*(NNOS(2) - 1))];

xpos = [xpos (xpos(end) + deltax(2)/2 + deltax(3)/2):deltax(3):...
    (xpos(end) + deltax(2)/2 + deltax(3)/2 + deltax(3)*(NNOS(3) - 1))];

xpos = [xpos (xpos(end) + deltax(3)/2 + deltax(4)/2):deltax(4):...
    (xpos(end) + deltax(3)/2 + deltax(4)/2 + deltax(4)*(NNOS(4) - 1))];

Cdx = p.*c; % Eq. 47
G = wb.*pb.*cb.*deltax; % Eq. 48
IF = (wb.*pb.*cb.*Tb + qmet).*deltax; % Eq. 49
R = deltax./k; % Eq. 50
Zx = deltat./(Cdx.*deltax); % Eq. 54
ZI = Zx./(Zx.*G + 2); % Eq. 53

% matriz de espalhamento
for meio = 1:4
    S{meio} = 1/(Zx(meio)*G(meio) + 2)*...
        [-Zx(meio)*G(meio) 2;
            2       -Zx(meio)*G(meio)];
end
%coeficientes de reflexão e transmissao
% ro(no, 1 ou 2). 1: no: numero do no; 1 porta da esquerda, 2 porta da direita
n = 1;
ro(1,1) = R(1)/(R(1) + 2*Zx(1)); % Eq. 59
tau(1,1) = 2*Zx(1)/(R(1) + 2*Zx(1)); % Eq. 60
    
ro(1,2) = R(1)/(R(1) + 2*Zx(1)); % Eq. 59
tau(1,2) = 2*Zx(1)/(R(1) + 2*Zx(1)); % Eq. 60

for n = 2:(NNOSTOTAL-1)
    
    %porta esquerda
    ro(n,1) = (R(meios(n)) + R(meios(n-1)) - 2*(Zx(meios(n)) - Zx(meios(n-1))))/...
        (R(meios(n)) + R(meios(n-1)) + 2*(Zx(meios(n)) + Zx(meios(n-1))));
    % Eq. 57
    tau(n,1) = 4*Zx(meios(n))/...
        (R(meios(n)) + R(meios(n-1)) + 2*(Zx(meios(n)) + Zx(meios(n-1))));
    % Eq. 58
    
    %porta direita
    ro(n,2) = (R(meios(n)) + R(meios(n+1)) - 2*(Zx(meios(n)) - Zx(meios(n+1))))/...
        (R(meios(n)) + R(meios(n+1)) + 2*(Zx(meios(n)) + Zx(meios(n+1))));
    % Eq. 57    
    tau(n,2) = 4*Zx(meios(n))/...
        (R(meios(n)) + R(meios(n+1)) + 2*(Zx(meios(n)) + Zx(meios(n+1))));
    % Eq. 58
end

ro(end,1) = R(end)/(R(end) + 2*Zx(end)); % Eq. 59
tau(end,1) = 2*Zx(end)/(R(end) + 2*Zx(end)); % Eq. 60


Vrp = ((2 - Zx(1)*G(1))*TC + Zx(1)*IF(1))/4; %tensão refletida no pseudo-no Eq. 62
for n = 1:4
    Vi(n) = (Zx(n)*G(n) + 2)*TC/4 - Zx(n)/4*IF(n); %tensão incidente inicial % Eq. 67
end
Vim = (Zx(4)*G(4) + 2)*TC/4 - Zx(4)/2*(IF(4) - qL); %tensão incidente inicial no meio nó Eq. 68

xt = zeros(NNOSTOTAL,2,2,ktotal); %matriz que contem os valores das tensões incidente (x,1,p,t) e refletidas (x,2,p,t) ao longo do tempo
% na porta da esquerda (x,ir,1,t) e da direita (x,ir,2,t)
T = zeros(NNOSTOTAL,ktotal); %matriz que contem os valores das temperaturas a longo dos nos e do tempo

% primeiro passo-de-tempo

% no 1
% cálculo das tensões incidentes
% porta 1
xt(1,1,1,1) = Vi(1);
% porta 2
xt(1,1,2,1) = Vi(1);

% cálculo das tensões refletidas
% porta 1
xt(1,2,1,1) = S{meios(1)}(1,1)*xt(1,1,1,1) + S{meios(1)}(1,2)*xt(1,1,2,1) + ZI(1)*IF(1);
% porta 2
xt(1,2,2,1) = S{meios(1)}(2,1)*xt(1,1,1,1) + S{meios(1)}(2,2)*xt(1,1,2,1) + ZI(1)*IF(1);

% cálculo da temperatura no nó
T(1,1) = 2*(xt(1,1,1,1) + xt(1,1,2,1))/(Zx(1)*G(1) + 2) + ...
    ZI(1)*IF(1);


% no 2 até NNOSTOTAL -1
for x = 2:(NNOSTOTAL - 1)
    % cálculo das tensões incidentes
    % porta 1
    xt(x,1,1,1) = Vi(meios(x));
    % porta 2
    xt(x,1,2,1) = Vi(meios(x));

    % cálculo das tensões refletidas
    % porta 1
    xt(x,2,1,1) = S{meios(x)}(1,1)*xt(x,1,1,1) + S{meios(x)}(1,2)*xt(x,1,2,1) + ZI(meios(x))*IF(meios(x));
    % porta 2
    xt(x,2,2,1) = S{meios(x)}(2,1)*xt(x,1,1,1) + S{meios(x)}(2,2)*xt(x,1,2,1) + ZI(meios(x))*IF(meios(x));
    
    % cálculo da temperatura no nó
    T(x,1) = 2*(xt(x,1,1,1) + xt(x,1,2,1))/(Zx(meios(x))*G(meios(x)) + 2) + ...
     + ZI(meios(x))*IF(meios(x));
end

% meio nó
% cálculo da tensão incidente
xt(end,1,1,1) = Vim;

% cálculo da tensão refletida
xt(end,2,1,1) = xt(end,1,1,1)*(2 - Zx(end)*G(end))/(Zx(end)*G(end) + 2) ...
    + 2*ZI(end)*(IF(end)/2 - qL);

% cálculo da temperatura
T(end,1) = 4*xt(end,1,1,1)/(Zx(end)*G(meios(end)) + 2) + ...
    2*ZI(end)*(IF(end)/2 - qL);

% passo-de-tempo 2 até ktotal
for t = 2:ktotal
    % nó 1
    % cálculo da tensão incidente
    % porta 1
    xt(1,1,1,t) = xt(1,2,1,t-1)*ro(1,1) + Vrp*tau(1,1);
    %porta 2
    xt(1,1,2,t) = xt(1,2,2,t-1)*ro(1,2) + xt(2,2,1,t-1)*tau(1,2);
    
    % cálculo das tensões refletidas
    % porta 1
    xt(1,2,1,t) = S{meios(1)}(1,1)*xt(1,1,1,t) + S{meios(1)}(1,2)*xt(1,1,2,t) + ZI(1)*IF(1);
    % porta 2
    xt(1,2,2,t) = S{meios(1)}(2,1)*xt(1,1,1,t) + S{meios(1)}(2,2)*xt(1,1,2,t) + ZI(1)*IF(1);
    % pseudo-nó
    Vrp = TC - (Vrp*ro(1,1) + xt(1,2,1,t-1)*tau(1,1));
    
    % cálculo da temperatura no nó
    T(1,t) = 2*(xt(1,1,1,t) + xt(1,1,2,t))/(Zx(1)*G(1) + 2) + ...
        Zx(1)/(Zx(1)*G(1) + 2)*IF(1);

    for x = 2:(NNOSTOTAL - 1)
        % cálculo da tensão incidente
        % porta 1
        xt(x,1,1,t) = xt(x,2,1,t-1)*ro(x,1) + xt(x-1,2,2,t-1)*tau(x,1);
        % porta 2
        xt(x,1,2,t) = xt(x,2,2,t-1)*ro(x,2) + xt(x+1,2,1,t-1)*tau(x,2);
        
        % cálculo das tensões refletidas no passo-de-tempo t
        % porta 1
        xt(x,2,1,t) = S{meios(x)}(1,1)*xt(x,1,1,t) + S{meios(x)}(1,2)*xt(x,1,2,t) + ZI(meios(x))*IF(meios(x));
        % porta 2
        xt(x,2,2,t) = S{meios(x)}(2,1)*xt(x,1,1,t) + S{meios(x)}(2,2)*xt(x,1,2,t) + ZI(meios(x))*IF(meios(x));  
        
         % cálculo da temperatura no nó
         T(x,t) = 2*(xt(x,1,1,t) + xt(x,1,2,t))/(Zx(meios(x))*G(meios(x)) + 2) + ...
             Zx(meios(x))/(Zx(meios(x))*G(meios(x)) + 2)*IF(meios(x));
    end
    
     % tensões incidentes na
     % porta 1
     xt(end,1,1,t) = xt(end,2,1,t-1)*ro(end,1) + xt(end-1,2,2,t-1)*tau(end,1);
    
     % tensões refletidas
     xt(end,2,1,t) = xt(end,1,1,t)*(2 - Zx(end)*G(end))/(Zx(end)*G(end) + 2) ...
         + 2*ZI(end)*(IF(end)/2 - qL);
     
     % cálculo da temperatura
    T(end,t) = 4*xt(end,1,1,t)/(Zx(end)*G(end) + 2) + ...
        2*ZI(end)*(IF(end)/2 - qL);

end

% for n = 1:ktotal
%     plot([0 xpos],[TC; T(:,n)])
%     pause(0.0001)
% end


%Encontrar a temperatura de superfície com o modelo analítico utilizando
%toda a dedução
%1 músculo
%2 Gordura
L1=Ltotal(1);L2=Ltotal(2);%comprimento do tecido
k1=k(1);k2=k(3);%condutividade térmica do tecido
wb1=wb(1);wb2=wb(3);%perfusão sanguínea
qmet1=qmet(1);qmet2=qmet(3);%produção de calor metabólico

TC1=TC;%temperatura na fronteira interna
qL2=qL;%fluxo de calor para o meio externo

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
%% teste
x1 = [deltax(1):deltax(1):(deltax(1)*NNOS(1))];
x1 = [x1 (x1(end) + deltax(1)/2 + deltax(2)/2):deltax(2):...
    (x1(end) + deltax(1)/2 + deltax(2)/2 + deltax(2)*(NNOS(2) - 1))];

x2 = [(deltax(3)/2):deltax(3):(deltax(3)/2 + deltax(3)*(NNOS(3) - 1))];

x2 = [x2 (x2(end) + deltax(3)/2 + deltax(4)/2):deltax(4):...
    (x2(end) + deltax(3)/2 + deltax(4)/2 + deltax(4)*(NNOS(4) - 1))];


T1 = TL1*sinh(l1*x1)/s1 ...
    + TC1*sinh(l1*(L1 - x1))/s1 ...
    + n1*(1 - (sinh(l1*x1) + sinh(l1*(L1 - x1)))/s1);

T2 = TL2*sinh(l2*x2)/s2 ...
    + TC2*sinh(l2*(L2 - x2))/s2 ...
    + n2*(1 - (sinh(l2*x2) + sinh(l2*(L2 - x2)))/s2);
Ta = [T1 T2];
xa = [x1 (x2 + x1(end) + deltax(2)/2)];
    
% plot(xa,Ta);

erro = (T(:,end)' - Ta)./Ta*100;

plot(xpos,erro)
% plot(xa,Ta);
plot([0 xpos],[TC T(:,end)'],'*',[0 xa],[TC Ta])