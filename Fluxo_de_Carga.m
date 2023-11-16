%%  Leitura e criacao dos vetores da rede
clear, clc;
% Leitura dos dados
[B Tipo V O bk Pg Qg Pc Qc] = textread('BARRAS_com_Centrais.dat', '%d %d %f %f %f %f %f %f %f', 'headerlines',1);
[Bi Bf Rkm Xkm bkmsh a fi] = textread('RAMOS_com_Centrais.dat', '%d %d %f %f %f %f %f', 'headerlines',1);

% Vetores e matrizes auxiliares
nb = size(B,1);
nr = size(Bi, 1);
bk = zeros(4, 1);
ykm = zeros(nr, 1);
Ykm = zeros(6, 6);
tol = 10^-4;

% Formacao dos elementos serie da rede em estudo
zkm = Rkm + j*Xkm;
for n = 1:nr
    ykm(n,1) = inv(zkm(n,1));
end

% Matriz admitancia
for n = 1:nr
    k = Bi(n, 1);
    m = Bf(n, 1);
    
    Ykm(k, m) = -a(n, 1)*exp(-j*fi(n, 1))*ykm(n, 1);
    Ykm(m, k) = -a(n, 1)*exp(-j*fi(n, 1))*ykm(n, 1);
    Ykm(k, k) = Ykm(k, k) + a(n, 1)^2*ykm(n, 1) + j*bkmsh(n, 1);
    Ykm(m, m) = Ykm(m, m) + a(n, 1)^2*ykm(n, 1) + j*bkmsh(n, 1);
end

% 
Ykm(1, 1) = Ykm(1, 1) + j*bk(1, 1);
Ykm(2, 2) = Ykm(2, 2) + j*bk(2, 1);
Ykm(3, 3) = Ykm(3, 3) + j*bk(3, 1);
Ykm(4, 4) = Ykm(4, 4) + j*bk(4, 1);
 %Ykm(5, 5) = Ykm(5, 5) + j*bk(5, 1);
 %Ykm(6, 6) = Ykm(6, 6) + j*bk(6, 1);


% Matrizes separadas (Condutancia e Susceptancia)
Gkm = real(Ykm);
Bkm = imag(Ykm);

%%  Aplicacao do metodo de Newton
Pk = zeros(nb, 1);
Qk = zeros(nb, 1);
DPk = zeros(nb, 1);
DQk = zeros(nb, 1);
H = zeros(nb, nb);
M = zeros(nb, nb);
N = zeros(nb, nb);
L = zeros(nb, nb);

Pesp = Pg-Pc;
Qesp = Qg-Qc;

% ov = [O; V];
for interacoes = 1:100
    % Reset de Pk e Qk
    Pk = zeros(nb, 1);
    Qk = zeros(nb, 1);

    % Calculo dos fluxos de potencia
    for k = 1:nb
        for m = 1:nb
            if Tipo(k, 1) ~= 1 % Check de barra PV e PQ para Pk
                Okm = O(k, 1) - O(m, 1);
                Pk(k, 1) = Pk(k, 1) + V(k, 1)*V(m, 1)*(Gkm(k, m)*cos(Okm) + Bkm(k, m)*sin(Okm));
                DPk(k, 1) = Pesp(k, 1) - Pk(k, 1);
            end
            if Tipo(k, 1) == 3 % Check de barra PQ para Qk
                Okm = O(k, 1) - O(m, 1);
                Qk(k, 1) = Qk(k, 1) + V(k, 1)*V(m, 1)*(Gkm(k, m)*sin(Okm) - Bkm(k, m)*cos(Okm));
                DQk(k, 1) = Qesp(k, 1) - Qk(k, 1);
            end
        end
    end
    
    % Teste de convergencia
    maxerroDPk = max(abs(DPk));
    maxerroDQk = max(abs(DQk));
    if maxerroDPk < tol && maxerroDQk < tol
        break
    end

    % Jacobiana
    for n = 1:nr
        k = Bi(n, 1);
        m = Bf(n, 1);
        % Calculo das Jacobianas das extremidades
        Okm = O(k, 1) - O(m, 1);
        H(k, m) = V(k, 1)*V(m, 1)*(Gkm(k, m)*sin(Okm) - Bkm(k, m)*cos(Okm));
        H(m, k) = V(k, 1)*V(m, 1)*(-Gkm(m, k)*sin(Okm) - Bkm(m, k)*cos(Okm));

        N(k, m) = V(k, 1)*(Gkm(k, m)*cos(Okm) + Bkm(k, m)*sin(Okm));
        N(m, k) = V(m, 1)*(Gkm(m, k)*cos(Okm) - Bkm(m, k)*sin(Okm));

        M(k, m) = -V(k, 1)*V(m, 1)*(Gkm(k, m)*cos(Okm) + Bkm(k, m)*sin(Okm));
        M(m, k) = -V(k, 1)*V(m, 1)*(Gkm(m, k)*cos(Okm) - Bkm(m, k)*sin(Okm));
        
        L(k, m) = V(k, 1)*(Gkm(k, m)*sin(Okm) - Bkm(k, m)*cos(Okm));
        L(m, k) = -V(m, 1)*(Gkm(m, k)*sin(Okm) + Bkm(m, k)*cos(Okm));
    end
    for n=1:nb
        % Calculo das Jacobianas na diagonal principal
        H(n, n) = -(V(n, 1)^2)*Bkm(n, n) - Qk(n, 1);
        N(n, n) = (1/V(n, 1))*(Pk(n, 1) + Gkm(n, n)*(V(n, 1)^2));
        M(n, n) = -(V(n, 1)^2)*Gkm(n, n) + Pk(n, 1);
        L(n, n) = (1/V(n, 1))*(Qk(n, 1) - Bkm(n, n)*(V(n, 1)^2));
        % Check para tipo de barra
        if Tipo(n, 1) == 1 
            H(n, n) = 10^10;
            L(n, n) = 10^10;
        end
        if Tipo(n, 1) == 2
            L(n, n) = 10^10;
        end
    end

    % Monta a jacobiana
    J = [H N; M L];
    dov = inv(J)*[DPk; DQk];    % Calcula o deltaV e deltaTheta
    for n = 1:nb                % Atualiza os valores de Tensao e defazagem tirados da tabela
        O(n, 1) = O(n, 1) + dov(n, 1);
        V(n, 1) = V(n, 1) + dov(n+nb, 1);
    end
end
%% Fluxos de potencia
Pkm = zeros(nb, nb);
Qkm = zeros(nb, nb);
Pmk = zeros(nb, nb);
Potencia_ativa_saida = zeros(nb, 1);
Potencia_ativa_entrada = zeros(nb, 1);
Potencia_reativa_saida = zeros(nb, 1);
Potencia_reativa_entrada = zeros(nb, 1);
Qmk = zeros(nb, nb);
pkm = zeros(nb, 1);
pmk = zeros(nb, 1);
qkm = zeros(nb, 1);
qmk = zeros(nb, 1);

% forma vetorial
for n = 1:nb
    for l = 1:nr
        k = Bi(l, 1);
        m = Bf(l, 1);
        Okm = O(k, 1) - O(m, 1);
        Omk = O(m, 1) - O(k, 1);

        if Bi(l, 1) == n

            pkm(n, 1) = pkm(n, 1) + a(l, 1)*V(k, 1)*V(m, 1) * (a(l, 1)*Gkm(k, m)*V(k, 1)/V(m, 1) - Gkm(k, m)*cos(Okm + fi(l, 1)) - Bkm(k, m)*sin(Okm + fi(l, 1)));
            qkm(n, 1) = qkm(n, 1) + a(l, 1)*V(k, 1)*V(m, 1) * (-a(l, 1)*(Bkm(k, m) + bkmsh(l, 1))*V(k, 1)/V(m, 1) - Gkm(k, m)*sin(Okm + fi(l, 1)) + Bkm(k, m)*cos(Okm + fi(l, 1)));
            pmk(n, 1) = pmk(n, 1) + V(m, 1)*V(k, 1) * (Gkm(k, m)*V(m, 1)/V(k, 1) - a(l, 1)*Gkm(k, m)*cos(Okm + fi(l, 1)) + a(l, 1)*Bkm(k, m)*sin(Okm + fi(l, 1)));
            qmk(n, 1) = qmk(n, 1) + V(m, 1)*V(k, 1) * (-(Bkm(k, m) + bkmsh(l, 1))*V(m, 1)/V(k, 1) + a(l, 1)*Gkm(k, m)*sin(Okm + fi(l, 1)) + a(l, 1)*Bkm(k, m)*cos(Okm + fi(l, 1)));
        elseif Bf(l, 1) == n

            pkm(n, 1) = pkm(n, 1) + a(l, 1)*V(m, 1)*V(k, 1) * (a(l, 1)*Gkm(m, k)*V(k, 1)/V(m, 1) - Gkm(m, k)*cos(Omk + fi(l, 1)) - Bkm(m, k)*sin(Omk + fi(l, 1)));
            qkm(n, 1) = qkm(n, 1) + a(l, 1)*V(m, 1)*V(k, 1) * (-a(l, 1)*(Bkm(m, k) + bkmsh(l, 1))*V(m, 1)/V(k, 1) - Gkm(m, k)*sin(Omk + fi(l, 1)) + Bkm(m, k)*cos(Omk + fi(l, 1)));
            pmk(n, 1) = pmk(n, 1) + V(k, 1)*V(m, 1) * (Gkm(m, k)*V(k, 1)/V(m, 1) - a(l, 1)*Gkm(m, k)*cos(Omk + fi(l, 1)) + a(l, 1)*Bkm(m, k)*sin(Omk + fi(l, 1)));
            qmk(n, 1) = qmk(n, 1) + V(k, 1)*V(m, 1) * (-(Bkm(m, k) + bkmsh(l, 1))*V(k, 1)/V(m, 1) + a(l, 1)*Gkm(m, k)*sin(Omk + fi(l, 1)) + a(l, 1)*Bkm(m, k)*cos(Omk + fi(l, 1)));
        end
    end   
end

% forma matrixial
for n = 1:nr
    k = Bi(n, 1);
    m = Bf(n, 1);
    Okm = O(k, 1) - O(m, 1);
    Omk = O(m, 1) - O(k, 1);
    % Calculo de Pkm e Pmk
    %Pkm(k, m) = Gkm(k, m)*V(k, 1)^2 - Gkm(k, m)*V(k, 1)*V(m, 1)*cos(Okm) - Bkm(k, m)*V(k, 1)*V(m, 1)*sin(Okm);
    %Pkm(m, k) = Gkm(m, k)*V(m, 1)^2 - Gkm(m, k)*V(m, 1)*V(k, 1)*cos(Okm) - Bkm(m, k)*V(m, 1)*V(k, 1)*sin(Okm);
    Pkm(k, m) = a(n, 1)*V(k, 1)*V(m, 1) * (a(n, 1)*Gkm(k, m)*V(k, 1)/V(m, 1) - Gkm(k, m)*cos(Okm + fi(n, 1)) - Bkm(k, m)*sin(Okm + fi(n, 1)));
    Pkm(m, k) = a(n, 1)*V(k, 1)*V(m, 1) * (a(n, 1)*Gkm(m, k)*V(k, 1)/V(m, 1) - Gkm(k, m)*cos(Okm + fi(n, 1)) - Bkm(k, m)*sin(Okm + fi(n, 1)));
    %Pkm(k, k) = Pkm(k, k) + a(n, 1)*V(k, 1)*V(m, 1) * (a(n, 1)*Gkm(k, k)*V(k, 1)/V(m, 1) - Gkm(k, m)*cos(Okm + fi(n, 1)) - Bkm(k, m)*sin(Okm + fi(n, 1)));
    %Pkm(m, m) = Pkm(m, m) + a(n, 1)*V(k, 1)*V(m, 1) * (a(n, 1)*Gkm(m, m)*V(k, 1)/V(m, 1) - Gkm(k, m)*cos(Okm + fi(n, 1)) - Bkm(k, m)*sin(Okm + fi(n, 1)));

    Pmk(k, m) = V(m, 1)*V(k, 1) * (Gkm(k, m)*V(m, 1)/V(k, 1) - a(n, 1)*Gkm(k, m)*cos(Okm + fi(n, 1)) + a(n, 1)*Bkm(k, m)*sin(Okm + fi(n, 1)));
    Pmk(m, k) = V(m, 1)*V(k, 1) * (Gkm(m, k)*V(m, 1)/V(k, 1) - a(n, 1)*Gkm(m, k)*cos(Okm + fi(n, 1)) + a(n, 1)*Bkm(m, k)*sin(Okm + fi(n, 1)));
    %Pmk(k, k) = Pmk(k, k) + V(m, 1)*V(k, 1) * (Gkm(k, k)*V(m, 1)/V(k, 1) - a(n, 1)*Gkm(k, k)*cos(Okm + fi(n, 1)) + a(n, 1)*Bkm(k, k)*sin(Okm + fi(n, 1)));
    %Pmk(m, m) = Pmk(m, m) + V(m, 1)*V(k, 1) * (Gkm(m, m)*V(m, 1)/V(k, 1) - a(n, 1)*Gkm(m, m)*cos(Okm + fi(n, 1)) + a(n, 1)*Bkm(m, m)*sin(Okm + fi(n, 1)));

    % Calculo de Qkm e Qmk
    Qkm(k, m) = a(n, 1)*V(k, 1)*V(m, 1) * (-a(n, 1)*(Bkm(k, m) + bkmsh(n, 1))*V(k, 1)/V(m, 1) - Gkm(k, m)*sin(Okm + fi(n, 1)) + Bkm(k, m)*cos(Okm + fi(n, 1)));
    Qkm(m, k) = a(n, 1)*V(k, 1)*V(m, 1) * (-a(n, 1)*(Bkm(m, k) + bkmsh(n, 1))*V(k, 1)/V(m, 1) - Gkm(m, k)*sin(Okm + fi(n, 1)) + Bkm(m, k)*cos(Okm + fi(n, 1)));
    %Qkm(k, k) = Qkm(k, k) + a(n, 1)*V(k, 1)*V(m, 1) * (-a(n, 1)*(Bkm(k, k) + bkmsh(n, 1))*V(k, 1)/V(m, 1) - Gkm(k, k)*sin(Okm + fi(n, 1)) + Bkm(k, k)*cos(Okm + fi(n, 1)));
    %Qkm(m, m) = Qkm(m, m) + a(n, 1)*V(k, 1)*V(m, 1) * (-a(n, 1)*(Bkm(m, m) + bkmsh(n, 1))*V(k, 1)/V(m, 1) - Gkm(m, m)*sin(Okm + fi(n, 1)) + Bkm(m, m)*cos(Okm + fi(n, 1)));
    
    Qmk(k, m) = V(m, 1)*V(k, 1) * (-(Bkm(k, m) + bkmsh(n, 1))*V(m, 1)/V(k, 1) + a(n, 1)*Gkm(k, m)*sin(Okm + fi(n, 1)) + a(n, 1)*Bkm(k, m)*cos(Okm + fi(n, 1)));
    Qmk(m, k) = V(m, 1)*V(k, 1) * (-(Bkm(m, k) + bkmsh(n, 1))*V(m, 1)/V(k, 1) + a(n, 1)*Gkm(m, k)*sin(Okm + fi(n, 1)) + a(n, 1)*Bkm(m, k)*cos(Okm + fi(n, 1)));
    %Qmk(k, m) = Qmk(k, k) + V(m, 1)*V(k, 1) * (-(Bkm(k, k) + bkmsh(n, 1))*V(m, 1)/V(k, 1) + a(n, 1)*Gkm(k, k)*sin(Okm + fi(n, 1)) + a(n, 1)*Bkm(k, k)*cos(Okm + fi(n, 1)));
    %Qmk(k, m) = Qmk(k, k) + V(m, 1)*V(k, 1) * (-(Bkm(m, m) + bkmsh(n, 1))*V(m, 1)/V(k, 1) + a(n, 1)*Gkm(m, m)*sin(Okm + fi(n, 1)) + a(n, 1)*Bkm(m, m)*cos(Okm + fi(n, 1)));
end

% Matriz de perdas
Perdas_Ativas = Pkm - Pmk;
Perdas_Reativas = Qkm - Qmk;
for i=1:nb
    for j=1:nb
        Potencia_ativa_saida(i, 1) = Potencia_ativa_saida(i, 1) + Pkm(i, j);
        Potencia_ativa_entrada(i, 1) = Potencia_ativa_entrada(i, 1) + Pmk(i, j);
        Potencia_reativa_saida(i, 1) = Potencia_reativa_saida(i, 1) + Qkm(i, j);
        Potencia_reativa_entrada(i, 1) = Potencia_reativa_entrada(i, 1) + Qmk(i, j);
    end
end
%% Relatorio
fidt = fopen('RELATORIO.txt', 'wt');
fprintf(fidt, 'RELATORIO DO FLUXO DE POTENCIA\n\n');
fprintf(fidt, 'Dados das Barras\n');
fprintf(fidt, '\t Barra \t MagV_a(pu) \t OV_a(*)\n');

for n = 1:nb
    fprintf(fidt, '\t %d \t      %.2f \t          %.2f\n', B(n, 1), V(n, 1), O(n, 1));
end


filename = 'RELATORIO.xls';
TabelaRelatorio = table(pkm, qkm, pmk, qmk, V, O);
writetable(TabelaRelatorio,filename,'Sheet','Relatorio','WriteVariableNames',true);
%% Plotagem de Graficos
% Plotando as tensões e correntes no diagrama de fasores
figure
V2_Incial= 0.9410 * exp(1j * (-0.11*pi/180));
V2_Final= 0.9934 * exp(1j * (0.20*pi/180));
V3_Inicial = 0.9646 * exp(1j * (-0.11*pi/180));
V3_Final = 0.9981 * exp(1j * (0.21*pi/180));
polarplot([0 angle(V2_Incial)], [0 abs(V2_Incial)], 'r', 'LineWidth', 2);
hold on
%polarplot([0 angle(V3_Inicial)], [0 abs(V3_Inicial)], 'g', 'LineWidth', 2);
polarplot([0 angle(V2_Final)], [0 abs(V2_Final)], 'g', 'LineWidth', 2);
%polarplot([0 angle(V3_Final)], [0 abs(V3_Final)], 'g--', 'LineWidth', 2);
rlim([0 1.05]); 
% Adicionando o fator de potência global ao título
title('Diagrama de Fasores');
legend('V2_I','V2_F', 'V3_F', 'Location', 'best');

%% Faltas

%% blank
% clc, clear;
% Y4 = [-3j 1j 0 1j; 1j -2j 1j 0; 0 1j -2j 1j; 1j 0 1j -2j];
% Z4 = inv(Y4);
% E = 1;
% Icc = E / Z4(4, 4);
% Ifalta = [0; 0; 0; -Icc];
% Vf = Z4 * Ifalta;
% Vpf = [1; 1; 1; 1];
% V = Vpf + Vf;
% 
% I14 = (V(1, 1) - V(4, 1)) / Z4(1,1)
% I24 = (V(2, 1) - V(4, 1)) / Z4(2,2)
% I34 = (V(3, 1) - V(4, 1)) / Z4(3,3)
% I44 = (V(4, 1) - V(4 ,1)) / Z4(4,4)
% I23 = (V(2, 1) - V(3, 1)) / 1j;

% %%

% In2 = 1*10^6/220
% Ic2 = 0.75*10^6/220;
% RTC2 = In2/5 
% In1 = 1*10^6/3800
% Ic1 = 0.75*10^6/3800;
% RTC1 = In1/5
% 
% RTC2 = 1000;
% RTC1 = 60;
% 
% Is1 = Ic1/RTC1;
% Is2 = Ic2/RTC2;
% 
% tap = Is2/Is1
% 
% Itap1 = 3.2;
% Itap2 = 2.9;
% Im = (Is1 + Is2)/2
% DeltaI = abs(Is1 - Is2)
% Adj = (DeltaI / Im)*100

% %%

% %%
% clc, clear;
% VBaixa = 69e3;
% 
% VAlta = 220e3;
% 
% S = 150e6;
% 
% InBaixa = S / (sqrt(3)*VBaixa);
% InAltaB = S / (sqrt(3)*(VAlta-20e3));
% InAltaM = S / (sqrt(3)*VAlta);
% InAltaA = S / (sqrt(3)*(VAlta+10e3));
% 
% RTCBaixa = InBaixa/(5/sqrt(3));
% RTCAlta = InAltaB/5;
% 
% RTCBaixa = 500;
% RTCAlta = 100;
% 
% %Corrente nos TCS
% 
% IsBaixa = InBaixa*sqrt(3)/RTCBaixa;
% IsAltaB = InAltaB/RTCAlta;
% IsAltaM = InAltaM/RTCAlta;
% IsAltaA = InAltaA/RTCAlta;
% 
% Em = 0;
% 
% Eb = abs(IsAltaM - IsAltaB)/IsAltaB*100;
% Ea = abs(IsAltaM - IsAltaA)/IsAltaA*100;
% 
% % Adj_Declividade
% Eajuste = 10;
% Evazio = 2;
% Etc = Eb;
% Er = Eajuste + Evazio + Etc
% Adj_Declividade = 30;


% %%
% clc, clear;
% RTG = 220 / 69;
% RTCAlta = 60/5;
% RTCBaixa = 200/5;
% 
% InBaixa = 7.39e6/(69e3);
% IsBaixa = InBaixa / RTCBaixa;
% 
% InAlta = InBaixa / RTG;
% IsAlta = InAlta / RTCAlta;
% 
% deltaI = abs(IsBaixa-IsAlta);
% Im = (IsAlta + IsBaixa)/2
% 
% ad= deltaI/Im*100;


% %%
% clc, clear;
% 
% S = 85e6;
% V = 138e3;
% kf = 1.35;
% Icc = 10e3;
% 
% In = S/(sqrt(3)*V)
% If = Icc/20
% RTC = 100
% Icm = kf*In;
% Is = Icm/RTC


%  %%
%  clear, clc;
%  p = 0.194;
%  km = 50;
%  zrede = 0.8 * km * p;
%  zrele = zrede * ((400/5)/(69e3/(115/sqrt(3))));


% %%
% clear, clc;
% 
% S = 125e6;
% Icc = 11.5e3;
% V = 138e3;
% kf = 1.2;
% t = 0.5;
% RTC = 120;
% 
% In = S/(sqrt(3)*V);
% 
% 
% 
% Itap = In * kf/RTC;
% Iacio = Itap*RTC;
% 
% M = Icc/Iacio;

