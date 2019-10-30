% --------------------------------------------------- %
% --------------------- INÍCIO ---------------------- %

clear %comando para limpar
clc %comando para limpar

menu1 = menu('MÉTODO DOS ELEMENTOS FINITOS','VIGA BI-ENGASTADA COM SEÇÃO VARIÁVEL');


% --------------------------------------------- %
% ------------ VALORES PARA ALTERAR ----------- %
% CARACTERÍSTICAS DO MATERIAL %                 
dados = inputdlg({'Módulo de Young do material (GPa):','Massa específica do material (kg/m³):','Base da seção transversal (cm):','Altura da viga na seção constante (cm):','Altura da viga no engaste (cm):','Comprimento da viga (m):','Comprimento da seção variável (m):','Quantidade de elementos:','Quantidade de formas modais:','Carregamento uniformemente distribuído (N/m):'});
E = (str2num(dados{1})*10^9); %módulo de Young: Pa                  
rho = (str2num(dados{2})); %massa por volume: kg/m³

% DIMENSÕES DA VIGA %
b = (str2num(dados{3})/100); %base da viga: m
h = (str2num(dados{4})/100); %altura da seção constante: m
H = (str2num(dados{5})/100); %altura da seção variável: m
L = (str2num(dados{6})); %comprimento total da viga: m
L1 = (str2num(dados{7})); %comprimento da mísula: m

% DIMENSÕES DO ELEMENTO %
n = (str2num(dados{8})); %quantidade elementos
modal = (str2num(dados{9})); %quantidade de formas modais

% CARREGAMENTO PERMANENTE DISTRIBUÍDO %
q = (str2num(dados{10})); %carregamento em N/m
% --------------------------------------------- %
% --------------------------------------------- %


% --------------------------------------------- %
% ---------------- NÃO ALTERAR ---------------- %

% CONDIÇÕES MÍNIMAS E MÁXIMAS %
if E<=0
    E=10
end
if rho<=0
    rho=100
end
if n<=0
    n=3;
end
if L<=0
    L=1;
end
if H<=0
    H=0.2;
end
if h<=0
    h=0.1;
end
if b<=0
    b=0.1;
end
if mod(n,3) ~= 0
    n = n + 1;
    if mod(n,3) ~=0
        n = n - 2;
    end
end
if L1>(L/2)
    L1 = L/2
end

% ELEMENTOS E GRAUS DE LIBERDADE %
delta = L/n; %tamanho do elemento
a = L/2/n; %tamanho do elemento em coordenada xi
syms xi
vet = [0:delta:n]; %vetor que caminha em cada elemento de 0 a L
Nn = n+1; %número de nós
Ngl= Nn*2; %número de graus de liberdade

% MATRIZES M, K E VETOR F %
M = zeros(Ngl,Ngl);    
K = zeros(Ngl,Ngl);    
F = zeros(Ngl,1);

% LAÇO DE REPETIÇÃO PARA ALTURA DE CADA ELEMENTO %
i = 0;
for x=(vet)
    i = i+1;
    x1 = x + delta/2;
    if x < L1
        h1(1,i) = ((h-H)/L1)*x1 + H;
    elseif x < (L-L1)
        h1(1,i) = h;
    elseif x < L
        h1(1,i) = (H-h)*x1/L1 + H - (H-h)*L/L1;
    end
end

% LAÇO DE REPETIÇÃO PARA A MATRIZ GLOBAL %
for i = (1:n)
    A = b*h1(1,i);
    I = b*(h1(1,i))^3/12;
    py = - q - (rho*A*10);
    N = [1/4*(2-3*xi+xi^3) a/4*(1-xi-xi^2+xi^3) 1/4*(2+3*xi-xi^3) a/4*(-1-xi+xi^2+xi^3)];
    Me = eval((rho*A*a*int(N'*N,xi,-1,1)));
    Ke = eval((E*I/a^3)*int(diff(N,2)'*diff(N,2),xi,-1,1));
    fe = eval((a*int(py*N',xi,-1,1)));
    
    Mee=zeros(Ngl,Ngl); 
    Kee=zeros(Ngl,Ngl);
    Fee=zeros(Ngl,1);
    
    Mee((2*i-1):(2*i+2),(2*i-1):(2*i+2))= Me;
    Kee((2*i-1):(2*i+2),(2*i-1):(2*i+2))= Ke;
    Fee((2*i-1):(2*i+2),1)= fe;
    
    M = M+Mee;
    K = K+Kee;
    F = F+Fee;
end

% CONDIÇÃO DE CONTORNO %
bc = [1 2   Ngl-1 Ngl];

K(bc,:) = [];    
K(:,bc) = [];    
K;

M(bc,:) = [];
M(:,bc) = [];
M;

F(bc,:) = [];
F;

% VETOR DESLOCAMENTO %
u = K\F;

DOF = (1:Ngl)';
DOF(bc,:) = []; %retira-se os gdl restringidos
DOF; %vetor com as direções não impedidas

uyz = zeros(Ngl,1); %vetor nulo com dimensão: Ngl x 1
    
for i=1:length(DOF)    
    uyz(DOF(i),1) = u(i,1); %relaciona a posição livre (gdl) com o deslocamento calculado
end

for i=1:Ngl/2
    uy(i,:) = uyz(2*i-1,:); %somente translação
end
results.uy = uy;

% AUTOVETORES E AUTOVALORES %
[Avet,Aval] = eig(K,M);        % Avet: eigenvectors
                               % Aval: eigenvalues

wn   = sqrt(diag(Aval)); %frequencias do sistema (rad/s)
fn   = wn/(2*pi); %frequencias do sistema (Hz)
T    = fn.^-1; %período de vibração natural (s)

results.fn = fn; %frequencias do sistema (Hz)
results.avet = Avet(1:2:Ngl-4,:);

% RESULTADOS %
format long

results.x = 0:L/n:L; %vetor com os nós para plotar em x
figure(1)
plot(results.x,results.uy,'--r');
title('Deslocamento da Viga','Interpreter','latex','fontsize',14)
xlabel('Comprimento da Viga [m]','Interpreter','latex','fontsize',14)
ylabel('Deslocamento [m]','Interpreter','latex','fontsize',14)
%legend(['n = ',num2str(n)])
set (gcf,'Color','w')
grid on
Deslocamento_lim_inferior = min(results.uy)*1000

figure(2)
bar(results.fn,'stacked')
title('Frequencias de Ressonancia [Hz]','Interpreter','latex','fontsize',14)
xlabel('Elemento','Interpreter','latex','fontsize',14)
ylabel('Frequencia [Hz]','Interpreter','latex','fontsize',14)
set (gcf,'Color','w')
grid on
%legend(['Discretização: n = ',num2str(n)])

% CONDIÇÃO MÁXIMA PARA FORMAS MODAIS %
if modal <=(Ngl-6)
    MODAL = modal+2;
else
    MODAL = Ngl-6;
end
%

for i = (3:MODAL)
    figure(i) %separar 1 figura por modo de vibraçao
    plot(results.avet(:,i-2))
    title('Modos de Vibracao','Interpreter','latex','fontsize',14)
    xlabel('Elemento','Interpreter','latex','fontsize',14)
    ylabel('Vibracao [m]','Interpreter','latex','fontsize',14)
    set (gcf,'Color','w')
    grid on
    %legend(['Discretização: n = ',num2str(n)])
end

% ----------------------- FIM ------------------------ %