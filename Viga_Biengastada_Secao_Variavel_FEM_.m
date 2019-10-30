% --------------------------------------------------- %
% --------------------- IN�CIO ---------------------- %

clear %comando para limpar
clc %comando para limpar

menu1 = menu('M�TODO DOS ELEMENTOS FINITOS','VIGA BI-ENGASTADA COM SE��O VARI�VEL');


% --------------------------------------------- %
% ------------ VALORES PARA ALTERAR ----------- %
% CARACTER�STICAS DO MATERIAL %                 
dados = inputdlg({'M�dulo de Young do material (GPa):','Massa espec�fica do material (kg/m�):','Base da se��o transversal (cm):','Altura da viga na se��o constante (cm):','Altura da viga no engaste (cm):','Comprimento da viga (m):','Comprimento da se��o vari�vel (m):','Quantidade de elementos:','Quantidade de formas modais:','Carregamento uniformemente distribu�do (N/m):'});
E = (str2num(dados{1})*10^9); %m�dulo de Young: Pa                  
rho = (str2num(dados{2})); %massa por volume: kg/m�

% DIMENS�ES DA VIGA %
b = (str2num(dados{3})/100); %base da viga: m
h = (str2num(dados{4})/100); %altura da se��o constante: m
H = (str2num(dados{5})/100); %altura da se��o vari�vel: m
L = (str2num(dados{6})); %comprimento total da viga: m
L1 = (str2num(dados{7})); %comprimento da m�sula: m

% DIMENS�ES DO ELEMENTO %
n = (str2num(dados{8})); %quantidade elementos
modal = (str2num(dados{9})); %quantidade de formas modais

% CARREGAMENTO PERMANENTE DISTRIBU�DO %
q = (str2num(dados{10})); %carregamento em N/m
% --------------------------------------------- %
% --------------------------------------------- %


% --------------------------------------------- %
% ---------------- N�O ALTERAR ---------------- %

% CONDI��ES M�NIMAS E M�XIMAS %
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
Nn = n+1; %n�mero de n�s
Ngl= Nn*2; %n�mero de graus de liberdade

% MATRIZES M, K E VETOR F %
M = zeros(Ngl,Ngl);    
K = zeros(Ngl,Ngl);    
F = zeros(Ngl,1);

% LA�O DE REPETI��O PARA ALTURA DE CADA ELEMENTO %
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

% LA�O DE REPETI��O PARA A MATRIZ GLOBAL %
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

% CONDI��O DE CONTORNO %
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
DOF; %vetor com as dire��es n�o impedidas

uyz = zeros(Ngl,1); %vetor nulo com dimens�o: Ngl x 1
    
for i=1:length(DOF)    
    uyz(DOF(i),1) = u(i,1); %relaciona a posi��o livre (gdl) com o deslocamento calculado
end

for i=1:Ngl/2
    uy(i,:) = uyz(2*i-1,:); %somente transla��o
end
results.uy = uy;

% AUTOVETORES E AUTOVALORES %
[Avet,Aval] = eig(K,M);        % Avet: eigenvectors
                               % Aval: eigenvalues

wn   = sqrt(diag(Aval)); %frequencias do sistema (rad/s)
fn   = wn/(2*pi); %frequencias do sistema (Hz)
T    = fn.^-1; %per�odo de vibra��o natural (s)

results.fn = fn; %frequencias do sistema (Hz)
results.avet = Avet(1:2:Ngl-4,:);

% RESULTADOS %
format long

results.x = 0:L/n:L; %vetor com os n�s para plotar em x
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
%legend(['Discretiza��o: n = ',num2str(n)])

% CONDI��O M�XIMA PARA FORMAS MODAIS %
if modal <=(Ngl-6)
    MODAL = modal+2;
else
    MODAL = Ngl-6;
end
%

for i = (3:MODAL)
    figure(i) %separar 1 figura por modo de vibra�ao
    plot(results.avet(:,i-2))
    title('Modos de Vibracao','Interpreter','latex','fontsize',14)
    xlabel('Elemento','Interpreter','latex','fontsize',14)
    ylabel('Vibracao [m]','Interpreter','latex','fontsize',14)
    set (gcf,'Color','w')
    grid on
    %legend(['Discretiza��o: n = ',num2str(n)])
end

% ----------------------- FIM ------------------------ %