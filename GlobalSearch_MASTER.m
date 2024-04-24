function Evidente = Master

fprintf('É um caso de Comando de Voo? Se sim, responda 1, caso não, digite 0:\n');
ComandoVoo = input(sprintf('1 ou 0? \n'));

tic
%%EVIDENTES
format long g
M = xlsread('PEEmestrado.xlsx')
[row_M, col_M] = size(M);
aux_M = M;

k = 1;
u =1;
for i = 1:row_M
    for j = 1:col_M
        %contar a qtd de Evidentes
        if M(i,j) > 10000
            e_aux(k) = M(i,j);
            k = k+1;
        end
        %contar a qtd de Latentes
        if M(i,j) > 1 && M(i,j)<10000
            L_aux(u) = M(i,j);
            u = u+1;
        end
    end
end

%cada numero maior que 10000 representa um flag de evidente
e = unique(e_aux)

%quantidade de eventos evidentes distintos
qtd_e = length(e)


%Usuario insere a probabilidade de cada evento evidente
fprintf('Insira a probabilidade relativa a cada evidente\n');

for count = 1:qtd_e
    ProbE(count)=input(sprintf('Insira a probabilidade do evidente %d:',count));   
end

%atribuir aos flags de evidentes as suas respectivas probabilidades
aux = e(1:end);
for m = 1:row_M
    for n = 1:col_M
        for z = 1:qtd_e
            if M(m,n)==aux(z)
                M(m,n) = ProbE(z);
            end
        end
    end
end

%cada numero maior que 1 e menor que 10000 representa um flag de latente
L = unique(L_aux)
edges = L;

%total de variaveis = qtd de eventos latentes
total_var = length(L)

%vetor de custo
fprintf('Insira os custos relativos a cada latente (a partir de a1, sendo a1 refente ao latente 1, identificado com número 2 na tabela excel):\n');

for count = 1:total_var
    a(count)=input(sprintf('Insira o custo referente ao latente %d:',count));   
end

%define tempos de inspeção iniciais (não nulo)
X0 = 50*ones(1,total_var);
A = 0*eye(total_var);
b = zeros(total_var,1);
Aeq = zeros(total_var);
beq = zeros(total_var,1);

fprintf('Insira os valores de lambda relativos a cada latente (a partir de lambda1, sendo lambda1 refente ao latente 1, identificado com número 2 na tabela excel):\n');

for count1 = 1:total_var
    taxa(count1)=input(sprintf('Insira o lambda referente ao latente %d:',count1));
end
lambda = taxa';

%Restricao Specific Risk (ii): A prob de ocorrencia de um latente é <= 1e-3
%lambda = 4.7e-9;%podem ser lambdas diferentes
ub_SR2 = (1e-3./(lambda)).*ones(total_var,1); %Caso com Specific Risk

%Restrição de hard-time
fprintf('Insira o Hard-Time de cada latente (a partir de T1, sendo T1 refente ao latente 1, identificado com o número 2 na tabela excel):\n');

for count2 = 1:total_var
    t(count2)=input(sprintf('Insira o Hard-Time referente ao latente %d:',count2));
end
HT = t';

%Comparando ambos os limites superiores para pegar o mais restritivo
ub = min(ub_SR2,HT)
%ub = [50000;50000;50000];
%Limite inferior
lb = ones(1,total_var);

options = optimset('Algorithm','active-set','TolCon',1e-60,'TolX',1e-100,'MaxFunEvals',1e5,'MaxIter',1e5);
x = fmincon(@(x)myfun(x,a),X0,A,b,Aeq,beq,lb,ub,@(x)mycon(M,aux_M,x,total_var,edges,lambda,e,qtd_e,ComandoVoo),options);
opts = optimoptions(@fmincon,'Algorithm','sqp');
problem = createOptimProblem('fmincon','objective',@(x)myfun(x,a),'x0',X0,'lb',lb,'ub',ub,'options',opts);
%gs.NumTrialPoints = 1e15;
gs = GlobalSearch;
x = run(gs,problem);

out = x

toc
format long
f = sum(a(:)./x(:))
end

%Aplicacao do requisito do §25.1309 e SR(i)
function [probTotal,ceq] = mycon(M,aux_M,x,total_var,edges,lambda,e,qtd_e,ComandoVoo)

%§25.1309
[new_row,new_col] = size(M);
i = 1;
aux_matrix_R = M;
aux = edges(1:end);
vector_aux = zeros(1,total_var);
for m = 1:new_row
    for n = 1:new_col
        for z = 1:total_var
            if M(m,n)==aux(z) && vector_aux(z)==0
                M(m,n) = x(i)*lambda(i);
                i = i+1;
                vector_aux(z)=vector_aux(z)+1;
            elseif M(m,n)==aux(z) && vector_aux(z)>0
                M(m,n) = x(z)*lambda(z);
                vector_aux(z)=vector_aux(z)+1;
            end
        end
    end
end

 prod_row = prod(M,2);
 probTotal = sum(prod_row) + M(1,1) - 1.4e-9;
ceq = sum(prod_row) + M(1,1) - 1e-9;

% %Restricao Specific Risk (i): Se um latente ocorreu, o residual <= 1e-5
for p = 1:length(edges)
     for m = 1:new_row
         for n = 1:new_col
             if aux_matrix_R(m,n)==edges(p)    
                 M(m,n) = 1;
             end
         end
     end
     prod_row_SR = prod(M,2);
     SR1(p) = sum(prod_row_SR) + M(1,1) - 1.4e-5;
     probTotal = [probTotal SR1(p)];
 end

% %-------------------------------------------------------------------------------
%CASOS COM COMANDO DE VOO (APLICAÇÃO RESIDUAL RISK) 
if ComandoVoo == 1
    for p = 1:qtd_e
        for m = 1:new_row
            for n = 1:new_col
                if aux_M(m,n)==e(p)
                    aux_M(m,n) = 1;
                end
            end
        end
        prod_row_RR = prod(aux_M,2);
        RR(p) = sum(prod_row_RR) + aux_M(1,1) - 1.4e-3;
        %probTotal = [probTotal RR(p)];
        ceq = RR(p);
    end
end
if ComandoVoo == 0
    %restrição de equacionamento não linear
    ceq=[];
end
end
%-------------------------------------------------------------------------------

%Definição da função objetivo de custo
function f = myfun(x,a)
f = sum(a(:)./x(:));
end
