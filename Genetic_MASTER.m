function Evidente = MestradoEstagioDocencia

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

%cada numero maior que 1000 representa um flag de evidente
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
HT = [50000;50000;50000];

fprintf('Insira o Hard-Time de cada latente (a partir de T1, sendo T1 refente ao latente 1, identificado com o número 2 na tabela excel):\n');

for count2 = 1:total_var
    t(count2)=input(sprintf('Insira o Hard-Time referente ao latente %d:',count2));
end
HT = t';

%Comparando ambos os limites superiores para pegar o mais restritivo
ub = min(ub_SR2,HT)

%Limite inferior
lb = ones(1,total_var);

%define active-set algorithm como padrão
%'active-set', 'interior-point', 'sqp', 'trust-region-reflective', or 'sqp-legacy'.
%options = optimset('Algorithm','active-set','TolCon',1e-60,'MaxFunEvals',1e3,'MaxIter',1e3);
options = optimset('Algorithm','active-set','TolCon',1e-60);
%options = optimset('Algorithm','sqp','TolCon',1e-60);
%options = optimoptions('ga','TolCon',1e-60,'TolFun',1e-60);

%You can set separate options for the hybrid function. Use optimset for fminsearch, or optimoptions for fmincon, patternsearch, or fminunc. For example:
%{
ga Hybrid Function
A hybrid function is another minimization function that runs after the genetic algorithm terminates. You can specify a hybrid function in the HybridFcn option. Do not use with integer problems. The choices are

[] — No hybrid function.

'fminsearch' — Uses the MATLAB® function fminsearch to perform unconstrained minimization.

'patternsearch' — Uses a pattern search to perform constrained or unconstrained minimization.

'fminunc' — Uses the Optimization Toolbox™ function fminunc to perform unconstrained minimization.

'fmincon' — Uses the Optimization Toolbox function fmincon to perform
constrained minimization.
%}

%hybridopts = optimoptions('fmincon','OptimalityTolerance',1e-60);
%options = optimoptions('ga','HybridFcn',{'fmincon',hybridopts},'TolCon',1e-60,'TolFun',1e-60);

%hybridopts = optimoptions('patternsearch');
%options = optimoptions('ga','HybridFcn',{'patternsearch',hybridopts},'TolCon',1e-60,'TolFun',1e-60);

%options = optimoptions('ga','CrossoverFraction',0.8,'CrossoverFcn',@crossoverscattered,'MaxGenerations',1e3,'Tolcon',1e-60,'MutationFcn',{@mutationadaptfeasible},'CreationFcn',{@gacreationnonlinearfeasible,'UseParallel',true,'NumStartPts',20})

%options = gaoptimset('TolCon',1e-20,'CreationFcn',@gacreationnonlinearfeasible,'NonlinConAlgorithm','auglag')

%hybridopts = optimoptions('fmincon');
%options = optimoptions('ga','PopulationSize',50,'Generations',300,'StallTimeLimit',3600,'TimeLimit',18000,'CrossoverFraction',0.8,'EliteCount',2,'HybridFcn',{'fmincon',hybridopts},'MutationFcn',@mutationadaptfeasible,'NonlinearConstraintAlgorithm','penalty');

%x = ga(fitnessfcn,nvars,A,b,Aeq,beq,LB,UB,nonlcon,options)
x = ga(@(x)myfun(x,a),3,A,b,Aeq,beq,lb,ub,@(x)mycon(M,aux_M,x,total_var,edges,lambda,e,qtd_e,ComandoVoo),options);

%The solver should be one of these four: fmincon, fminunc, lsqnonlin, lsqcurvefit.
%problem = createOptimProblem('fmincon','objective',@(x)myfun(x,a),'x0',X0,'lb',lb,'ub',ub,'options',options);
%pt = paretosearch('FunctionTolerance',1e-60);
%x = run(pt,problem,24);
%paretosearch(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
%x = patternsearch(@(x)myfun(x,a),X0,A,b,Aeq,beq,lb,ub,@(x)mycon(new_R,x,total_var,edges,lambda),options);
%x = surrogateopt(@(x)myfun(x,a),X0,A,b,Aeq,beq,lb,ub,@(x)mycon(new_R,x,total_var,edges,lambda),options);

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

%Restricao Specific Risk (i): Se um latente ocorreu, o residual <= 1e-5
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

%-------------------------------------------------------------------------------
%%CASOS COM COMANDO DE VOO (APLICAÇÃO RESIDUAL RISK)

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
%f = a(1)/x(1) + a(2)/x(2)+ a(3)/x(3);
f = sum(a(:)./x(:));
end
