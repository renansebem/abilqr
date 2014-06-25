format long
clear all
clc 

tic
m = 40;
inds = 4; %number of variables for fitness function
outs = 2;
nbit = 24;
P = fi(10*rand(m,inds), 0, nbit,16);
system_order = 3;
system_input_n = 1;
target = [0.005 48];
eps = 0.05*target;
error_final = [1 1000];
iterations = 0;

% prealocating Q & R
Q = zeros(system_order,system_order);
R = zeros(system_input_n,system_input_n);

while max(error_final > eps) && (iterations < 150)    
    for n = 1:inds
        if n>1
%         S1 = S
%         C1 = C
        end
%         P1 = P
        %select the best values for the fitness function
        [S,g] = selection(P,m,n,outs,system_order,system_input_n);
        %randomize number of crossovers and the cut points
        cuts = randi([2 nbit/4]);
%         S2 = S
%         P2 = P
        
        for k = 1:(cuts-1)
            posCuts = randi([ceil(nbit*(k-1)/cuts)+1 ceil(nbit*(k)/cuts)]);
            C = crossover(S(:,n),m,n,posCuts,nbit);
            S(:,n) = C;
        end
%         S3 = S
%         P3 = P
%         C3 = C
        %randomize the number of individual mutated and the number of bit
        %mutated
        mutsIn = 1;
        mutsEx = randi([1 m/2]);
        for k = 1:mutsEx
            indM = randi([1 m/2]);
            for v = 1:mutsIn 
                posM = randi([1 nbit*n]);
                M = mutation(C,posM,indM,nbit,n);
                C = M;
            end
        end
%         C4 = C
        
        %divide the vector bit into 2 vectors and rebuild the integers
%         for l=1:m/2
%             C(l,n) = double(Cbin(l,n));
%         end
        for l=1:m/2
            P(l,:) = (S(l,:));
        end
%         P5 = P
        for l=(m/2+1):m
            P(l,n) = (C(l-m/2,:));
        end
%         P6 = P

        for k=1:m
            for index_o = 1:system_order
                Q(index_o,index_o) = double(P(k,index_o));
            end
            [~, p] = chol(Q);
            if p
                Q = eye(system_order);
            end
            for index_o = 1:system_input_n
                R(index_o,index_o) = double(P(k,system_order + index_o));
            end
            if ~R
                R = eye(system_input_n);
            end
            g(k,:) = desired_function(Q,R);
            error = (g(k,:)-target);
            if ~max(error > eps)
                P_final = P(k,:);
                error_final = error;
                g_final = g(k,:);
            end
        end
    end
    iterations = iterations+1
    %[f,x] = min(abs(g-ones(length(g),1)*target));
end
% FinalP = P
% Finalg = g
iterations
Xminimum = P_final
Yminimum = g_final
toc