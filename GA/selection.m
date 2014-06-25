function [S,g] = selection(P,m,n,outs,system_order,system_input_n)

%evaluate the fitness function for the N variables

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
    f(k,:) = desired_function(Q,R);
end

%calculate the probability of fitness for each f(k)
num = 0;
den = 0;
for k=1:m
    for j = 1:outs
        prob_indvidual(k,j) = abs(f(k,j))/sum(f(:,j));
    end
    probF(k,1) = sum(prob_indvidual(k,:));
    probF(k,2) = k;
end


%buble sort to ordering the probabilities (ascendent)
l = m;
while l > 1;
    for k=1:(m-1)
        if probF(k,1) > probF(k+1,1)
            tmp = probF(k,:);
            probF(k,:) = probF(k+1,:);
            probF(k+1,:) = tmp;
        end
    end
    l = l-1;
end
for k=1:m/2
    a = probF(k,2);
    S(k,:) = P(a,:);
    g(k,:) = f(k,:);
end


