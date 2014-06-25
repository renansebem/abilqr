function [C] = crossover(S,m,n,posBin,nbit)
%transform integer vectors into bit vectors
%Sbin = dec2bin(m/2,0,nbit);
% S_fi = fi(S);
% for l = 1:m/2
%     Sbin(l,:) = bin(S_fi(l,n));
% end

C = S;
Sbin = bin(S);
Cbin = bin(C);

%crossover for bit vector
ix = randperm(m/2,m/2);
for l=1:(m/2-1)
    for k=1:posBin
       Cbin(l,k) = Sbin(ix(l),k);
    end
    for k=posBin:nbit
       Cbin(l,k) = Sbin(ix(l+1),k);
    end
end
for k=1:posBin
   Cbin(m/2,k) = Sbin(ix(l+1),k);
end
for k=posBin:nbit
   Cbin(m/2,k) = Sbin(ix(1),k);
end

C.bin = Cbin;

