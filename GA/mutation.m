function M = mutation(C,posM,indM,nbit,n)

%Mutation: permuting one or two bits of the bit Vector.
M = C;
Cbin = bin(C);
l = indM;    
for k = 1:nbit
    ix = randi([1 nbit]);
    if (k == posM) 
        tmp = Cbin(l,k);
        Cbin(l,k) = Cbin(l,ix);
        Cbin(l,ix) = tmp;
    end
end

%Mutation: toggling one bit of the bit vector
ix = randi([1 nbit]);
if Cbin(l,ix) == dec2bin(0)
    Cbin(l,ix) = dec2bin(1);
else
    Cbin(l,ix) = dec2bin(0);
end


M.bin = Cbin;
