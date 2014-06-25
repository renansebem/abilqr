function [out] = desired_function(Q, R)
format long e

% A =   [0.8949 -0.0099; 0.0947 0.9995  ];
% 
% B = [0.0947 ; 0.0048];
% 
% C = [0 0.0689 ];

A = [0.8925 -0.0123; 0.0945 0.9994];

B = [0.0945; 0.0048];

C = [0 0.0883];


D = 0;
[slA, scA] = size(A);
[slB, scB] = size(B);
[slC, scC] = size(C);

Aqsi = [ A zeros(slA,1); -C 1 ];
Bqsi = [B; 0];
Cqsi = [C 0];
Dqsi = [0];

%Q = [1 0; 0 0.1];
%R = 25;

%K = dlqr(A,B,Q,R);

K = dlqr(Aqsi,Bqsi,Q,R);

Kx = K(1:slA);

Kw = K(slA+1);
%[K, S, E] = dlqr(A,B,Q,R);

N = 150; % ms

QSI = zeros(slA+1,N);
u =  zeros(1,N);
%r =  zeros(1,N);
r =  0.05*32767*ones(N);
x =  zeros(length(A),N);
w =  zeros(1,N);
y =  zeros(1,N);
flag = 1;
settlingtime = 10000;

for k = 1:N
    
    u(k) = -K * QSI(:,k);%-((Kx * x(:,k)) + (Kw * w(:,k))) ;
    if u(k) > 0.2*32767 
        u(k) = 0.2*32767;
    elseif u(k) < -0.8*32767
        u(k) = -0.8*32767;
    end
    y(:,k) = Cqsi*QSI(:,k);
    
    w(:,k+1) = w(:,k) + r(k) - y(:,k); 
    x(:,k+1) = A*x(:,k) + B*u(k);
    QSI(:,k+1) = [x(:,k+1); w(:,k+1)];
    
    %QSI(:,k + 1) = Aqsi*QSI(:,k) + Bqsi*u(k) + [ zeros(slA,1); 1]*r(k) - [0;1]*y(:,k);
    n = 50;
    if k>n && abs(sum(y(k-n:k))/n - 0.05*32767) < 0.03*0.05*32767 && flag && y(k-n)> (0.05*32767 -0.03*0.05*32676)
        settlingtime = k-n;
        flag = 0;
    end

end

overshoot = max(y)/(y(end)) - 1;

out = [overshoot, settlingtime];


% subplot(311)
% plot(y)
% title('Sa?da')
% subplot(312)
% plot(u)
% title('A??o de controle')
% subplot(313)
% plot(w)
% title('Estado estendido')