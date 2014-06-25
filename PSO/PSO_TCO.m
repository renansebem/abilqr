clc, clear all, close all%, format long E
tic
%field x vs y
bottom = [0, 0, 0, 0];
top = [100, 100, 100, 100];
goal = [0.005 53]; % overshoot, settling time
eps = 0.055*goal;

system_order = 3;
system_input_n = 1;

inertia = 0.4;
initial_inertia = 90;
self_confidence = 0.55;
toward_best = 0.05;
speed = 1;


%swarm
particle_number = 400;
particle_parameters = length(bottom);
for k=1:particle_parameters
    particles(k,:) = bottom(k)*ones(1,particle_number)+(top(k)-bottom(k))*rand(1, particle_number); 
    
end

for k = 1:particle_number
    particles_dpos(:,k) = initial_inertia*(bottom' + (top-bottom)'.*rand(particle_parameters,1).*ones(particle_parameters,1));
end

best_particle = 1;
best_z = 100*goal;
z_best = 100*goal;
particles_best = 10*ones(particle_parameters, particle_number);



m=0;
error = [1 500];

while max(error > eps)
    %m = m+1
    %error_acum = 0;
    

        
        
    for i=1:1:particle_number
        
        



        for index_o = 1:system_order
            Q(index_o,index_o) = particles(index_o,i);
        end
        [~, p] = chol(Q);
        if p
            Q = eye(system_order);
            particles(1:system_order,i) = ones(system_order,1);
        end
        for index_o = 1:system_input_n
            R(index_o,index_o) = particles(system_order + index_o,i);
            if R <= 0
                R = eye(system_input_n);
                particles(system_order + index_o,i) = eye(system_input_n);
            end
        end
        



        %best_z = desired_function(best_Q, best_R);
        z = desired_function(Q,R);
        if sum((abs(z-goal)./goal).^2) < sum((abs(z_best-goal)./goal).^2) || min(abs(z - goal) < eps)
            Q_best = Q;
            R_best = R;
            particles_best(:,i) = particles(:,i);
            z_best = z;%desired_function(Q_best, R_best);
            if sum((abs(z-goal)./goal).^2) < sum((abs(best_z-goal)./goal).^2) || min(abs(z - goal) < eps)
                best_Q = Q;
                best_R = R;
                best_particle = i;
                best_z = z;%desired_function(Q, R);
            end
        end

        
        if particles(:,i) ~= particles(:,best_particle)
            for j=1:1:particle_parameters

                particles_dpos(j,i) = speed*(toward_best*rand*(particles(j,best_particle) - particles(j,i)) + self_confidence*rand*(particles_best(j,i) - particles(j,i)) + inertia*particles_dpos(j,i));
                particles(j,i) = particles_dpos(j,i) + particles(j,i);
                if particles(j,i) > top(j) 
                    particles(j,i) = top(j)-abs(particles_dpos(j,i)/2);
                end
                if particles(j,i) < bottom(j)
                    particles(j,i) = bottom(j)+abs(particles_dpos(j,i)/2);
                end
            end
        end
    end
    %error = error_acum/i;
    error = abs((best_z) - goal)
    
    %plot(particles(1,:), particles(2,:), 'o',particles(1,best_particle),particles(2,best_particle),'rx')
    
    %axis([-5,5,-5,5])
    %axis equal
    %movie1(m) = getframe;
end

%movie(movie1,1,2)
best_over_settling = best_z
particles(:,best_particle)
toc
%z = x^2 + y^2
