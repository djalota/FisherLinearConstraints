%% Tatonnement Implementation Fisher Market

%Consider N people and M goods
N = 20;
M = 20;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 0.1;

x_iter = 1;

p = p_init;
supply_demand = [];
while x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x'  );%maximize( B(i)*log(V(i, :)*x') - p*x');
        subject to
        x >= zeros(1, M);
        cvx_end
        opt_xi(i, :) = x;
    end
    p = p + 0.05*(sum(opt_xi, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand = [supply_demand, abs(my_sum)];
    x_iter = x_iter + 1;
end

semilogy(1:x_iter-1, supply_demand)
hold on


%% ADMM Implementation Fisher Market

%Consider N people and M goods
N = 10;
M = 10;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand2 = [];

beta = 1;

y_in = (1/M)*ones(N, M);

while  x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x' - (beta/2)* square_pos(norm(x-y_in(i, :), 2)) );
        subject to
        x >= zeros(1, M);
        cvx_end
        opt_xi(i, :) = x;
    end
    
    %Find optimal y values
    cvx_begin
    variable y(N, M);
    maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* square_pos(norm(sum(y, 1) - C, 2)) );
    cvx_end
    y_in = y;
    
    
    p = p + beta*(sum(y_in, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand2 = [supply_demand2, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand2)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% ADMM Implementation Fisher Market with supply inequality constraint

%Consider N people and M goods
N = 5;
M = 5;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand2 = [];

beta = 0.3;

y_in = (1/M)*ones(N, M);

while  x_iter <= 50 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x' - (beta/2)* square_pos(norm(x-y_in(i, :), 2)) );
        subject to
        x >= zeros(1, M);
        cvx_end
        opt_xi(i, :) = x;
    end
    
    %Find optimal y values
    cvx_begin
    variable y(N, M);
    maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* sum_square_pos( sum(y, 1) - C, 2) );
    cvx_end
    y_in = y;
    
    for j = 1:M
        if sum(y_in(:, j), 1) - C >= 0
            p(j) = p(j) + beta*(sum(y_in(:, j), 1) - C(j));
            
        else
            p(j) = 0;
        end
    end
    %p = p + beta*(sum(y_in, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand2 = [supply_demand2, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand2)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% ADMM Implementation Fisher Market with x=y inequality constraint

%Consider N people and M goods
N = 10;
M = 10;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand2 = [];

beta = 1;

y_in = (1/M)*ones(N, M);

while  x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*max(0, x-y_in(i, :))' - (beta/2)* sum_square_pos(x-y_in(i, :)) );
        subject to
        x >= zeros(1, M);
        cvx_end
        opt_xi(i, :) = x;
    end
    
    %Find optimal y values
    cvx_begin
    variable y(N, M);
    maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* square_pos(norm(sum(y, 1) - C, 2)) );
    cvx_end
    y_in = y;
    
%     for j = 1:M
%         if sum(y_in(:, j), 1) - C >= 0
%             p(j) = p(j) + beta*(sum(y_in(:, j), 1) - C(j));
%             
%         else
%             p(j) = 0;
%         end
%     end
    p = p + beta*(sum(y_in, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand2 = [supply_demand2, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand2)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% AMA Implementation Fisher Market

%Consider N people and M goods
N = 20;
M = 20;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = 0.1*ones(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand3 = [];

beta = 0.05;

y_in = (1/M)*ones(N, M);
mu = 0.0001;

while x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x' - mu/2*sum(x.^2)  );%- mu/2*sum(x.^2) 
        subject to
        x >= zeros(1, M);
        cvx_end
        opt_xi(i, :) = x;
    end
    
    %Find optimal y values
    cvx_begin
    variable y(N, M);
    maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* square_pos(norm(sum(y, 1) - C, 2)) );
    cvx_end
    y_in = y;
    
    
    p = p + beta*(sum(y_in, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand3 = [supply_demand3, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand3, 'LineWidth', 3)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'ADMM')
hold off

% %% AAMA Implementation Fisher Market
% 
% %Consider N people and M goods
% N = 5;
% M = 5;
% 
% %The budget B is a N vector containing the total budget available for all
% %people
% B = ones(1, N)';
% 
% %The capacity C is a M vector containing the maximum number of people that
% %can use a given good
% C = (N/M)*ones(1, M);
% 
% %The valuation V is a N*M matrix
% %V = rand(N, M);
% 
% %Initialize the Prices in the market
% p_init = rand(1, M);
% 
% %Create matrix to store the optimal xi
% opt_xi = 0.1*ones(N, M);
% 
% %Convergence criteria
% epsilon = 0.001;
% 
% x_iter = 1;
% 
% p = p_init;
% supply_demand4 = [];
% 
% beta = 0.05;
% alpha_old = 1;
% 
% y_in = (1/M)*ones(N, M);
% 
% while x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
%     %Find optimal x values
%     for i = 1:N
%         cvx_begin
%         variable x(1, M);
%         maximize( B(i)*log(V(i, :)*x') - p*x'  );
%         subject to
%         x >= zeros(1, M);
%         cvx_end
%         opt_xi(i, :) = x;
%     end
%     
%     %Find optimal y values
%     cvx_begin
%     variable y(N, M);
%     maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* square_pos(norm(sum(y, 1) - C, 2)) );
%     cvx_end
%     y_in = y;
%     
%     
%     p_new = p + beta*(sum(y_in, 1) - C);
%     
%     alpha_new = (1+sqrt(1+4*alpha_old^2))/2;
%     p = p_new + (alpha_old-1)/(alpha_new)*(p_new - p_init);
%     
%     alpha_old = alpha_new;
%     p_init = p_new;
%     
%     sup_dem = sum(opt_xi, 1) - C;
%     my_sum = 0;
%     for j = 1:M
%         my_sum = my_sum + sup_dem(j)^2;
%     end
%     supply_demand4 = [supply_demand4, abs(my_sum)];
%     x_iter = x_iter + 1;
% end
% semilogy(1:x_iter-1, supply_demand4)
% %ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
% %xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
% %legend('Dual Ascent', 'ADMM')
% hold off

%% Augmented Lagrangian Implementation Fisher Market

%Consider N people and M goods
N = 5;
M = 5;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = 0.1*ones(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand5 = [];

beta = 1;

y_in = (1/M)*ones(N, M);

while x_iter <= 50 %max(abs(sum(opt_xi, 1) - C)) > epsilon & 
    %Find optimal x values
    
    cvx_begin
    variable x(N, M);
    maximize( B'*log(sum(V.*x, 2)) - p*sum(x, 1)' - (beta/2)* square_pos(norm(sum(x, 1) - C, 2)) );
    subject to
    x >= zeros(N, M);
    cvx_end
    x_in = x;
    opt_xi = x_in;
    
    p = p + beta*(sum(x_in, 1) - C);
    
    sup_dem = sum(x_in, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand5 = [supply_demand5, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand5)
ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% Accelerated Augmented Lagrangian Implementation Fisher Market

%Consider N people and M goods
N = 5;
M = 5;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
%V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = 0.1*ones(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand6 = [];

alpha_old = 1;
beta = 1;

y_in = (1/M)*ones(N, M);

while x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon & 
    %Find optimal x values
    
    cvx_begin
    variable x(N, M);
    maximize( B'*log(sum(V.*x, 2)) - p*sum(x, 1)' - (beta/2)* square_pos(norm(sum(x, 1) - C, 2)) );
    subject to
    x >= zeros(N, M);
    cvx_end
    x_in = x;
    opt_xi = x_in;
    
    p_new = p + beta*(sum(x_in, 1) - C);
    
    alpha_new = (1+sqrt(1+4*alpha_old^2))/2;
    p = p_new + (alpha_old-1)/(alpha_new)*(p_new - p_init) + (alpha_old/alpha_new)*(p_new - p);
    
    alpha_old = alpha_new;
    p_init = p_new;
    
    sup_dem = sum(x_in, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand6 = [supply_demand6, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand6)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Augmented Lagrangian', 'Accelerated Augmented Lagrangian')
hold off

%% ADMM Different Norm Implementation Fisher Market - epsilon constraint

%Consider N people and M goods
N = 5;
M = 5;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 1;

x_iter = 1;

p = p_init;
supply_demand2 = [];

beta = 1;

y_in = (1/M)*ones(N, M);

while  x_iter <= 50 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x'  ); %beta*norm(x-y_in(i, :), 1)
        subject to
        x >= zeros(1, M);
        x - y_in(i, :) <= epsilon/x_iter;
        y_in(i, :) - x <= epsilon/x_iter;
        cvx_end
        opt_xi(i, :) = x;
    end
    opt_xi = opt_xi + 0.5*beta*(y_in-opt_xi);
    %Find optimal y values
    cvx_begin
    variable y(N, M);
    maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* square_pos(norm(sum(y, 1) - C, 2)) );
    subject to
    y >= zeros(N, M);
    cvx_end
    y_in = y;
    
    
    p = p + beta*(sum(y_in, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand2 = [supply_demand2, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand2)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% Benchmark Devanur Implementation Fisher Market

%Consider N people and M goods
N = 5;
M = 5;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 0.1;

x_iter = 1;

p = p_init;
supply_demand = [];
while x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x'  );%maximize( B(i)*log(V(i, :)*x') - p*x');
        subject to
        x >= zeros(1, M);
        cvx_end
        opt_xi(i, :) = x;
    end
    p = p.*exp((sum(opt_xi, 1) - C)/20); %+ 0.05*(sum(opt_xi, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand = [supply_demand, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% Benchmark Devanur Implementation Fisher Market

%Consider N people and M goods
N = 5;
M = 5;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 0.1;

x_iter = 1;

p = p_init;
supply_demand = [];
while x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x'  );%maximize( B(i)*log(V(i, :)*x') - p*x');
        subject to
        x >= zeros(1, M);
        cvx_end
        opt_xi(i, :) = x;
    end
    p = p.*exp((sum(opt_xi, 1) - C)/20); %+ 0.05*(sum(opt_xi, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand = [supply_demand, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% ADMM Different Norm Implementation Fisher Market - epsilon constraint in objective

%Consider N people and M goods
N = 5;
M = 5;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 1;

x_iter = 1;

p = p_init;
supply_demand2 = [];

beta = 0.2;

gamma_plus = 0.1*ones(N, M);
gamma_minus = -0.1*ones(N, M);

y_in = (1/M)*ones(N, M);

while  x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x' - gamma_plus(i, :)*(max(x - y_in(i, :)-epsilon/x_iter, 0)')  - gamma_minus(i, :)*(min(x - y_in(i, :)+epsilon/x_iter, 0)') ); %beta*norm(x-y_in(i, :), 1)
        subject to
        x >= zeros(1, M);%x - y_in(i, :) <= epsilon/x_iter;%y_in(i, :) - x <= epsilon/x_iter;
        cvx_end
        opt_xi(i, :) = x;
    end
    opt_xi = opt_xi + 0.5*beta*(y_in-opt_xi);
    %Find optimal y values
    cvx_begin
    variable y(N, M);
    maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* square_pos(norm(sum(y, 1) - C, 2)) );
    subject to
    y >= zeros(N, M);
    cvx_end
    y_in = y;
    
    
    p = p + beta*(sum(y_in, 1) - C);
    gamma_plus = max(gamma_plus + beta*(opt_xi - y_in -epsilon/x_iter), 0);
    gamma_minus = min(gamma_minus - beta*(opt_xi - y_in+epsilon/x_iter), 0);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand2 = [supply_demand2, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand2)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% ADMM Different Norm Implementation Fisher Market - epsilon constraint for both x and y

%Consider N people and M goods
N = 5;
M = 5;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 1;

x_iter = 1;

p = p_init;
supply_demand2 = [];

beta = 1;

y_in = (1/M)*ones(N, M);

while  x_iter <= 100 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x'  ); %beta*norm(x-y_in(i, :), 1)
        subject to
        x >= zeros(1, M);
        x - y_in(i, :) <= epsilon/x_iter;
        y_in(i, :) - x <= epsilon/x_iter;
        cvx_end
        opt_xi(i, :) = x;
    end
    opt_xi = opt_xi + 0.5*beta*(y_in-opt_xi);
    %Find optimal y values
    cvx_begin
    variable y(N, M);
    maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* square_pos(norm(sum(y, 1) - C, 2)) );
    subject to
    y >= zeros(N, M);
    opt_xi - y <= epsilon/x_iter;
    y - opt_xi <= epsilon/x_iter;
    cvx_end
    y_in = y;
    
    
    p = p + beta*(sum(y_in, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand2 = [supply_demand2, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand2)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% ADMM Implementation Fisher Market - agent maximizes IOP with epsilon constraint

%Consider N people and M goods
N = 10;
M = 10;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
V = rand(N, M);

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 1;

x_iter = 1;

p = p_init;
supply_demand2 = [];

beta = 0.5;

%gamma_plus = 0.1*ones(N, M);
%gamma_minus = -0.1*ones(N, M);

y_in = (1/M)*ones(N, M);

while  x_iter <= 50 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( V(i, :)*x' ); %beta*norm(x-y_in(i, :), 1)
        subject to
        p*x'<= B(i);
        x >= zeros(1, M);%x - y_in(i, :) <= epsilon/x_iter;%y_in(i, :) - x <= epsilon/x_iter;
        x - y_in(i, :) <= epsilon/x_iter;
        y_in(i, :) - x <= epsilon/x_iter;
        cvx_end
        opt_xi(i, :) = x;
    end
    opt_xi = opt_xi + 0.5*beta*(y_in-opt_xi);
    %Find optimal y values
    cvx_begin
    variable y(N, M);
    maximize( -(beta/2)*square_pos(norm(opt_xi-y, 2)) - (beta/2)* square_pos(norm(sum(y, 1) - C, 2)) );
    subject to
    y >= zeros(N, M);
    p*y'<= B';
    opt_xi - y <= epsilon/x_iter;
    y - opt_xi <= epsilon/x_iter;
    cvx_end
    y_in = y;
    
    
    p = p + beta*(sum(y_in, 1) - C);
    %gamma_plus = max(gamma_plus + beta*(opt_xi - y_in -epsilon/x_iter), 0);
    %gamma_minus = min(gamma_minus - beta*(opt_xi - y_in+epsilon/x_iter), 0);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand2 = [supply_demand2, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand2)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off