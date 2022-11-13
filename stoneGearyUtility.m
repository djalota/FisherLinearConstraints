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

%Vector of W values
W = rand(N, M)/10;

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = zeros(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand2 = [];

beta = 5;

y_in = (1/M)*ones(N, M);

while  x_iter <= 200 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*(log(x - W(i, :)) *V(i, :)' ) - p*x' - (beta/2)* square_pos(norm(x-y_in(i, :), 2)) );
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
    
    disp(x_iter)
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
N = 10;
M = 10;

%The budget B is a N vector containing the total budget available for all
%people
B = ones(1, N)';

%The capacity C is a M vector containing the maximum number of people that
%can use a given good
C = (N/M)*ones(1, M);

%The valuation V is a N*M matrix
%V = rand(N, M);

%Vector of W values
%W = rand(N, M)/10;

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = 0.1*ones(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand3 = [];

beta = 1;

y_in = (1/M)*ones(N, M);

while x_iter <= 200 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*(log(x - W(i, :)) *V(i, :)' ) - p*x'  );
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
    
    disp(x_iter)
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

%% ADMM Implementation Fisher Market - Linear Utilities

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
supply_demand3 = [];

beta = 1;

y_in = (1/M)*ones(N, M);

while  x_iter <= 200 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
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
    supply_demand3 = [supply_demand3, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand3)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'Two-Block ADMM', 'Two-Block AMA', 'Two-Block AAMA', 'Augmented Lagrangian')
hold off

%% AMA Implementation Fisher Market - Linear Utilities

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
%V = rand(N, M);

%Vector of W values
%W = rand(N, M)/10;

%Initialize the Prices in the market
p_init = rand(1, M);

%Create matrix to store the optimal xi
opt_xi = 0.1*ones(N, M);

%Convergence criteria
epsilon = 0.001;

x_iter = 1;

p = p_init;
supply_demand4 = [];

beta = 0.1;

y_in = (1/M)*ones(N, M);

while x_iter <= 200 %max(abs(sum(opt_xi, 1) - C)) > epsilon &
    %Find optimal x values
    for i = 1:N
        cvx_begin
        variable x(1, M);
        maximize( B(i)*log(V(i, :)*x') - p*x'  );
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
    
    disp(x_iter)
    p = p + beta*(sum(y_in, 1) - C);
    
    sup_dem = sum(opt_xi, 1) - C;
    my_sum = 0;
    for j = 1:M
        my_sum = my_sum + sup_dem(j)^2;
    end
    supply_demand4 = [supply_demand4, abs(my_sum)];
    x_iter = x_iter + 1;
end
hold on
semilogy(1:x_iter-1, supply_demand4, 'LineWidth', 3)
%ylabel('Difference in Supply and Demand', 'FontSize', 20, 'Interpreter','latex')
%xlabel('Iterations', 'FontSize', 20, 'Interpreter','latex')
%legend('Dual Ascent', 'ADMM')
hold off