n = input('Enter the dimension of your observables: ');
iterations = input('Enter the number of iterations: ');
M = input('Enter the Matrix coefficient: ');
P=input('what are the eigenvalues of your operator? (in ascending order) ');
prompt="how many times do you want to run the program?";
Q=input(prompt);

L=[];
function l=quantum_optimization(n,iterations,M,P)
    % Get user inputs
    
    % Initialize variables
    [Ma, Mb] = size(M);  % Measurement settings for Alice and Bob
    [A, B, lambda] = initialize_operators(n, Ma, Mb,P);
    
    % Storage for optimization values
    optimization_values = zeros(1, 3*iterations);
    
    % Main optimization loop
    for k = 1:iterations
        % Optimize Alice's operators
        [A, current_val] = optimize_alice(A, B, M, lambda, n, Ma, Mb,P);
        optimization_values(3*k-2) = current_val;
        
        % Optimize Bob's operators
        [B, current_val] = optimize_bob(A, B, M, lambda, n, Ma, Mb,P);
        optimization_values(3*k-1) = current_val;
        
        % Optimize state vector
        [lambda, current_val] = optimize_state(A, B, M, n, Ma, Mb);
        optimization_values(3*k) = current_val;
    end
    
    % Final phase adjustment
    [A, B, lambda] = adjust_phases(A, B, lambda, Ma, Mb, n);
    
    % Calculate final state
    rho = calculate_final_state(lambda, n);
    
    % Plot results
    %plot(1:3*iterations, optimization_values);
    %title('Optimization Progress');
    %xlabel('Iteration Step');
    %ylabel('Objective Value');
    l=optimization_values(end);

end

%% Helper Functions

function [A, B, lambda] = initialize_operators(n, Ma, Mb,P)
    % Initialize Alice's operators
    A = zeros(n, n, Ma);
    A(:,:,1) = eye(n);
    
    % Initialize Bob's operators with random unitaries
    B = zeros(n, n, Mb);
    B(:,:,1) = eye(n);
    
    for i = 2:Mb
        % Create random diagonal matrix with ±1 eigenvalues
        diag_vals = randi([1,2], [1,n]);
        for t=1:n
            B(t,t,i)=P(diag_vals(1,t));
        end
        % Apply random unitary transformations
        for q = 1:3
            G = RandomUnitary(n, 0);
            B(:,:,i) = G * B(:,:,i) * G';
        end
    end
    
    % Initialize Schmidt coefficients
    lambda = diag(normalize(rand(n,1), "norm"));
end

function [A, current_val] = optimize_alice(A, B, M, lambda, n, Ma, Mb,P)
    % Compute X matrices
    X = zeros(n, n, Ma);
    
    for mu = 1:Ma
        X_coef = zeros(n, n);
        for nu = 1:Mb
            X_coef = X_coef + M(mu, nu) * B(:,:,nu);
        end
        X(:,:,mu) = (lambda * X_coef * lambda').';
    end
    
    % Diagonalize X matrices and update Alice's operators
    for mu = 2:Ma
        a=zeros(n);
        [V, D] = eig(X(:,:,mu));
        for i=1:n
            if real(D(i,i))>0
                a(i,i)=P(2);
            else
                a(i,i)=P(1);
            end
        

        end
        % Diagonal matrix with sign of eigenvalues
        A(:,:,mu) = V * a * V';
    end
    
    % Compute current optimization value
    current_val = compute_objective(A, X, Ma, n);
end

function [B, current_val] = optimize_bob(A, B, M, lambda, n, Ma, Mb,P)
    % Compute Y matrices
    Y = zeros(n, n, Mb);
    for nu = 1:Mb
        Y_coef = zeros(n, n);
        for mu = 1:Ma
            Y_coef = Y_coef + M(mu, nu) * A(:,:,mu);
        end
        Y(:,:,nu) = (lambda * Y_coef * lambda').';
    end
    
    % Diagonalize Y matrices and update Bob's operators
    for nu = 2:Mb
        b=zeros(n);
        [V, D] = eig(Y(:,:,nu));
        for i=1:n
            if real(D(i,i))>0
                b(i,i)=P(2);
            else
                b(i,i)=P(1);
            end
        end  
        % Diagonal matrix with sign of eigenvalues
        B(:,:,nu) = V * b * V';
    end
    
    % Compute current optimization value
    current_val = compute_objective(B, Y, Mb, n);
end

function [lambda, current_val] = optimize_state(A, B, M, n, Ma, Mb)
    % Compute inequality matrix
    IE = zeros(n, n);
    for mu = 1:Ma
        for nu = 1:Mb
            IE = IE + M(mu, nu) * (A(:,:,mu) .* B(:,:,nu));
        end
    end
    [m1,e]=eig(IE);
    [~, linear_index] = max(e(:));
    [row, ~] = ind2sub(size(e), linear_index);
    raw_schmidt=m1(:,row);
    % Normalize and create diagonal matrix
    lambda = diag(normalize(raw_schmidt, "norm"));
    
    % Compute current optimization value
    Y = zeros(n, n, Mb);
    Y_coef = zeros(n, n, Mb);
    for nu = 1:Mb
        for mu = 1:Ma
            Y_coef(:,:,nu) = Y_coef(:,:,nu) + M(mu, nu) * A(:,:,mu);
        end
        Y(:,:,nu) = (lambda * Y_coef(:,:,nu) * lambda').';
    end
    
    N = zeros(n, n);
    for nu = 1:Mb
        N = N + Y(:,:,nu) * B(:,:,nu);
    end
    current_val = abs(trace(N));
end

function val = compute_objective(operators, matrices, count, n)
    N = zeros(n, n);
    for i = 1:count
        N = N + operators(:,:,i) * matrices(:,:,i);
    end
    val = abs(trace(N));
end

function [A, B, lambda] = adjust_phases(A, B, lambda, Ma, Mb, n)
    % Remove complex phases from the state vector
    phases = exp(1i * angle(diag(lambda)));
    U = diag(conj(phases));
    
    for mu = 1:Ma
        A(:,:,mu) = U' * A(:,:,mu) * U;
    end
    
    for nu = 1:Mb
        B(:,:,nu) = U' * B(:,:,nu) * U;
    end
    
    lambda = abs(lambda);
end

function rho = calculate_final_state(lambda, n)
    % Create Schmidt vector
    schmidt_vector = zeros(n*n, 1);
    for c = 1:n
        schmidt_vector(n*(c-1)+c, 1) = lambda(c,c);
    end
    
    % Create density matrix
    rho = schmidt_vector * schmidt_vector';
end
for l=1:Q
    L(end+1)=quantum_optimization(n,iterations,M,P);

end
J=1:Q;
plot(J,L);
