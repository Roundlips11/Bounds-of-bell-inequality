prompt1 = "Tell us the dimension of your observables? ";
prompt2 = "Tell us the number of iterations? ";
prompt3 = "Tell us the Matrix coefficient? ";

n = input(prompt1);
itr = input(prompt2);
M = input(prompt3);


N = zeros(n,n);
Mb = size(M,2); % number of measurement settings for Bob
Ma = size(M,1); % number of measurement settings for Alice

Xv = zeros(n,n,Ma); % Eigen vectors while optimizing Alice's operators
Xd = zeros(n,n,Ma); % Eigen values while optimizing Alice's operators ?

Yv = zeros(n,n,Mb); % Eigen vectors while optimizing Bob's operators
Yd = zeros(n,n,Mb); % Eigen values while optimizing Bob's operators ?

X_coef = zeros(n,n,Ma);
X = zeros(n,n,Ma);

Y_coef = zeros(n,n,Mb);
Y = zeros(n,n,Mb);

IE = zeros(n,n);

b = zeros(n,n,Mb); % Diagonal matrices for Y with 0 or -1 for negative and 1 for positive values
B = zeros(n, n, Mb); % Diagonal matrices for Bob's operators
B(:,:,1) = eye(n);
b(:,:,1)=eye(n);

a = zeros(n,n,Ma); % Diagonal matrices for X with 0 or -1 for negative and 1 for positive values
A = zeros(n,n,Ma); % Alice's operators
A(:,:,1) = eye(n);
a(:,:,1)=eye(n);




% Initializing random diagonal matrices for Bob's operators
for i = 2:Mb
    pick=[-1,1];
    initial_diag1 = randi([1,2],[1,n]);
    for t=1:n
        B(t,t,i)=pick(initial_diag1(1,t));
    end
end


% Performing random unitary transformations for 100 iterations
for q = 1:3
    for nu = 2:Mb
        P = RandomUnitary(n, 0);
        B(:,:,nu) = P * B(:,:,nu) * P';
    end
end


% Schmidt coefficient
rawschmi = rand(n,1);
eleschmi = normalize(rawschmi, "norm");
lamda = diag(eleschmi);
lamda_conj= lamda'; % conjugate leaving ther complex
B_val = [];

for k = 1:itr

    % Optimizing Alice's measurement operator while fixing Bob's and state vector

    % Step-1: Compute matrix X
    for mu = 1:Ma
        X_coef(:,:,mu) = zeros(n, n); % Reset X_coef for each mu
        for nu = 1:Mb
            X_coef(:,:,mu) = X_coef(:,:,mu) + M(mu, nu) * B(:,:,nu);
        end
        X(:,:,mu) = (lamda * (X_coef(:,:,mu)) * lamda_conj).';
    end

    % Diagonalizing X matrices(problematic for X=0 cases)
    for mu = 2:Ma
        [Xv(:,:,mu), Xd(:,:,mu)] = eig(X(:,:,mu));
    end

    % Step-2: Create diagonal matrices for X
    for mu = 2:Ma
        for i = 1:n
            if real(Xd(i,i,mu)) >= 0
                a(i,i,mu) = 1;
            else
                a(i,i,mu) = -1;
            end
        end
        A(:,:,mu) = Xv(:,:,mu) * a(:,:,mu) * (Xv(:,:,mu)');
    end

    % Compute Alice's operators
    N = zeros(n, n); % Reset N for each iteration
    for mu = 1:Ma
        N = N +  A(:,:,mu)*X(:,:,mu)  ;
    end
    val = abs(trace(N))
    if val>3
        Weirdlamda=lamda;
        WeirdA=A;
        WeirdB=B;
    end
    B_val(end+1) = val;

    % Optimizing Bob's measurement operator while fixing Alice's and state vector

    % Step-1: Compute matrix Y
    for nu = 1:Mb
        Y_coef(:,:,nu) = zeros(n, n); % Reset Y_coef for each nu
        for mu = 1:Ma
            Y_coef(:,:,nu) = Y_coef(:,:,nu) + M(mu, nu) * A(:,:,mu);
        end
        Y(:,:,nu) = (lamda * (Y_coef(:,:,nu)) * lamda_conj).';
    end

    % Diagonalizing Y matrices
    for nu = 1:Mb
        [Yv(:,:,nu), Yd(:,:,nu)] = eig(Y(:,:,nu));
    end

    % Step-2: Create diagonal matrices for Y
    for nu = 2:Mb
        for i = 1:n
            if real(Yd(i, i, nu)) >= 0
                b(i, i, nu) = 1;
            else
                b(i, i, nu) = -1;
            end
        end
        B(:,:,nu) = Yv(:,:,nu) * b(:,:,nu) * (Yv(:,:,nu)');
    end

    % Compute Bob's operators
    N = zeros(n, n); % Reset N for each iteration
    for nu = 1:Mb
        
        N = N + Y(:,:,nu) * B(:,:,nu);
    end
    val = abs(trace(N))
    if val>3
        Weirdlamda=lamda;
        WeirdA=A;
        WeirdB=B;
    end
    B_val(end+1) = val;

    % Optimizing the state vector while fixing Alice's and Bob's operators

    % Step-1: Define the Inequality matrix
    IE = zeros(n, n); % Reset IE for each iteration
    for mu = 1:Ma
        for nu = 1:Mb
            IE = IE + M(mu, nu) * (A(:,:,mu) .* B(:,:,nu));
        end
    end

    
    [IEv, IEd] = eig(IE);
    IE_max = max(eig(IE));


    for r = 1:n
        if IE_max == IEd(r, r)
            rawschmi = IEv(:, r);
        end
    end

    lamda = diag(rawschmi);
    lamda_conj=lamda';
    N = zeros(n, n);
    for nu = 1:Mb
        Y(:,:,nu) = (lamda * (Y_coef(:,:,nu)) * lamda_conj).';
        N = N + Y(:,:,nu) * B(:,:,nu);
    end
    val = abs(trace(N))
    if val>3
        Weirdlamda=lamda;
        WeirdA=A;
        WeirdB=B;
    end
    B_val(end+1) = val;
end
K=cos(angle(rawschmi))+i*sin(angle(rawschmi));
U=diag(K);
for mu=1:Ma
    A(:,:,mu)=U'*A(:,:,mu)*U;
end
for nu=1:Mb
    B(:,:,nu)=U'*B(:,:,nu)*U;
end
lamda=abs(lamda);
rho_a=lamda*lamda;
rho_b=lamda*lamda;
badschmi=zeros(n*n,1);
for c=1:n
    badschmi(n*(c-1)+c,1)=lamda(c,c);
end
rho=badschmi*badschmi';
H = 1:3*itr;
plot(H, B_val);
