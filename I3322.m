prompt1="Tell us the dimension of your observables? ";
prompt2="Tell us the number of iterations?";
n=input(prompt1);
h=input(prompt2);

%random initialization
B0=eye(n);
b1=zeros(n,n);
b2=zeros(n,n);
b3=zeros(n,n);
while all(all(b1==zeros(n,n))) || all(all(b1==eye(n)))
    initial_diag1 = randi([0, 1], [1, n]);
    b1=diag(initial_diag1);
end
while all(all(b2==zeros(n,n))) || all(all(b2==eye(n)))
    initial_diag2 = randi([0, 1], [1, n]);
    b2=diag(initial_diag2);
end
while all(all(b3==zeros(n,n))) || all(all(b3==eye(n)))
    initial_diag3 = randi([0, 1], [1, n]);
    b3=diag(initial_diag3);
end
B1=b1;
B2=b2;
B3=b3;
for q=1:r
    P=RandomUnitary(n,0);
    B1=P*B1*P';
    B2=P*B2*P';
    B3=P*B3*P';
end
rawschmi = rand(n,1);
eleschmi= normalize(rawschmi,"norm");
fatschmi= diag(eleschmi);
O=[];
K=[];

for k=1:h
    X0=(fatschmi*(-1*B1-2*B2)*fatschmi).';
    X1=(fatschmi*(B1+B2-B3)*fatschmi).';
    X2=(fatschmi*(B1+B2+B3-B0)*fatschmi).';
    X3=(fatschmi*(-B1+B2)*fatschmi).';
    % diagonalising Xs
    [x01,x02]=eig(X0);
    [x11,x12]=eig(X1);
    [x21,x22]=eig(X2);
    [x31,x32]=eig(X3);
    %step 1
    a1=zeros(n);
    a2=zeros(n);
    a3=zeros(n);
    for i=1:n
        if real(x12(i,i))>0
            a1(i,i)=1;
        end
        if real(x22(i,i))>0
            a2(i,i)=1;
        end
        if real(x32(i,i))>0
            a3(i,i)=1;
        end
    end
    A0=eye(n);
    A1=x11*a1*(x11)';
    A2=x21*a2*(x21)';
    A3=x31*a3*(x31)';
    N=A0*X0+A1*X1+A2*X2+A3*X3;
    o=trace(N);
    o=real(o);
    O(end+1)=o;
    %step 2
    Y0=(fatschmi*(-A2)*fatschmi).';
    Y1=(fatschmi*(A1+A2-A3-A0)*fatschmi).';
    Y2=(fatschmi*(A1+A2+A3-2*A0)*fatschmi).';
    Y3=(fatschmi*(-A1+A2)*fatschmi).';
    [y01,y02]=eig(Y0);
    [y11,y12]=eig(Y1);
    [y21,y22]=eig(Y2);
    [y31,y32]=eig(Y3);
    b1=zeros(n);
    b2=zeros(n);
    b3=zeros(n);
    for j=1:n
        if real(y12(j,j))>0
            b1(j,j)=1;
        end
        if real(y22(j,j))>0
            b2(j,j)=1;
        end
        if real(y32(j,j))>0
            b3(j,j)=1;
        end
    end
    B1=y11*b1*(y11');
    B2=y21*b2*(y21');
    B3=y31*b3*(y31');
    N=Y0*B0+Y1*B1+Y2*B2+Y3*B3;
    o=trace(N);
    o=real(o);
    O(end+1)=o;
    %step 3 (can be optimized further)
    M=-(A0.*B1)-2*(A0.*B2)+(A1.*B1)+(A1.*B2)-(A1.*B3)+(A2.*B1)+(A2.*B2)+(A2.*B3)-(A2.*B0)-(A3.*B1)+(A3.*B2);
    m2=max(eig(M));
    [m1,e]=eig(M);
    for u=1:n
        if m2==e(u,u)
            rawschmi=m1(:,u);
        end
    end

    eleschmi= normalize(rawschmi,"norm");
    fatschmi= diag(rawschmi);
    X0=(fatschmi*(-B1-2*B2)*fatschmi).';
    X1=(fatschmi*(B1+B2-B3)*fatschmi).';
    X2=(fatschmi*(B1+B2+B3-B0)*fatschmi).';
    X3=(fatschmi*(-B1+B2)*fatschmi).';
    N=A0*X0+A1*X1+A2*X2+A3*X3;
    o=abs(trace(N));
    m2;
    O(end+1)=o;
    K(end+1)=m2;
end
H=1:3*h;
plot(H,O)

























