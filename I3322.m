prompt1="Tell us the dimension of your observables? ";
prompt2="Tell us the number of iterations?";
n=input(prompt1);
h=input(prompt2);

%random initialization
B0=eye(n);

initial_diag1 = randi([0, 1], [1, n]);
b1=diag(initial_diag1);
P=RandomUnitary(n,0);
B1=P*b1*P';

initial_diag2 = randi([0, 1], [1, n]);
b2=diag(initial_diag2);
B2=P*b2*P';

initial_diag3 = randi([0, 1], [1, n]);
b3=diag(initial_diag3);
B3=P*b3*P';

rawschmi = rand(n,1);
eleschmi= normalize(rawschmi,"norm");
fatschmi= diag(eleschmi);
O=[];

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
        if x12(i,i)>0
            a1(i,i)=1;
        end
        if x22(i,i)>0
            a2(i,i)=1;
        end
        if x32(i,i)>0
            a3(i,i)=1;
        end
    end
    A0=eye(n);
    A1=x11*a1*(x11)';
    A2=x21*a2*(x21)';
    A3=x31*a3*(x31)';
    N=A0*X0+A1*X1+A2*X2+A3*X3;
    o=trace(N);
    o=real(o)
    O(end+1)=o;
    %step 2
    Y0=(fatschmi*(-1*A2)*fatschmi).';
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
    for i=1:n
        if y12(i,i)>0
            b1(i,i)=1;
        end
        if y22(i,i)>0
            b2(i,i)=1;
        end
        if y32(i,i)>0
            b3(i,i)=1;
        end
    end
    B1=y11*b1*(y11');
    B2=y21*b2*(y21');
    B3=y31*b3*(y31');
    N=Y0*B0+Y1*B1+Y2*B2+Y3*B3;
    o=trace(N);
    o=real(o)
    O(end+1)=o;
    %step 3
    M=-(A0.*B2)-(A0.*B1)-2*(A0.*B2)+(A1.*B1)+(A1.*B2)+(A2.*B1)+(A2.*B2)-(A1.*B3)+(A2.*B3)-(A3.*B1)+(A3.*B2);
    rawschmi=eig(M);
    eleschmi= normalize(rawschmi,"norm");
    fatschmi= diag(eleschmi);
    Y0=(fatschmi*(-1*A2)*fatschmi).';
    Y1=(fatschmi*(A1+A2-A3-A0)*fatschmi).';
    Y2=(fatschmi*(A1+A2+A3-2*A0)*fatschmi).';
    Y3=(fatschmi*(-A1+A2)*fatschmi).';
    N=Y0*B0+Y1*B1+Y2*B2+Y3*B3;
    o=trace(N);
    o=real(o)
    O(end+1)=o;
end
H=1:3*h;
plot(H,O)

























