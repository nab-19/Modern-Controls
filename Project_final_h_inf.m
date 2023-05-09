%% H_infinity_controller
g = 9.8;
m0 = 3;
m1 = 0.5;
m2 = 1;
l1 = 1.5;
l2 = 0.75;
d0 = (m1^2)*(l1^2) - (m1^2)*(l2^2) - m0*m1*(l2^2)+m0*m2*(l1^2)-m0*m2*(l2^2)+m1*m2*(l1^2)-m1*m2*(l2^2);
c = -1*d0 + (m2^2)*(l2^2) + m1*m2*(l2^2) - m1*m2*(l1^2);

j_inv_n = zeros(3,3);
j_inv_n(1,2) = g*m1*(m1+m2)*(l1^2)/d0;
j_inv_n(1,3) = -g*m2*(m1*(l1^2)-m1*(l2^2)+m2*(l1^2)-m2*(l2^2))/d0;
j_inv_n(2,2) = -g*(m0+m1)*(m1+m2)*l1/d0;
j_inv_n(2,3) = g*m0*m2*l1/d0;
j_inv_n(3,2) = g*m0*(m1+m2)*(l1^2)/d0;
j_inv_n(3,3) = -g*c/(l2*d0);

j_inv_w = zeros(3,1);
j_inv_w(1,1) = (-m1*(l2^2) + m2*(l1^2) - m2*l2)/d0;
j_inv_w(2,1) = m1*l1/d0;
j_inv_w(3,1) = (m1*(l1^2) - m1*(l2^2) + m2*(l1^2) - m2*(l2^2))/(-l2*d0);

B = zeros(6,1);
B(4:6,1) = j_inv_w(1:3,1);
A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4:6,1:3) = j_inv_n(1:3,1:3);
C = eye(6);

% controller synthesis

B2 = B;
B1 = zeros(6,1);
C1 = C;
C2 = C;
D11 = zeros(6,1);
D12 = zeros(6,1);
D21 = zeros(6,1);
D22 = zeros(6,1);

cvx_begin sdp
    variable rho3
    variable X1(6,6) symmetric
    variable Y1(6,6) symmetric
    variables An(6,6) Bn(6,6) Cn(1,6) Dn(1,6)
    minimize rho3
    subject to
        [Y1     eye(6);
         eye(6)     X1] >= 0;

        [A*Y1 + Y1*A' + B2*Cn + Cn'*B2'        (A' + An + (B2*Dn*C2)')'            B1 + B2*Dn*D21       (C1*Y1 + D12*Cn)';
             A' + An + (B2*Dn*C2)'        X1*A + A'*X1 + Bn*C2 + C2'*Bn'         X1*B1+Bn*D21           (C1+D12*Dn*C2)';
             (B1 + B2*Dn*D21)'                   (X1*B1+Bn*D21)'                    -rho3           (D11 + D12*Dn*D21)';
             C1*Y1 + D12*Cn                       C1 + D12*Dn*C2                 D11 + D12*Dn*D21         -rho3*eye(6)] <=0;
cvx_end

Y2 = eye(6);
X2 = eye(6) - X1*Y1;

K21 = zeros(7,7);
K21(1:6,1:6) = X2;
K21(1:6,7) = X1*B2;
K21(7,7) = 1;

K22 = zeros(7,12);
K22(1:6,1:6) = An-X1*A*Y1;
K22(1:6,7:12) = Bn;
K22(7,1:6) = Cn;
K22(7,7:12) = Dn;

K23 = zeros(12,12);
K23(1:6,1:6) = Y2';
K23(7:12,1:6) = C2*Y1;
K23(7:12,7:12) = eye(6);

K2 = inv(K21)*K22*inv(K23);
disp('K2 matrix');
disp(K2);

AK2 = K2(1:6,1:6);
BK2 = K2(1:6,7:12);
CK2 = K2(7,1:6);
DK2 = K2(7,7:12);

DK = (eye(1)+DK2*D22)\DK2;
BK = BK2*(eye(6)-D22*DK);
CK = (eye(1)-DK*D22)*CK2;
AK = AK2 - BK/(eye(6)-D22*DK)*D22*CK;

disp('AK');
disp(AK);
disp('BK');
disp(BK);
disp('CK');
disp(CK);
disp('DK');
disp(DK);

Acl = zeros(12,12);
Acl(1:6,1:6) = A + B2*DK2*C2;
Acl(1:6,7:12) = B2*CK2;
Acl(7:12,1:6) = BK2*C2;
Acl(7:12,7:12) = AK2;

Bcl = zeros(12,1);
Bcl(1:6,1) = B1+B2*DK2*D21;
Bcl(7:12,1) = BK2*D21;

Ccl = zeros(6,12);
Ccl(1:6,1:6) = C1+D12*DK2*C2;
Ccl(1:6,7:12) = D12*CK2;

Dcl = D11+D12*DK2*D21;

disp('Acl');
disp(Acl);
disp('Bcl');
disp(Bcl);
disp('Ccl');
disp(Ccl);
disp('Dcl');
disp(Dcl);


%% Simulations with H Infinity Controller
sys = ss(Acl,[],Ccl,0);
sys.OutputName = 'θ';

t1 = 0:0.05:15;
x1 = [0 0.5 -0.5 0 0 0 0 0.5 -0.5 0 0 0];
figure(1);
w1 = zeros(length(t1),1);
w(t1>1)=1;
lsim(sys,w1,t1,x1);
title('Small Initial Conditions with Controller');
grid on;

t2 = 0:0.05:30;
x2 = [0 pi()/2 pi()/2 0 0 0 0 pi()/2 pi()/2 0 0 0];
figure(2);
w2 = zeros(length(t2),1);
lsim(sys,w2,t2,x2);
title('Half-Pi Initial Conditions with Controller');
grid on;

t3 = 0:0.05:45;
x3 = [0 pi() pi() 0 0 0 0 pi() pi() 0 0 0];
figure(3);
w3 = zeros(length(t3),1);
lsim(sys,w3,t3,x3);
title('Pi Initial Conditions with Controller');
grid on;


    