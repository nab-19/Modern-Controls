%% Setting constant values and finding matrices

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

%disp("J_inv_n");
%disp(j_inv_n);

j_inv_w = zeros(3,1);
j_inv_w(1,1) = (-m1*(l2^2) + m2*(l1^2) - m2*l2)/d0;
j_inv_w(2,1) = m1*l1/d0;
j_inv_w(3,1) = (m1*(l1^2) - m1*(l2^2) + m2*(l1^2) - m2*(l2^2))/(-l2*d0);

%disp("j_inv_w");
%disp(j_inv_w);

c_1 = j_inv_w;
c_2 = j_inv_n*j_inv_w;
c_3 = j_inv_n*j_inv_n*j_inv_w;

cont = zeros(3,3);
cont(1:3,1) = c_1(1:3,1);
cont(1:3,2) = c_2(1:3,1);
cont(1:3,3) = c_3(1:3,1);

%disp("cont");
%disp(cont);
%disp(rank(cont));

I3 = [1 0 0; 0 1 0; 0 0 1];
B = zeros(6,1);
B(4:6,1) = j_inv_w(1:3,1);
A = zeros(6,6);
A(1:3,4:6) = I3;
A(4:6,1:3) = j_inv_n(1:3,1:3);
C = [1 0 0 0 0 0 ; 0 1 0 0 0 0; 0 0 1 0 0 0 ; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];


%% LQR Controller, Continuous Time, Infinite Horizon


N = 500;
Qf = C; % cost of each state is unity 
Q = [1 0 0 0 0 0;
    0 100 0 0 0 0;
    0 0 100 0 0 0;
    0 0 0 4 0 0;
    0 0 0 0 400 0;
    0 0 0 0 0 400]; % strict penalty for large angles or large angular velocities
R = 1; % cost of input is equal to magnitude of input
G = -B*B';
[P1,K1,L1] = icare(A,B,Q,R,[],[],G);
disp("P");
disp(P1);
disp("K");
disp(K1);

sys2 = ss(A-B*K1,[],C,0);
sys2.OutputName = 'θ';

t3 = 0:0.05:15;
x3 = [0 0.05 -0.05 0 0 0];
figure(1);
u3 = zeros(length(t3),1);
lsim(sys2,u3,t3,x3);
title('Small Initial Conditions with Controller');
grid on;

t5 = 0:0.05:20;
x5 = [0 0.2 0.2 0 0 0];
figure(2);
u5 = zeros(length(t5),1);
lsim(sys2,u5,t5,x5);
title('Larger Initial Conditions with Controller');
grid on;

t6 = 0:0.05:30;
x6 = [0 pi()/2 pi()/2 0 0 0];
figure(3);
u6 = zeros(length(t6),1);
lsim(sys2,u6,t6,x6);
title('Half-Pi Initial Conditions with Controller');
grid on;

t7 = 0:0.05:45;
x7 = [0 pi() pi() 0 0 0];
figure(4);
u7 = zeros(length(t7),1);
lsim(sys2,u7,t7,x7);
title('Pi Initial Conditions with Controller');
grid on;


    