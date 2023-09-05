clc;
clear;
g   =   9.81;           %m/s^2
b   =   3.59e-5;        %thrust coefficient        
d   =   2.081e-6;       %drag(torque) coefficient
m   =   0.82;           %kg
Ixx =   0.00963;        %kg*m^2
Iyy =   0.00963;        %kg*m^2
Izz =   0.019;          %kg*m^2
Jr  =   5.225e-5;       %kg*m^2
l   =   0.24;           %m
Ax  =   0.0;            %kg/s  
Ay  =   0.0;            %kg/s
Az  =   0.0;            %kg/s

%motor simülasyonu için katsayılar
Jr_motor =   6.5e-7;         %rotor moment of inertia from data sheet
b_motor  =   2.415e-6;       %damping coef
K_emf    =   0.00255;        %Nm/A
K_t      =   0.00255;        %Nm/A
R_m      =   0.117;          %ohm
L        =   1.17e-4;        %henry
omega_max=(90/100)*9000*(2*pi/(60)); % 90% of 9000 rpm to rad/s
% aerodynamic force and moments constant
Kf=3.13e-5;
KM=7.5e-7;

%Linearization
%x=[x y z xdot ydot zdot phi theta psi phidot thetadot psidot]'
A= [0 0 0 1 0 0 0 0 0 0 0 0;...   
    0 0 0 0 1 0 0 0 0 0 0 0;...
    0 0 0 0 0 1 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 -g 0 0 0 0;...
    0 0 0 0 0 0 g 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 1 0 0;...
    0 0 0 0 0 0 0 0 0 0 1 0;...
    0 0 0 0 0 0 0 0 0 0 0 1;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0];
 
 B=[0 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    1/m 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    0 1/Ixx 0 0;...
    0 0 1/Iyy 0;...
    0 0 0 1/Izz;];
 
 C =[1 0 0 0 0 0 0 0 0 0 0 0;...
     0 1 0 0 0 0 0 0 0 0 0 0;...
     0 0 1 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 1 0 0 0 0 0;...
     0 0 0 0 0 0 0 1 0 0 0 0;...
     0 0 0 0 0 0 0 0 1 0 0 0;];

%controllabilty matrix
Cc=[B A*B A^2*B A^3*B A^4*B A^5*B A^6*B A^7*B A^8*B A^9*B A^10*B A^11*B A^12*B];
rank(Cc);

Q_1=[45 45 45 50 50 50 1000 1000 1000 1000 1000 1000];
Q = diag(Q_1);

D = zeros(6,4);

R = eye(4);
sys=ss(A,B,C,zeros(6,4));
K = lqr(A,B,Q,R);

Q_1=[45 45 900 50 50 100 1000 1000 1000 1000 1000 1000];
Q = diag(Q_1);           
        
sys=ss(A,B,C,zeros(6,4));
K_1 = lqr(A,B,Q,R);
[D D]
G = tf(sys);
p = size(C,1);
[n,m] = size(B);
R = eye(m);
Bnoise = eye(n);
W = eye(12);
V = 0.01*eye(4);
Estss = ss(A,[B B],C,[D D]);
[Kess,Ll] = kalman(Estss,W,V);
