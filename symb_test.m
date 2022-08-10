clearvars
close all
clc;
tic;

syms beta alpha
syms l1 l2 l3 l4 l5 l6 l7
syms Z1 Z2 Z3 Z4 Z5 Z6 Z7
syms mh2 mH2 mA2 mHp2 mbar2 m122


% mbar2 == 2*m122/sin(2*beta)

% assuming Z2-symmetric potential (l6=l7=0)

% l = [l1,l2,l3,l4,l5];
% 
% Z = [ 
%         Z1 == l1*cos(beta)^4 + l2*sin(beta)^4 + (1/2)*(l3+l4+l5)*sin(2*beta)^2,
%         Z2 == l1*sin(beta)^4 + l2*cos(beta)^4 + (1/2)*(l3+l4+l5)*sin(2*beta)^2,
%         Z3 == (1/4)*sin(2*beta)^2*(l1+l2-2*(l3+l4+l5))+l3,
%         % Z4 == (1/4)*sin(2*beta)^2*(l1+l2-2*(l3+l4+l5))+l4,
%         Z5 == (1/4)*sin(2*beta)^2*(l1+l2-2*(l3+l4+l5))+l5,
%         Z6 == (-1/2)*sin(2*beta)*(l1*cos(beta)^2-l2*sin(beta)^2-(l3+l4+l5)*cos(2*beta)),
%         % Z7 == (-1/2)*sin(2*beta)*(l1*sin(beta)^2-l2*cos(beta)^2+(l3+l4+l5)*cos(2*beta)),
% ];

% assuming CP-conserving potential

l = [l1,l2,l3,l4,l5,l6,l7];

Z = [ Z1 == l1*cos(beta)^4 + l2*sin(beta)^4 + (1/2)*(l3+l4+l5)*sin(2*beta)^2 + 2*sin(2*beta)*(l6*cos(beta)^2+l7*sin(beta)^2),
        Z2 == l1*sin(beta)^4 + l2*cos(beta)^4 + (1/2)*(l3+l4+l5)*sin(2*beta)^2 + 2*sin(2*beta)*(l6*cos(beta)^2+l7*sin(beta)^2),
        Z3 == (1/4)*sin(2*beta)^2*(l1+l2-2*(l3+l4+l5))+l3+(l6-l7)*sin(2*beta)*cos(2*beta),
        Z4 == (1/4)*sin(2*beta)^2*(l1+l2-2*(l3+l4+l5))+l4+(l6-l7)*sin(2*beta)*cos(2*beta),
        Z5 == (1/4)*sin(2*beta)^2*(l1+l2-2*(l3+l4+l5))+l5+(l6-l7)*sin(2*beta)*cos(2*beta),
        Z6 == (-1/2)*sin(2*beta)*(l1*cos(beta)^2-l2*sin(beta)^2-(l3+l4+l5)*cos(2*beta))+(l6*cos(beta)*cos(3*beta)+l7*sin(beta)*sin(3*beta)),
        Z7 == (-1/2)*sin(2*beta)*(l1*sin(beta)^2-l2*cos(beta)^2+(l3+l4+l5)*cos(2*beta))+(l6*cos(beta)*cos(3*beta)+l7*sin(beta)*sin(3*beta)),
];

Y = solve(Z,l);

Y.l1 = simplify(Y.l1,'Steps',100);
Y.l2 = simplify(Y.l2,'Steps',100);
Y.l3 = simplify(Y.l3,'Steps',100);
Y.l4 = simplify(Y.l4,'Steps',100);
Y.l5 = simplify(Y.l5,'Steps',100);

Y.l6 = simplify(Y.l6,'Steps',100);
Y.l7 = simplify(Y.l7,'Steps',100);


Y.l1
Y.l2
Y.l3
Y.l4
Y.l5
Y.l6
Y.l7

% mass -> Higgs



