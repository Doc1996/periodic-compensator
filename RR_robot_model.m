%% SYMBOLIC PARAMETERS OF THE ROBOT

clear all;
clc;
close all;

% a vertical RR robot is observed whose initial position is such that its first and second articles are directed to
% the right, with its axes x0, x1 and x2 directed to the right, and z0, z1 and z2 from the paper

syms theta_k d_k a_k alpha_k;
syms q1 q2 der_q1 der_q2 der2_q1 der2_q2;
syms L1 L2 L1_cm L2_cm m1 m2 I1 I2 g;

DH_10 = [q1, 0, L1, 0];
DH_21 = [q2, 0, L2, 0];


%% DIRECT KINEMATICS OF THE ROBOT

% T_kBkA represents the homogeneous transformation matrix of system k-1 into system k

T_kBkA = [cos(theta_k), -cos(alpha_k)*sin(theta_k), sin(alpha_k)*sin(theta_k), a_k*cos(theta_k);
          sin(theta_k), cos(alpha_k)*cos(theta_k), -sin(alpha_k)*cos(theta_k), a_k*sin(theta_k);
          0, sin(alpha_k), cos(alpha_k), d_k;
          0, 0, 0, 1];

T_10 = subs(T_kBkA, [theta_k, d_k, a_k, alpha_k], DH_10);
T_21 = subs(T_kBkA, [theta_k, d_k, a_k, alpha_k], DH_21);

T_20 = simplify(T_10*T_21);

p_10 = T_10(1:3, 4);
p_20 = T_20(1:3, 4);

p_10_cm = subs(p_10, L1, L1_cm);
p_20_cm = subs(p_20, L2, L2_cm);

r3_10 = T_10(1:3, 3);
r3_20 = T_20(1:3, 3);

r3_10_cm = subs(r3_10, L1, L1_cm);
r3_20_cm = subs(r3_20, L2, L2_cm);


%% DYNAMICS OF THE ROBOT

% the positive direction of the y0 axis of robot represents the vector along which acts -g

% since the derivative of the Lagrangian must, and the Jacobians must not, be calculated with the time variable, some
% expressions will have to be substituted before calculating the variable

% Jacobians have the same number of columns as the number of joints

% the Jacobian of the translation of a coordinate system towards the base is defined as the partial derivative of the
% center of mass coordinates for that transformation by all joint variables

% The Jacobian of the rotation of a coordinate system towards the base is defined as the matrix of all the ordered
% vectors of approaching the centers of mass for the transformations to the selected joint (eg J_10 is [0; 0; 1], and
% J_30 is [[0; 0; 1], r3_10, r3_30_cm ], only the last vector refers to the center of mass), where the first vector is
% always [0; 0; 1], and the column whose joint is translational instead of an approach vector has a zero-vector of the
% same dimension, while the columns for which approach vectors cannot be defined are equal to zero-vectors

% for the potential energy of an individual link in relation to the base of the robot, the y coordinate of that link
% is observed


syms q1_t(t) q2_t(t);

q = [q1; q2];
der_q = [der_q1; der_q2];
der2_q = [der2_q1; der2_q2];

q_t = [q1_t; q2_t];
der_q_t = diff(q_t, t);
der2_q_t = diff(der_q_t, t);


J_p_10_cm = jacobian(p_10_cm, q);
der_p_10_cm = J_p_10_cm * der_q;

J_r3_10_cm = [[0; 0; 1], [0; 0; 0]];
der_phi_10_cm = J_r3_10_cm * der_q;

J_p_20_cm = jacobian(p_20_cm, q);
der_p_20_cm = J_p_20_cm * der_q;

J_r3_20_cm = [[0; 0; 1], r3_10_cm];
der_phi_20_cm = J_r3_20_cm * der_q;


der_q_10 = der_q1;
der_q_20 = der_q1 + der_q2;

K_10 = m1/2*(der_p_10_cm.')*der_p_10_cm + I1/2*(der_phi_10_cm.')*der_phi_10_cm;
K_20 = m2/2*(der_p_20_cm.')*der_p_20_cm + I2/2*(der_phi_20_cm.')*der_phi_20_cm;

K = K_10 + K_20;

U_10 = - m1*([0; -g; 0].')*p_10_cm;
U_20 = - m2*([0; -g; 0].')*p_20_cm;

U = U_10 + U_20;


Lagr = simplify(K - U);

der_Lagr_by_q = simplify(jacobian(Lagr, q));

der_Lagr_by_der_q = simplify(jacobian(Lagr, der_q));

der_Lagr_t_by_der_q = subs(der_Lagr_by_der_q, [q, der_q, der2_q], [q_t, der_q_t, der2_q_t]);

der_Lagr_t_by_der_q_and_t = diff(der_Lagr_t_by_der_q, t);

der_Lagr_by_der_q_and_t = simplify(subs(der_Lagr_t_by_der_q_and_t, [q_t(t), diff(q_t(t), t), diff(q_t(t), t, t)], ...
    [q, der_q, der2_q]));

tau = simplify(der_Lagr_by_der_q_and_t - der_Lagr_by_q).';


% matrix M can be calculated as SUM_BY_i[m_i*(J_p_i_cm.')*J_p_i_cm + I_i*(J_r3_i_cm.')*J_r3_i_cm], where i represents
% the transformation from the base to the coordinate system of one joint

% matrix C can be calculated as der_M_by_t*der_q - 1/2*MATRIX_WITH_NEW_ROWS_FOR_i[(der_q.')*der_M_by_q_i*der_q]

% Jacobian of the matrix is calculated as a matrix composed of partial derivatives of original matrix per one variable,
% where new rows are used for each variable

% matrix G can be calculated as a partial derivative of the total potential energy (the sum of potential energies from
% the base to the coordinate system of each joint) by joint variables


M = simplify(m1*(J_p_10_cm.')*J_p_10_cm + I1*(J_r3_10_cm.')*J_r3_10_cm + m2*(J_p_20_cm.')*J_p_20_cm + ...
    I2*(J_r3_20_cm.')*J_r3_20_cm);

G = simplify(jacobian(U, q).');

C_with_der_q = simplify(tau - M*der2_q - G);


% der2_q1 and der2_q2 expressed through the control laws of individual joints (located only next to matrix M) can be
% written through the following system of equations, where control variables u1 and u2 are used instead of tau):

% M11*der2_q1 + M21*der2_q2 = N1
% M21*der2_q1 + M22*der2_q2 = N2


syms u1 u2;

M11 = M(1, 1); M12 = M(1, 2); M21 = M(2, 1); M22 = M(2, 2);
N1 = u1 - C_with_der_q(1) - G(1);
N2 = u2 - C_with_der_q(2) - G(2);

der2_q_expressed = simplify(M\([u1; u2] - C_with_der_q - G));
der2_q1_expressed = der2_q_expressed(1);
der2_q2_expressed = der2_q_expressed(2);


% der2_q1_expressed:
% (2*I2*u1 - 2*I2*u2 + 2*L2_cm^2*m2*u1 - 2*L2_cm^2*m2*u2 - L1*L2_cm^2*g*m2^2*cos(q1) + L1^2*L2_cm^2*der_q1^2*m2^2*sin(2*q2) + L1*L2_cm^2*g*m2^2*cos(q1 + 2*q2) - 2*I2*L1*g*m2*cos(q1) - 2*I2*L1_cm*g*m1*cos(q1) - 2*L1*L2_cm*m2*u2*cos(q2) + 2*L1*L2_cm^3*der_q1^2*m2^2*sin(q2) + 2*L1*L2_cm^3*der_q2^2*m2^2*sin(q2) + 2*I2*L1*L2_cm*der_q1^2*m2*sin(q2) + 2*I2*L1*L2_cm*der_q2^2*m2*sin(q2) - 2*L1_cm*L2_cm^2*g*m1*m2*cos(q1) + 4*L1*L2_cm^3*der_q1*der_q2*m2^2*sin(q2) + 4*I2*L1*L2_cm*der_q1*der_q2*m2*sin(q2))/(2*I1*I2 + L1^2*L2_cm^2*m2^2 + 2*I2*L1^2*m2 + 2*I2*L1_cm^2*m1 + 2*I1*L2_cm^2*m2 + 2*L1_cm^2*L2_cm^2*m1*m2 - L1^2*L2_cm^2*m2^2*cos(2*q2))
% der2_q2_expressed:
% -(I2*u1 - I1*u2 - I2*u2 - L1^2*m2*u2 - L1_cm^2*m1*u2 + L2_cm^2*m2*u1 - L2_cm^2*m2*u2 + L1^2*L2_cm*g*m2^2*cos(q1 + q2) - L1*L2_cm^2*g*m2^2*cos(q1) + I1*L2_cm*g*m2*cos(q1 + q2) + L1^2*L2_cm^2*der_q1^2*m2^2*sin(2*q2) + (L1^2*L2_cm^2*der_q2^2*m2^2*sin(2*q2))/2 - I2*L1*g*m2*cos(q1) - I2*L1_cm*g*m1*cos(q1) + L1*L2_cm*m2*u1*cos(q2) - 2*L1*L2_cm*m2*u2*cos(q2) + L1*L2_cm^3*der_q1^2*m2^2*sin(q2) + L1^3*L2_cm*der_q1^2*m2^2*sin(q2) + L1*L2_cm^3*der_q2^2*m2^2*sin(q2) + L1*L2_cm^2*g*m2^2*cos(q1 + q2)*cos(q2) - L1^2*L2_cm*g*m2^2*cos(q1)*cos(q2) + L1_cm^2*L2_cm*g*m1*m2*cos(q1 + q2) + I1*L1*L2_cm*der_q1^2*m2*sin(q2) + I2*L1*L2_cm*der_q1^2*m2*sin(q2) + I2*L1*L2_cm*der_q2^2*m2*sin(q2) - L1_cm*L2_cm^2*g*m1*m2*cos(q1) + L1^2*L2_cm^2*der_q1*der_q2*m2^2*sin(2*q2) + 2*L1*L2_cm^3*der_q1*der_q2*m2^2*sin(q2) + L1*L1_cm^2*L2_cm*der_q1^2*m1*m2*sin(q2) + 2*I2*L1*L2_cm*der_q1*der_q2*m2*sin(q2) - L1*L1_cm*L2_cm*g*m1*m2*cos(q1)*cos(q2))/(I1*I2 + L1^2*L2_cm^2*m2^2 + I2*L1^2*m2 + I2*L1_cm^2*m1 + I1*L2_cm^2*m2 + L1_cm^2*L2_cm^2*m1*m2 - L1^2*L2_cm^2*m2^2*cos(q2)^2)
