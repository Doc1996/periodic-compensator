function der_depend_vars = diff_function(t, depend_vars)

global K_P K_D K_I K_nD alpha Q_k gama N_osc bool_for_omega_est bool_for_PInD_control bool_for_RC_control;
global q1_ref_bias q2_ref_bias q1_ref_phase q2_ref_phase q1_ref_amp q2_ref_amp q1_ref_omega q2_ref_omega;
global N_harms_of_ref_signal T_sim bool_for_25 bool_for_50 bool_for_75 bool_for_100;


% constants of the robot model

[L1, L2, L1_cm, L2_cm, m1, m2, I1, I2, g] = deal(0.35, 0.3, 0.15, 0.12, 4.5, 3.5, 0.15, 0.08, 9.81);


% reference signals

q1_ref = q1_ref_bias;
q2_ref = q2_ref_bias;

der_q1_ref = 0;
der_q2_ref = 0;

for i = 1:N_harms_of_ref_signal
    q1_ref = q1_ref + q1_ref_amp(i) * sin(q1_ref_omega(i) * t + q1_ref_phase);
    q2_ref = q2_ref + q2_ref_amp(i) * sin(q2_ref_omega(i) * t + q2_ref_phase);
    
    der_q1_ref = der_q1_ref + q1_ref_amp(i) * q1_ref_omega(i) * cos(q1_ref_omega(i) * t + q1_ref_phase);
    der_q2_ref = der_q2_ref + q2_ref_amp(i) * q2_ref_omega(i) * cos(q2_ref_omega(i) * t + q2_ref_phase);
end


% variables of the system

q1 = depend_vars(1);
q2 = depend_vars(2);
der_q1 = depend_vars(3);
der_q2 = depend_vars(4);

q1_err = q1 - q1_ref;
q2_err = q2 - q2_ref;

der_q1_err = der_q1 - der_q1_ref;
der_q2_err = der_q2 - der_q2_ref;

first_harm_omega_1_est = depend_vars(7);
first_harm_omega_2_est = depend_vars(8);


% variables of the robot controller

der_depend_vars = zeros(8 + 4 * N_osc + 2, 1);

z0_1 = depend_vars(8 + 1);
z0_2 = depend_vars(8 + 2);

der_depend_vars(8 + 1) = K_I(1) * (der_q1_err + alpha(1) * q1_err);
der_depend_vars(8 + 2) = K_I(2) * (der_q2_err + alpha(2) * q2_err);

sum_of_products_of_Q_k_and_z2_1 = 0;
sum_of_products_of_Q_and_z2_2 = 0;

factor_for_first_harm_omega = 1;

for i = (8 + 2 + 1):2:(8 + 2 * N_osc + 2)
    z1_1 = depend_vars(i);
    z1_2 = depend_vars(i + 1);
    
    z2_1 = depend_vars(2 * N_osc + i);
    z2_2 = depend_vars(2 * N_osc + i + 1);
    
    sum_of_products_of_Q_k_and_z2_1 = sum_of_products_of_Q_k_and_z2_1 + Q_k(1) * z2_1;
    sum_of_products_of_Q_and_z2_2 = sum_of_products_of_Q_and_z2_2 + Q_k(2) * z2_2;
    
    der_depend_vars(i) = z2_1;
    der_depend_vars(i + 1) = z2_2;
    
    der_depend_vars(2 * N_osc + i) = - first_harm_omega_1_est^2 * factor_for_first_harm_omega^2 * z1_1 + ...
        Q_k(1) * (der_q1_err + alpha(1) * q1_err);
    der_depend_vars(2 * N_osc + i + 1) = - first_harm_omega_2_est^2 * factor_for_first_harm_omega^2 * z1_2 + ...
        Q_k(2) * (der_q2_err + alpha(2) * q2_err);
    
    factor_for_first_harm_omega = factor_for_first_harm_omega + 1;
end

if bool_for_PInD_control == true && bool_for_RC_control == true
	u1 = - K_P(1) * q1_err - K_D(1) * der_q1_err - K_nD(1) * abs(q1_err) * der_q1_err - K_I(1) * z0_1 - ...
        sum_of_products_of_Q_k_and_z2_1;
	u2 = - K_P(2) * q2_err - K_D(2) * der_q2_err - K_nD(2) * abs(q2_err) * der_q2_err - K_I(2) * z0_2 - ...
        sum_of_products_of_Q_and_z2_2;
elseif bool_for_PInD_control == true && bool_for_RC_control == false
    u1 = - K_P(1) * q1_err - K_D(1) * der_q1_err - K_nD(1) * abs(q1_err) * der_q1_err - K_I(1) * z0_1;
	u2 = - K_P(2) * q2_err - K_D(2) * der_q2_err - K_nD(2) * abs(q2_err) * der_q2_err - K_I(2) * z0_2;
else
    u1 = 0;
    u2 = 0;
end


% variables of the robot model

M11 = m2*L1^2 + 2*m2*cos(q2)*L1*L2_cm + m1*L1_cm^2 + m2*L2_cm^2 + I1 + I2;
M12 = m2*L2_cm^2 + L1*m2*cos(q2)*L2_cm + I2;
M21 = m2*L2_cm^2 + L1*m2*cos(q2)*L2_cm + I2;
M22 = m2*L2_cm^2 + I2;

C_with_der_q_1 = -L1*L2_cm*der_q2*m2*sin(q2)*(2*der_q1 + der_q2);
C_with_der_q_2 = L1*L2_cm*der_q1^2*m2*sin(q2);

% trigonometric functions at G_1 and G_2 depend on the reference pose of the robot
G_1 = g*(L2_cm*m2*cos(q1 + q2) + L1*m2*cos(q1) + L1_cm*m1*cos(q1));
G_2 = L2_cm*g*m2*cos(q1 + q2);

M_det = M11 * M22 - M12 * M21;

M11_inv = M22 / M_det;
M12_inv = - M12 / M_det;
M21_inv = - M21 /M_det;
M22_inv = M11 / M_det;

M_inv_with_C_1 = M11_inv * C_with_der_q_1 + M12_inv * C_with_der_q_2;
M_inv_with_C_2 = M21_inv * C_with_der_q_1 + M22_inv * C_with_der_q_2;

M_inv_with_G_1 = M11_inv * G_1 + M12_inv * G_2;
M_inv_with_G_2 = M21_inv * G_1 + M22_inv * G_2;

M_inv_with_u_1 = M11_inv * u1 + M12_inv * u2;
M_inv_with_u_2 = M21_inv * u1 + M22_inv * u2;

der2_q1 = - M_inv_with_C_1 - M_inv_with_G_1 + M_inv_with_u_1;
der2_q2 = - M_inv_with_C_2 - M_inv_with_G_2 + M_inv_with_u_2;


% variables of the circular frequency estimator

if bool_for_PInD_control == true && bool_for_RC_control == true && bool_for_omega_est == true && N_osc > 0
    z1_1_for_first_harm = depend_vars(8 + 2 + 1);
    z1_2_for_first_harm = depend_vars(8 + 2 + 2);
    
    % circular frequency estimation law
    
    if first_harm_omega_1_est > 0
        der_first_harm_omega_1_est = - gama(1) * z1_1_for_first_harm * (der_q1_err + alpha(1) * q1_err);
    else
        der_first_harm_omega_1_est = 0;
    end
    if first_harm_omega_2_est > 0
        der_first_harm_omega_2_est = - gama(2) * z1_2_for_first_harm * (der_q2_err + alpha(2) * q2_err);
    else
        der_first_harm_omega_2_est = 0;
    end
else
    der_first_harm_omega_1_est = 0;
    der_first_harm_omega_2_est = 0;
end

der_depend_vars(1) = der_q1;
der_depend_vars(2) = der_q2;
der_depend_vars(3) = der2_q1;
der_depend_vars(4) = der2_q2;
der_depend_vars(5) = u1;
der_depend_vars(6) = u2;
der_depend_vars(7) = der_first_harm_omega_1_est;
der_depend_vars(8) = der_first_harm_omega_2_est;


% printing data about the simulation progress

if (t >= T_sim / 4) && (bool_for_25 == false); fprintf('    25%%'); bool_for_25 = true; end
if (t >= T_sim / 2) && (bool_for_50 == false); fprintf('    50%%'); bool_for_50 = true; end
if (t >= T_sim *3 / 4) && (bool_for_75 == false); fprintf('    75%%'); bool_for_75 = true; end
if (t == T_sim) && (bool_for_100 == false); fprintf('    100%%\n'); bool_for_100 = true; end

end