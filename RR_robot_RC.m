%% INITIAL CONDITIONS AND PARAMETERS OF THE ROBOT, REFERENCES AND CONTROLLERS FOR SIMULATION

clear all;
clc;
format compact;

global K_P K_D K_I K_nD alpha Q_k gama N_osc bool_for_omega_est bool_for_PInD_control bool_for_RC_control;
global q1_ref_bias q2_ref_bias q1_ref_phase q2_ref_phase q1_ref_amp q2_ref_amp q1_ref_omega q2_ref_omega;
global N_harms_of_ref_signal T_sim bool_for_25 bool_for_50 bool_for_75 bool_for_100;


% initial conditions of the robot and controller
q1_0 = 0.5;  %% ADJUSTABLE
q2_0 = 0.5;  %% ADJUSTABLE
der_q1_0 = 0;  %% FIXED
der_q2_0 = 0;  %% FIXED
u1_0 = 0;  %% FIXED
u2_0 = 0;  %% FIXED

% parameters of the reference
N_harms_of_ref_signal = 1;  %% ADJUSTABLE  % the number of reference harmonics
q1_ref_first_harm_omega = 6;  %% ADJUSTABLE
q2_ref_first_harm_omega = 6;  %% ADJUSTABLE
limits_for_bias = [0, 0];  %% ADJUSTABLE
limits_for_phase_multiplied_by_45_deg = [0, 0];  %% ADJUSTABLE
limits_for_amp = [1, 1];  %% ADJUSTABLE

% simulation time and step for plotting response
T_sim = 30;  %% ADJUSTABLE
T_step = 0.05;  %% ADJUSTABLE

% tolerance of the differential equation solver
bool_for_finer_solver_tol = true;  %% ADJUSTABLE  % true if a finer tolerance is requested, otherwise it is false

% use of the PInD part of controller
bool_for_PInD_control = true;  %% ADJUSTABLE  % true if a PInD part of controller is used, otherwise it is false

% use of the part of controller for periodic compensation
bool_for_RC_control = true;  %% ADJUSTABLE  % true if a part of controller for periodic comp. is used, otherwise it is false

% estimation of the circular frequency of reference
bool_for_omega_est = true;  %% ADJUSTABLE  % true if the base omega is estimated, otherwise it is false


% conditions for adjusting the controller parameters:
% k_g = (L1 + L2) * (m1 + m2) * g = 15
% M_max = (L1 + L2)^2 * (m1 + m2) = 0.75
% K_P > k_g = 15  -->  empiric: K_P = 50
% K_D such that an aperiodic response is obtained  -->  K_D > 10  -->  empiric: K_D = 80
% K_I < sqrt((K_P - k_g) * K_D / M_max) = 60  -->  K_I < 60  -->  empiric: K_I = 8
% K_nD < K_I / (K_P - k_g) * (M_max + k_g) = 3.6
% alfa > K_I / (K_P - k_g) = 0.25  -->  empiric: alfa = 0.5
% checkup: K_D > alfa * M_max  -->  80 > 0.4  --> satisfies

% parameters of the controller, set so that each cell element defines other parameters
K_P_cell = {[50, 50]};  %% FIXED
K_D_cell = {[80, 80]};  %% FIXED
K_I_cell = {[8, 8]};  %% FIXED
K_nD_cell = {[2, 2]};  %% FIXED
alpha_cell = {[0.5, 0.5]};  %% FIXED
Q_k_cell = {[10, 10]};  %% FIXED
N_osc_cell = {12};  %% ADJUSTABLE  % number of oscillators (should be greater than the number of reference harmonics)
first_harm_omega_est_0_cell = {[8, 8]};  %% ADJUSTABLE
gama_cell = {[10, 180]}; %% ADJUSTABLE  % gamma_2 should be about 18 times larger than gamma_1


% known properties of circular frequency estimation:
% - reference values of z1 and z2 are equal to 0
% - z1 is an ideally symmetrical sinusoid of constant amplitude mostly less than 1 (it differs for different reference
% signals), with a period equal to the period of the reference signal
% - z2 is ideally a constant of value 0
% - due to the mentioned reference values, the error of z1, i.e. z2 is equal to z1, i.e. z2, so it will continue to be
% used that way


%% SIMULATION AND SYSTEM RESPONSES

t = (0:T_step:T_sim).';

q_ref_first_harm_omega = [q1_ref_first_harm_omega, q2_ref_first_harm_omega];

q1_ref_bias = randi(limits_for_bias);
q2_ref_bias = randi(limits_for_bias);
q1_ref_phase = pi/4 * randi(limits_for_phase_multiplied_by_45_deg);
q2_ref_phase = pi/4 * randi(limits_for_phase_multiplied_by_45_deg);

q1_ref_amp = zeros(1, N_harms_of_ref_signal);
q2_ref_amp = zeros(1, N_harms_of_ref_signal);
q1_ref_omega = zeros(1, N_harms_of_ref_signal);
q2_ref_omega = zeros(1, N_harms_of_ref_signal);

for i = 1:N_harms_of_ref_signal
    q1_ref_amp(i) = randi(limits_for_amp);
    q2_ref_amp(i) = randi(limits_for_amp);
    q1_ref_omega(i) = q1_ref_first_harm_omega * i;
    q2_ref_omega(i) = q2_ref_first_harm_omega * i;
end

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

q_ref = [q1_ref, q2_ref];
der_q_ref = [der_q1_ref, der_q2_ref];


[freq_ref_1, one_sided_trans_signal_ref_1] = calc_fft_of_resulting_signal(t, q_ref(:, 1), q1_ref_first_harm_omega);
[freq_ref_2, one_sided_trans_signal_ref_2] = calc_fft_of_resulting_signal(t, q_ref(:, 2), q2_ref_first_harm_omega);


figure('Name', 'Responses of the system', 'NumberTitle', 'off');
set(gcf, 'Color', 'w');

legend_str_cell = cell(1, length(K_P_cell) * length(K_D_cell) * length(K_I_cell) * length(alpha_cell) * ...
    length(K_nD_cell) * length(Q_k_cell) * length(N_osc_cell) * length(first_harm_omega_est_0_cell) * ...
    length(gama_cell) + 1);
legend_str_cell{1} = 'reference';

cell_index = 1;

for i = 1:length(K_P_cell), for j = 1:length(K_D_cell), for k = 1:length(K_I_cell), for m = 1:length(alpha_cell), ...
for n = 1:length(K_nD_cell), for p = 1:length(Q_k_cell), for r = 1:length(N_osc_cell), ...
for s = 1:length(first_harm_omega_est_0_cell), for v = 1:length(gama_cell)
    
    K_P = K_P_cell{i};
    K_D = K_D_cell{j};
    K_I = K_I_cell{k};
    alpha = alpha_cell{m};
    K_nD = K_nD_cell{n};
    Q_k = Q_k_cell{p};
    N_osc = N_osc_cell{r};
    first_harm_omega_est_0 = first_harm_omega_est_0_cell{s};
    gama = gama_cell{v};
    
    if bool_for_PInD_control == true && bool_for_RC_control == true && bool_for_omega_est == true && N_osc > 0
        legend_str_cell{cell_index + 1} = sprintf(['controlled with PID+RC and omega estim., N_h_a_r_m = %d, ', ...
            'K_P = %d, K_D = %d, K_I = %d, alpha = %0.1f, K_n_D = %d, Q_k = %d, N_o_s_c = %d, init. omega = ', ...
            '%0.1f / %0.1f, gama = %d / %d'], N_harms_of_ref_signal, K_P(1), K_D(1), K_I(1), alpha(1), K_nD(1), ...
            Q_k(1), N_osc, first_harm_omega_est_0(1), first_harm_omega_est_0(2), gama(1), gama(2));
    elseif bool_for_PInD_control == true && bool_for_RC_control == true && bool_for_omega_est == false
        legend_str_cell{cell_index + 1} = sprintf(['controlled with PInD+RC and no omega estim., N_h_a_r_m = %d, ', ...
            'K_P = %d, K_D = %d, K_I = %d, alpha = %0.1f, K_n_D = %d, Q_k = %d, N_o_s_c = %d, init. omega = ', ...
            '%0.1f / %0.1f'], N_harms_of_ref_signal, K_P(1), K_D(1), K_I(1), alpha(1), K_nD(1), Q_k(1), N_osc, ...
            first_harm_omega_est_0(1), first_harm_omega_est_0(2));
    elseif bool_for_PInD_control == true && bool_for_RC_control == false && bool_for_omega_est == false
        legend_str_cell{cell_index + 1} = sprintf(['controlled with PInD and no omega estim., N_h_a_r_m = %d, K_P = ', ...
            '%d, K_D = %d, K_I = %d, alpha = %0.1f, K_n_D = %d, init. omega = %0.1f / %0.1f'], N_harms_of_ref_signal, ...
            K_P(1), K_D(1), K_I(1), alpha(1), K_nD(1), first_harm_omega_est_0(1), first_harm_omega_est_0(2));
    else
        legend_str_cell{cell_index + 1} = sprintf(['uncontrolled and without omega estim., N_h_a_r_m = %d, init. ', ...
            'omega = %0.1f / %0.1f'], N_harms_of_ref_signal, first_harm_omega_est_0(1), first_harm_omega_est_0(2));
    end
    
    if bool_for_finer_solver_tol == true
        options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    else
        options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
    end
    
    bool_for_25 = false; bool_for_50 = false; bool_for_75 = false; bool_for_100 = false;
    z_all_0 = zeros(4 * N_osc + 2, 1);
    % first_harm_omega_sqr_est_0 = first_harm_omega_est_0.^2;
    % depend_vars_0 = [q1_0; q2_0; der_q1_0; der_q2_0; u1_0; u2_0; first_harm_omega_sqr_est_0(1); ...
    %     first_harm_omega_sqr_est_0(2); z_all_0];
    depend_vars_0 = [q1_0; q2_0; der_q1_0; der_q2_0; u1_0; u2_0; first_harm_omega_est_0(1); ...
        first_harm_omega_est_0(2); z_all_0];
    
    tic;
    
    % performing simulation of the system
    [t, depend_vars] = ode15s('diff_function', t, depend_vars_0, options);
    
    q = [depend_vars(:, 1), depend_vars(:, 2)];
	q_err = q - q_ref;
    der_q = [depend_vars(:, 3), depend_vars(:, 4)];
	der_q_err = der_q - der_q_ref;
    u = [diff(depend_vars(:, 5)) ./ diff(t), diff(depend_vars(:, 6)) ./ diff(t)];
    u = [u; u(end, :)];
    first_harm_omega_ref = ones(size(t)) * q_ref_first_harm_omega;
    % first_harm_omega_est = sqrt([depend_vars(:, 7), depend_vars(:, 8)]);
    first_harm_omega_est = [depend_vars(:, 7), depend_vars(:, 8)];
    first_harm_omega_err = first_harm_omega_est - first_harm_omega_ref;
    
    [freq_1, one_sided_trans_signal_1] = calc_fft_of_resulting_signal(t, q(:, 1), q1_ref_first_harm_omega);
    [freq_2, one_sided_trans_signal_2] = calc_fft_of_resulting_signal(t, q(:, 2), q2_ref_first_harm_omega);
    
    if cell_index == 1
        % plotting the system reference
        invisible_subplot_handle = plot_results_depending_on_gain(t, q_ref, q_err, der_q_ref, der_q_err, u, ...
            first_harm_omega_ref, first_harm_omega_err, freq_ref_1, one_sided_trans_signal_ref_1, freq_ref_2, ...
            one_sided_trans_signal_ref_2, true);
    end
    
    % plotting the system response
	plot_results_depending_on_gain(t, q, q_err, der_q, der_q_err, u, first_harm_omega_est, first_harm_omega_err, ...
        freq_1, one_sided_trans_signal_1, freq_2, one_sided_trans_signal_2, false);
    
    delta_time = toc;
    fprintf('plotted graph %d of %d in %0.1f seconds\n\n', cell_index, length(K_P_cell) * length(K_D_cell) * ...
        length(K_I_cell) * length(alpha_cell) * length(K_nD_cell) * length(Q_k_cell) * length(N_osc_cell) * ...
        length(first_harm_omega_est_0_cell) * length(gama_cell), delta_time);
    
    cell_index = cell_index + 1;
end, end, end, end, end, end, end, end, end

set(invisible_subplot_handle, 'Visible', 'off');
legend(invisible_subplot_handle, legend_str_cell, 'Location', 'northwestoutside');
