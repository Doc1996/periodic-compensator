function invisible_subplot_handle = plot_results_depending_on_gain(t, q_can_be_ref, q_err, der_q_can_be_ref, ...
    der_q_err, u, first_harm_omega_can_be_ref, first_harm_omega_err, freq_1_can_be_ref, one_sided_trans_signal_1_can_be_ref, ...
    freq_2_can_be_ref, one_sided_trans_signal_2_can_be_ref, bool_for_ref)

if bool_for_ref == true
    plot_results(1, t, q_can_be_ref(:, 1), 'joint 1 position', 't [s]', 'q_1(t)', false, true);
    plot_results(3, t, q_can_be_ref(:, 2), 'joint 2 position', 't [s]', 'q_2(t)', false, true);
    plot_results(5, t, der_q_can_be_ref(:, 1), 'joint 1 velocity', 't [s]', 'der q_1(t)', false, true);
    plot_results(7, t, der_q_can_be_ref(:, 2), 'joint 2 velocity', 't [s]', 'der q_2(t)', false, true);
    plot_results(11, t, first_harm_omega_can_be_ref(:, 1), 'primary \omega_1', 't [s]', '\omega_1(t)', false, true);
    plot_results(12, t, first_harm_omega_can_be_ref(:, 2), 'primary \omega_2', 't [s]', '\omega_2(t)', false, true);
    if length(freq_1_can_be_ref) > 1 && length(freq_2_can_be_ref) > 1
        plot_results(15, freq_1_can_be_ref, one_sided_trans_signal_1_can_be_ref, 'sampled spectrum f_1', ...
            'f [Hz]', 'amplitude 1', false, true);
        plot_results(16, freq_2_can_be_ref, one_sided_trans_signal_2_can_be_ref, 'sampled spectrum f_2', ...
            'f [Hz]', 'amplitude 2', false, true);
    end
    invisible_subplot_handle = plot_results(17, 1, NaN, 'invisible subplot for legend', '-', '-', false, true);
else
    plot_results(1, t, q_can_be_ref(:, 1), 'joint 1 position', 't [s]', 'q_1(t)', false, false);
    plot_results(2, t, abs(q_err(:, 1)), 'joint 1 position error', 't [s]', '| q_1_,_e_r_r(t) |', true, false);
    plot_results(3, t, q_can_be_ref(:, 2), 'joint 2 position', 't [s]', 'q_2(t)', false, false);
    plot_results(4, t, abs(q_err(:, 2)), 'joint 2 position error', 't [s]', '| q_2_,_e_r_r(t) |', true, false);
    plot_results(5, t, der_q_can_be_ref(:, 1), 'joint 1 velocity', 't [s]', 'der q_1(t)', false, false);
    plot_results(6, t, abs(der_q_err(:, 1)), 'joint 1 velocity error', 't [s]', '| der q_1_,_e_r_r(t) |', true, false);
    plot_results(7, t, der_q_can_be_ref(:, 2), 'joint 2 velocity', 't [s]', 'der q_2(t)', false, false);
    plot_results(8, t, abs(der_q_err(:, 2)), 'joint 2 velocity error', 't [s]', '| der q_2_,_e_r_r(t) |', true, false);
    plot_results(9, t, u(:, 1), 'control signal 1', 't [s]', 'u_1(t)', false, false);
    plot_results(10, t, u(:, 2), 'control signal 2', 't [s]', 'u_2(t)', false, false);
    plot_results(11, t, first_harm_omega_can_be_ref(:, 1), 'primary \omega_1', 't [s]', '\omega_1(t)', false, false);
    plot_results(12, t, first_harm_omega_can_be_ref(:, 2), 'primary \omega_2', 't [s]', '\omega_2(t)', false, false);
    plot_results(13, t, abs(first_harm_omega_err(:, 1)), 'primary \omega_1 error', 't [s]', ...
        '| \omega_1_,_e_r_r(t) |', true, false);
    plot_results(14, t, abs(first_harm_omega_err(:, 2)), 'primary \omega_2 error', 't [s]', ...
        '| \omega_2_,_e_r_r(t) |', true, false);
    if length(freq_1_can_be_ref) > 1 && length(freq_2_can_be_ref) > 1
        plot_results(15, freq_1_can_be_ref, one_sided_trans_signal_1_can_be_ref, 'sampled spectrum f_1', ...
            'f [Hz]', 'amplitude 1', false, false);
        plot_results(16, freq_2_can_be_ref, one_sided_trans_signal_2_can_be_ref, 'sampled spectrum f_2', ...
            'f [Hz]', 'amplitude 2', false, false);
    end
    invisible_subplot_handle = plot_results(17, 1, NaN, 'invisible subplot for legend', '-', '-', false, false);
end


function subplot_handle = plot_results(num_of_order, x_var, y_var, title_str, x_label_str, y_label_str, ...
        bool_for_y_log, bool_for_ref)

subplot_handle = subplot(5, 4, num_of_order);

if bool_for_ref == true
    plot(x_var, y_var, 'k--', 'LineWidth', 1);
else
    plot(x_var, y_var, 'LineWidth', 1);
end

hold on;

if bool_for_y_log == true
	set(gca,'YScale', 'log');
end

title(title_str); xlabel(x_label_str); ylabel(y_label_str);
xlim([0 x_var(end)]); grid on;

end

end