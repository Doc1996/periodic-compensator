function [freq, one_sided_trans_signal] = calc_fft_of_resulting_signal(time, signal, q1_ref_first_harm_omega)

global N_harms_of_ref_signal;

freq_max_needed = floor(q1_ref_first_harm_omega / 2 / pi) * N_harms_of_ref_signal * 4;

% the precision of the spectrum increases with a longer sampling time, which is one of the properties of the Fourier transform
T_sim = time(end);
resoulution_of_sampled_time = 1 / freq_max_needed;
sampled_time = 0:resoulution_of_sampled_time:T_sim;
num_of_samples = length(sampled_time);

% due to sampling, a continuous spectrum appears instead of a discrete one
sampled_signal = interp1(time, signal, sampled_time);

two_sided_trans_signal = fft(sampled_signal);
two_sided_trans_signal = abs(two_sided_trans_signal / num_of_samples);

one_sided_trans_signal = 2 * (two_sided_trans_signal(1:(floor(num_of_samples / 2 + 1)))).';

freq = (0:(floor(num_of_samples / 2))).' / T_sim;

end