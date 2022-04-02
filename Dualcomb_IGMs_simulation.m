% This code demonstrates the simulation method for dual-comb IGMs primarily based on the model in [1].
% This algorithm can be used for the parameter optimization of dual-comb interferometer and the verification of error correction algorithm.
% References
% [1] H. Yu, K. Ni, Q. Zhou, X. Li, X. Wang, G. Wu, "Digital error correction of dual-comb interferometer without external optical referencing information," Optics Express, 2019, 27(20): 29425-29438.
% [2] H. Yu, Q. Zhou, X. Li, X. Wang, K. Ni, “Mode-resolved dual-comb spectroscopy using error correction based on single optical intermedium,” Optics Express, 2021, 29(4): 6271-6281.
% Copyright, Haoyang Yu(yu.haoyang@csu.edu.cn) & Kai Ni(ni.kai@sz.tsinghua.edu.cn)
clear all;
clc;
close all;
%% Basic parameters
c = 3e8;							
lambda_c = 1560e-9;             	 
lambda_bw = 1.8e-9;                  
fr = [57.199e6, 57.200e6];          
fceo = [2e6, 0];	                 
Fs = fr(2);                     
tau0_LO = 0;
tau0_mea = 8e-9;
phi0_LO = 0;phi0_mea = 0;
I_IGMs = 1;
t_count = 1;                
act_n = Fs / abs(fr(1) - fr(2));                   
t_act_n = act_n / Fs;
%% Derived parameters
dfr = abs(fr(1) - fr(2));
T_RF = 1 / (fr(1) - fr(2));
miu_c = c / lambda_c;
miu_mea_c = round((miu_c - fceo(1)) / fr(1)) * fr(1) + fceo(1);
miu_LO_c = round((miu_mea_c - fceo(2)) / fr(2)) * fr(2) + fceo(2);
fc_RF = miu_mea_c - miu_LO_c;
f_BW_opt = c * lambda_bw / (lambda_c ^ 2);
f_BW_RF = f_BW_opt / fr(1) * dfr;
tau_BW_RF = 4 * log(2) / pi / f_BW_RF;            
tau0_RF = (tau0_mea * fr(1) - tau0_LO * fr(2)) * T_RF;
phic_RF = phi0_mea + 2 * pi * tau0_mea * miu_mea_c - phi0_LO - 2 * pi * tau0_LO * miu_LO_c;
envelope = @(x)  exp(- ((pi * f_BW_RF * x) .^ 2) / 4 / log(2));
if (abs(fc_RF) - f_BW_RF * 3 / 2 <= 0) | (abs(fc_RF) + f_BW_RF * 3 / 2 >= min(Fs,fr) / 2)
    warning('Spectrum Aliasing') 
end
n_count = round(t_count * Fs);
t_lab = (1:n_count) / Fs;
N1 = round(tau0_RF / T_RF);
N2 = round((t_lab(end) + tau0_RF) / T_RF);
N = N1:(N2-N1)/abs(N2-N1):N2;
peak_ideal = round((N * T_RF - tau0_RF) * Fs);
%% Noise simulation
n_peak = ceil(n_count / Fs * dfr);
dfr_sim = zeros(1,n_peak);dfr_sim_max = 0.02;
dfc_sim = zeros(1,n_peak);dfc_sim_max = 0.5e6;
std1 = 0.01;
std2 = 0.1e6;
random_mat1 = randn(1,n_peak - 1);
random_mat2 = randn(1,n_peak - 1);
for i = 2:n_peak
  dfr_sim(i)=median([-dfr_sim_max, dfr_sim(i-1) + random_mat1(i-1)*std1, dfr_sim_max]);
  dfc_sim(i)=median([-dfc_sim_max, dfc_sim(i-1) + random_mat2(i-1)*std2, dfc_sim_max]);
end
n_frame = Fs/dfr;
dfr_sim = repmat(dfr_sim,n_frame,1);dfc_sim = repmat(dfc_sim,n_frame,1);
dfr_sim = dfr_sim(:)';dfc_sim = dfc_sim(:)';
dfr_sim = dfr_sim(1:n_count);dfc_sim = dfc_sim(1:n_count);
f1 = [0:n_count - 1] * Fs / n_count;
f2 = fftshift(f1); zi = find(f2==0); f2(1:zi-1) = f2(1:zi-1) - Fs;
Ffr = fftshift(fft(dfr_sim)); Ffr(abs(f2) > dfr / 2) = 0;
Ffr(abs(f2) <= dfr / 2) = Ffr(abs(f2) <= dfr / 2) .* 0.5 .* (cos(f2(abs(f2) <= dfr / 2) / (dfr / 2) * pi) + 1);
Ffc = fftshift(fft(dfc_sim)); Ffc(abs(f2) > dfr / 2) = 0;
Ffc(abs(f2) <= dfr / 2) = Ffc(abs(f2) <= dfr / 2) .* 0.5 .* (cos(f2(abs(f2) <= dfr / 2) / (dfr / 2) * pi) + 1);
delta_fr_RF = real(ifft(ifftshift(Ffr)));
delta_fc_RF = real(ifft(ifftshift(Ffc)));
delta_tau0_RF = cumsum(2 * pi * delta_fr_RF / Fs) / 2 / pi / dfr;
delta_phic_RF = cumsum(2 * pi * delta_fc_RF / Fs);
%% IGMs synthesis
IGMs = zeros(1,n_count);
IGMs_ideal = zeros(1,n_count);
for i = 1:numel(peak_ideal)
    flag = peak_ideal(i) - floor(act_n / 2):peak_ideal(i) + floor(act_n / 2 - 0.5);
    flag(flag < 1 | flag > n_count) = [];
    N_temp = N(i);
    temp = t_lab(flag) - N_temp *  T_RF + tau0_RF + delta_tau0_RF(flag);
    IGMs(flag) = IGMs(flag) + I_IGMs * envelope(temp) .* cos(2 * pi * fc_RF * t_lab(flag) + phic_RF + delta_phic_RF(flag)); 
    IGMs_ideal(flag) = IGMs_ideal(flag) + I_IGMs * envelope(t_lab(flag) - N_temp *  T_RF + tau0_RF) .* cos(2 * pi * fc_RF * t_lab(flag) + phic_RF);    
end
std_noise = 0.01;                              
E_noise = std_noise * randn (1,n_count);
IGMs = IGMs + E_noise;
save('SimuDCI', 'IGMs', 'IGMs_ideal','delta_tau0_RF','delta_fc_RF','delta_phic_RF','fc_RF');
%% Display
flag = peak_ideal(2) - floor(act_n / 2):peak_ideal(2) + floor(act_n / 2 - 0.5);
figure;
plot(t_lab(flag) * 1000, IGMs(flag),'k')
title(strcat('Interferogram'),'fontsize',20)
xlabel('Time/ms','fontsize',20)
ylabel('Signal voltage/a.u.','fontsize',20);
