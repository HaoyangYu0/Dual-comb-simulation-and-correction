% This code demonstrates the self-referencing error correction method for dual-comb IGMs referred to [1,2].
% This algorithm can be used for the frame-to-frame alignment and coherent averaging of dual-comb IGMs in the time domain.
% The explaination of the algorithm can be found in Ref [1]. The requirement of the algorithm can be found in Ref [2,3].
% References
% [1] H. Yu, K. Ni, Q. Zhou, X. Li, X. Wang, G. Wu, "Digital error correction of dual-comb interferometer without external optical referencing information," Optics Express, 2019, 27(20): 29425-29438.
% [2] H. Yu, Q. Zhou, X. Li, X. Wang, K. Ni, “Mode-resolved dual-comb spectroscopy using error correction based on single optical intermedium,” Optics Express, 2021, 29(4): 6271-6281.
% [3] H. Yu, K. Ni, Q. Zhou, X. Li, X. Wang, G. Wu, “Continuous Jitter Reconstruction of Dual-comb interferometer for self-referencing error correction,” in 2020 Conference on Lasers and Electro-Optics Pacific Rim, CLEO-PR 2020. 
% Copyright, Haoyang Yu(yu.haoyang@csu.edu.cn) & Kai Ni(ni.kai@sz.tsinghua.edu.cn)
clear all;
close all;
clc;
load('SimuDCI', 'IGMs', 'IGMs_ideal','delta_tau0_RF','delta_fc_RF','delta_phic_RF','fc_RF');
%% Parameters
fr = [57.199e6,57.2e6];dfr = abs(diff(fr));
Fs = fr(2);
n_dfr = Fs / dfr;
t_count = 1;
act_n = 100;nn = 4096;left = 200;right = 200;
n_count = round(t_count * Fs / n_dfr) * n_dfr;
IGMs = IGMs(1:n_count);
%% Findpeaks
Env_I = abs(hilbert(IGMs)); 
[~, peak_I,~,~] = findpeaks(Env_I,'MinPeakDistance',1000,'MinPeakHeight',0.1 * max(Env_I));
n_peak = length(peak_I);
%% Error extraction
f = (0:nn/2-1) * Fs / nn;
cv_init = IGMs(peak_I(1) - floor(act_n / 2):peak_I(1) + floor(act_n / 2 - 0.5));
spec_cv_init = fft(cv_init,nn);
mag_cv_init = abs(spec_cv_init(1:nn/2));
idx_cv_init = mag_cv_init > 0.5 * max(mag_cv_init);
fc_cv_init = sum(f(idx_cv_init) .* mag_cv_init(idx_cv_init)) / sum(mag_cv_init(idx_cv_init));
phi_cv_init = angle(spec_cv_init(1:nn/2));
TJ = zeros(1,n_peak);FJ = zeros(1,n_peak);PJ = zeros(1,n_peak);
% fc_IGMs = round(fc_cv_init / dfr) * dfr;FJ(1) = fc_cv_init - fc_IGMs;
for i = 1 : n_peak
    flag = peak_I(i) - floor(act_n / 2):peak_I(i) + floor(act_n / 2 - 0.5);
    cv = IGMs(flag);
    spec_cv = fft(cv,nn);
    mag_cv = abs(spec_cv(1:nn/2));
    phi_cv = angle(spec_cv(1:nn/2));
    idx_cv = mag_cv > 0.5 * max(mag_cv);
    fc_cv = sum(f(idx_cv) .* mag_cv(idx_cv)) / sum(mag_cv(idx_cv));
    dphi = unwrap(phi_cv - phi_cv_init);
    temp = polyfit(2 * pi * f(idx_cv & idx_cv_init),dphi(idx_cv & idx_cv_init),1);
    TJ(i) = temp(1);
    PJ(i) = temp(2);
    FJ(i) = fc_cv - fc_cv_init;
end
true_peak = peak_I - TJ * Fs;
TJ = (0:n_peak - 1) / dfr - true_peak / Fs;
TJ = TJ - TJ(1);
PJ = PJ - PJ(1);
%% Error comparison
t_lab = true_peak / Fs;
delta_phi0_RF = delta_phic_RF -  2 * pi * fc_RF * delta_tau0_RF;
TJ_ideal = delta_tau0_RF(round(true_peak));
TJ_ideal = TJ_ideal - TJ_ideal(1);
PJ_ideal = delta_phi0_RF(round(true_peak));
PJ_ideal = PJ_ideal - PJ_ideal(1);
FJ_ideal = delta_fc_RF(round(true_peak));
FJ_ideal = FJ_ideal - FJ_ideal(1);
fontsize1 = 14;
set(0,'DefaultAxesFontSize',fontsize1);
set(0,'DefaultTextFontSize',fontsize1);
figure;set(gcf,'position',[54 85 1400 600]);
subplot 231;plot(t_lab,TJ_ideal*1e9,'r',t_lab,TJ*1e9,'b');
xlabel('Time(s)');ylabel('Time jitter(ns)'); title('(a)','FontName','Times New Roman');
subplot 234;plot(t_lab,(TJ-TJ_ideal-mean(TJ-TJ_ideal))*1e9,'g');
xlabel('Time(s)');ylabel('Time jitter residuals(ns)'); title('(d)','FontName','Times New Roman');
subplot 232;plot(t_lab,FJ_ideal/1e6,'r',t_lab,FJ/1e6,'b');
xlabel('Time(s)');ylabel('Frequency jitter(MHz)'); title('(b)','FontName','Times New Roman');
subplot 235;plot(t_lab,(FJ-FJ_ideal-mean(FJ-FJ_ideal))/1e6,'g');
xlabel('Time(s)');ylabel('Frequency jitter residuals(MHz)'); title('(e)','FontName','Times New Roman');
subplot 233;plot(t_lab,wrapToPi(PJ_ideal),'r',t_lab,wrapToPi(PJ),'b');
xlabel('Time(s)');ylabel('Phase jitter(rad)'); title('(c)','FontName','Times New Roman');
yticks([-pi,-pi/2,0,pi/2,pi]);yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
subplot 236;plot(t_lab,wrapToPi(PJ-PJ_ideal)-mean(wrapToPi(PJ-PJ_ideal)),'g');
xlabel('Time(s)');ylabel('Phase jitter residuals(rad)'); title('(f)','FontName','Times New Roman');
%% Error correction
flag = peak_I(1) - floor(act_n / 2):peak_I(1) + floor(act_n / 2 - 0.5);
flag_long = peak_I(1) - floor(150 / 2):peak_I(1) + floor(300 / 2 - 0.5);
t_N = true_peak / Fs - (0:n_peak - 1) / dfr;
t = (1:n_dfr)' / Fs;
IGMs_reshape = reshape(IGMs,n_dfr,n_peak);
MAG = [];
for i = 1:n_peak
IGMs_reshape1(:,i) = interp1(t + TJ(i),IGMs_reshape(:,i),t,'spline',0);
IGMs_reshape1(:,i) = hilbert(IGMs_reshape1(:,i));
[~,n_temp] = max(abs(IGMs_reshape(:,i)));
flag_temp = n_temp - floor(act_n / 2):n_temp + floor(act_n / 2 - 0.5);
spec = fft(IGMs_reshape(flag_temp,i),nn);
mag = abs(spec(1:nn/2))/max(abs(spec(1:nn/2)));
MAG = [MAG mag];
end
SPEC1 = fft(real(IGMs_reshape1(flag,:)),nn);
MAG1 = abs(SPEC1(1:nn/2,:)) ./ max(abs(SPEC1(1:nn/2,:)));
phase1 = 2 * pi * repmat(FJ,n_dfr,1) .* (repmat(t,1,n_peak) - repmat(t_N(1),n_dfr,n_peak));
phase2 = repmat(PJ,n_dfr,1);
IGMs_reshape2 = IGMs_reshape1 .* exp(-1i * phase1);
SPEC2 = fft(real(IGMs_reshape2(flag,:)),nn);
MAG2 = abs(SPEC2(1:nn/2,:))./ max(abs(SPEC2(1:nn/2,:)));
IGMs_reshape3 = IGMs_reshape2 .* exp(-1i * phase2);
SPEC3 = fft(real(IGMs_reshape3(flag,:)),nn);
MAG3 = abs(SPEC3(1:nn/2,:))./ max(abs(SPEC3(1:nn/2,:)));
figure;set(gcf,'position',[54 85 1600 600]);
t_temp_long = t(flag_long) - t(flag_long(1));
t_temp = t(flag) - t(flag_long(1));
subplot 241;plot(t_temp_long * 1e6,real(IGMs_reshape(flag_long,:)),'k','Linewidth',1);
ylim([-1.1 1.1]);xlabel('Time(μs)');ylabel('Amp(a.u.)'); title('(a)','FontName','Times New Roman');
subplot 242;plot(t_temp * 1e6,real(IGMs_reshape1(flag,:)),'k','Linewidth',1);
xlim([min(t_temp)*1e6 max(t_temp)*1e6]);ylim([-1.1 1.1]);xlabel('Time(μs)');ylabel('Amp(a.u.)'); title('(b)','FontName','Times New Roman');
subplot 243;plot(t_temp * 1e6,real(IGMs_reshape2(flag,:)),'k','Linewidth',1);
xlim([min(t_temp)*1e6 max(t_temp)*1e6]);ylim([-1.1 1.1]);xlabel('Time(μs)');ylabel('Amp(a.u.)'); title('(c)','FontName','Times New Roman');
subplot 244;plot(t_temp * 1e6,real(IGMs_reshape3(flag,:)),'k','Linewidth',1);
xlim([min(t_temp)*1e6 max(t_temp)*1e6]);ylim([-1 1]);xlabel('Time(μs)');ylabel('Amp(a.u.)'); title('(d)','FontName','Times New Roman');
subplot 245;plot(f/1e6,MAG,'k','Linewidth',1);
xlim([min(f)/1e6 max(f)/1e6]);ylim([0 1.1]);xlabel('RF frequency(MHz)');ylabel('Amp(a.u.)'); title('(e)','FontName','Times New Roman');
subplot 246;plot(f/1e6,MAG1,'k','Linewidth',1);
xlim([min(f)/1e6 max(f)/1e6]);ylim([0 1.1]);xlabel('RF frequency(MHz)');ylabel('Amp(a.u.)'); title('(f)','FontName','Times New Roman');
subplot 247;plot(f/1e6,MAG2,'k','Linewidth',1);
xlim([min(f)/1e6 max(f)/1e6]);ylim([0 1.1]);xlabel('RF frequency(MHz)');ylabel('Amp(a.u.)'); title('(g)','FontName','Times New Roman');
subplot 248;plot(f/1e6,MAG3,'k','Linewidth',1);
xlim([min(f)/1e6 max(f)/1e6]);ylim([0 1.1]);xlabel('RF frequency(MHz)');ylabel('Amp(a.u.)'); title('(h)','FontName','Times New Roman');

