%==========================================================================
%           A controllable Impedance Spectroscopy Device
%
%         E. Sadeghi, M. H. Zand, M. Hamzeh, and S.M.M. Alavi++
% 
%
%++ corresponding author: 
%          Department of Electrical Engineering and
%          Shahid Beheshti University, Tehran, Iran.
%          e-mail: m_alavi@sbu.ac.ir
%          website: arg.sbu.ac.ir 
%          June 2019
%==========================================================================
% This Matlab file runs test #2 of the above manuscript 
%==========================================================================
% The computation of battery impedance spectra by using FFT, and
% system identification (sysid) method in cEIS device. 
%
% Two equivalent circuit models (ECMs) are estimated to fit the impedance spectra 
% obtained from FFT. 
%
% The first ECM model is a first-order Randles circuit as shown below:
%                             R1
%                    |-----/\/\/\/------|
%          Rinf      |                  |
%  -------/\/\/\-----|                  |----
%                    |      1/(C1s)     |
%                    |-------||---------|
%
%                        v      m1*s + m0
%                   Z = --- = --------------
%                        i        s + n0
%
% where Rinf=m1, R1=m0/n0-m1, and C1=1/(n0*R1). 
%
%
%
% The second identified ECM model is a second-order Randles circuit as shown below:
%
%                             R1                      R2
%                    |-----/\/\/\/------|    |-----/\/\/\/-----|
%          Rinf      |                  |    |                 |
%  -------/\/\/\-----|                  |----|                 |--------
%                    |      1/(C1s)     |    |      1/(C2s)    |
%                    |-------||---------|    |--------||--------
%
%
%                               v      m2*s^2 + m1*s + m0
%                          Z = --- = ----------------------
%                               i        s^2 + n1*s + n0
% where, 
% tau1=R1*C1; tau2=R2*C2; a1=1/tau1; a2=1/tau2; b1=1/C1; b2=1/C2;
% m2=Rinf; m1=Rinf*(a1+a2)+b1+b2; m0=Rinf*a1*a2+b1*a2+b2*a1;
% n1=a1+a2; n0=a1*a2;
% 
%

close all
clear all
clc

%% Load Data
load 18signal-from-1hz-upto-2khz
%load 21signal-from-1hz-upto-2khz
%% Setup Sample Time
L=20;  %put the signal time duration in seconds
Ts=L/length(current);
fs=1/Ts;
t = Ts:(1/fs):L;

%% coherency of input(current) and output data(voltage)
[Cxy,F] = mscohere(current,voltage,[],[],[],fs); % this will go up to half fs
F_good = find((Cxy > .998) & (F < 2000));  % find indices of frequencies where coherence is good, up to fcut

figure(3);
plot(F,Cxy);
hold on
plot(F(F_good),Cxy(F_good),'or')
grid on
ax2=gca;
ylabel('Data Correlaion')
xlabel('Frequency (Hz)')
ax2.FontName = 'Times New Roman';
ax2.FontSize = 14;
xlim([0 2300]);
ylim([.9 1.001]);
saveas(gcf,sprintf('datacorr-18signal-from-1hz-upto-2khz.fig'))
saveas(gcf,sprintf('datacorr-18signal-from-1hz-upto-2khz.png'))
saveas(gcf,sprintf('datacorr-18signal-from-1hz-upto-2khz.eps'),'epsc')

%% Computation of Impedance Spectra Using FFT
NFFT=2^nextpow2(length(current));   % Number of FFT
m1_1st_order_ecm=fft(current,NFFT);   %fast fourier transform of current in NFFT point
f2_2nd_order_ecm=fft(voltage,NFFT);    %fast fourier transform of voltage in NFFT point
f3=m1_1st_order_ecm(1:NFFT/2) ;           
f4=f2_2nd_order_ecm(1:NFFT/2) ; 
xfft1=fs*(0:(NFFT/2)-1)/NFFT ;   %convert x axis to frequency scale
%%  peak detection of FFT curves
[y,locs]=findpeaks(abs(f3),'MinPeakDistance',50,'MinPeakHeight',350);  %should be set for each data set differently
n=1; 
for i=1:1:length(locs)
     current_peak(n)=f3(locs(i)) ;
     voltage_peak(n)=f4(locs(i)) ;
     f_good(n)=xfft1(locs(i));
     n=n+1;
end
%%  Nyquist Plot of impedance
w=logspace(-1,6,5000);
figure(4);
impedance=voltage_peak./current_peak ;
plot(real(impedance),-imag(impedance),'--*k','MarkerSize',7,'LineWidth',1);

%% Estimation of first order ECM and Impedance Spectra Using SYSID

% Make data ready for system identification
ddata = iddata(voltage, current, Ts); 
% Filter data for sysid
ddata_filt = idfilt(ddata,2*pi*[.01 1000]);% filter the impednace data below the real axis 
% Define the structure of the model to be identified 
Rinf_1st_order_ecm=0.9*(0.1+(0.2-0.1)*rand);%Ohm
R1_1st_order_ecm=0.1*(.1+(.15-.1)*rand);%Ohm
C1_1st_order_ecm=100*(.01+(.05-.01)*rand);%F
m1_1st_order_ecm=Rinf_1st_order_ecm;
m0_1st_order_ecm=(Rinf_1st_order_ecm+R1_1st_order_ecm)/(R1_1st_order_ecm*C1_1st_order_ecm);
n0_1st_order_ecm=1/(R1_1st_order_ecm*C1_1st_order_ecm);

model_struc=idtf([m1_1st_order_ecm m0_1st_order_ecm],[1 n0_1st_order_ecm]);%initial model
model_struc.Structure.num.Minimum = eps;
model_struc.Structure.den.Minimum = eps;
opt = tfestOptions('SearchMethod', 'auto');

% Estimate the model
estimated_tf_1st_order_ecm = tfest(ddata_filt,model_struc,opt) % estimated model

% Get the parameters of the estimated model
params=getpvec(estimated_tf_1st_order_ecm);
m1hat_1st_order_ecm=params(1);
m0hat_1st_order_ecm=params(2);
n0hat_1st_order_ecm=params(3);

% Compute the estimation of parameters
Rinfest_1st_order_ecm=m1hat_1st_order_ecm
R1est_1st_order_ecm=m0hat_1st_order_ecm/n0hat_1st_order_ecm-m1hat_1st_order_ecm
C1est_1st_order_ecm=1/(n0hat_1st_order_ecm*R1est_1st_order_ecm)

% Plot impedance Spectra of the estimated model
Zest_1st_order_ecm=Rinfest_1st_order_ecm+R1est_1st_order_ecm./(R1est_1st_order_ecm*C1est_1st_order_ecm*1i*w+1);
figure(4)
hold on
plot(real(Zest_1st_order_ecm),-imag(Zest_1st_order_ecm),'b','LineWidth',2);

%% Estimation of first order ECM and Impedance Spectra Using SYSID
% Model Structure
Rinf_2nd_order_ecm=0.9*(0.1+(0.2-0.1)*rand);%Ohm

R1_2nd_order_ecm=0.1*(.1+(.15-.1)*rand);%Ohm
C1_2nd_order_ecm=100*(.01+(.05-.01)*rand);%F

R2_2nd_order_ecm=.5*R1_2nd_order_ecm;%1*(0.12+(.18-0.12)*rand);
C2_2nd_order_ecm=100*(.001+(.005-.001)*rand);

tau1_2nd_order_ecm=R1_2nd_order_ecm*C1_2nd_order_ecm;
tau2_2nd_order_ecm=R2_2nd_order_ecm*C2_2nd_order_ecm;
a1=1/tau1_2nd_order_ecm;
a2=1/tau2_2nd_order_ecm;
b1=1/C1_2nd_order_ecm;
b2=1/C2_2nd_order_ecm;
f2_2nd_order_ecm=Rinf_2nd_order_ecm;
f1_2nd_order_ecm=Rinf_2nd_order_ecm*(a1+a2)+b1+b2;
f0_2nd_order_ecm=Rinf_2nd_order_ecm*a1*a2+b1*a2+b2*a1;
g1_2nd_order_ecm=a1+a2;
g0_2nd_order_ecm=a1*a2;

model_struc=idtf([f2_2nd_order_ecm f1_2nd_order_ecm f0_2nd_order_ecm],[1 g1_2nd_order_ecm g0_2nd_order_ecm]);%initial model

model_struc.Structure.num.Minimum = eps;
model_struc.Structure.den.Minimum = eps;
opt = tfestOptions('SearchMethod', 'auto');

% Estimate the model
estimated_tf_2nd_order_ecm = tfest(ddata_filt,model_struc,opt) % estimated model

% Compute parameters from the estimated model
params=getpvec(estimated_tf_2nd_order_ecm);
f2hat_2nd_order_ecm=params(1);
f1hat_2nd_order_ecm=params(2);
f0hat_2nd_order_ecm=params(3);
g1hat_2nd_order_ecm=params(4);
g0hat_2nd_order_ecm=params(5);


Rinfest_2nd_order_ecm=f2hat_2nd_order_ecm
ROOTS=-roots([1 g1hat_2nd_order_ecm g0hat_2nd_order_ecm]);
a1hat=max(ROOTS);
a2hat=min(ROOTS);
tau1hat_2nd_order_ecm=1/a1hat;
tau2hat_2nd_order_ecm=1/a2hat;
B=[f1hat_2nd_order_ecm-(a1hat+a2hat)*f2hat_2nd_order_ecm;f0hat_2nd_order_ecm-a1hat*a2hat*f2hat_2nd_order_ecm];
A=[1 1;a2hat a1hat];
X=inv(A)*B;
C1est_2nd_order_ecm=1/X(1)
C2est_2nd_order_ecm=1/X(2)
R1est_2nd_order_ecm=tau1hat_2nd_order_ecm/C1est_2nd_order_ecm
R2est_2nd_order_ecm=tau2hat_2nd_order_ecm/C2est_2nd_order_ecm

% plot corresponding impedance spectra
Zjwhat=(f2hat_2nd_order_ecm*(j*w).^2+f1hat_2nd_order_ecm*j*w+f0hat_2nd_order_ecm)./((j*w).^2+g1hat_2nd_order_ecm*j*w+g0hat_2nd_order_ecm);
h1=figure(4)
hold on
plot(real(Zjwhat),-imag(Zjwhat),'r','LineWidth',2)


%% Improve Plot 
figure(4)
ax4=gca;
ax4.FontName = 'Times New Roman';
ax4.FontSize = 14;
strmax = ['0.2 C-rate'];
text(.164,8.5e-3,strmax,'HorizontalAlignment','left','FontName', 'Times New Roman', 'FontSize', 14);
set(gca,'TickLabelInterpreter','latex');
xlabel('Real\{Z\}');
ylabel('-Imag\{Z\}');
grid on
legend({'FFT','1st-order estimated ECM','2nd-order estimated ECM'},'box','off','location','best','FontSize',14,'Interpreter','latex');
saveas(gcf,sprintf('eis-18signal-from-1hz-upto-2khz.fig'))
saveas(gcf,sprintf('eis-18signal-from-1hz-upto-2khz.png'))
saveas(gcf,sprintf('eis-18signal-from-1hz-upto-2khz.eps'),'epsc')

%% current plot in time domain
figure(1);
plot(t,current,'b');
xlabel('Time (sec)')
ylabel('Current (A)')
ax1=gca;
ax1.FontName = 'Times New Roman';
ax1.FontSize = 14;
saveas(gcf,sprintf('current-18signal-from-1hz-upto-2khz.fig'))
saveas(gcf,sprintf('current-18signal-from-1hz-upto-2khz.png'))
saveas(gcf,sprintf('current-18signal-from-1hz-upto-2khz.eps'),'epsc')

%% Plot the voltage signals
h1=figure(2)

% plot of the measured voltage
hold on
plot(t,voltage,'--k','LineWidth',1);

% plot response of the estimated model to the injected current
y_1st_order_ecm=lsim(estimated_tf_1st_order_ecm,current,t);
figure(2)
hold on
plot(t,y_1st_order_ecm,'b','LineWidth',1)

% plot response of the estimated model to the injected current
y_2nd_order_ecm=lsim(estimated_tf_2nd_order_ecm,current,t);
hold on
plot(t,y_2nd_order_ecm,'r','LineWidth',1)

xlabel('Time (sec)')
ylabel('voltage (V)')
ax2=gca;
ax2.FontName = 'Times New Roman';
ax2.FontSize = 14;
box on;

legend({'Measured','1st-order estimated ECM','2nd-order estimated ECM'},'box','off','location','southeast','FontSize',14,'Interpreter','latex');


ylim([-.06 .08]);
MagInset(h1, -1, [1 1.0025 -0.01 .005], [5 15 .035 .075], {'NW','NW';'SE','SE'});% real
ax3=gca;
ax3.FontName = 'Times New Roman';
ax3.FontSize = 10;

saveas(gcf,sprintf('voltage-18signal-from-1hz-upto-2khz.fig'))
saveas(gcf,sprintf('voltage-18signal-from-1hz-upto-2khz.png'))
saveas(gcf,sprintf('voltage-18signal-from-1hz-upto-2khz.eps'),'epsc')
