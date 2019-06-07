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
% This Matlab file runs test #1 of the above manuscript 
%==========================================================================
% The Matlab codes compute the impedance spectra by using FFT, 
% and by using system identification methods. 
%
% Two first-order Randles equivalent circuit model (ECM) is estimated to
% fit the impedance spectra obtained from FFT. 
%
% The first-order Randles circuit is shown below:
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


%% setup verification using a multisine signal consist of two freq. , 5Hz and 250Hz
close all
load setup-verification-by-2-sinusoidal-5Hz-250Hz

%% Setup Sample Time
L=20;  %put the signal time duration in seconds
Ts=L/length(current);
fs=1/Ts;
t = Ts:(1/fs):L;
w=logspace(-2,6,1000);
Rinf_true=0.2199;  %thess values have been measured by LCR meter device
R1_true=0.1108;
C1_true=0.0239;
Z=Rinf_true+R1_true./(R1_true*C1_true*1i*w+1);
plot(real(Z),-imag(Z),'k','LineWidth',2);
hold on
%% Estimation of first order ECM and Impedance Spectra Using SYSID

% Make data ready for system identification
ddata = iddata(voltage, current, Ts); 
% Filter data for sysid
ddata_filt = idfilt(ddata,2*pi*[.01 1000]);% filter the high frequency ripples on data
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
%%
% Plot impedance Spectra of the estimated model
Zest_1st_order_ecm=Rinfest_1st_order_ecm+R1est_1st_order_ecm./(R1est_1st_order_ecm*C1est_1st_order_ecm*1i*w+1);

plot(real(Zest_1st_order_ecm),-imag(Zest_1st_order_ecm),'--r','LineWidth',2);
legend({'Original Randles Circuit','Estimated Circuit'},'location','best','FontSize',14,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',14);
xlabel('Real\{Z\}','FontSize',14);
ylabel('-Imag.\{Z\}','FontSize',14);
axis equal
saveas(gcf,sprintf('EIS-2signal-test.fig'))
saveas(gcf,sprintf('EIS-2signal-test.png'))
saveas(gcf,sprintf('EIS-2signal-test.pdf'))
saveas(gcf,sprintf('EIS-2signal-test.eps'),'epsc')