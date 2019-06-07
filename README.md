
A controllable Electrochemical Impedance Spectroscopy (cEIS) Device

By E. Sadeghi, M. H. Zand, M. Hamzeh, and S.M.M. Alavi
 
Department of Electrical Engineering, Shahid Beheshti University, Tehran, Iran. m_alavi@sbu.ac.ir, http://arg.sbu.ac.ir
      

The computation of battery impedance spectra by using FFT, and system identification (sysid) method in cEIS device. 
Two equivalent circuit models (ECMs) are estimated to fit the impedance spectra obtained from FFT. The first ECM model is a first-order Randles circuit as shown below:

                               R1
                      |-----/\/\/\/------|
            Rinf      |                  |
    -------/\/\/\-----|                  |----
                      |      1/(C1s)     |
                      |-------||---------|

                        v      m1*s + m0
                   Z = --- = --------------
                        i        s + n0

where Rinf=m1, R1=m0/n0-m1, and C1=1/(n0*R1). 

The second identified ECM model is a second-order Randles circuit as shown below:

                               R1                      R2
                       |-----/\/\/\/------|    |-----/\/\/\/-----|
             Rinf      |                  |    |                 |
     -------/\/\/\-----|                  |----|                 |--------
                       |      1/(C1s)     |    |      1/(C2s)    |
                       |-------||---------|    |--------||-------|


                              v      m2*s^2 + m1*s + m0
                          Z = --- = ----------------------
                               i        s^2 + n1*s + n0
where, tau1=R1*C1; tau2=R2*C2; a1=1/tau1; a2=1/tau2; b1=1/C1; b2=1/C2; m2=Rinf; m1=Rinf*(a1+a2)+b1+b2; m0=Rinf*a1*a2+b1*a2+b2*a1; 
n1=a1+a2; n0=a1*a2;
