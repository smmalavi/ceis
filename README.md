# A controllable Electrochemical Impedance Spectroscopy (cEIS) Device

## By E. Sadeghi, M. H. Zand, M. Hamzeh, and S.M.M. Alavi
## Department of Electrical Engineering, Shahid Beheshti University, Tehran, Iran.
      
This repository provides you with the Matlab codes for test #1 and test #2 of the above manuscript. The current and voltage data is loaded in the beginning of each test. I would refer you to the manuscript for the information about tests, and data. The Matlab codes compute the impedance spectra by using FFT, and system identification methods. 

In test #1, a first-order Randles equivalent circuit model (ECM) is estimated. In test #2, the first- and second-order Randles ECMs are estimated and their accuracies are compared. 

The first-order Randles ECM is shown below:

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

The second-order Randles ECM is shown below:

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
n1=a1+a2; n0=a1*a2

Feel free to contact me for any further information you have

Mahdi Alavi

