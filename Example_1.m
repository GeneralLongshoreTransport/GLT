% Example 1
%
% GIVEN: Offshore waves with significant wave height Hs_o= 2 (m), wave period 
% Tp = 6 (s) and offshore wave angle Theta_o = 20 (deg) with respect to the 
% shore-normal direction. The median grain size Dn50 is equal to 12.5 mm. 
% FIND: Calculate the losgshore transport rate by GLT procedure.
% SOLUTION: The longshore transport rate is 0.0230 [m^3/s]

Hs_o= 2;    % significant breaking wave height (m)
Theta_o= 20;% breaking wave angle (deg)
Tp= 6;      % wave period (s)
d50= 12.5;  % nominal diameter of the units (mm)

Q_GLT= GLT(Hs_o,Tp,Theta_o,d50,1) % longshore transport rate [m^3/s]