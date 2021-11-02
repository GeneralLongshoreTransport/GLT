% Example 2
%
% GIVEN: Breaking waves with significant wave height Hs_b= 2 (m), wave period 
% Tp = 6 (s) and breaking wave angle Theta_b = 20 (deg) with respect to the 
% shore-normal direction. The median grain size D50 is equal to 12.5 mm. 
% FIND: Calculate the losgshore transport rate by GLT procedure.
% SOLUTION: The longshore transport rate is 0.0626 [m^3/s]

Hs_b= 2;    % significant breaking wave height (m)
Theta_b= 20;% breaking wave angle (deg)
Tp= 6;      % wave period (s)
d50= 12.5;  % nominal diameter of the units (mm)

Q_GLT= GLT(Hs_b,Tp,Theta_b,d50,2) % longshore transport rate [m^3/s]