function [QGlt] = GLT(Hs,Tp,theta,d50,dataloc)
% DESCRIPTION:
% Matlab function for the longshore transport rate calculation by the general 
% longshore transport model
% 
% [QGlt] = GLT(Hs,Tp,theta,d50,dataloc)
% where QGlt is the longshore transport rate in volume per unit time [m^3/s]
%  
% INPUT VARIABLES:
% Hs         % significant wave height (m)
% Tp         % peak wave period (s)
% theta      % wave angle (deg)
% d50        % nominal diameter of the units (mm)
% dataloc    % wave data input location (1 - offshore; 2 - breaking point)
%
% OUTPUT VARIABLES:
% QGlt       % LT rate in volume per unit time [m^3/s]
if dataloc == 1                                                         % offshore wave data input location
  Hso=Hs;                                                               % offshore significant wave height [m]
  theta_o=theta;                                                        % offshore wave angle [deg]
  g=9.806;                                                              % gravity acceleration [m/s^2]
  gamma=1.42;                                                           % breaking index [-]     
  n=0.4;                                                                % porosity index [-] 
  deltaglt=1.65;                                                        % relative mass density of the unit [-]
  Tm=0.78.*Tp;                                                          % mean wave period [s]
  Lo=(g.*(Tm.^2))/(2.*pi);                                              % wavelength in offshore condition [m] 
  depth_o=Lo/2;                                                         % depth in offshore condition [m] 
  Hko=1.55.*Hso;                                                        % characteristic wave height [m]
  co=Lo./Tm;                                                            % offshore wave celerity [m/s]  
  ko=2.*pi./Lo;                                                         % offshore wave number [m^-1] 
  cg=(co./2).*(1.+((2.*ko.*depth_o)./(sinh(2.*ko.*depth_o))));          % wave group celerity [m/s]
  Hkb=((Hko.^2).*cg.*(cosd(theta_o)).*(gamma./g).^(0.5)).^(0.4);        % characteristic wave height at breaking [m] (Eq.8)
  ckb=sqrt(g.*Hkb/gamma);                                               % characteristic wave celerity at breaking point [m/s] (Eq.9)
  sintkb=(ckb.*sind(theta_o))./co;                                      % sin of characteristic wave angle at breaking (Eq.10)
  Ns=(0.89*Hkb./(1.55*deltaglt*(d50*10^-3)));                           % modified stability number [-] (Eq.3)
  Hsb=Hkb./1.55;                                                        % significant wave height at breking point [m]
  depth_b=Hsb./gamma;                                                   % depth at wave breaking point [m]
  Lb=zeros(size(depth_b));                                              
  for i=1:size(depth_b);                                                
   guess = Lo;                                                                
   Lb = (g.*Tm.^2)/(2.*pi).*tanh((2.*pi).*(depth_b./guess));            % ...iterative solution procedure of...
   diff = abs(Lb-guess);                                                % ...linear dispersion relationship 
    while diff>0.01;
     diff = abs(Lb-guess);
     guess = Lb+(0.5.*diff);
     Lb = (g.*Tm.^2)/(2.*pi).*tanh((2.*pi).*(depth_b./guess));  
   end
  end
  kb=2.*pi./Lb;                                                         % wave number at breaking point [m^-1] 
  ld=(((1.4.*Ns)-1.3).*d50.*10.^-3)./((tanh(depth_b.*kb)).^(2));        % displacement length [m] (Eq.2)
  Nod_cal=(Ns.^6.87).*exp(0.051.*(log(Ns).^3)-0.87.*(log(Ns).^2)-4.92); % displaced particles at the end of 1000 waves attack [-] (Eq.4)
  Sn_cal=ld.*Nod_cal./(d50*10^-3)./1000.*sin(sintkb);                   % longshore transport as number of displaced units per wave [-] (Eq.1)
  Q_cal_sec=Sn_cal.*(d50.*10.^-3).^3./(Tm.*(1-n));                      % longshore transport rate as volume per unit time [m^3/s] (Eq.11)
  QGlt = Q_cal_sec;                                                     % QGlt MatLab function output [m^3/s]
    
elseif dataloc == 2                                                     % offshore wave data input location
  Hsb=Hs;                                                               % significant wave height at breaking point [m]
  theta_b=theta;                                                        % wave angle at breaking point [deg]
  g=9.806;                                                              % gravity acceleration [m/s^2]
  gamma=1.42;                                                           % breaking index [-] 
  depth_b=Hsb./gamma;                                                   % depth at wave breaking point [m] 
  n=0.4;                                                                % porosity [-]
  deltaglt=1.65;                                                        % relative mass density of the unit [-]
  Tm=0.78.*Tp;                                                          % mean wave period [s]
  for i=1:size(depth_b);
   guess = (g.*(Tm.^2))/(2.*pi);
   Lb = (g.*Tm.^2)/(2.*pi).*tanh((2.*pi).*(depth_b./guess));            % wavelenght at breaking [m]: ...
   diff = abs(Lb-guess);                                                % ... iterative solution procedure ...
   while diff>0.01                                                      % ... of linear dispersion relationship
    diff = abs(Lb-guess);
    guess = Lb+(0.5.*diff);
    Lb = (g.*Tm.^2)/(2.*pi).*tanh((2.*pi).*(depth_b./guess));  
   end
  end  
  Hkb=1.55.*Hsb;                                                        % characteristic wave height at breaking [m] (Eq.5) 
  cb=Lb./Tm;                                                            % wave celerity at breaking point [m/s] 
  kb=2.*pi./Lb;                                                         % wave number at breaking point [m^-1]
  ckb=sqrt(g.*Hkb/gamma);                                               % characteristic wave celerity at breaking point [m/s] (Eq.6)
  sintkb=(ckb.*sind(theta_b))./cb;                                      % sin of characteristic wave angle at breaking (Eq.7)
  Ns=(0.89*Hkb./(1.55*deltaglt*(d50*10^-3)));                           % modified stability number [-] (Eq.3)
  ld=(((1.4.*Ns)-1.3).*d50.*10.^-3)./((tanh(depth_b.*kb)).^(2));        % displacement length [m] (Eq.2)
  Nod_cal=(Ns.^6.87).*exp(0.051.*(log(Ns).^3)-0.87.*(log(Ns).^2)-4.92); % displaced particles at the end of 1000 waves attack [-] (Eq.4)   
  Sn_cal=ld.*Nod_cal./(d50*10^-3)./1000.*sin(sintkb);                   % longshore transport as number of displaced units per wave [-] (Eq.1)
  Q_cal_sec=Sn_cal.*(d50.*10.^-3).^3./(Tm.*(1-n));                      % longshore transport rate as volume per unit time [m^3/s] (Eq.11)
  QGlt = Q_cal_sec;                                                     % QGlt MatLab function output [m^3/s]
 end
end