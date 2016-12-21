function dy = func(~,y)

global  Fao R catdens U a Ta Cn2o Vdot P Tfeed epsi thetaN2 Tref ThetaB

dy = zeros(2,1);    % Column vector

%% Flow Rates


Fco2= Fao*(1-y(1))/((1+epsi*y(1))*y(2)/Tfeed);              %  A
Fh2=  Fao*(ThetaB-4*y(1))/((1+epsi*y(1))*y(2)/Tfeed);       %  B
Fch4= Fao*(y(1))/((1+epsi*y(1))*y(2)/Tfeed);                %  C
Fh2o= Fao*(2*y(1))/((1+epsi*y(1))*y(2)/Tfeed);              %  D
Fn2=  Cn2o*Vdot;                                            %  I  Constant

%% Partial Pressures

Ft=Fch4+Fco2+Fh2+Fh2o+Fn2;

Pch4=(Fch4/Ft)*P; %[atm]
Pco2=(Fco2/Ft)*P; %[atm]
Ph2=(Fh2/Ft)*P;   %[atm]
Ph2o=(Fh2o/Ft)*P; %[atm]

%% Cp Data

CpData=[5.457, 1.045E-3, 0,-1.157E5;3.249, 0.422E-3, 0,0.083E5;1.702, 9.081E-3, -2.164E-6,0;3.470, 1.450E-3, 0,0.121E5];
CpDataN2=[3.28,.593E-3,0,.04E5];



%% Specific Heat

CpA =R*(CpData(1)+CpData(5)*y(2)+CpData(13)*y(2)^-2) ;		    	% Heat capacity of A
CpB =R*(CpData(2)+CpData(6)*y(2)+CpData(14)*y(2)^-2) ;		    	% Heat capacity of B
CpC =R*(CpData(3)+CpData(7)*y(2)+CpData(11)*y(2)^2) ;               % Heat capacity of C
CpD =R*(CpData(4)+CpData(8)*y(2)+CpData(16)*y(2)^-2);               % Heat capacity of D
CpN2=R*(CpDataN2(1)+CpDataN2(2)*y(2)+CpDataN2(4)*y(2)^-2);          % Heat capacity of I

dCp = 2*CpD + CpC - CpA - 4*CpB;

%% Heat of reaction

DeltaHr=-1*164.732*1000;         % Heat of Reaction at Reference Temperature

dHR=DeltaHr+(dCp*(y(2)-Tref));

dHR=0.25*dHR;                    % On Stoichiometric Ratio Basis

%% Rate Constants
Keq=1/(2.06199512*10^11*exp(-22430/y(2)));         % Converted to atm
k=catdens*(1.0635*10^11)*exp(-113497.4/(R*y(2)));  % [mol/m3*hr]
Kco2 = (9.099*10^-7)*exp(69691.8/(R*y(2)));        % unitless
Kh2 = (9.6104*10^-4)*exp(39942/(R*y(2)));          % unitless

%% Reaction Rates
Beta = (1/Keq)*(Pch4*Ph2o^2)/(Pco2*Ph2^4);
rco2 = (1-Beta)*(k*Kco2*Kh2^4*Pco2*Ph2^4)/((1+Kco2*Pco2+Kh2*Ph2)^5);

%% Differntial Equations

dy(1)= rco2/Fao;
dy(2)= ((U*a*(Ta-y(2)))-(rco2*dHR))/(Fao*((CpA+(4.5*CpB)+(thetaN2*CpN2))+y(1)*dCp));

