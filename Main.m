%% Header

% Project
% Reaction of Intrest

% C02 + 4 H2 --> CH4 + 2 H20


clc
close all
clear all

global  Fao R catdens U a Ta Cn2o Vdot P Tfeed epsi thetaN2 Tref ThetaB

%% Given Variables

disp('Welcome');
crew=input('Please enter the number of crew members ---> ');

while crew<0
    disp(' ');
    disp('Error: Please enter a positive interger');
    disp(' ');
    crew=input('Please enter the number of crew members ---> ');
    disp(' ');
end

disp(' ');
disp(' ');
disp(' ');
disp('Be advised the reactor size is 10 mL and is designed for 10 people');

Tfeed=573.15;       %    [K]  Preheated

P=0.9;              %   [atm] Lower than atmospheric
R=8.314;            % [m3*Pa/K/mol]

Cao=(P*101325)/(R*Tfeed*6.05); % [mol/m3] Initial Concentration from Ideal Gas Law
Cbo=4.5*Cao;                   % [mol/m3] Initial Concentration
Cn2o=.1*(Cao+Cbo);             % [mol/m3] Initial Concentration
HumanCO2offgas=15; %[L/hr]

IdealFao=crew*HumanCO2offgas/1000; %[m3/hr] Co2 off gas

Vdot=IdealFao/Cao; % Solving for appropriate Vdot

Fao=Cao*Vdot;                  % [mol/h]
Fbo=Cbo*Vdot;                  % [mol/h]
Fn2o=Cn2o*Vdot;                % [mol/h]

ThetaB=Cbo/Cao;

Tref=298.15; %[K]

delta=-2;
Fto=(Fao+Fbo+Fn2o);
epsi=(Fao/Fto)*delta;
thetaN2=Cn2o/Cao;
Ta=294.15; %Ambient Air in K
catdens=1180000; %[g/m3] Density of Catalyst

U=7.9;  % Calculated by hand
a=80;  % Calculated by hand a=4/D D = 0.05 m


xo=0; %Inital conversion

Volume=.00001; %[m3]

%% Functions and Solvers


Initial= [xo, Tfeed];

%Solve Differential
[V,y] = ode45(@Projectfunc,[0,Volume], Initial);

%Plotting the results
subplot(4,1,1)
plot(V*1000000, y(:,1)) %Convert Volume to mL from m^3
ylabel('Conversion','FontSize',14)
xlabel('Volume (mL)','FontSize',14)
%title('Reactor Data')
subplot(4,1,2)
plot(V*1000000, y(:,2)) %Convert Volume to mL from m^3
xlabel('Volume (mL)','FontSize',14)
ylabel('Temperature, T(K)','FontSize',14)


%% Varying Reactor Inlet Temperature

Tfeed=300;       %    [K]
n=1;

while Tfeed<600

Initial= [xo, Tfeed];

%Solve Differential
[V,y] = ode45(@func,[0,Volume], Initial);

Conversion(n)=y(end,1);
Temperature(n)=y(end,2);

Tin(n)=Tfeed;

Tfeed=Tfeed+10;
n=n+1;

end



subplot(4,1,3)
plot(Tin,Temperature) %Convert Volume to mL from m^3
xlabel('Temperature In (K)','FontSize',14)
ylabel('Outlet Temperature (K)','FontSize',14)

subplot(4,1,4)
plot(Tin, Conversion) %Convert Volume to mL from m^3
xlabel('Temperature In (K)','FontSize',14)
ylabel('Conversion','FontSize',14)

disp(' ');
disp('Optimium Operating Parameters');
disp(' ');

[Conv,i]=max(Conversion);

Max=Tin(i)-273;

fprintf('Please preheat the feed in degrees C to %3.4f\n',Max) %Display the answer

disp(' ');
disp('Operating Conditions');
disp(' ');


fprintf('Max conversion achievable %3.4f\n',Conv) %Display the answer



ReactT=Temperature(i)-273;

fprintf('Reactor Temperature in degrees C %3.4f\n',ReactT) %Display the answer