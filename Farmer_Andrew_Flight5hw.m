%Andrew Farmer Flight Homework 5
clc
clear all
close all
format long
delAC = 17; %distance between horizontial tail and aerodynamic center
h = 2000; %alt
To = 288.15; %Temp
bw = 30; %wing span
sw = 300; %surface area of wing
aw = bw^2/sw; %aspect ratio wing
cw = 11.32; %chord of wing
bh = 18.29; %span of horizontial tail
sh = 63.7; %surface area of horizontial tail
ch = 3.48;  %chord of horizontial tail
ah = bh^2/sh; %aspect ratio of horizontial tail
cmu = 0; 
xacw = .25; %delXacW
xach = (delAC - xacw*cw)*(1/cw);
W = 20500; %weight
g = 32.17;
m = (20500/g);
Iyy = 55814;
uo = 600;
wo = 12.5;
ao = 1.2*(pi/180); %alpha naught
bo = 0*(pi/180); %sideslip angle naught
phio = 0*(pi/180); 
de = -.8*(pi/180); % elevator control stiffness
claw  = 5.7;
clah = 5.7;
cdo = .0175;
iw = 0;
aow  = 0;
ih = 0;
cmu = 0; 
rho = .00237;
M = 0.00198469182644; %air molar mass
R = 6.1324; %universial gas constant
L = .0019812; %standard temperature laps rate
rhoinf = rho*exp(-(((g*M*h)/(R*To))-((L*h)/(To))));
qinf = (rhoinf*uo^2)*.5; %calculates rhoinf
cla = claw + clah *(sh/sw);
clw = claw*(ao-aow);
%xach = 14.17/cw;
clh = clah*(ao+ih+de);
cda = (2/pi)*(clw/(aw)*claw + clh/(ah)*clah*(sh/sw));
cma = claw*xacw - clah*xach*(sh/sw);
cm = claw*(ao+iw-aow)*xacw - clah*(ao+aow-iw+ih+de)*xach*(sh/sw);
%Part a
fprintf('Cla %f \n',cla)
fprintf('Cda %f \n',cda)
fprintf('Cma %f \n',cma)
%PARTB
SM = -cma/cla;
disp('The aircraft is stable in pitch as the static margin is positive')
%PARTC
distance = -SM;
fprintf('The static Margin is %f\n',SM)
fprintf('The Center of Gravity can move AFT %f percent of the mean aerodynamic chord before the aircraft becomes unstable in pitch\n',SM*100)
%part D
cl =  clw +clh*(sh/sw);
cd = cdo + (1/(pi*aw))*(clw^2 +clh^2 *(aw/ah)*(sh/sw));
fprintf('CL %f \n',cl)
fprintf('CD %f \n',cd)
%Part e
xu = -cos(bo)*sw*rhoinf*uo*cd;
zu = -sw*rhoinf*uo*cl;
mu = qinf*sw*bw*cmu + cmu*sw*bw*rhoinf*uo; %Notes say Cm at x but hw solution says cm dot idk which is correct
xa = cda*cos(bo)*qinf*sw + cl*qinf*sw;
za = -cda*cos(bo)*qinf*sw - cla*qinf*sw;
ma = (1/uo)*cma*qinf*sw*cw;
xq = 0;
zq = -clah*qinf*sh*((xach*cw)/(uo));
mq = -clah*qinf*sh*((xach*cw)^2)/(uo);
format bank
partEVals = [xu zu mu xa za ma xq zq mq]';
disp(partEVals)
alon = [xu/m, xa/m, (xq/m-wo),-g*1;
        zu/(uo*m), za/(uo*m), zq/(uo*m)+1,0;
        (mu/Iyy) + (ma*zu)/(m*uo*Iyy),(ma/Iyy) + (ma*za)/(m*uo*Iyy),(mq/Iyy) + (ma*(zq+uo*m))/(m*uo*Iyy),0;
        0 0 1 0];
disp(alon)
