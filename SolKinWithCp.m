%Author: Alka Potdar
% Run using MATLAB R2009b
% Computational modeling and analysis of macrophage iron release (April 2014)

% Fit solution kinetics experimental data in the presence of both cp and Tf, that is compute k2 and k3 and use
% k1,k11 and k4from prior simulation without  cp

function SolKinWithCp

clear 
clc

xdata=importdata('x_new_withcp.mat');
ydata=importdata('y_new_withcp.mat');

% ydata is in millimolar, convert to micromolar


ydata=transpose(ydata*1000);


% Initial values of k2 and k3 found by series of runs of with lsqcurvefit

[x,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit(@allfun,[793.1700    1.6601    5.6005],xdata,ydata,[0.000001 0.000001 0.0000001],[1000 1000 1000])



%--------------------------------------------------------------------------

% function, f for iron transferrin (ferric bound to transferrin, monoferric and holoferric forms)

%--------------------------------------------------------------------------
% computing the values of function f using ode15s


function f=allfun(x,xdata)


y = ode15s(@odefun,[0 120],[120 0 0 0 55 0 0.9],[],x(1),x(2),x(3));

% Fe+2: 120 uM; Tf: 55 uM; Cp2+ : 0.9 uM 

xfull= 0:5:120;

sxint = deval(y,xdata);

sxint_full = deval(y,xfull);

f=(sxint(3,:)+2*sxint(4,:));

f_full=(sxint_full(3,:)+2*sxint_full(4,:));

g=transpose(f);



dlmwrite('ysim.csv',f_full);

dlmwrite('xsim.csv',xfull);

dlmwrite('ferrous_withcp.csv',sxint_full(1,:));

dlmwrite('species_withcp_decay.csv',sxint_full);


% amount of non-transferrin bound fe+3

decay=sxint_full(2,:);

dlmwrite('decay_withcp.csv',decay);

 

%--------------------------------------------------------------------------
% defining the set of ordinary differential equations for the reaction system

function df=odefun(t,c,k2,k3,k6)


df = zeros(7,1);    % a column vector


    co2= 0.13; %concentration of oxygen (uM) for 1% oxygen
    
    k1 = 3.09; % (per (uM.s)
    
    k11= 153.5; % (per s)
    
    k4 = 0.04; % (per (uM.s)
    
    

df(1)=-k1*c(1)*co2+k11*c(2)-k3*c(1)*c(7);% eqn for fe+2

df(2)=k1*c(1)*co2-k11*c(2)+k3*c(1)*c(7)-k4*c(2)*c(5)-k4*c(2)*c(3)-k6*c(2);% eqn for fe+3,k6 is  first order decay rate constant for fe3+ ions in solution

df(3)=k4*c(2)*c(5)-k4*c(2)*c(3);% eqn for monoferric transferrin

df(4)=k4*c(2)*c(3);% eqn for holo transferrin

df(5)=-k4*c(2)*c(5);% eqn for transferrin

df(6)=-k2*c(6)*co2+k3*c(1)*c(7);%eqn for cp+1

df(7)=k2*c(6)*co2-k3*c(1)*c(7);%eqn for cp+2

% oxygen constant at co2;  

tdata=t;
oxid=df(1,1);
tfbind=df(3,1);
%fprintf('%6f\n ', t)
%fprintf('%6f\n', tfbind);


%--------------------------------------------------------------------------



