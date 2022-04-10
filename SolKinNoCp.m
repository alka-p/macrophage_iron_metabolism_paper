%Author: Alka Potdar
% Run using MATLAB R2009b
% Computational modeling and analysis of macrophage iron release (April 2014)

%Fit solution kinetics experimental data without cp in presence of Tf, that is compute k1,k11, k4 and k6

function SolKinNoCp
clear 
clc



xdata=importdata('xdata_nocp.mat');
ydata=importdata('ydata_nocp.mat');


% ydata is in millimolar, convert to micromolar
ydata=transpose(ydata*1000);



% Initial values of k1,k11, k4 and k6 found by series of runs with lsqcurvefit


[x,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit(@allfun,[3.0925  153.6074    0.0409   10.7243],xdata,ydata,[0.00001 0.001 0.001 0.001],[1*1000 1*1000 1*1000 1*1000])


%--------------------------------------------------------------------------

% function, f for iron transferrin (ferric bound to transferrin, monoferric and holoferric forms)

%--------------------------------------------------------------------------
% computing the values of function f using ode15s

function f=allfun(x,xdata)

y = ode15s(@odefun,[0 120],[120 0 0 0 55],[],x(1),x(2),x(3),x(4));

% Fe+2: 120 uM; Tf: 55 uM

xfull= 0:5:120;
sxint_full = deval(y,xfull);
f_full=(sxint_full(3,:)+2*sxint_full(4,:));


sxint = deval(y,xdata);

f=(sxint(3,:)+2*sxint(4,:));

g=transpose(f);


dlmwrite('xsim_nocp.csv',xfull);

dlmwrite('ysim_nocp.csv',f_full);

% amount of non-transfrrin bound fe+3

decay=sxint_full(2,:);

dlmwrite('species_nocp_nodecay.csv',sxint_full)



%--------------------------------------------------------------------------
% defining the set of ordinary differential equations for the reaction system

function df=odefun(t,c,k1,k11,k4,k6)


df = zeros(5,1);    % a column vector


co2= 0.13; %concentration of oxygen (uM) for 1% oxygen

df(1)=-k1*c(1)*co2+k11*c(2); % eqn for fe+2

df(2)=k1*c(1)*co2-k11*c(2)-k4*c(2)*c(5)-k4*c(2)*c(3)-k6*c(2); % eqn for fe+3

df(3)=k4*c(2)*c(5)-k4*c(2)*c(3); % eqn for monoferric transferrin

df(4)=k4*c(2)*c(3); % eqn for holo transferrin

df(5)=-k4*c(2)*c(5); % eqn for transferrin


% oxygen constant at co2;   

% k6 is the removal of fe+3 ions by ppt etc

tdata=t;
oxid=df(1,1);
tfbind=df(3,1);
%fprintf('%6f\n ', t);
%fprintf('%6f\n', oxid);

%--------------------------------------------------------------------------



