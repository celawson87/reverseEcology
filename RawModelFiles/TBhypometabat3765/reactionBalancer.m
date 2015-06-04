%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear all
clc; clear all;

%%% Import the file
model = readCbModel('TBhypometabat3765Balanced');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check the mass- and charge-balancing
[massImbalance,imBalancedMass,imBalancedCharge,imBalancedBool] = checkMassChargeBalance(model);