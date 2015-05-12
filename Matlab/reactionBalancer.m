%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear all
clc; clear all;

%%% Import the file
model = readCbModel('TBepimetabat2973/TBepimetabat2973Balanced.xml');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check the mass- and charge-balancing
[massImbalance,imBalancedMass,imBalancedCharge,imBalancedBool] = checkMassChargeBalance(model);