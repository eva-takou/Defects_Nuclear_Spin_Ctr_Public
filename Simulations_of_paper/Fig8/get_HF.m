function [A,B,wL]=get_HF
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 24, 2023
%--------------------------------------------------------------------------

%HF parameters from: 
%Abobeih, M.H., Randall, J., Bradley, C.E. et al. 
%Atomic-scale imaging of a 27-nuclear-spin cluster using a quantum sensor. 
%Nature 576, 411–415 (2019).

A=[-20.72,...
   -23.22,...
   -31.25,...
   -14.07,...
   -11.346,...
   -48.58,...
   -8.32,...
   -9.79,...
   213.154,...
    17.643,...
    14.548,...
    20.569,...
    8.029,...
    -19.815,...
    -13.961,...
    -4.66,...
    -5.62,...
    -36.308,...
    24.399,...
    2.690,...
    1.212,...
    7.683,...
    -3.177,...
    -4.225,...
    -3.873,...
    -3.618,...
    -4.039];
B=[12,...
    13,...
    8,...
    13,...
    59.21,...
    9,...
    3,...
    5,...
    3,...
    8.6,...
    10,...
    41.51,...
    21,...
    5.3,...
    9,...
    7,...
    5,...
    26.62,...
    24.81,...
    11,...
    13,...
    4,...
    2,...
    0,...
    0,...
    0,...
    0];

B0 = 403; %in Gauss
B0 = B0*1e-4;

gamma_C13 = 6.728284 * 1e7; %rad/T*s
wL        = gamma_C13 * B0 ; 
wL        = wL*1e-6; %2.7 MHz
wL        = wL*1e3; %in KHz
wL        = wL/(2*pi);  %431.5484 kHz

end
