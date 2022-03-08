clc
clear
close all

% refractive index of the medium 1
n1 = 1.45;

% refractive index of the dielectric layer, medium 2
n2 = 2.25;

% --- d1 = thickness of medium 1 
% --- d2 = thickness of medium 2
d1 = 259*10^-9; % in meters
d2 = 167*10^-9; % in meters

% ******************
% Number of segments
N = 10;

% ** how many frequency values do you want?
freqRange = 100000;

% ********* Constants *********

% speed of light in vacuum
c0 = 3*10^8; % in m/s

% ********* Derivations *********


% Grating period
Lambda = d1 + d2;

% average refractive index
n_avg = (n1*d1 + n2*d2)/(Lambda);

% Bragg frequency  - p.275
% mu_B = c0/(n_avg*(2*Lambda));
mu_B = c0/(n1*(1.5*10^-6));


% frequency range;
freqValues = linspace(-2*mu_B,2*mu_B,freqRange);


% Interface
I12 = 1/(2*n2) * [(n2+n1),(n2-n1);
        (n2-n1),(n2+n1)];
I21 = 1/(2*n1) * [(n1+n2),(n1-n2);
        (n1-n2),(n1+n2)];
    

R_Q2 = zeros(1,freqRange);
for count = 1:1:freqRange
    
    % phase, Phi
    Phi_1 = (2*pi*n1*freqValues(count)/c0)*(n1 * d1);
    Phi_2 = (2*pi*n2*freqValues(count)/c0)*(n2 * d2);
    

    e_1 = [exp(-1i*(Phi_1)), 0; 
        0, exp(1i*(Phi_1))];
    
    e_2 = [exp(-1i*(Phi_1)), 0; 
        0, exp(1i*(Phi_1))];
    
    
        % Unit cell Wave matrix
    M_temp = [I21]*[e_2]*[I12]*[e_1];
    
   
    % Wave Transfer Matrix for each frequency
    M = (M_temp^N) * e_1 * (M_temp^N);

    % Reflectance
    r_small_Q2 = M(1,2)/M(2,2);
    R_Q2(count) = (abs(r_small_Q2)).^2;
end

figure(1)
plot(freqValues/mu_B,R_Q2,'-b');
xlabel('Frequency, \nu');
ylabel('Reflectance intensity, \Re');
xticklabels({'-2\nu_B','','-\nu_B','',0,'', '\nu_B','', '2\nu_B'});
title('Question 2:')
