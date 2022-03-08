clc
clear
close all

% ENEE 692 - Introduction to Photonics
% Assignment #4
% Problem #1


% ********* Variables *********
% Question 1
% refractive index of the medium 1
n1 = 1.5;

% refractive index of the dielectric layer, medium 2
n2 = 3.5;

% Grating period: 

% --- d1 = thickness of medium 1 
% --- d2 = thickness of medium 2
d1 = 1.5*10^-6; % in meters
d2 = d1;


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


% ********* Variables *********
% % Question 1
% Bragg frequency  - p.275
mu_B = c0/(n_avg*(2*Lambda));

% frequency range;
freqValues = linspace(0,4*mu_B,freqRange);

% ******************

% eta
eta = (n1*d1 - n2*d2)/(n1*d1 + n2*d2);

% denominator
denom = 4 * n1 * n2;

% phi_1 + phi_2
phi1_p_phi2 = pi * (freqValues/mu_B);

% phi_1 - phi_2
phi1_m_phi2 = eta * phi1_p_phi2;

% ********* Wave Transfer Matrix *********

% attempt 2: M0 = [I21][e2][I12][e1]
A = 1/denom * (((n1+n2)^2)*exp(-1i*(phi1_p_phi2)) - ((n1-n2)^2)*exp(-1i*(phi1_m_phi2)));
B = 1/denom * (((n1^2) - (n2^2))*exp(1i*(phi1_p_phi2)) - ((n1^2) - (n2^2))*exp(1i*(phi1_m_phi2)));
C = conj(B);
D = conj(A);

% reflection 
r_small = B./D;

% transmission
t_small = 1./D;


% **************************
% ********* Part a *********

% ********* Finding Phi from eq. 7.1-46 *********
% cos(Phi) = Re{1/t} = Re(D) [from 7.1-43 and 7.1-46]
Phi = acos(real(1./t_small));


% ********* Finding Psi from eq. 7.1-46 *********
Psi_N = sin(N*Phi) ./ sin(Phi);

% ********* Finding R_one fom equation 6.2-14 *********
% R_one = |r_small|^2
% R_one = intensity of reflection for one frequency
R_one = (abs(r_small)).^2;


% ********* Finding R_N fom equation 7.1-51 *********
R_N = ((Psi_N.^2).*R_one) ./ (1 - R_one + ((Psi_N.^2).*R_one));


% **************************
% ********* Part b *********

% Step 1: loops through each values of A, B, C, and D to generate the wave
% transfer matrix (WTM)
% Step 2: finds the Nth power of the WTM
% Step 3: finds A_b, B_b, C_b, and D_b;
% Step 4: finds r
% Step 5: finds R_N

R_Q1b = zeros(1,freqRange);
for counter = 1:1:freqRange
    % Unit cell Wave matrix
    M_temp = [A(counter), B(counter); 
        C(counter), D(counter)];

    % Wave Transfer Matrix for each frequency
    M_Q1b = M_temp^10;


    % Reflectance
    r_small_Q2 = M_Q1b(1,2)/M_Q1b(2,2);
    R_Q1b(counter) = (abs(r_small_Q2)).^2;
end


% ********* Generating R_N vs. frequency *********
figure(1)

plot(freqValues/mu_B,R_N,'-b');
xlabel('Frequency, \nu');
ylabel('Reflectance intensity, \Re');
xticklabels({0,'','\nu_B','','2\nu_B','', '3\nu_B','', '4\nu_B'});
title('Question 1: using \Re=((\Psi_N^2)*R)/(1-R+((\Psi_N^2)*R))')

figure(2)
% part 1b
plot(freqValues/mu_B, R_Q1b,'-r');
xlabel('Frequency, \nu');
ylabel('Reflectance intensity, \Re');
xticklabels({0,'','\nu_B','','2\nu_B','', '3\nu_B','', '4\nu_B'});
title('Question 1b: using M_0^N')

figure(3)
plot(freqValues/mu_B,R_N,'-b');
hold on
plot(freqValues/mu_B, R_Q1b,'.r','MarkerSize',0.1);
hold off
xlabel('Frequency, \nu');
ylabel('Reflectance intensity, \Re');
xticklabels({0,'','\nu_B','','2\nu_B','', '3\nu_B','', '4\nu_B'});
title('Supoerimposing 1a and 1b')
legend('1a','1b');