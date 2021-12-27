%% EQE Analysis

clear all;
close all;
clc

%CONSTANTS

k = 1.38064852e-23;                     % [(m^2)kg(s^-2)(K^-1)], Boltzmann constant
T = 300;                                % [K], temperature
q = 1.602 * 10^-19;
VT= 25.8e-3; %(k*T)/q;                             % [V], 25.8mV thermal voltage at 300K
h = 6.62607004 * 1e-34 ; 
c = 3 * 10^8 ;                  

Vocs =[];
Jsc = [];
Urbach_energy = [];
bandgap_EQE = [];
J0s = [];

%% Analyse EQE

% plot EQE
names = '';
file_names = '';

headerlinesIn = 8;

Folder_name = 'EQE_tandem';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Folder = dir(strcat(Folder_name,'\','*.xls'));

%%

for i1 = 1:numel(Folder)
    
filename = Folder(i1).name;
names = [names;{filename(1:end-4)}];
delimiterIn = '\t';

data = xlsread (strcat(Folder_name,'\',filename));
EQE = data(:,2);
wavelength = data (:,1);

%
irr = importdata('Irradiance.txt');
irr_new = irr (:,2);

wavelength_new = irr(:,1).*10e8;
EQE_new_pero = interp1(wavelength, EQE, wavelength_new,'linear','extrap');
Jsc_pero = (q/(h*c))*trapz(wavelength_new,(EQE_new_pero.*wavelength_new.*irr_new))./(10*10e17);
Jsc = [Jsc; Jsc_pero];

% filenames
kerem = figure (1);

plot (wavelength, EQE,'MarkerSize',3,'Marker','o','LineWidth',3,'Color',[i1/numel(Folder) 0 0])

hold all

xlim ([320,850])
ylim ([0,1])

xlabel ('Wavelength (nm)')
ylabel ('EQE')
saveas(kerem, 'EQE_figure_single_junction.png');

% Jdark

E = flip(1240./wavelength);
L = wavelength_new;
BB_1 = 2*pi*c*c*h./(L.*L.*L.*L.*L.*L);
BB_2 = 1./(exp((h*c/k*T)./L-1));%black body spectrum BB (s^-1cm^-2(eV)^-1)(in photon flux!)
BB = BB_1.*BB_2;

%% from CMWs code

bb2 = q*2.*pi./(4.1356e-15.^3*3e8.^2).*E.^2./(exp(E./0.0259)-1);
j_0_rad = trapz(E, flip(EQE).*bb2);

%%

Voc_rad=(VT)*log((Jsc./j_0_rad)+1); %open-circuit voltage Voc (V)
voc=Voc_rad/VT;  
FF_rad=(voc-log(voc+0.72))./(voc+1);                          %fill factor FF
eff_rad=100*(Jsc.*Voc_rad.*FF_rad)./(100);                 %efficieny in (%)

%J0 = (q/(h*c))*trapz(wavelength_new,(EQE_new_pero.*wavelength_new.*BB))./(100);
J0s = [J0s; j_0_rad];

%calculate radiative voltage and estimate deficit with a known Vocs

% Bandgap calculation

derivative_EQE = diff(EQE)/10; % step size is 10 nm
energy = (1240./wavelength);

index = find(derivative_EQE == min(derivative_EQE));
bandgap_EQE = [bandgap_EQE; energy(index+1)]; % eV calculate

%figure (2)

%plot((energy(2:length(derivative_EQE)+1)),diff(EQE)/10)

%hold all

Urbach_energy = [Urbach_energy; 15]; % mV calculate

%writetable(x, strcat(Folder_name,'\','all_analysed.txt'));

end

x = table(names, Jsc, Urbach_energy,bandgap_EQE, J0s)
