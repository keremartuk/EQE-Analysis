clear all;
close all;

% EQE Analysis

Jsc = [];
Urbach_energy = [];
bandgap_pero = [];

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
EQE_pero = data(:,2);
EQE_si = data (:,3);
wavelength = data (:,1);
%Jsc loss

q = 1.602 * 10^-19;
h = 6.62607004 * 1e-34 ; 
c = 3 * 10^8 ;                  

irr = importdata('Irradiance.txt');
irr_new = irr (:,2);

wavelength_new = irr(:,1).*10e8;
EQE_new_pero = interp1(wavelength, EQE_pero, wavelength_new,'linear','extrap');
EQE_new_si = interp1(wavelength, EQE_si, wavelength_new,'linear','extrap');

Jsc_pero = (q/(h*c))*trapz(wavelength_new,(EQE_new_pero.*wavelength_new.*irr_new))./(10*10e17)
Jsc_si = (q/(h*c))*trapz(wavelength_new,(EQE_new_si.*wavelength_new.*irr_new))./(10*10e17)+0.59

J_mismatch = Jsc_pero - Jsc_si % negative in top limited, positive in bottom limited

% extract artificial zeros

for i=1:1:length(EQE_si)

if EQE_si(i) == 0

    EQE_si(i) = nan;

end

if EQE_pero(i) == 0

    EQE_pero(i) = nan;

end

end

% filenames

kerem = figure(1);

plot (wavelength, EQE_pero.*100,'MarkerSize',3,'Marker','o','LineWidth',3,'Color',[0 0 0])

hold all

plot (wavelength, EQE_si.*100,'MarkerSize',3,'Marker','o','LineWidth',3,'Color',[1 0 0])

xlim ([310,1180])
ylim ([0,100])

xlabel ('Wavelength (nm)')
ylabel ('EQE(%)')

saveas(kerem, 'EQE_figure_tandem.png');

% urbach pero

Urbach_energy = [Urbach_energy; 15]; % mV calculate

% Bandgap calculation

derivative_EQE = diff(EQE_pero)/10; % step size is 10 nm
energy = (1240./wavelength);  

index = find(derivative_EQE == min(derivative_EQE));
bandgap_pero = [bandgap_pero; energy(index+1)]; % eV calculate

%figure (2)

%plot((energy(2:length(derivative_EQE)+1)),diff(EQE_pero)/10)

%hold all


%writetable(x, strcat(Folder_name,'\','all_analysed.txt'));

end

x = table(names, Jsc_pero, Jsc_si, J_mismatch, bandgap_pero)
