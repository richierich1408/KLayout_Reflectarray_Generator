% Script Name: Set_Biasing.m
% Description: 
%     Sets the voltage values in a 1D array of BST capacitors to achieve beam steering.
%     The script calculates the required phase shifts for each unit cell based on the 
%     desired beam-steering angle and outputs the corresponding phase values.

% Author: David Hardy
% Date: 2025-03-27

clc; clear all; close all; 

title_font_size = 18;
set(groot, DefaultTextInterpreter="latex", DefaultTextFontSize=18);
set(groot, DefaultLegendInterpreter="latex", DefaultLegendFontSize=14);
set(groot, DefaultAxesTickLabelInterpreter="latex", DefaultAxesFontSize=14);

theta_r = 20;           % Beam-steering angle [degrees] 

% Set the path of the phase curve data and the output file path
raw_data = "12.5um";
phase_curve_data_path = sprintf("src/%s.csv",raw_data);  
output_lookup_table_path = sprintf("%s_%sdeg.csv",raw_data,num2str(theta_r)); 

%phase_curve_data_path = sprintf("Rectangular Patch/src/%s.csv",raw_data);  
%output_lookup_table_path = sprintf("Rectangular Patch/%s_%sdeg.csv",raw_data,num2str(theta_r)); 

% Constants
lambda0 = 10;           % Free-space wavelength [um]
k0 = 2*pi/lambda0;      % Free-space wavenumber [rad/um]

% Array Geometry
p_x = 5;                % Period of unit cell (x-direction) [um]
p_y = 2*p_x/sqrt(3);    % Period of unit cell (y-direction) [um]

% Phase gradient, in [deg]
phase_gradient = rad2deg(-k0*p_x*sind(theta_r));
disp("Phase Gradient: " + num2str(phase_gradient) + " [degrees]"); 

% Lookup table 
data = readmatrix(phase_curve_data_path); 
length_col = 1; 
phase_col = 2; 
length_vect = data(:,length_col);
phase_vect = data(:, phase_col);

% Normalize to start at 360, and decrease to 0
starting_phase = 360; 
phase_differential = starting_phase - phase_vect(1); 
for i=1:length(phase_vect)
    phase_vect(i) = phase_vect(i) + phase_differential; 
end

% "pchip" was found to be better for interpolating, but "spline" can be
%       used too
length_vect_interp = linspace(min(length_vect), max(length_vect), 1000); 
phase_vect_interp = interp1(length_vect, phase_vect, length_vect_interp, 'pchip'); 

% Plot
figure; 
plot(length_vect, phase_vect, LineWidth=2, DisplayName="Raw Data");
hold on; 
plot(length_vect_interp, phase_vect_interp, LineWidth=2, DisplayName="Interpolated", LineStyle="--");
hold off; 
xlabel("$\ell_\mathrm{rod}$ [$\mu$m]");
ylabel("Phase [$^\circ$]"); 
xlim([min(length_vect_interp), max(length_vect_interp)]); 
legend(Location="northeastoutside");
% Enable minor ticks
ax = gca; 
ax.XMinorTick = "on"; % Enable minor ticks on x-axis
ax.YMinorTick = "on"; % Enable minor ticks on y-axis
% Set major tick length (same for x and y in 2D)
ax.TickLength = [0.03, 0]; % [2D length, 3D length];
ax.Box = 'on'; 

% If we have a negative phase gradient, set the starting phase to 360
%       degrees. If positive, starting phase is 0.
if phase_gradient < 0 
    starting_phase = 360;
    starting_length = length_vect(1);
else
    starting_phase = phase_vect(end); 
    starting_length = length_vect(end); 
end

% We would like to find N, the number of cells which we can fit in a super
%       cell 
N_cells = floor((360 - min(phase_vect)) / abs(phase_gradient)) + 1;
disp("Number of cells per supercell: " + num2str(N_cells));

% Initialize 
l_rod = zeros(N_cells,1); 
ideal_phase = zeros(N_cells,1); 

l_rod(1) = starting_length; 
ideal_phase(1) = starting_phase; 

% Calculate the ideal phase of the n-th unit cell 
% Interpolate from the UC's phase curve to set the proper l_rod 
for n=2:N_cells 
    ideal_phase(n) = ideal_phase(n-1) + phase_gradient; 
    ideal_phase(n) = round(ideal_phase(n), 2); 

    % Interpolate to set the proper l_rod value 
    l_rod(n) = interp1(phase_vect_interp, length_vect_interp, ideal_phase(n), 'linear', 'extrap');  
end

% % Save V and Gamma to a CSV 
output_data = [l_rod, ideal_phase];
% csv_header = "n,Voltage,Phase"; 
fid = fopen(output_lookup_table_path, 'w');
fprintf(fid,"Length (um), Phase (deg)");
fclose(fid); 

writematrix(output_data, output_lookup_table_path, 'WriteMode', 'append');  

% Display confirmation
disp(['Data saved to ', output_lookup_table_path]);
disp(output_data);
