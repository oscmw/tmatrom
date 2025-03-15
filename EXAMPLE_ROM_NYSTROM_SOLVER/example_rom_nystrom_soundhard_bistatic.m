
% Copyright 2014, 2015, 2016, 2017, 2018, 2022, 2023, 2024 Stuart C. Hawkins and M. Ganesh.
% 	
% This file is part of TMATROM.
% 
% TMATROM is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TMATROM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TMATROM.  If not, see <http://www.gnu.org/licenses/>.

clear all
close all

%-----------------------------------------
% set main parameters
%-----------------------------------------

% wavenumber
kwave = 4*pi;

% origin
center = 0;

figure(1)
% Robin BC parameter (Sound hard case == 0)
robin_parameter = 0;

% set up scatterer object
scatterer = obstacleCircle();

%-----------------------------------------
% setup the solver
%-----------------------------------------

% setup solver with empty incident field for now
solver = solverNystromRobin(kwave,[],scatterer,robin_parameter);
solver.setup(60)

%-----------------------------------------
% derived parameters
%-----------------------------------------

% setup incident plane wave
p = plane_wave(0,kwave);

% get radius of scatterer
radius = solver.getRadius();

% get suggested order for the wavefunction expansion
nmax = suggestedorder(kwave,radius);

%-----------------------------------------
% setup the T-matrix
%-----------------------------------------

% setup T-matrix
tmat = ghtmatrix(nmax,kwave,solver,center);

% display T-matrix error check... based on Symmetry condition
fprintf('T-matrix error estimate %0.2e \n',tmat.error());

disp('Computed: Input/output independent  characterization of a sound-hard pinched ball')
%-----------------------------------------
% solve scattering problem
%-----------------------------------------

% create wave function expansion of plane wave
a = regularwavefunctionexpansion(nmax,center,p);

% compute wave function expansion of scattered wave using T-matrix
b = tmat * a;

disp('Computed: Output independent ROM object for the acoustic model')

disp('Evaluating and Visualizing (Fig. 1):')
disp ('Output total field at 25,000 grid points (Fig. 1)')
%-----------------------------------------
% visualize the total field
%-----------------------------------------
% visualize the total field
figure(1)
% setup a grid
mesh = load('mesh.txt');  
%t = linspace(-2,2,100);
%[x, y] = meshgrid(t, t);
x = mesh(:,1);
y = mesh(:,2);
z = x + y*1i;

% get a mask for the scatterer
mask = abs(z-center) > 0.999999*radius;

% compute the total field
total_field = real(b.evaluate(z,mask) + p.evaluate(z,mask));
allField = [x(:), y(:), total_field(:)];

fid = fopen('meshField.txt','w');
for i = 1:size(allField, 1)
    fprintf(fid,'%8.6f  %8.6f  %12.8f\n', allField(i,:));
end
fclose(fid);


% save data to a .mat file
%save('total_field_data.mat', 'x', 'y', 'total_field');
%text dosyası olarak kaydet


% plot
surf(x, y, total_field);
view([0 90]);
shading interp;
colorbar
title('Total field (in the plane [-10,10]x[-10,10]) exterior to a pinched-ball sound-hard scatterer')


% add the scatterer
hold on
solver.visualize()
hold off


%-----------------------------------------
% visualize the far field
%-----------------------------------------

% setup points on the circle
theta = linspace(0, 2*pi, 1000);
z = exp(1i*theta);

% compute bistatic ACS
bistatic_acs = 10*log10(2*pi*abs(b.evaluateFarField(z)).^2);

% save data to a .mat file
save('bistatic_acs_data.mat', 'theta', 'bistatic_acs');

% plot the cross section
figure(2)
plot(theta, bistatic_acs, 'r-')
xlabel('Receiver direction angles')
ylabel('Bistatic ACS (dB)')



%--------------------------------------
% Noktasal veriyi elde etmek
%--------------------------------------


