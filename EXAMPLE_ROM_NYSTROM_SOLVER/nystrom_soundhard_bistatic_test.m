clear all
close all

%-----------------------------------------
% Set main parameters
%-----------------------------------------
kwave = 4*pi;           % Wavenumber
center = 0;             % Origin
robin_parameter = 0;    % Sound-hard boundary condition

%-----------------------------------------
% Set up scatterer and solver
%-----------------------------------------
scatterer = obstacleKite();
%qx = scatterer(:,3)
solver = solverNystromRobin(kwave, [], scatterer, robin_parameter);
solver.setup(60);

%-----------------------------------------
% Incident field and T-matrix setup
%-----------------------------------------
p = plane_wave(0, kwave);               % Incident plane wave
radius = solver.getRadius();            % Get scatterer radius
nmax = suggestedorder(kwave, radius);  % Expansion order
tmat = ghtmatrix(nmax, kwave, solver, center);  % T-matrix

%-----------------------------------------
% Compute scattered field coefficients
%-----------------------------------------
a = regularwavefunctionexpansion(nmax, center, p);
b = tmat * a;

%-----------------------------------------
% Create grid and kite-shaped mask
%-----------------------------------------

%setup a grid
mesh = load('meshPython.txt')
x = mesh(:,1);
y = mesh(:,2);

t = linspace(-5, 5, 500);  % Higher resolution
[x, y] = meshgrid(t, t);
z = x + 1i*y;

% Get high-resolution kite boundary
t_kite = linspace(0, 2*pi, 5000);
[~, ~, qx, qy] = scatterer.geom(t_kite);


% Create precise kite mask with small buffer
kite_mask = ~inpolygon(x, y, 0.9995*qx, 0.9995*qy);  

%-----------------------------------------
% Field calculation using ONLY kite mask
%-----------------------------------------

total_field = real(b.evaluate(z, kite_mask) + p.evaluate(z, kite_mask));
total_field(~kite_mask) = NaN;  % Explicitly mask interior points

save('total_field_matlab.mat', 'total_field');

allField = [x(:), y(:), total_field(:)];

fid = fopen('meshField_kite.txt','w');
for i = 1:size(allField, 1)
    fprintf(fid,'%8.6f  %8.6f  %12.8f\n', allField(i,:));
end
fclose(fid);

%-----------------------------------------
% Plot results
%-----------------------------------------
figure(1)
surf(x, y, total_field, 'EdgeColor', 'none');
view(0, 90);
axis equal tight;
shading interp;
colorbar;
caxis([-2 2]);
title('Total Field with Kite-Shaped Scatterer');

% Overlay exact kite boundary
hold on;
plot(qx, qy, 'k-', 'LineWidth', 1.5);  % Direct plot of boundary points
hold off;

%-----------------------------------------
% Far-field pattern
%-----------------------------------------
theta = linspace(0, 2*pi, 1000);
z_ff = exp(1i*theta);
bistatic_acs = 10*log10(2*pi*abs(b.evaluateFarField(z_ff)).^2);

figure(2)
plot(theta, bistatic_acs, 'r-');
xlabel('Angle (rad)');
ylabel('Bistatic ACS (dB)');
title('Far-Field Pattern');
