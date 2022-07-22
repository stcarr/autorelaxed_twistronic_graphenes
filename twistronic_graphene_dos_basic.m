params = struct();

params.structure_type = 0; % 0: TBG, 1: Mono-on-bilayer, 2: alternating-twist trilayer
params.relax_strength = 1; % scaling of atomic relaxations (0.0 turns off relaxations)

params.AA_tunnel_strength = 1; % controls A-to-A interlayer tunneling strength (e.g. scales w0)
params.AB_tunnel_strength = 1; % rescales A-to-B interlayer tunneling strength (e.g. scales w1)
params.E_field = 0; % DOES NOT DO ANYTHING HERE   % vertical displacement field in eV (total potential energy across the layers)


params.theta = 1.2;

% truncation of momentum basis
params.k_cutoff = 0.15;         % momentum cutoff, in units of inverse Angstroms, translates to an energy cutoff of 0.15*v_D ~
params.grid_search = 30;       % momentum < cutoff found within [-grid_search,grid_search]^2 

% interlayer coupling sampling mesh
params.r_max = 8; % maximum radius, in Angstroms
params.dr = 0.2; % spacing of r mesh, in Angstroms

% interlayer fourier transform control
params.inter_q_cut_scale = 5; % maximum scattered momentum, in units of K0
params.inner_k_rad_scale = 5; % radius of each "island", in units of b12
params.dk_scale = 0.5; %  spacing of k mesh, in units of b12

% band structure line cut settings
params.q_cut_type = 2;         % what kind of line-cut we do in momentum space (2 for supplied k-points)
params.nq = 30;  

% create layer data structures
theta_list = params.theta*[-0.5, 0.5];  % twist of each layer
for t = 1:2
   layers(t) = Layer(t,deg2rad(theta_list(t)));
end

 % reciprocal geometry
G1 = layers(1).G; 
G2 = layers(2).G;
b12 = G2 - G1; % moire reciprocal vectors

Gm1 =  [b12(1,1);-b12(2,1)];
Gm2 = -b12(:,2);

q_spacing = 0.1;
q_arr = 0:q_spacing:(1-q_spacing);
[mesh_i,mesh_j] = meshgrid(q_arr,q_arr);
qmesh_x = mesh_i*Gm1(1) + mesh_j*Gm2(1);
qmesh_y = mesh_i*Gm1(2) + mesh_j*Gm2(2);

params.qx_list = qmesh_x(:);
params.qy_list = qmesh_y(:);

A = 2*pi*inv(b12)';
A1 = A(:,1);
A2 = A(:,2);

locs = [0,1/3,1/2].*(A1+A2);
% AA, AB, and DW (SP)

params.ldos_locations = locs;
[bands, ~, weights] = twistronic_graphene_continuum_relax_bandcalc(params);

%%
clf

ax_m = 150; % axis max, in meV
plot(1e3*(bands),'-k')
axis([-inf inf -ax_m ax_m])
%%
kpts = [qmesh_x(:),qmesh_y(:);];

nb = size(bands,2);
b_size = 30; % include this many bands above and below the Fermi energy
b_size = min(b_size, nb/2 - 1); % b_size needs to be smaller than half the # of bands
max_E = 0.6; % maximum energy in the DOS window, in eV
dE = 1e-3; % DOS window spacing, in eV
%[dos, idos, E_list] = interp_kp_dos(params.theta, bands, kpts, b_size, max_E, dE);
[dos_gauss, idos_gauss, E_list_gauss] = interp_kp_dos_gaussian(params.theta, bands, kpts, b_size, max_E, dE);
[dos, ldos, E_list] = interp_kp_ldos(params.theta, bands, weights, kpts, b_size, max_E, dE);
%%
clf
hold on
box on
plot(E_list,log10(dos),'k','DisplayName','Interp method')
plot(E_list_gauss,log10(dos_gauss),'r','DisplayName','Gauss method')
xlabel('Energy (eV)')
ylabel('DOS (log_{10}, states per eV per A^2')
axis([-max_E max_E -2 2])
legend()
%%
clf
hold on
plot(E_list,log10(ldos(:,1)),'-r','DisplayName','AA LDOS')
plot(E_list,log10(ldos(:,2)),'-b','DisplayName','AB LDOS')
plot(E_list,log10(ldos(:,3)),'-g','DisplayName','SP LDOS')
plot(E_list,log10(dos),'-k','DisplayName','DOS')

box on
axis([-max_E max_E -1 2.5])
xlabel('Energy (eV)')
ylabel('LDOS (log_{10}, states per eV per A^2)')
legend()

