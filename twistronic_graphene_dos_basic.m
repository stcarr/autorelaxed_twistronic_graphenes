params = struct();

params.structure_type = 2; % 0: TBG, 1: Mono-on-bilayer, 2: alternating-twist trilayer
params.relax_strength = 0; % scaling of atomic relaxations (0.0 turns off relaxations)

% For A-T-A (alternating twist trilayer), the AA stacking (r = 0) structure corresponds to:
% L3: o-o--o-o 
% L2: o-o--o-o (<-- motion along A1+A2 direction corresponds to this layer moving LEFT)
% L1: o-o--o-o


% For mono-on-bilayer, the AA-stacking (r = 0) correpsonds to the structure:
% L3: o-o--o-o (<-- motion along A1+A2 direction corresponds to this layer moving LEFT)
% L2: o-o--o-o
% L1: --o-o--o
% i.e. the A orbital of L1 is coupled to the B orbital of L2.

params.AA_tunnel_strength = 1; % controls A-to-A interlayer tunneling strength (e.g. scales w0)
params.AB_tunnel_strength = 1; % rescales A-to-B interlayer tunneling strength (e.g. scales w1)
params.E_field = 0; % DOES NOT DO ANYTHING HERE   % vertical displacement field in eV (total potential energy across the layers)


params.theta = 2.0;

% truncation of momentum basis
params.k_cutoff = 0.25;         % momentum cutoff, in units of inverse Angstroms, translates to an energy cutoff of 0.15*v_D ~
params.grid_search = 30;       % momentum < cutoff found within [-grid_search,grid_search]^2 

% interlayer coupling sampling mesh
params.r_max = 6; % maximum radius, in Angstroms
params.dr = 0.2; % spacing of r mesh, in Angstroms

% interlayer fourier transform control
params.inter_q_cut_scale = 10; % maximum scattered momentum, in units of K0
params.inner_k_rad_scale = 10; % radius of each "island", in units of b12
params.dk_scale = 2; %  spacing of k mesh, in units of b12

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

q_spacing = 0.025; % make sure 0.5 is divisible by this number, or we will miss the K-point!
q_arr = 0.5 - [0:q_spacing:(1-q_spacing)]; % K centered mesh
[mesh_i,mesh_j] = meshgrid(q_arr,q_arr);
qmesh_x = mesh_i*Gm1(1) + mesh_j*Gm2(1);
qmesh_y = mesh_i*Gm1(2) + mesh_j*Gm2(2);

params.qx_list = qmesh_x(:);
params.qy_list = qmesh_y(:);

A = 2*pi*inv(b12)';
A1 = A(:,1);
A2 = A(:,2);

locs = [0,1/3,1/2,2/3].*(A1+A2);
% AA, AB, DW (SP), BA
% for monolayer on bilayer, AB = ABC, BA = ABA

params.ldos_locations = locs;
[bands, ~, weights] = twistronic_graphene_continuum_relax_bandcalc(params);

%%
clf

ax_m = 500; % axis max, in meV
plot(1e3*(bands),'-k')
axis([-inf inf -ax_m ax_m])
%%
kpts = [qmesh_x(:),qmesh_y(:);];

nb = size(bands,2);
b_size = 30; % include this many bands above and below the Fermi energy
b_size = min(b_size, nb/2 - 1); % b_size needs to be smaller than half the # of bands
b_center = nb/2 - 0;
max_E = 0.5; % maximum energy in the DOS window, in eV
dE = 1e-3; % DOS window spacing, in eV
sig = 3*dE; % gaussian smearing

%[dos, idos, E_list] = interp_kp_dos(params.theta, bands, kpts, b_size, max_E, dE);
[dos_gauss, idos_gauss, E_list_gauss] = interp_kp_dos_gaussian(params.theta, bands, kpts, b_size, max_E, dE, sig);
[dos, ldos, E_list] = interp_kp_ldos(params.theta, bands, weights, kpts, b_center, b_size, max_E, dE);
%%
clf
hold on
box on
plot(E_list,dos,'k','DisplayName','Interp method')
plot(E_list_gauss,dos_gauss,'-r','DisplayName','Gauss method')
xlabel('Energy (eV)')
ylabel('DOS, states per eV per nm^2')
axis([-max_E max_E 0 7])
legend()
%%
clf
hold on
% ldos is indexed as [Energy, orbital index, location index]
% there are two orbitals per layer, e.g.
% ldos(:,1,:) corresponds to the bottom layer's A orbital
% ldos(:,2,:) corresponds to the bottom layer's B orbital
% ldos(:,3,:) correpsonds to the second layer's A orbital, etc.

plot(E_list,(ldos(:,1,1)),'-r','DisplayName','AA L1A LDOS')
plot(E_list,(ldos(:,3,1)),'-b','DisplayName','AA L2A LDOS')
plot(E_list,(ldos(:,1,2)),'-g','DisplayName','AB L1A LDOS')
plot(E_list,(ldos(:,3,2)),'color',[0.5,0.5,0],'DisplayName','AB L2A LDOS')

plot(E_list,(dos),'-k','DisplayName','DOS')

box on
axis([-max_E max_E 0 1.2*max(ldos(:))])
xlabel('Energy (eV)')
ylabel('LDOS (states per eV per A^2)')
legend()

