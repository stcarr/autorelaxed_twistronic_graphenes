params = struct();

params.structure_type = 2; % 0: TBG, 1: Mono-on-bilayer, 2: alternating-twist trilayer
params.relax_strength = 1; % scaling of atomic relaxations (0.0 turns off relaxations)

params.AA_tunnel_strength = 1; % controls A-to-A interlayer tunneling strength (e.g. scales w0)
params.AB_tunnel_strength = 1; % rescales A-to-B interlayer tunneling strength (e.g. scales w1)
params.E_field = 0; % DOES NOT DO ANYTHING HERE   % vertical displacement field in eV (total potential energy across the layers)


params.theta = 1.53;

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
params.q_cut_type = 1;         % what kind of line-cut we do in momentum space (2 for supplied k-points)
params.nq = 30;  

params.qx_list = [0];
params.qy_list = [0];
params.ldos_locations = [];

[bands, q_scale] = twistronic_graphene_continuum_relax_bandcalc(params);

%%
clf
nb = size(bands,2);
Ef = bands(1,nb/2);
ax_m = 150; % axis max, in meV
plot(q_scale,1e3*(bands-Ef),'-k')
axis([-inf inf -ax_m ax_m])