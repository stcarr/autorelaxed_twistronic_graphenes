% authors: ziyan (zoe) zhu, stephen carr 
% email: zzhu1@g.harvard.edu
% Example calculation of the TBG band structure
clear all
f_size = 12;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)

addpath('..')

% !! Settings

% turn on/off saving
savedata = false;           % save useful variables to folder ./data/
fname_suffix = 'test';
data_folder = 'test_data/';

relax_str = 1.0; % scaling of atomic relaxations (0.0 turns off relaxations)
E_field = 0.0;   % vertical displacement field in eV (total potential energy across the three layers)
AA_str = 1; % controls AA interlayer tunneling strength (e.g. scales w0)
AB_str = 1; % rescales AB interlayer tunneling strength (e.g. scales w1)

theta = 0.3;

% truncation of momentum basis
k_cutoff = 4;           % momentum cutoff, in units of norm(G1)
grid_search = 20;       % [G1,G2] are in [-grid_search,grid_search]^2 

% interlayer coupling sampling mesh
r_max = 8; % maximum radius, in Angstroms
dr = .25; % spacing of r mesh, in Angstroms
% interlayer fourier transform control
inter_q_cut_scale = 5; % maximum scattered momentum, in units of K0
inner_k_rad_scale = 5; % radius of each "island", in units of b12
dk_scale = 1/2; %  spacing of k mesh, in units of b12

% band structure line cut settings
q_cut_type = 1;         % what kind of line-cut we do in momentum space                        % 1: high symmetry line in the L12 bilayer moire Brillouin zone (single valley)
nq = 50;                % number of k points to sample on each high symmetry line segment 

theta_list = theta*[-0.5, 0.5];  % twist of each layer
% create layer data structures
for t = 1:2
   layers(t) = Layer(t,deg2rad(theta_list(t)));
end

% generate orbital positions
nsheets = 2;
norbs = 2;
orb_pos = zeros(nsheets,norbs,2);
A1 = layers(1).A;
A2 = layers(2).A;
orb_pos(1,:,:) = layers(1).orbPos;
orb_pos(2,:,:) = layers(2).orbPos;

% reciprocal geometry
G1 = layers(1).G; 
G2 = layers(2).G;
b12 = G2 - G1; % moire reciprocal vectors

G11 = G1(:, 1);
G12 = G1(:, 2);
G13 = G11 + G12; 
K0 = 1/3 * (G11 + G13); % a K point of Layer 1

inter_q_cut = inter_q_cut_scale*norm(K0); % maximum scattered momentum
inner_k_rad = inner_k_rad_scale*norm(b12(:,1)); % size of scattering "island"
dk = dk_scale*norm(b12(:,1)); % spacing of k-mesh  


% get the K-point of each monolayer
for t = 1:2
   th = layers(t).theta - layers(1).theta;
   K(t,:) = [cos(th) -sin(th); sin(th) cos(th)]*K0;
end

% define the k-point sampling
samp = linspace(0,1,nq)';
samp = samp(1:end-1);

d = (1.0/2.0)*sqrt(sum(K(2,:) - K(1,:).^2)); 

K_1 = K(1,:);
K_2 = K(2,:);
M = (K(1, :) + K(2, :))/2; % M point of supercell

rot120 = [cos(2*pi/3), sin(2*pi/3); -sin(2*pi/3), cos(2*pi/3)];

% definition of q vectors (NN coupling separations in k space)
% these are used in the BMD model, not in the DFT model.
q1_12 = K_2'-K_1';
q2_12 = rot120 * q1_12;
q3_12 = rot120 * q2_12; 

% defining high symmetry points
switch q_cut_type 
    case 1 % K-Gamma-M of the L12 supercell
        k_sc = q3_12;
        m_sc = 0.5*(q3_12-q2_12);
        gamma_sc = [0, 0];

        pt(1,:) = k_sc;
        pt(2,:) = gamma_sc;
        pt(3,:) = m_sc;
        pt(4,:) = k_sc;
        type = 'L12 supercell';
        xt_labels = {'$K_{12}$', '$\Gamma_{12}$', '$M_{12}$', '$K_{12}$'};
end

max_seg = size(pt,1)-1;

% making the line segments through the selected points
for seg = 1:max_seg
    if seg == 1
        q_list_x = pt(1,1)*(1-samp) + pt(2,1)*samp;
        q_list_y = pt(1,2)*(1-samp) + pt(2,2)*samp;
    else
        q_list_x = [q_list_x; pt(seg,1)*(1-samp) + pt(seg+1,1)*samp];
        q_list_y = [q_list_y; pt(seg,2)*(1-samp) + pt(seg+1,2)*samp];
    end  
end

q_list = [q_list_x, q_list_y]; % the list of center sites

q_list(end+1, :) = [pt(max_seg+1,1), pt(max_seg+1,2)];

q_list2 = -q_list - k_sc';
q_list = q_list - k_sc';
% calculate the path length at q point 
ni(1) = 1; 
for i = 2:max_seg+1 
    ni(i) = (i-1)*size(samp,1)+1;
end 
%ni(max_seg+1) = ni(max_seg+1)+1;

for p_idx = 1:max_seg
    dis_here = norm(pt(p_idx+1,:)-pt(p_idx,:));
    if p_idx == 1
        qarr = linspace(0,dis_here,ni(p_idx+1)-ni(p_idx)+1);
    else 
        qarr(ni(p_idx)+1:ni(p_idx+1)) = linspace(qarr(length(qarr))+dis_here/(ni(p_idx+1)-ni(p_idx)), ...
            qarr(length(qarr))+dis_here,ni(p_idx+1)-ni(p_idx));
    end 
end 

for i = 1:length(ni)
    xt(i) = qarr(ni(i));
end 

% Generate momentum lattice
% setup kp model (output the list of scattered k's)
dof_list = kDoF_bi(layers,k_cutoff,grid_search);

% generate structure
dof_list.gen_dof()

% get the kpoints
k_list = dof_list.k_list();
ndof = size(k_list,1);
fprintf("%d total k points \n",ndof)

% find the indices corresponding to the center site (of K lattice, not H)
layer_list = k_list(:, 5);
k1 = k_list(layer_list == 1, :);
k2 = k_list(layer_list == 2, :);
for k_idx = 1:size(k_list,1)
    for tar_sheet = 1:2
        if min( k_list(k_idx,5:end) == [tar_sheet 0 0 0 0])
           tar_dofs(tar_sheet, :) = k_idx;
        end
    end 
end


% generate interlayercouplings object (used for interlayer Hamiltonian)

% create R and K meshes
intercoupling = InterCouplings(layers, 1, 2, r_max, dr, inter_q_cut, dk, inner_k_rad, K0, theta, orb_pos);
intercoupling.gen_tR();

% generate interlayercouplings object, but for full kmesh
intercoupling_plot = InterCouplings(layers, 1, 2, r_max, dr, inter_q_cut, dk, inner_k_rad, K0, theta, orb_pos);
intercoupling_plot.tR = intercoupling.tR;
dk_grid = linspace(-8,8,100);
[Kx_grid, Ky_grid] = meshgrid(dk_grid,dk_grid);
intercoupling_plot.K_meshx = Kx_grid(:);
intercoupling_plot.K_meshy = Ky_grid(:);
intercoupling_plot.gen_tK(AA_str,AB_str);
tK_plots_unrelax = intercoupling_plot.tK;
Rx = intercoupling_plot.R_meshx;
Ry = intercoupling_plot.R_meshy;
tR_plots_unrelax = intercoupling_plot.tR;

% apply relaxations
tic
fprintf("Applying atomic relaxations... \n")
intercoupling.relax_configs(relax_str);
fprintf("Relaxations done: ")
toc

intercoupling.gen_tR();
% generate interlayercouplings object, but for full kmesh
intercoupling_plot = InterCouplings(layers, 1, 2, r_max, dr, inter_q_cut, dk, inner_k_rad, K0, theta, orb_pos);
intercoupling_plot.tR = intercoupling.tR;
dk_grid = linspace(-8,8,100);
[Kx_grid, Ky_grid] = meshgrid(dk_grid,dk_grid);
intercoupling_plot.K_meshx = Kx_grid(:);
intercoupling_plot.K_meshy = Ky_grid(:);
intercoupling_plot.gen_tK(AA_str,AB_str);
tK_plots_relax = intercoupling_plot.tK;
Rx = intercoupling_plot.R_meshx;
Ry = intercoupling_plot.R_meshy;
tR_plots_relax = intercoupling_plot.tR;

%% plot the interlayer coupling, in both configuration (b) and momentum (q) space
close all;
fig = figure('Position',[0 200 800 400]);

surf_on = true;
text_on = true;%~surf_on;

T_b(1,:) = [.66 .66 1];
T_b(2,:) = [1 1 1];
T_b(3,:) = [1 0 0];

T_q(1,:) = [1 1 1];
T_q(2,:) = [1 .2 0];
T_q(3,:) = [1 1 0];

x = 0.1 + [-0.1
            0.0 
            0.3];
 
cmap_b = interp1(x/x(end),T_b,linspace(0,1,255));
cmap_q = interp1(x/x(end),T_q,linspace(0,1,255));

o1 = 1;
splot_idx = 1;

unfreezeColors()
for o2 = 1:2
    if o2 == 1
       o2_str = 'A'; 
    else
       o2_str = 'B';
    end

    tR_relax = squeeze(tR_plots_relax(o1,o2,:));
    tR_unrelax = squeeze(tR_plots_unrelax(o1,o2,:));

    x = intercoupling.R_meshx;
    y = intercoupling.R_meshy;

    ax_m = 5;
    dx = linspace(-ax_m,ax_m,100);
    [mx,my] = meshgrid(dx,dx);
    m_tR_relax = griddata(x,y,tR_relax,mx,my);
    m_tR_unrelax = griddata(x,y,tR_unrelax,mx,my);

    m_tK_relax = reshape(tK_plots_relax(o1,o2,:),100,100);
    m_tK_unrelax = reshape(tK_plots_unrelax(o1,o2,:),100,100);

    splot(splot_idx) = subplot(2,4,splot_idx);
    splot_idx = splot_idx+1;
    hold on
    if surf_on
        surf(mx,my,m_tR_unrelax)
    end
    shading interp
    axis equal
    axis([-5 5 -5 5 -inf inf])
    view(2)
    caxis([-.1 .3])
    set(gca,'Xtick',[])
    str_here = ['${h}_{12}^{A' o2_str '}(b)$ unrelax'];
    if text_on
        text(4.6,4.4,1,str_here,'color','k','HorizontalAlignment','right')
    end
    colormap(cmap_b)
    if o2 == 1
       cbar_b = colorbar('TickLabelInterpreter', 'latex');
    end
    freezeColors()
    
    splot(splot_idx) = subplot(2,4,splot_idx);
    splot_idx = splot_idx+1;
    hold on
    if surf_on
        surf(Kx_grid,Ky_grid,abs(m_tK_unrelax))
    end
    shading interp
    axis equal
    axis([-7 7 -7 7 -inf inf])
    view(2)
    caxis([0 .3])    
    set(gca,'Ytick',[])
    set(gca,'Xtick',[])
    str_here = ['$\tilde{h}_{12}^{A' o2_str '}(q)$ unrelax'];
    if text_on
        text(6.4,6.1,1,str_here,'color','k','HorizontalAlignment','right')
    end
    colormap(cmap_q)
    if o2 == 1
       cbar_q = colorbar('TickLabelInterpreter', 'latex');
    end
    
    freezeColors()
    
    
    splot(splot_idx) = subplot(2,4,splot_idx);
    splot_idx = splot_idx+1;
    hold on
    if surf_on
        surf(mx,my,m_tR_relax)
    end
    shading interp
    axis equal
    axis([-5 5 -5 5 -inf inf])
    view(2)
    caxis([-.1 .3])
    xlabel('$b_x$ (\AA)')
    str_here = ['${h}_{12}^{A' o2_str '}(b)$ relax'];
    if text_on
        text(4.6,4.4,1,str_here,'color','k','HorizontalAlignment','right')
    end
    colormap(cmap_b)
    freezeColors()
    
    splot(splot_idx) = subplot(2,4,splot_idx);
    splot_idx = splot_idx+1;
    hold on
    if surf_on
        surf(Kx_grid,Ky_grid,abs(m_tK_relax))
    end
    shading interp
    axis equal
    box on
    axis([-7 7 -7 7 -inf inf])
    view(2)
    caxis([0 .3])
    set(gca,'Ytick',[])
    xlabel('$q_x$ (\AA$^{-1}$)')
    str_here = ['$\tilde{h}_{12}^{A' o2_str '}(q)$ relax'];
    if text_on
        text(6.4,6.1,1,str_here,'color','k','HorizontalAlignment','right')
    end
    colormap(cmap_q)
    freezeColors()
end

%unfreezeColors()

set(splot(1),'pos',[.05 .58 .2 .4])
set(splot(2),'pos',[.265 .58 .2 .4])
set(splot(3),'pos',[.05 .15 .2 .4])
set(splot(4),'pos',[.265 .15 .2 .4])

set(splot(5),'pos',[.55 .58 .2 .4])
set(splot(6),'pos',[.765 .58 .2 .4])
set(splot(7),'pos',[.55 .15 .2 .4])
set(splot(8),'pos',[.765 .15 .2 .4])

set(cbar_b,'pos',[0.475,0.58,.01,.35])
set(cbar_q,'pos',[0.475,0.15,.01,.35])

colormap(cmap_q)
