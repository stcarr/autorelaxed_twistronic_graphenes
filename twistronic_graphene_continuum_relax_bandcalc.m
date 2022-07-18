function [bands, qarr] = twistronic_graphene_continuum_relax_bandcalc(params)

    % authors: ziyan (zoe) zhu, stephen carr 
    % email: zzhu1@g.harvard.edu, stephen_carr1@brown.edu
    % calculation of an alternating-twist trilayer graphene band structure

    % !! Settings

    structure_type = params.structure_type; % 0: TBG, 1: Mono-on-bilayer, 2: alternating-twist trilayer
    relax_str = params.relax_strength; % scaling of atomic relaxations (0.0 turns off relaxations)
    E_field = params.E_field;   % vertical displacement field in eV (total potential energy across the three layers)
    AA_str = params.AA_tunnel_strength; % controls A-to-A interlayer tunneling strength (e.g. scales w0)
    AB_str = params.AB_tunnel_strength; % rescales A-to-B interlayer tunneling strength (e.g. scales w1)


    theta = params.theta;

    % truncation of momentum basis
    k_cutoff = params.k_cutoff;           % momentum cutoff, in units of norm(G1)
    grid_search = params.grid_search;       % [G1,G2] are in [-grid_search,grid_search]^2 

    % interlayer coupling sampling mesh
    r_max = params.r_max; % maximum radius, in Angstroms
    dr = params.dr; % spacing of r mesh, in Angstroms

    % interlayer fourier transform control
    inter_q_cut_scale = params.inter_q_cut_scale; % maximum scattered momentum, in units of K0
    inner_k_rad_scale = params.inner_k_rad_scale; % radius of each "island", in units of b12
    dk_scale = params.dk_scale; %  spacing of k mesh, in units of b12

    % band structure line cut settings
    q_cut_type = params.q_cut_type;         % what kind of line-cut we do in momentum space                        % 1: high symmetry line in the L12 bilayer moire Brillouin zone (single valley)
    nq = params.nq;                % number of k points to sample on each high symmetry line segment 

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
            qcut_name = 'L12 supercell';
            xt_labels = {'$K_{12}$', '$\Gamma_{12}$', '$M_{12}$', '$K_{12}$'};

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

        case 2 % supplied by params object
            q_list_x = params.qx_list;
            q_list_y = params.qy_list;
            qarr = zeros(size(q_list_x));

    end

    q_list = [q_list_x, q_list_y]; % the list of center sites

    if (q_cut_type == 1)

        q_list(end+1, :) = [pt(max_seg+1,1), pt(max_seg+1,2)];

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

    % apply relaxations
    if (relax_str ~= 0)
        tic
        fprintf("Applying atomic relaxations... \n")
        intercoupling.relax_configs(relax_str);
        fprintf("Relaxations done: ")
        toc
    else
        fprintf("No relaxation strength, skipping relax step... \n")        
    end

    % compute realspace couplings
    intercoupling.gen_tR();

    % apply fourier transform
    tic
    fprintf("FT of interlayer coupling starting... \n")
    intercoupling.gen_tK(AA_str,AB_str);
    fprintf("FT done: ")
    toc

    % check average AA vs AB coupling on first shell
    %H_inter_K0 = gen_interlayer_terms_dft(k_list(tar_dofs,:),layers,intercoupling,K0,inter_q_cut);
    %w0 = abs(H_inter_K0(3,1));
    %w1 = abs(H_inter_K0(4,1));
    %fprintf("w0 (inter_AA) = %s meV \n",num2str(1000*w0,'%.0f'));
    %fprintf("w1 (inter_AB) = %s meV \n",num2str(1000*w1,'%.0f'));


    % get H for each k point, and compute the band structure
    for q_idx = 1:size(q_list,1)

        tar_q = q_list(q_idx,:)+K0';
        H_inter = gen_interlayer_terms_dft(k_list,layers,intercoupling,tar_q,inter_q_cut);
        H_intra = gen_intralayer_terms_dft(k_list,layers,tar_q,E_field);

        H_blg = H_intra+H_inter;
        
        norb_l1 = 2*size(k1,1);
        norb_l2 = 2*size(k2,1);
        
        if (structure_type == 0)
            H = H_blg;
        elseif (structure_type == 1)
            % get size of mono-on-bilayer Hamiltonian
            size_fullH = 2*norb_l1 + norb_l2;
            l1_idxs = 1:norb_l1;
            l2_idxs = norb_l1+[1:norb_l1];
            l3_idxs = 2*norb_l1+[1:norb_l2];
            
            H = zeros(size_fullH, size_fullH);
            % first make the last 2/3'rds block identical to BLG Ham.
            H(norb_l1+1:end,norb_l1+1:end) = H_blg;
            % now assign diagonal of first 1/3 to that of the 2nd 2/3
            % (identical Dirac Hamiltonians)
            H(l1_idxs,l1_idxs) = H(l2_idxs,l2_idxs);
            % now add AB type coupling on identical k points
            t_AB = 0.31;
            for idx=1:2:norb_l1
                l1_A_orb = idx;
                l2_B_orb = norb_l1 + (idx+1); 
                H(l1_A_orb, l2_B_orb) = t_AB;
                H(l2_B_orb, l1_A_orb) = t_AB;
            end

        elseif (structure_type == 2)
            size_fullH = 2*norb_l1 + norb_l2;
            l1_idxs = 1:norb_l1;
            l2_idxs = norb_l1+[1:norb_l2];
            l3_idxs = norb_l1+norb_l2+[1:norb_l1];

            H = zeros(size_fullH, size_fullH);
            % first make the first 2/3'rds block identical to BLG Ham.
            H(1:norb_l1+norb_l2,1:norb_l1+norb_l2) = H_blg;
            % now assign diagonal of last 1/3 to that of the first 2/3
            % (identical Dirac Hamiltonians)
            H(l3_idxs,l3_idxs) = H(l1_idxs,l1_idxs);
            % now add twisted coupling between 2nd and 3rd layer
            H(l2_idxs,l3_idxs) = H(l2_idxs,l1_idxs);
            H(l3_idxs,l2_idxs) = H(l1_idxs,l2_idxs);
        end
        
        num_eigs = size(H,1);
        [~, raw_vals] = eigs(H,num_eigs); 
        [vals(q_idx,:), ~] = sort(real(diag(raw_vals))); % sort eigenvalues from smallest to largest 

        if (mod(q_idx,10) == 0 || q_idx == 1)
            fprintf("H Diag: %d / %d \n",q_idx,size(q_list,1));
        end
    end

    bands = vals;

end
