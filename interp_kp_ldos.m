function [dos, ldos, E_list] = interp_kp_ldos(theta, sweep_vals, sweep_weights, sweep_kpts, ...
                                                b_size, max_E, dE)

    % number of extra bands to include
    if ~exist('b_size','var')
        b_size = 10;
    end
    % size of dos window
    if ~exist('max_E','var')
        max_E = 0.3;
    end
    % spacing of dos window
    if ~exist('dE','var')
        dE = 1e-3;
    end

    E_list = [-max_E:dE:max_E];


    fprintf("Starting LDOS calculation (via interp. method) \n");
    tic

    tri_idx = 1;

    all_kpts1 = sweep_kpts;
    allbands1 = sweep_vals;
    allweights1 = sweep_weights; % norbs x nr x tot_dim x nk
    allweights1 = squeeze(sum(allweights1,1)); % sum over all layers and orbitals
    nr = size(allweights1,1);

    kpoints = all_kpts1(:,[1 2]);
    eig_vals = allbands1;

    nk = sqrt(size(eig_vals,1));
    nb = size(eig_vals,2);

    n_tri = nk*nk*(2*b_size + 2);

    triangles1 = zeros(n_tri,3);
    triangles2 = zeros(n_tri,3);
    weights1 = zeros(n_tri,nr);
    weights2 = zeros(n_tri,nr);
    k_tris1 = zeros(n_tri,3,2);
    k_tris2 = zeros(n_tri,3,2);



    for tar_b = (nb/2)-b_size:(nb/2+1)+b_size
    %for tar_b = 1;
    %fprintf("on band %d / %d \n",tar_b - ((nb/2)-b_size) + 1, 2*b_size + 2);

        for i = 1:nk
            for j = 1:nk
                tar_band(i,j) = eig_vals((i-1)*nk + j,tar_b);
                tar_weights(i,j,:) = allweights1(:,tar_b,(i-1)*nk + j);
                k_mesh(i,j,:) = kpoints((i-1)*nk + j,:);
                %k_mesh(i,j,1) = 0.1*(j-nk/2);
                %k_mesh(i,j,2) = 0.1*(i-nk/2);

                %tar_band(i,j) = norm( squeeze(k_mesh(i,j,:)) ).^2;
            end
        end
        %nk = sqrt(size(k_mesh,1)*size(k_mesh,2));

        % make triangles!
        for i = 1:nk
            for j = 1:nk

                ip = i+1;
                im = i-1;
                jp = j+1;
                jm = j-1;

                if (ip > nk)
                   ip = 1; 
                end
                if (im < 1)
                   im = nk; 
                end
                if (jp > nk)

                   jp = 1; 
                end
                if (jm < 1)
                   jm = nk; 
                end

                triangles1(tri_idx,1) = tar_band(i,j);
                triangles1(tri_idx,2) = tar_band(ip,j);
                triangles1(tri_idx,3) = tar_band(i,jp);
                triangles1(tri_idx,4) = tar_band(i,jp);
                weights1(tri_idx,:) = tar_weights(i,j,:);

                k_tris1(tri_idx,1,:) = k_mesh(1,1,:);
                k_tris1(tri_idx,2,:) = k_mesh(2,1,:);
                k_tris1(tri_idx,3,:) = k_mesh(1,2,:);
                k_0s1(tri_idx,:) = k_mesh(i,j,:);
                %k_tris1(tri_idx,1,:) = k_tris1(tri_idx,1,:) + k_mesh(i,j,:);
                %k_tris1(tri_idx,2,:) = k_tris1(tri_idx,2,:) + k_mesh(i,j,:);
                %k_tris1(tri_idx,3,:) = k_tris1(tri_idx,3,:) + k_mesh(i,j,:);

                triangles2(tri_idx,1) = tar_band(i,j);
                triangles2(tri_idx,2) = tar_band(im,j);
                triangles2(tri_idx,3) = tar_band(i,jm);
                weights2(tri_idx,:) = tar_weights(i,j,:);

                k_tris2(tri_idx,1,:) = k_mesh(end,end,:);
                k_tris2(tri_idx,2,:) = k_mesh(end-1,end,:);
                k_tris2(tri_idx,3,:) = k_mesh(end,end-1,:);
                %k_tris2(tri_idx,1,:) = k_tris2(tri_idx,1,:) + k_mesh(i,j,:);
                %k_tris2(tri_idx,2,:) = k_tris2(tri_idx,2,:) + k_mesh(i,j,:);
                %k_tris2(tri_idx,3,:) = k_tris2(tri_idx,3,:) + k_mesh(i,j,:);


                tri_idx = tri_idx+1;            

            end
        end
    end

    dos = zeros(length(E_list),1);
    ldos = zeros(length(E_list),nr);

    for E_idx = 1:length(E_list)
        E = E_list(E_idx);
        for t = 1:length(triangles1)
           if E > min(triangles1(t,:)) && E < max(triangles1(t,:))  
            %dos(E_idx) = dos(E_idx)+1;
            % replace this with proper gradient computation function!
            % need slope |b| and cross sectional area (length) f!
            % v is dk1
            v(1) = k_tris1(t,2,1) - k_tris1(t,1,1);
            v(2) = k_tris1(t,2,2) - k_tris1(t,1,2);
            % w is dk2
            w(1) = k_tris1(t,3,1) - k_tris1(t,1,1);
            w(2) = k_tris1(t,3,2) - k_tris1(t,1,2);
            % NOTE: here w(2) == 0!

            E0 = E - triangles1(t,1);
            Ev = triangles1(t,2) - triangles1(t,1);
            Ew = triangles1(t,3) - triangles1(t,1);

            if (abs(w(2)) > 1e-12)
               fprintf('WARNING: dk2 is not purely along x dir! \n');
               pause(10) 
            end

            % b is the vector such that E(k) = E0 + dot(k,b)
            b(1) = Ew/w(1);
            b(2) = (Ev - b(1)*v(1)) / v(2);

            % For w(2) not equal to zero:
            %{ 
            b(1) = Ew*v(1) - Ev*w(1);
            b(2) = Ev*w(2) - Ew*v(2);
            b = b/( w(2)*v(1) - w(1)*v(2));
            %}

            % ? dot(b,v) = Ev, dot(b,w) = Ew
            t1 = E0/dot(b,v); % E0/Ev
            t2 = E0/dot(b,w); % E0/Ew
            t3 = (E0 - dot(b,v)) / dot(b,w-v); % (E0-Ev)/(Ew-Ev)

            f = 0;

            if (t3 < 0 || t3 > 1)
                p1 = t1*v;
                p2 = t2*w;
                f = sqrt( sum((p1 - p2).^2) );
            end
            if (t1 < 0 || t1 > 1)
                p2 = t2*w;
                p3 = v + t3*(w - v);
                f = sqrt( sum((p2 - p3).^2) );
            end
            if (t2 < 0 || t2 > 1)
                p1 = t1*v;
                p3 = v + t3*(w - v);
                f = sqrt( sum((p1 - p3).^2) );
            end

            if(norm(b) == 0)
                fprintf('WARNING: b = 0! \n');
                pause(10) 
            end
            if (f/norm(b) > 10)
                fprintf('WARNING!! \n');
                pause(10)                
            end

            dos(E_idx) = dos(E_idx) + f/norm(b);
            ldos(E_idx,:) = ldos(E_idx,:) + weights1(t,:)*f/norm(b);

           end
           if E > min(triangles2(t,:)) && E < max(triangles2(t,:))  
            % replace this with proper gradient computation function!
            % need slope |b| and cross sectional area (length) 
            % v is dk1
            v(1) = k_tris2(t,2,1) - k_tris2(t,1,1);
            v(2) = k_tris2(t,2,2) - k_tris2(t,1,2);
            % w is dk2

            w(1) = k_tris2(t,3,1) - k_tris2(t,1,1);
            w(2) = k_tris2(t,3,2) - k_tris2(t,1,2);

            E0 = E - triangles2(t,1);
            Ev = triangles2(t,2) - triangles2(t,1);
            Ew = triangles2(t,3) - triangles2(t,1);

            if (abs(w(2)) > 1e-12)
               fprintf('WARNING: dk2 is not purely along x dir! \n');
               pause(10) 
            end

            b(1) = Ew/w(1);
            b(2) = (Ev - b(1)*v(1)) / v(2);

            % For w(2) not equal to zero:
            %{ 
            b(1) = Ew*v(1) - Ev*w(1);
            b(2) = Ev*w(2) - Ew*v(2);
            b = b/( w(2)*v(1) - w(1)*v(2));
            %}

            t1 = E0/dot(b,v);
            t2 = E0/dot(b,w);
            t3 = (E0 - dot(b,v)) / dot(b,w-v);

            f = 0;

            if (t3 < 0 || t3 > 1)
                p1 = t1*v;
                p2 = t2*w;
                f = sqrt( sum((p1 - p2).^2) );
            end
            if (t1 < 0 || t1 > 1)
                p2 = t2*w;
                p3 = v + t3*(w - v);
                f = sqrt( sum((p2 - p3).^2) );
            end
            if (t2 < 0 || t2 > 1)
                p1 = t1*v;
                p3 = v + t3*(w - v);
                f = sqrt( sum((p1 - p3).^2) );
            end

            if(norm(b) == 0)
                fprintf('WARNING: b = 0! \n');
                pause(10) 
            end
            if (f/norm(b) > 10)
                fprintf('WARNING!! \n');
                pause(10)                
            end

            dos(E_idx) = dos(E_idx) + f/norm(b);
            ldos(E_idx,:) = ldos(E_idx,:) + weights2(t,:)*f/norm(b);

           end
        end
        %E_idx/length(E_list)

    end

    %dos_sweep{t_idx} = dos;

    idos = zeros(size(dos));
    for x = 1:length(E_list)
        if (x > 1)
            idos(x) = trapz(E_list(1:x),dos(1:x));
        end
    end

    %tot_bands = 2*(b_size+1);
    %tot_bands = 4*tot_bands; % 2 for valley, 2 for spin
    %idos_rescale = tot_bands/idos(end);
    alpha = 2.47;
    sc_alpha = alpha/(2*sind(theta/2));
    sc_area = sc_alpha^2*sind(60)*1e-2; %area in nm^2
    %n0 = 1/sc_area;
    %dos_rescale = idos_rescale*n0;

    dos_rescale = 100*4/(2*pi)^2;
    idos_rescale = dos_rescale*sc_area;

    idos(:) = idos_rescale*(idos(:) - 0.5*idos(end));

    dos = dos_rescale*dos;
    ldos = dos_rescale*ldos; % ldos in same units as dos

    toc
end