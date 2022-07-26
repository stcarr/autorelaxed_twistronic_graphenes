function [dos, idos, E_list] = interp_kp_dos_gaussian(theta, sweep_vals, sweep_kpts, ...
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


    fprintf("Starting DOS calculation (via Gaussian method) \n");
    tic

    tri_idx = 1;

    all_kpts1 = sweep_kpts;
    allbands1 = sweep_vals;

    kpoints = all_kpts1(:,[2 1]);
    eig_vals = allbands1;

    nk = sqrt(size(eig_vals,1));
    nb = size(eig_vals,2);

    E_list = [-max_E:dE:max_E];
    dos = zeros(length(E_list),1);

    sig = .005;

    for tar_b = (nb/2)-b_size:(nb/2+1)+b_size

        for i = 1:nk
            for j = 1:nk
                Eh = eig_vals((i-1)*nk + j,tar_b);
                dos = dos + exp(-(E_list'-Eh).^2/(2*sig^2))/(sqrt(2*pi*sig^2))/nk^2;
            end
        end

    end

    %dos_sweep{t_idx} = dos;

    idos = zeros(size(dos));
    for x = 1:length(E_list)
        if (x > 1)
            idos(x) = trapz(E_list(1:x),dos(1:x));
        end
    end

    alpha = 2.47;
    sc_alpha = alpha/(2*sind(theta/2));
    sc_area = sc_alpha^2*sind(60)*1e-2; %area in nm^2
    
    dos_rescale = 4/sc_area;
    idos_rescale = 4;

    idos(:) = idos_rescale*(idos(:) - 0.5*idos(end));

    dos = dos_rescale*dos;

    toc

    
end