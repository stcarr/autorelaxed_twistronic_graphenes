function [dos, E_list] = interp_kp_dos_gaussian(theta, sweep_vals, sweep_kpts, weights, ...
                                                b_size, E_list, sig)

    fprintf("Starting DOS calculation (via Gaussian method) \n");
    tic

    tri_idx = 1;

    all_kpts1 = sweep_kpts;
    allbands1 = sweep_vals;

    kpoints = all_kpts1(:,[2 1]);
    eig_vals = allbands1;

    nk = sqrt(size(eig_vals,1));
    nb = size(eig_vals,2);

    norbs = size(weights,1);
    dos = zeros(norbs,length(E_list));

    %sig = .005;

    for tar_b = (nb/2)-b_size:(nb/2+1)+b_size

        for i = 1:nk
            for j = 1:nk
                idx_h = (i-1)*nk + j;
                w_h = weights(:,tar_b,idx_h);
                Eh = eig_vals(idx_h,tar_b);
                dos = dos + w_h*(exp(-(E_list-Eh).^2/(2*sig^2))/(sqrt(2*pi*sig^2))/nk^2);
            end
        end

    end


    %idos = zeros(size(dos));
    %for x = 1:length(E_list)
    %    if (x > 1)
    %        idos(:,x) = trapz(E_list(1:x),dos(1:x));
    %    end
    %end

    alpha = 2.47;
    sc_alpha = alpha/(2*sind(theta/2));
    sc_area = sc_alpha^2*sind(60)*1e-2; %area in nm^2
    
    dos_rescale = 4/sc_area;
    %idos_rescale = 4;

    %idos(:) = idos_rescale*(idos(:) - 0.5*idos(end));

    dos = dos_rescale*dos'; % swap dimensions to match that returned by the ldos function

    toc

    
end