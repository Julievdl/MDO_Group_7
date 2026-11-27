function [] = EMWET_write_init(design_vars, derived_vars, fixed)

    % DEFINITIONS FOR EMWET INIT FILE
    % writes .init file with the name EMWET_iteration

    Airfoil     = 'e553';
    section_num = 3;     % root, kink, tip
    airfoil_num = 3;     % define airfoil at root, kink, tip just in case?

    fid = fopen( 'EMWET_iteration.init','wt');
    fprintf(fid, '%g %g \n', derived_vars.m_TO, derived_vars.m_ZF);
    fprintf(fid, '%g \n', fixed.n_max);
    fprintf(fid, '%g %g %g %g \n', derived_vars.S, derived_vars.b, section_num, airfoil_num);
    
    %Airfoil Section definition - one per location (root, kink, tip)
    fprintf(fid, '0 %s \n', Airfoil);
    fprintf(fid, '%g %s \n', fixed.b_k/design_vars.s, Airfoil);
    fprintf(fid, '1 %s \n', Airfoil);
    
    %Wing section definition (root, kink, spar)
    fprintf(fid, '%g %g %g %g %g %g \n', design_vars.c_r, 0, 0, 0, fixed.xi_f, fixed.xi_a);
    fprintf(fid, '%g %g %g %g %g %g \n', derived_vars.c_k, fixed.b_k*tand(design_vars.Lambda_LE), fixed.b_k, fixed.b_k * tand(fixed.Gamma_dihedral), fixed.xi_f, fixed.xi_a);
    fprintf(fid, '%g %g %g %g %g %g \n', derived_vars.c_t, design_vars.s*tand(design_vars.Lambda_LE), design_vars.s, design_vars.s * tand(fixed.Gamma_dihedral), fixed.xi_f, fixed.xi_a);
    
    fprintf(fid, '%g %g \n', fixed.b_tank_start/design_vars.s, fixed.eta_tank_end);
    fprintf(fid, '%g \n', fixed.N_eng/2);

    %Engine position and weight definition
    fprintf(fid, '%g  %g \n', fixed.b_eng/design_vars.s, fixed.m_eng);
    
    %Material property definition
    for k = 1:4   % upper panel, lower pnel, front spar, rear spar
        fprintf(fid, '%g %g %g %g \n', fixed.E_al*1e6, fixed.rho_al, fixed.sigma_yield_tens*1e6, fixed.sigma_yield_comp*1e6);
    end

    fprintf(fid, '%g %g \n', fixed.eta_stringer, fixed.p_rib);
    fprintf(fid, '1 \n');
    fclose(fid);
end
