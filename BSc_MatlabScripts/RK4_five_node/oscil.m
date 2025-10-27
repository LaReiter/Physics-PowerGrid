function B = oscil(A,F)
% returnerer delta_i,j*sin(theta_i*dtheta - theta_j*dtheta) matrix
% hvor delta_i,j = K_i,j*a_i,j  (dvs. koblingsstyrke gange "aktiv"-delta
% funktion)
B = (A.*sin(repmat(F,1,length(F))' - repmat(F,1,length(F))))*ones(length(F),1);
% B = sum(sin(meshgrid(F)-meshgrid(F)'),2);
end

