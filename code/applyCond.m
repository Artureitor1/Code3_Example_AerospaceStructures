function [vL,vR,uR] = applyCond(n_i,n_dof,fixNod)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_i      Number of DOFs per node
%                  n_dof    Total number of DOFs
%   - fixNod  Prescribed displacements data [Npresc x 3]
%              fixNod(k,1) - Node at which the some DOF is prescribed
%              fixNod(k,2) - DOF (direction) at which the prescription is applied
%              fixNod(k,3) - Prescribed displacement magnitude in the corresponding DOF
%--------------------------------------------------------------------------
% It must provide as output:
%   - vL      Free degree of freedom vector
%   - vR      Prescribed degree of freedom vector
%   - uR      Prescribed displacement vector
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering to
% determine at which DOF in the global system each displacement is prescribed.

    uR = zeros(size(fixNod,1),1);
    vR = zeros(size(fixNod,1),1);
    v = linspace(1,n_dof,n_dof);
        for i = 1:size(uR,1)
            if (fixNod(i,2))== 2
                vR(i,1) = 2*fixNod(i,1);
                uR(i,1) = fixNod(i,3);
            end
            if (fixNod(i,2))== 1
                vR(i,1) = 2*(fixNod(i,1))-1;
                uR(i,1) = fixNod(i,3);
            end
        end 
    vL = transpose (setdiff(v,vR));
end