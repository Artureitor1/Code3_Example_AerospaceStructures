function F_el = computeF(n_el,n_ne,n_i,Qy_vect,x_an,Tmat,mat,x,Tn,Me,g,L1)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_i         Number of DOFs per node
%                  n_dof       Total number of DOFs
%   - Fdata  External nodal forces [Nforces x 3]
%            Fdata(k,1) - Node at which the force is applied
%            Fdata(k,2) - DOF (direction) at which the force acts
%            Fdata(k,3) - Force magnitude in the corresponding DOF
%--------------------------------------------------------------------------
% It must provide as output:
%   - Fext  Global force vector [n_dof x 1]
%            Fext(I) - Total external force acting on DOF I
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering to
% determine at which DOF in the global system each force is applied.

F_el=zeros(n_ne*n_i,n_el); 

for e=1:n_el %internal forces thermal expansion
    x1e=x(Tn(e,1),1);
    x2e=x(Tn(e,2),1);
    le=sqrt((x2e-x1e)^2);
    
    cond1=0;
    cond2=0;

    for i=1:length(x_an)
        if (x_an(i)>=x1e) && (cond1==0)
            lim_ienf=i;
            cond1=1; 
        end
        if (x_an(i)>x2e) && (cond2==0)
            lim_sup=i-1;
            cond2=1;
        end
    end

    qe=mean(Qy_vect(lim_inf:lim_sup,1));
    if x1e==L1 % contribution of half Me in node
        F_e=(qe*le/2)*[ 
            1
            le/6
            1
            -le/6];
        F_e=F_e+[-g*Me/2;0;0;0];
    end
    if x2e==L1 % contribution of half Me in node
        F_e=(qe*le/2)*[ 
            1
            le/6
            1
            -le/6];
        F_e=F_e+[0;0;-g*Me/2;0];
    end
    if (x1e~=L1) && (x2e~=L1)
        F_e=(qe*le/2)*[ 
            1
            le/6
            1
            -le/6];
    end

    F_el(:,e)=F_e;

end

end