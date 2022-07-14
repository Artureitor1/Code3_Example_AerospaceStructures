function Kel = computeKelBar(n_ne,n_i,n_el,x,Tnod,mat,Tmat)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
%                  n_el       Total number of elements
%   - x     Nodal coordinates matrix [n x n_d]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - % Material properties matrix
%             mat(m,1) = Young modulus of material m
%             mat(m,2) = Section area of material m
%             mat(m,3) = Section inertia of material m
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - Kel   Elemental stiffness matrices [n_el_dof x n_el_dof x n_el]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------

Kel=zeros(n_ne*n_i,n_ne*n_i,n_el);
    for e=1:n_el
        x1e=x(Tnod(e,1),1);
        x2e=x(Tnod(e,2),1);
        le=sqrt((x2e-x1e)^2);
        
     Ke=[ 12 6*le -12 6*le
         6*le 4*le^2 -6*le 2*le^2
         -12 -6*le 12 -6*le
         6*le 2*le^2 -6*le 4*le^2
        ];
     Ke=Ke*mat(Tmat(e),3)*mat(Tmat(e),1)/le^3;
        
     Kel(:,:,e) = Ke(:,:); 
    end
end