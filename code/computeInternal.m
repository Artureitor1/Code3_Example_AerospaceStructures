function [pu,pt,Fy,Mz] = computeInternal(n_el,n_ne,n_i,x,u,Kel,Td,Tnod,mat,Tmat)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
%                  n_el       Total number of elements
%   - u     Global displacement vector [n_dof x 1]
%            u(I) - Total displacement on global DOF I
%   - Td    DOFs connectivities table [n_el x n_el_dof]
%            Td(e,i) - DOF i associated to element e
%   - x     Nodal coordinates matrix [n x n_d]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - eps   Strain vector [n_el x 1]
%            eps(e) - Strain of bar e
%   - sig   Stress vector [n_el x 1]
%            sig(e) - Stress of bar e
%--------------------------------------------------------------------------

Fy=zeros(n_el,2);
Mz=zeros(n_el,2);
pu=zeros(n_el,4);
pt=zeros(n_el,3);

for e=1:n_el
    x1e=x(Tnod(e,1),1);
    x2e=x(Tnod(e,2),1);
    le=sqrt((x2e-x1e)^2);

    ue=zeros(n_ne*n_i,1);
    for i=1:n_ne*n_i
        I=Td(e,i);
        ue(i,1)=u(I);
    end

    Finte=Kel(:,:,e)*ue;
        Fy(e,1)=-Finte(1);
        Fy(e,2)=Finte(3);
        Mz(e,1)=-Finte(2);
        Mz(e,2)=Finte(4);

    coef = (1/le^3)*[2 le -2 le
                  -3*le -2*le^2 3*le -le^2
                  0 le^3 0 0
                  le^3 0 0 0]*ue;
    pu(e,:)=[coef(1),coef(2),coef(3),coef(4)];
    pt(e,:)=[3*coef(1),2*coef(2),coef(3)];
end


end