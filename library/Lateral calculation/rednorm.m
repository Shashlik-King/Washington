function rnorm = rednorm(f,Ndof,BC)
% 
% 
%--------------------------------------------------------------------------
% PURPOSE
% Compute the reduced Euclidian norm of the vector, f. Constrained terms in
% f are identified from the boundary conditions matrix, BC.
% 
% INPUT:  f       = Force vector, [ Ndof x 1 ],
%         Ndof    = Total number of degrees of freedom.
%         BC      = Global prescribed nodal degrees of freedom matrix,
%                   BC = [dof value]
%
% OUTPUT: rnorm   : Reduced Euclidian norm of f, scalar
%
% CODE            : AHA
% APPROVED        : LBI, ML, LA

% LAST MODIFIED   : AHA    13.08.2007   Programming
% 
%--------------------------------------------------------------------------


% -------- Set up array 'fixed' marking constrained terms in f ------------

fixed             = zeros(Ndof,1);

for m = 1:length(BC(:,1))

    dof_no        = BC(m,1);
    fixed(dof_no) = 1;
    
end


% -------- Evaluate reduced norm ------------------------------------------

rnorm             = 0;

for m = 1:length(f(:,1))

    if fixed(m) ~= 1

        rnorm     = rnorm + f(m)^2;

    end

end

rnorm             = sqrt(rnorm);
