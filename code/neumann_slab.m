function [phi,phiL] = neumann_slab(L,SigT,SigS,Q,N,M)
% function phi = neumann_slab(L,SigT,SigS,Q,N,M)
%   This function solves the one-speed transport equation in slab geometry
%   for the scalar flux using a Neumann series solution.  The integrals
%   are computed using the trapezoid rule.
%  
%   The problem is a slab of width L, with uniform total cross-section
%   SigT, uniform (isotropic) scattering cross-section SigS, and
%   uniform isotropic source of volumetric strength Q.
%   
%   Inputs:
%       L       -- slab width [cm]
%       SigT    -- total cross-section [1/cm]
%       SigS    -- scattering cross-section [1/cm]
%       Q       -- uniform isotropic source [1/cm^3-s]
%       N       -- number of points for evaluating flux
%       M       -- number of Neumann terms to compute
%   Output
%       phi     -- scalar flux [N,1]
%       phiL    -- Neumann terms [N,M+1] (M=0 yields uncollided only)
phi  = zeros(N,1);      % scalar flux
phiL = zeros(N,M+1);    % the Lth collided fluxes
QL   = Q*ones(N,1);     % source vector
% uncollided flux
phiL(:,1) = ithcollided(L,N,QL,SigT);
% collided fluxes
for i = 1:M
   QL          = SigS*phiL(:,i);     	     % compute new source
   phiL(:,i+1) = ithcollided(L,N,QL,SigT);   % get ith collided flux term
end
% scalar flux is just the sum
for i = 1:M+1
    phi = phi + phiL(:,i);
end
end % function neumann_slab

function phi = ithcollided(L,N,Q,SigT)
% function phi = ithcollided(L,N,Q,SigT)
%   This function computes a Neumann term for the given source.  The
%   exponential integral is computed using Matlab's symbolic function Ei,
%   which is pretty slow!  The method of subtraction of singularities has
%   been used.  The flux is computed as follows:
%     phi(x)   = int( k(x,x')*Q(x')dx', x=0,L )
%     phi(x_i) = h * ( 0.5*k(x_i,x_0)*Q(x_0) + 1.0*k(x_i,x_1) + ... )
%   
%   Inputs:
%       L       -- slab width [cm]
%       N       -- number of points for evaluating flux
%       Q       -- source
%       SigT    -- total cross-section [1/cm]
%   Output
%       phi     -- the current Neumann term
phi = zeros(N,1);
h = L / (N-1);
% loop through all N
for i = 1:N
    xi = h*(i-1);
    phi(i) = 0;
    for j = 1:N % begin trapezoid integration
        if (i~=j) % if i=j, Q(i)=Q(j) and the term vanishes; avoids NaN!
            if (j==0 || j==N) % end points get coefficient of 0.5
                A = 0.5;
            else
                A = 1.0;
            end
            xj      = h*(j-1);
            phi(i)  = phi(i) + ...
                      A*h*mfun('Ei',1,(SigT*abs(xi-xj)))*(Q(j)-Q(i));
        end 
    end
    % extra part from subtraction of the singularity
    phi(i) = phi(i) + 0.5*Q(i)/SigT * ...
             (2-mfun('Ei',2,SigT*(L-xi))-mfun('Ei',2,SigT*(xi)));
end
end % function ithcollided


    
