function phi = cpm(Delta,SigT,SigS,Source)
% function phi = cpm(Delta,SigT,SigS,Source)
%   This function solves the one-speed transport equation in slab geometry
%   for the scalar flux using a the collision probability method.
%  
%   Inputs: (all must be column vectors, i.e dimension [n,1]
%       Delta   -- region widths [cm]
%       SigT    -- region total cross-section [1/cm]
%       SigS    -- region scattering cross-section [1/cm]
%       Source  -- region uniform isotropic source [1/cm^3-s]
%   Output
%       phi     -- scalar flux in each region
tau = comptau(Delta,SigT);  % compute tau
P = ffcp(Delta,SigT,tau);   % compute collision probabilities
H = zeros(length(Delta));   % set up H and right hand side
s = zeros(length(Delta),1); % right hand side
for i = 1:length(Delta)
    for j = 1:length(Delta)
        H(i,j) = -SigS(j)/SigT(j)*P(i,j);        % add diagonal 1's below
        s(i) = s(i) + Source(j)*Delta(j)*P(i,j); % right hand side
    end
end
H = H + eye(length(Delta));     % add 1's to the diagonal
f = ( H \ s);                   % solve for the collision rates
phi = f ./ ( SigT .* Delta );   % solve for the flux
end 

function tau = comptau(Delta,SigT)
% function tau = comptau(Delta,SigT)
%   This function computes the optical pathlengths between regions.
%   
%   Inputs:
%       Delta   -- region widths [cm]
%       SigT    -- region total cross-section [1/cm]
%   Output:
%       tau     -- optical path length
tau = zeros(length(Delta)); 
DeltaSigT = Delta.*SigT;    % length of regions in mfp's 
for i = 1:length(Delta)
    % tau(i,i) remains 0
    for ip = i+1:length(Delta) % only ip > i, since tau(x,x')=tau(x',x)
        % adding optical distance between x(i+1/2) and x(ip-1/2)
        tau(i,ip) = sum( DeltaSigT(i+1:ip-1) ); 
        tau(ip,i) = tau(i,ip);
    end
end
end

function P = ffcp(Delta,SigT,tau)
% function P = ffcp(Delta,SigT,tau)
%   This function computes the first-flight collision probabilities.  Note,
%   it uses MATLAB's symbolic function for E3, which is REALLY SLOW.  The
%   student is strongly encouraged to look up fast and accuracte numerical
%   approximations to the En functions.  See Hebert's textbook.
%
%   Inputs:
%       Delta   -- region widths [cm]
%       SigT    -- region total cross-section [1/cm]
%       tau     -- optical path length
%   Outputs:
%       P       -- first-flight collision probabilities
P = zeros(length(Delta));
for i = 1:length(Delta)
    P(i,i) = 1 - 0.5/SigT(i)/Delta(i) * ...
            ( 1 - 2*mfun('Ei', 3, SigT(i)*Delta(i)) );
    for j = 1:i-1 
        P(i,j) = 0.5/SigT(j)/Delta(j) * ...
            ( mfun('Ei', 3, tau(i,j)) - ...
              mfun('Ei', 3, tau(i,j)+Delta(i)*SigT(i)) - ...
              mfun('Ei', 3, tau(i,j)+Delta(j)*SigT(j)) + ...
              mfun('Ei', 3, tau(i,j)+Delta(j)*SigT(j)+Delta(i)*SigT(i) ));
        P(j,i) = P(i,j); % reciprocity in action!
    end
end
end

 