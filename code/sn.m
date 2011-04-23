function [phi,psi] = sn(EDGE,NFM,SigT,SigS,RegMat,Source)
% function [phi,psi] = sn(EDGE,NFM,Sigma,RegMat,Source)
%   This function solves the one-speed transport equation in slab geometry
%   using the discrete ordinates method.  The S4 and diamond difference 
%   approximations are implemented, and only vacuum boundaries are treated.
%  
%   Inputs: 
%       EDGE    -- region edges 
%       NFM     -- number of fine meshes per region
%       SigT    -- total cross-section for each material
%       SigS    -- scattering cross-section for each material
%       RegMat  -- which material goes in each region 
%       Source  -- region uniform isotropic source [1/cm^3-s]
%   Output
%       phi     -- fine mesh-centered scalar flux
%       psi     -- fine mesh-edge angular flux

% assign angles and weights
mu = [ -0.8611363115 -0.3399810435  0.3399810435  0.8611363115 ];
wt = [  0.3478548451  0.6521451549  0.6521451549  0.3478548451 ];
% allocate
totNFM = sum(NFM);            % total number of fine meshes
psi    = zeros(totNFM+1,4);   % current angular flux 
phi    = zeros(totNFM,1);     % current scalar flux
S      = zeros(totNFM,4);     % fine mesh source
Q      = S;                   % emission density
fmmid  = zeros(totNFM,1);     % fine mesh material id
% compute discretization
j = 0;
for i = 1:length(NFM)
    Delta( (j+1):(j+NFM(i)) )    = ( EDGE(i+1) - EDGE(i) ) / NFM(i);
    S( (j+1):(j+NFM(i)), : )     = Source(i)/2; 
    fmmid( (j+1):(j+NFM(i))   )  = RegMat(i);  
    j = sum(NFM(1:i));
end
% precompute coefficients, following Eqs. 20.16-20.20.
alpha = zeros(4,1); 
A     = zeros(totNFM,4);
B     = A;
for i = 1:totNFM
    m = fmmid(i); 
    for n = 1:4  % use the sign of mu to choose appropriate form
        smu    = sign(mu(n)); 
        denom  = 2*mu(n)+smu*(1+smu*alpha(n))*SigT(m)*Delta(i);
        A(i,n) = (2*mu(n)-smu*(1-smu*alpha(n))*SigT(m)*Delta(i)) / denom;
        B(i,n) = smu * 2 * Delta(i) / denom;
    end
end
% convergence parameters
eps_phi = 1e-5; max_it  = 200;
err_phi = 1;    it = 0;
% Begin source iterations
while (err_phi > eps_phi && it <= max_it )
    % Save old scalar flux
    phi0 = phi; 
    % Update source
    for i = 1:totNFM
        Q(i,:) = S(i,:) + 0.5*SigS(fmmid(i))*phi(i);
    end
    % Perform sweeps
    for i = totNFM:-1:1  % right-to-left
        psi(i,1:2) = A(i,1:2).*psi(i+1,1:2) + B(i,1:2).*Q(i,1:2);
    end
    for i = 1:totNFM     % left-to-right
        psi(i+1,3:4) = A(i,3:4).*psi(i,3:4) + B(i,3:4).*Q(i,3:4);
    end    
    % --- Insert Acceleration Here ---
    % Update phi (need cell-centered psi, so use Eq. 20.15)
    for i = 1:totNFM
         phi(i) = sum( wt(:) .* (0.5*( (1+alpha(:)).*psi(i+1,:)' + ...
                                       (1-alpha(:)).*psi(i,:)'    ) ));
    end
    % Update error and iteration counter
    err_phi =  max(  abs(phi-phi0)./phi0 );
    it = it + 1;
end
if (it <= max_it)
    disp(['Converged in ',num2str(it), ' iterations.'])  
else
    disp('Failed to converge.')
end

end 
