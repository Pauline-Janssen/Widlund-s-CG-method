function CGW_solution = CGW(L,b, limit)

M = (L+L.')/2; % symm. part of L
n = size(L,1);
u = zeros(n,1); % u_0
temp_u = u; % u_(-1)

w = 1; % w_1
r = b; % residual r_0
v = linsolve(M,r); % v_1 

rho = v'*M*v; % rho_1

u_neu = temp_u + w * (v + u - temp_u); % u_1
temp_u = u; % store u_0 for next iteration
u = u_neu; % store u_1 for next iteration
r = b - L*u; % residual r_1

residual = norm(r); % stopping criteria

% for loop starts with k = 2 
for k = 2:limit
    if residual < 1e-12
        break;
    end
    
    v = linsolve(M,r);
    
    temp_rho = rho;     % store current rho as rho_(k-1)
    rho = v'*M*v;    % calculate rho_k

    w = (1+(rho/temp_rho)/w)^(-1);     

    u_neu = temp_u + w * (v + u - temp_u); % approx. solution
    temp_u = u; % store u_(k-1)
    u = u_neu;
    
    r = b - L*u; % calculate the residual
    residual = norm(r); % stopping criteria
end 

CGW_solution = u;

end