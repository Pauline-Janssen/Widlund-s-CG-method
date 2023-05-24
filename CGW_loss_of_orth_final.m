function CGW_solution = CGW_loss_of_orth_final(L, x, n, limit, reorth)

M = (L+L.')/2; % symm. part of L
S = (L.'-L)/2; % skew symm. part of L

b = L * x; % calculate b

u = zeros(n,1); % u_0 
temp_u = u; % u_(-1)

w = 1; % w_1 
r = b; % residual r_0
v = linsolve(M,r); % v_1 

rho = v'*M*v;  % rho_1 

u_neu = temp_u + w * (v + u - temp_u); % u_1
temp_u = u; % store u_0 for next iteration
u = u_neu; % store u_1 for next iteration
r = b - L*u; % residual r_1

residual = norm(r); % stopping criteria

arr = []; % array for orthogonality or reorthogonalisation of vectors 
h = v / sqrt(v'*M*v);  
arr = [arr h]; % add normed v_1

orth_step = [];
orth_step_1 = norm(eye(1) - (arr' * M * arr), "fro");
orth_step = [orth_step orth_step_1];

help_var = 'without reorthogonalisation';

% for loop starts with k = 2
for k = 2:limit
    if residual < 1e-12
        break;
    end
    
    v = linsolve(M,r);
    
    % double reorthogonalisation with modfied Gram-Schmidt
    if reorth
        help_var = 'with reorthogonalisation';
        for  l = 1 :2
            for j = 1:k-1
                v = v - (arr(:,j)'*M*v)*arr(:,j); 
                % vectors in arr already normed, Gram-Schmidt simplifies
            end
        end
    end
    
    arr = [arr v/sqrt(v'*M*v)]; 

    temp_rho = rho; % store current rho as rho_(k-1)
    rho = v'*M*v; % calculate rho_k

    w = (1+(rho/temp_rho)/w)^(-1);    

    u_neu = temp_u + w * (v + u - temp_u); % calculate new approx. solution
    temp_u = u; % store u_(k-1)
    u = u_neu;
    
    r = b - L*u; % calculate the residual  
    residual = norm(r); % stopping criteria
    orth_step(k) = norm(eye(k) - (arr' * M * arr), "fro");
end 

CGW_solution = u;

% plotting of the orthogonality with or without reorthogonalisation
k = length(orth_step);
k_vector = 1:1:k;

figure(1);
txt = [num2str(help_var)];
legend('location', 'east');
semilogy(k_vector, orth_step,'-','LineWidth',1,'DisplayName',txt)

hold on

end