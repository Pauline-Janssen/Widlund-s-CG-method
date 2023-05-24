function CGW_solution = CGW_rank_defi_final(L, x, alpha, n, limit, reorth, tol)

M = (L+L.')/2; % symm. part of L
S = (L.'-L)/2; % skew symm. part of L

b = L * x; % calculate b

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

arr = []; % array for reorthogonalisation of vectors 
h = v / sqrt(v'*M*v);     
arr = [arr h]; % add normed v_1

dimension_krylovspace = [];
S = svd(arr); % singular values of arr
max_sv = S(1);
rel_S = S/max_sv;
number_of_non_zeros = sum(rel_S > tol); 
rank_k = number_of_non_zeros;
dimension_krylovspace = [dimension_krylovspace rank_k];

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

    temp_rho = rho;     % store current rho as rho_(k-1)
    rho = v'*M*v;    % calculate rho_k

    w = (1+(rho/temp_rho)/w)^(-1);     

    u_neu = temp_u + w * (v + u - temp_u); % calculate new approx. solution
    temp_u = u; % store u_(k-1)
    u = u_neu;
    
    r = b - L*u; % calculate the residual
    residual = norm(r); % stopping criteria

    S = svd(arr); % singular values of arr
    max_sv = S(1);
     
    for i = 1:length(S)
         rel_ratio = S(i)/max_sv;
         if rel_ratio < tol
             S(i) = 0;
         end
     end
     
     number_of_non_zeros = sum(S > 0) ;
     rank_k = number_of_non_zeros;
     dimension_krylovspace = [dimension_krylovspace rank_k];
  
end 

CGW_solution = u;

% plotting of the rank deficiency with or without reorthogonalisation
k = length(dimension_krylovspace);
k_vector = 1:1:k;

figure(3);
txt = ['Krylovsubspace dimension, ' newline 'alpha = ' num2str(alpha) ', ', help_var];
plot(k_vector, dimension_krylovspace,'LineWidth',1, 'DisplayName',txt)
legend('location', 'southeast');
hold on

end