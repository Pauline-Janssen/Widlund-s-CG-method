function CGW_solution = CGW_alphas(L, x, alpha, n, limit)

M = (L+L.')/2; % symm. part of L
S = (L.'-L)/2; % skew symm. part of L
K = M\S;
lambda = abs(eigs(K,1)); 

b = L * x; 

u = zeros(n,1); % u_0
temp_u = u; % u_(-1)
e0 = sqrt((u-x)'*M*(u-x))/(sqrt(x'*M*x)); 
b0 = 2*sqrt(1+lambda^2)*(sqrt(1+lambda^(-2))+lambda^(-1))^(-0);

w = 1; % w_1
r = b; % residual r_0
v = linsolve(M,r); % v_1 

rho = v'*M*v;  % rho_1 

u_neu = temp_u + w * (v + u - temp_u); % u_1
temp_u = u; % store u_0 for next iteration
u = u_neu; % store u_1 for next iteration
r = b - L*u; % residual r_1 
e1 = sqrt((u-x)'*M*(u-x))/(sqrt(x'*M*x)); 
b1 = 2*sqrt(1+lambda^2)*(sqrt(1+lambda^(-2))+lambda^(-1))^(-1);

residual = norm(r);

err_vector = [e0 e1];
bound_vector = [b0 b1];

% for loop starts with k = 2 
for k = 2:limit
    if residual < 1e-10
        break;
    end
    
    v = linsolve(M,r);
     
    temp_rho = rho; % store current rho as rho_(k-1)
    rho = v'*M*v; % calculate rho_k

    w = (1+(rho/temp_rho)/w)^(-1);     

    u_neu = temp_u + w * (v + u - temp_u); % calculate new approx. solution
    temp_u = u; % store u_(k-1) 
    u = u_neu;
    
    r = b - L*u; % calculate the residual 
    residual = norm(r); % stopping criteria
    error =  sqrt((u-x)'*M*(u-x))/(sqrt(x'*M*x)); % relative error u_k-x in the M-norm
    err_vector(k+1) = error; % store the error for plotting
    bound_vector(k+1) = 2*sqrt(1+lambda^2)*(sqrt(1+lambda^(-2))+lambda^(-1))^(-k); % store the bound for plotting
    
end 

CGW_solution = u;

% plotting of the error and bound
k = length(err_vector); 
k_vector = 0:1:k-1;

figure(1);
txt = ['error for alpha = ', num2str(alpha)];
txt2 = ['error bound for alpha = ', num2str(alpha)];

colorMap = dictionary(0.01, "#008B00", 0.1, "#9C661F", 0.5, "#9400D3", 0.9, "#EEB422", 3, "#00688B",50.0, "#9C661F", 100.0, "#008B00");

semilogy(k_vector, err_vector,'Color', colorMap(alpha), 'LineWidth',1,'DisplayName',txt)
legend('location', 'east');
hold on

semilogy(k_vector, bound_vector,'Color', colorMap(alpha), 'LineWidth',1, 'LineStyle', '--','DisplayName',txt2)
legend('location', 'east');
hold on

end