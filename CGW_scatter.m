function CGW_solution = CGW_scatter(L,line,alpha,n)
    M = (L+L.')/2; % symm. part of L
    S = (L.'-L)/2; % skew symm. part of L
    K = M\S;
    
    colorMap = dictionary(0.01, "#008B00", 0.1, "#9C661F", 0.5, "#9400D3",1,"#9C661F",0.8, "#008B00", 0.9, "#EEB422", 3, "#00688B",50.0, "#9C661F", 100.0, "#008B00" );
    
    subplot(2,3,line)
    scatter(line*ones(n,1),imag(eig(K)),50,'.','Color',colorMap(alpha))
    hold on;