function rhs = rhsAdvect1dF1(u,v,h,nx,tol)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coded by Manuel Diaz, NTU, 2015.08.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Artebrant, Robert, and H. Joachim Schroll. "Limiter-free third order
%     logarithmic reconstruction." SIAM Journal on Scientific Computing
%     28.1 (2006): 359-381. 
% [2] Čada, Miroslav, and Manuel Torrilhon. "Compact third-order limiter
%     functions for finite volume methods." Journal of Computational Physics
%     228.11 (2009): 4118-4145.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Arrays
rhs=zeros(size(u)); 
uL=zeros(size(u));
uR=zeros(size(u));

if v>0
    for i = 2:nx-1
        % Slopes
        %dL=(u(i)-u(i-1))/h; dR=(u(i+1)-u(i))/h;
    
        % Coefficients
        %a=u(i)-(dR-dL)/24;
        %b=(dR+dL)/(2*h)-(3/4)*d*dx^2;
        %c=(dR-dL)/(2*h^2);
        
        % Parameter d
        d=0;
            
        % Left Approx. for u_{i+1/2} 
        % uL(x) = a+b*x+c*x^2+d*x^3; % where x= -dx/2
        uL(i)=(2*u(i-1)+5*u(i)-u(i+1))/6+d*h^3/4;
    
        % rhs
        rhs(i) = (uL(i)-uL(i-1))/h;
    end
else
    for i = 2:nx-2
        % Slopes
        %dL=(u(i+1)-u(i))/h; dR=(u(i+2)-u(i+1))/h;
        
        % Coefficients
        %a=u(i)-(dR-dL)/24;
        %b=(dR+dL)/(2*h)-(3/4)*d*dx^2;
        %c=(dR-dL)/(2*h^2);
        
        % Parameter d
        d=0;
        
        % Right Approx. for u_{i+1/2}
        % uR(x) = a+b*x+c*x^2+d*x^3; % where x= dx/2
        uR(i)=(-u(i)+5*u(i+1)+2*u(i+2))/6-d*h^3/4;
        
        % rhs
        rhs(i) = v*(uR(i)-uR(i-1))/h;
    end
end

end