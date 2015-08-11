function rhs = rhsAdvect1dF3(u,v,h,nx,tol)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite Difference LDLR Formulation of the RHS advective term
% coded by Manuel Diaz, NTU, 2015.08.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Artebrant, Robert, and H. Joachim Schroll. "Limiter-free third order
%     logarithmic reconstruction." SIAM Journal on Scientific Computing
%     28.1 (2006): 359-381. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants parameters
q=1.4; 

% Initialize Arrays
rhs=zeros(size(u)); 
uL=zeros(size(u));
uR=zeros(size(u));

if v>0
    for i = 2:nx-1
        % Slopes
        d1=(u(i)-u(i-1))/h; d2=(u(i+1)-u(i))/h;
    
        % Coefficients
        a=(1-tol)*(1+tol-( (2*abs(d1).^q.*abs(d2).^q)+tol)./...
            (abs(d1).^(2*q)+abs(d2).^(2*q)+tol) ); % a(d1,d2)
        b=a./(a-1); c=(a-1).*(d2.*(1-b)-d1)./(b-a); d=d1-c;
    
        % Left Approx. for u_{i+1/2}
        uL(i)=u(i)+(c*h).*pEta(a)+(d*h).*pEta(b);
    
        % rhs
        rhs(i) = (uL(i)-uL(i-1))/h;
    end
else
    for i = 2:nx-2
        % Slopes
        d1=(u(i+1)-u(i))/h; d2=(u(i+2)-u(i+1))/h;
        
        % Coefficients
        a=(1-tol)*(1+tol-( (2*abs(d1).^q.*abs(d2).^q)+tol)./...
            (abs(d1).^(2*q)+abs(d2).^(2*q)+tol) ); % a(d1,d2)
        b=a./(a-1); c=(a-1).*(d2.*(1-b)-d1)./(b-a); d=d1-c;
        
        % Right Approx. for u_{i+1/2}
        uR(i)=u(i+1)+(c*h).*nEta(a)+(d*h).*nEta(b);
        
        % rhs
        rhs(i) = v*(uR(i)-uR(i-1))/h;
    end
end

end

% $\eta(t)$ functions:
function eta=pEta(t); eta=-(log(1-t)+t)./t.^2; eta(t==0)=0.5; end
function eta=nEta(t); eta=((t-1).*log(1-t)-t)./t.^2; eta(t==0)=-0.5; end