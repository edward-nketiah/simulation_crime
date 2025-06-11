function [space_time_data,N1,N_b] = fgenerate_data(S,T)
%------------------------------------------------------------------------------%
% Modified Fast Algorithm: Based on zhuang, ogata and vere-jones, 2004.        %
% Generate univariate multi-dimensional space-time self-exciting point process.%
%------------------------------------------------------------------------------%
% Generate G0 (G(1), background) based on the thinning methods
rng(055);
[lamstar,~] = lambda([0,0,0]);
N = poissrnd(lamstar*T*S(1)*S(2));
x = [rand(N,1)*S(1)-10,rand(N,1)*S(2)-10,rand(N,1)*T];
% x = sortrows(x,3);%There is no need to rearrange here and it will not affect the results: rearrange together at the end
[fun1,~] = lambda(x);
G = x(rand(N,1)<fun1/lamstar,:);
N_b = size(G,1);
G(:,4) = 1;
space_time_data = G;
% [~,fun2x,fun2y,fun2t] = lambda(G);
n = N_b;
% Generate offspring
theta = 0.6;
omega = 0.1;
sigmax = 0.01;
sigmay = 0.01;

while n>0
%---------------------------------------------------------------------------------------------------------%
%Expand the boundary to (-inf,+inf)*(-inf,+inf)*(0,+inf), so that the integration step in the loop        % 
%(corresponding to all percent signs in this program) can be omitted                                      %
%In the end, only points beyond the boundary need to be eliminated: because their cumulative intensity is % 
%limited, although some more points (outside the boundary) will be calculated, the integration part takes % 
%more time.                                                                                               %
%This part requires proof.                                                                                %                                                                                        %
lamoff = theta*ones(n,1);
%---------------------------------------------------------------------------------------------------------%
    n = poissrnd(lamoff);
    n1 = sum(n);
    %Position from the center point: In fact, fun2x, fun2y and fun2t should be used, here is a simplified calculation
    G2 = [randn(n1,2).*[sigmax,sigmay],exprnd(1/omega,n1,1),zeros(n1,1)];
    r = 1;
    G_off = zeros(n1,4);
        while n1>0
            n1 = sum(n>0);
            G_off(r:r+n1-1,1:3) = G2(r:r+n1-1,1:3)+G(n>0,1:3);
            
            n = n-1;               
            r = n1+r;        
        end
%collection point
    space_time_data = [space_time_data;G_off];
    G = G_off;
    n = size(G,1);
end
%Remove points outside the boundary
space_time_data(space_time_data(:,1)<-10|space_time_data(:,1)>10,:) = [];   
space_time_data(space_time_data(:,2)<-10|space_time_data(:,2)>10,:) = []; 
space_time_data(space_time_data(:,3)>T,:) = []; 
[space_time_data, N1] = sortrows(space_time_data,3);
end

function [fun1,fun2x,fun2y,fun2t] = lambda(x)
%This is a generalized parametrization intensity function 
mu = 5.71;
theta = 0.6;
omega = 0.1;
sigmax = 0.01;
sigmay = 0.01;
% fun1 = mu*(cos(x(:,3)/120)+2)./(2*pi*4.5^2).*exp(-x(:,1).^2/(2*4.5^2)-x(:,2).^2/(2*4.5^2));
fun1 = mu/(2*pi*4.5^2).*exp(-x(:,1).^2/(2*4.5^2)-x(:,2).^2/(2*4.5^2)); 
% fun1 = mu./(2*pi*4.5^2); 
fun2x = @(s)1/(sqrt(2*pi)*sigmax).*exp(-(s-x(:,1)).^2/(2*sigmax^2));
fun2y = @(s)1/(sqrt(2*pi)*sigmay).*exp(-(s-x(:,2)).^2/(2*sigmay^2));
fun2t = @(s)theta*omega*exp(-omega*(s-x(:,3)));
end