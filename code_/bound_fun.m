function [X1,Y1,T1,D_1,Xs,Ys,D_s,Tt,D_t] = bound_fun(result,r)
n1 = size(result,1);
x = result(:,1);
y = result(:,2);
t = result(:,3);
if r == 1
    D1 = max(result(:,6),0.02);% Nearest neighbor space-time distance (root)
    boundt = result(:,10);%Time boundary, used to approximate the kernel function
    [X1,Y1,T1,D_1] = Bound(x,y,t,n1,D1,boundt);
    Xs = 0;
    Ys = 0;
    D_s = 0;
    Tt = 0;
    D_t = 0;
else
    D1 =  max(result(:,6),0.02);%Nearest neighbor space distance (root)
    D2 = max(result(:,9),0.05);%nearest neighbor time distance
    bounds = result(:,7);%Space boundary, used to approximate the kernel function
    boundt = result(:,10);%Time boundary, used to approximate the kernel function
    [Xs,Ys,~,D_s] = Bound(x,y,[],n1,D1,bounds);
    [~,~,Tt,D_t] = Bound([],[],t,n1,D2,boundt);
    X1 = 0;
    Y1 = 0;
    T1 = 0;
    D_1 = 0;
end
end



function [X1,Y1,T1,D_1] = Bound(x,y,t,n1,D1,bound)
if isempty(t)&&~isempty(y)%Spatial boundary: Delimited by y coordinate,
    %and the boundary value is determined by the overall spatial distance.
    i = 100;
    y_sort = sort(y);
    if quantile(y_sort(i:end)-y_sort(1:end-i+1),0.05)>quantile(bound,0.9)
        i = 4;
        k = 2;
    else
        i = 150;
        k = 50;
    end
    while quantile(y_sort(i:end)-y_sort(1:end-i+1),0.05)<quantile(bound,0.9)
        i = i+k;
        if i>1000&&i<n1
            disp('Note: huge computation! Please change B (in fwkde.m) smaller.');
            break
        elseif i>n1
            i = n1;
            break
        end
    end
    if n1>2*i
        X = zeros(n1-2*i,2*i,'single');
        Y = zeros(n1-2*i,2*i,'single');
        D1_new = zeros(n1-2*i,2*i,'single');
        for j = 1:2*i+1
            X(:,j) = x(j:end-2*i+j-1);
            Y(:,j) = y(j:end-2*i+j-1);
            D1_new(:,j) = D1(j:end-2*i+j-1);
        end
        X1 = [repmat(X(1,:),i,1);X(1:end,:);repmat(X(end,:),i,1)];
        clear X
        Y1 = [repmat(Y(1,:),i,1);Y(1:end,:);repmat(Y(end,:),i,1)];
        clear Y
        D_1 = [repmat(D1_new(1,:),i,1);D1_new(1:end,:);repmat(D1_new(end,:),i,1)];
        clear D1
    else
        X1 = repmat(x',n1,1);
        Y1 = repmat(y',n1,1);
        D_1 = repmat(D1',n1,1);
    end
    T1 = 0;
elseif ~isempty(t)&&isempty(y)%time boundary
    i = 100;
    if quantile(t(i:end)-t(1:end-i+1),0.05)>quantile(bound,0.9)
        i = 4;
        k = 2;
    else
        i = 150;
        k = 50;
    end
    while quantile(t(i:end)-t(1:end-i+1),0.05)<quantile(bound,0.9)
        i = i+k;
        if i>1000&&i<n1
            disp('Note: huge computation! Please change B (in fwkde.m) smaller.');
            break
        elseif i>n1
            i = n1;
            break
        end
    end
    if n1>2*i
        T = zeros(n1-2*i,2*i,'single');
        D1_new = zeros(n1-2*i,2*i,'single');
        for j = 1:2*i+1
            T(:,j) = t(j:end-2*i+j-1);
            D1_new(:,j) = D1(j:end-2*i+j-1);
        end
        T1 = [repmat(T(1,:),i,1);T(1:end,:);repmat(T(end,:),i,1)];
        clear T
        D_1 = [repmat(D1_new(1,:),i,1);D1_new(1:end,:);repmat(D1_new(end,:),i,1)];
        clear D1
    else
        T1 = repmat(t',n1,1);
        D_1 = repmat(D1',n1,1);
    end
    X1 = 0;
    Y1 = 0;
elseif ~isempty(t)&&~isempty(y)%space-time boundary
    i = 100;
    if quantile(t(i:end)-t(1:end-i+1),0.05)>quantile(bound,0.9)
        i = 4;
        k = 2;
    else
        i = 150;
        k = 50;
    end
    while quantile(t(i:end)-t(1:end-i+1),0.05)<quantile(bound,0.9)
        i = i+k;
        if i>1000&&i<n1
            disp('Note: huge computation! Please change B (in fwkde.m) smaller.');
            break
        elseif i>n1
            i = n1;
            break
        end
    end
    if n1>2*i
        X = zeros(n1-2*i,2*i,'single');
        Y = zeros(n1-2*i,2*i,'single');
        T = zeros(n1-2*i,2*i,'single');
        D1_new = zeros(n1-2*i,2*i,'single');
        for j = 1:2*i+1
            X(:,j) = x(j:end-2*i+j-1);
            Y(:,j) = y(j:end-2*i+j-1);
            T(:,j) = t(j:end-2*i+j-1);
            D1_new(:,j) = D1(j:end-2*i+j-1);
        end
        X1 = [repmat(X(1,:),i,1);X(1:end,:);repmat(X(end,:),i,1)];
        clear X
        Y1 = [repmat(Y(1,:),i,1);Y(1:end,:);repmat(Y(end,:),i,1)];
        clear Y
        T1 = [repmat(T(1,:),i,1);T(1:end,:);repmat(T(end,:),i,1)];
        clear T
        D_1 = [repmat(D1_new(1,:),i,1);D1_new(1:end,:);repmat(D1_new(end,:),i,1)];
        clear D1
    else
        X1 = repmat(x',n1,1);
        Y1 = repmat(y',n1,1);
        T1 = repmat(t',n1,1);
        D_1 = repmat(D1',n1,1);
    end
end
end