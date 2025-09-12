function [b_data,bs_sort,bt_sort,o_data,o_sort,x1,y1,t1,n] = MC_data(P,space_time_data,n1)
%Monte Carlo method--separate data into background events and trigger events
%The storage format of %P is as follows: the first column is its diagonal element, followed 
%by the probability of triggering adjacent events, and the first few lines are padded with zeros.
%rng(999);
R = rand(size(P,1),1);
x = R < P(:,1);
b_data = space_time_data(x,:);
bt_sort = find(x==1);
[b_data,temps_sort]= sortrows(b_data,2);%Sort spatial distances according to y coordinate
b_data(:,3) = space_time_data(x,3);%Sorting the time distance requires that the time order remains unchanged, because the latter is calculated separately in time and space.
bs_sort=bt_sort(temps_sort);

n = size(space_time_data,1);
x = space_time_data(:,1);
y = space_time_data(:,2);
t = space_time_data(:,3);
x1 = zeros(n1-1,n-1);%This value means that when calculating g in lambda, the most recent n1 data at the time point are used for calculation.
y1 = x1;
t1 = x1;
for j = 1:n1-1
    x1(j,:) = x(2:end)-[Inf*ones(j-1,1);x(1:end-j)];
    y1(j,:) = y(2:end)-[Inf*ones(j-1,1);y(1:end-j)];
    t1(j,:) = t(2:end)-[Inf*ones(j-1,1);t(1:end-j)];
end
x1 = reshape(x1',[],1);
y1 = reshape(y1',[],1);
t1 = reshape(t1',[],1);
o_sort = 1:length(t1);
x = R<cumsum(P,2);
x = x(2:end,2:end)-x(2:end,1:end-1);
x = (x==1);
o_data = [x1(x),y1(x),t1(x)];
o_sort = o_sort(x);
[o_data,temp_sort]= sortrows(o_data,3);
o_sort=o_sort(temp_sort);