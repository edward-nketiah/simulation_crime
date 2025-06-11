function P1 = Initial_P(space_time_data,n1,p)
%The storage format of %P is as follows: the first column is its diagonal element, 
%followed by the probability of triggering adjacent events, and the first few lines are padded with zeros.
%In a given time range, the closer the point, the greater the probability
n = size(space_time_data,1);
x = space_time_data(:,1);
y = space_time_data(:,2);
x1 = zeros(n1-1,n-1);%This value means that when calculating g in lambda, the most recent n1 data at the time point are used for calculation.
y1 = x1;
for j = 1:n1-1
    x1(j,:) = x(2:end)-[Inf*ones(j-1,1);x(1:end-j)];
    y1(j,:) = y(2:end)-[Inf*ones(j-1,1);y(1:end-j)];
end
x1 = x1';
y1 = y1';
temp0 = ((x1+y1).^2);
temp = ((x1+y1).^2);
temp(temp0==0) = eps;
temp = 1./temp;
P1 = temp./sum(temp,2);
Temp = [zeros(1,n1-1);1./temp0];
P1 = [p.*ones(n,1),[zeros(1,n1-1);(1-p).*P1]];
P1(1,1) = 1;
P1(isnan(sum(Temp,2)-Inf),1) = 0;
P1(isnan(sum(Temp,2)-Inf),2:end) = P1(isnan(sum(Temp,2)-Inf),2:end)/(1-p);