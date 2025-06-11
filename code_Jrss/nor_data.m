function [result, V] = nor_data(space_time_data,S,T,k1,k2,k0,k3,r1,r2,r3)
%---------------------------------------------------------------------------%
%              ~~~ Split space S, for all time (ignoring order) data ~~~                %
%Calculates the variance of variance normalized data. %
% Calculate the distance between the points in each sub-interval block. The distance between the r3th point (1 dimension) and the r1 and 2nd points (higher than 1 dimension) is Di%
%Reference: G. O. MOHLER, M. B. SHORT, P. J. BRANTINGHAM, F. P. SCHOENBERG, %
%and G. E. TITA, 2011, Self-Exciting Point Process Modeling of Crime %
%k1 and k2 correspond to the number of divisions on the x and y axes respectively; k3 and k4 respectively correspond to the blocks divided by integers on the x and y axes, %
%k0 split time                                                                 %
%Consider the end of the block each time to determine the attributes of the entire segment to improve operating efficiency; %
%Special note: This algorithm has strict restrictions on space constraints, but is extremely loose on time constraints.                 %
%This algorithm calculates spatial distance, space-time distance, and time distance at the same time, but the space-time distance is not strict %
%Since the calculation and summation of the subsequent kernel density estimation takes too long, a boundary value is assigned to the points in each block, and its effect is %
%When calculating the kernel function, ignore values ??outside the boundary, and the boundary constraint is value<e-15. This approximation is only available if the data size is %
%When % is large (in the order of millions), the effect is good. When the amount of data is small (within 10,000), the exact solution can be obtained directly. However, when %
%The amount of data exceeds 10 million, and the notebook cannot handle it.                                     %
%---------------------------------------------------------------------------%
% Data variance normalization
V = std(space_time_data);
nor_data = space_time_data./V;
%Divide the area S into k1*k2 rectangles and take their center
if isempty(S)
    E = [min(nor_data)-10^-10;max(nor_data)+10^-10]';%Small program bug, make sure there are no points on the border
else
    E = [S(1) S(2); S(3) S(4);T(1)-10^-10 T(2)+10^-3]./V';%Boundaries of S, considering all boundaries
end
S_new(:,1) = E(1,1)+(E(1,2)-E(1,1))/k1*(0:k1);%The dividing point corresponding to the first spatial coordinate
S_new(:,2) = E(2,1)+(E(2,2)-E(2,1))/k2*(0:k2);%The dividing point corresponding to the second spatial coordinate
S_new1 = E(3,1)+(E(3,2)-E(3,1))/k0*(0:k0);%The dividing point corresponding to the time coordinate
% Find the radius corresponding to each sub-area
% Find the space-time radius and space radius of the 15th and 100th closest points respectively
n = length(nor_data);
nor_data1 = sortrows(nor_data,1);
if n>k3
    a = k3;
else
    k3 = n;
    a = k3;
end
x_begin = 1;
y_begin = 1;
i = 1;
j = 1;
clear nor_data
while i<=k1%The block where the x-axis of the judgment point is located
    p = S_new(i+1,1);
    while nor_data1(a,1)<p&&a<n
        a = min(a+k3,n);
    end
    nor_data2 = nor_data1(a-k3+1:a,1);
    x_end = a-k3+sum(nor_data2<p);
    nor_data1(x_begin:x_end,4) = i;
    x_begin = x_end+1;
    i = i+1;
end
nor_data_2 = sortrows(nor_data1,2);
if n>k3
    a = k3;
else
    k3 = n;
    a = k3;
end
C = cell(k1,k2);
D = zeros(k1,k2);
Num = zeros(k1,k2);
clear nor_data1 nor_data2
while j<=k2%The block where the y-axis of the judgment point is located
    p = S_new(j+1,2);
    while nor_data_2(a,2)<p&&a<n
        a = min(a+k3,n);
    end
    nor_data2 = nor_data_2(a-k3+1:a,2);
    y_end = a-k3+sum(nor_data2<p);
    nor_data_2(y_begin:y_end,5) = j;
    for i = 1:k1%Determine the number of points in each block
        temp = nor_data_2(y_begin -1 + find(nor_data_2(y_begin:y_end,4)==i),:);
        Num(i,j) = size(temp,1);
        if Num(i,j)>99%Determine whether the number of points in the block is greater than 99
            if ~isempty(r1)
                temp1 = sort(pdist2(temp,temp),2)';
                temp(:,6) = temp1(r1,:);%Calculate the space-time distance of the r1th point
            else
                temp_1 = sort(pdist2(temp(:,1:2),temp(:,1:2)),2)';
                temp(:,6) = temp_1(r2,:);%Calculate the spatial distance of the r2th point
            end
            D(i,j) = 1;%If the number of points in the block is not less than 100, assign a value of 1
        else
            D(i,j) = 0;
        end
        C{i,j} = temp;%The points in each block are placed into the corresponding cells            
    end
    y_begin = y_end+1;
    j = j+1;
end
clear temp temp1 temp_1 nor_data_2
nor_data = zeros(n,6); 
a1 = 1;
for i=1:k1
  for j=1:k2    
    k = 1;
    while D(i,j)==0&&Num(i,j)~=0%Don't count blocks without points to save time
        if sum(sum(Num))<100
            disp('Initial background number is too small.');
            break;
        end
        Temp1 = [];
        Temp_1 = i+(-1*k:k);
        Temp_2 = j+(-1*k:k);
        Temp_1(Temp_1<1|Temp_1>k1) = [];
        Temp_2(Temp_2<1|Temp_2>k2) = [];
        if sum(sum(Num(Temp_1,Temp_2)))>99
            for k_1 = 1:length(Temp_1)
                for k_2 = 1:length(Temp_2)
                    Temp1 = [Temp1; C{Temp_1(k_1),Temp_2(k_2)}(:,1:3)];
                end
            end
            D(i,j) = 1;
            Temps = C{i,j}(:,1:3);
            if ~isempty(r1)                          
                Temp2 = sort(pdist2(Temps,Temp1),2)';
                C{i,j}(:,6) = Temp2(r1,:);
            else
                Temp_2 = sort(pdist2(Temps(:,1:2),Temp1(:,1:2)),2)'; 
                C{i,j}(:,6) = Temp_2(r2,:);
            end
        else
            k = k+1;
        end
    end
    if ~isempty(C{i,j})%Eliminate these points to save time
        nor_data(a1:a1+size(C{i,j},1)-1,1:6) = C{i,j};
        if length(space_time_data)>10000%Spatial boundary constraints, please note that the boundary is restored to its original scale (the same below)
            nor_data(a1:a1+size(C{i,j},1)-1,7) = max(V(1:2))*C{i,j}(:,6)*sqrt(12*2);%The boundary constraint is value<e-10
        else
            nor_data(a1:a1+size(C{i,j},1)-1,7) = max(V(1:2))*C{i,j}(:,6)*sqrt(15*2);%The boundary constraint is value<e-15
        end
        a1 = a1+size(C{i,j},1);
    end
  end
end
%Find the time radius of the 100th closest point
clear C
nor_data = sortrows(nor_data,3);
if n>k3
    a = k3;
else
    k3 = n;
    a = k3;
end
t_begin = 1;
j = 1;
D = zeros(k0,1);
C1 = cell(k0,1);
Num = zeros(k0,1);
while j<=k0%The block where the judgment point t-axis is located
    p = S_new1(j+1);
    while nor_data(a,3)<p&&a<n
        a = min(a+k3,n);
    end
    nor_data3 = nor_data(a-k3+1:a,3);
    t_end = a-k3+sum(nor_data3<p);
    nor_data(t_begin:t_end,8) = j;
    temp = nor_data(t_begin:t_end,:);
    Num(j) = size(temp,1);
        if Num(j) > 199 %Determine whether the number of points in the block is greater than 99
            temp1 = sort(pdist2(temp(:,3),temp(:,3)),2)';
            temp(:,9) = temp1(r3,:);%Calculate the time distance of the r3th point
            D(j) = 1;%If the number of points in the block is not less than 100, assign a value of 1
        else
            D(j) = 0;
        end
    C1{j} = temp;%The points in each block are placed into the corresponding cells            
    t_begin = t_end+1;
    j = j+1;
end
nor_data = zeros(n,10); 
a1 = 1;
for i = 1:k0 
    k = 1;
    while D(i)==0&&Num(i)~=0
        Temp1 = [];
        Temp = i + (-1*k:k);
        Temp(Temp<1|Temp>k0) = [];    
        Temp(Num(Temp)==0) = [];
        if sum(Num(Temp))>299
            for k_1 = 1:length(Temp)
                Temp1 = [Temp1; C1{Temp(k_1)}(:,3)];
            end
            D(i) = 1;
            Tempt = C1{i}(:,3);
            Temp2 = sort(pdist2(Tempt,Temp1),2)';
            C1{i}(:,9) = Temp2(r3,:);
        else
            k = k+1;
        end
    end
    if ~isempty(C1{i})     
        nor_data(a1:a1+size(C1{i},1)-1,1:9) = C1{i};
        if length(space_time_data)>10000
            nor_data(a1:a1+size(C1{i},1)-1,10) = V(3)*C1{i}(:,9)*sqrt(4*2);% time boundary constraint is value<e-10
        else
            nor_data(a1:a1+size(C1{i},1)-1,10) = V(3)*C1{i}(:,9)*sqrt(5*2);%The boundary constraint is value<e-15
        end
        a1 = a1+size(C1{i},1);
    end       
end
clear C1
result = nor_data;
result = [space_time_data, result(:,4:end)];
clear nor_data