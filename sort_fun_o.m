function [temp_bs,temp_bt,temp_o] = sort_fun_o(t1,space_time_data,bs_sort,bt_sort,o_sort,n1,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matching descendant strength KDE: Because KDE uses truncation, the truncation 
%used for each estimation point is different (the length is the same) %
%Predetermine the KDE (approximately: select the nearest point) cutoff range 
%corresponding to each estimation point.               %
%Because the Monte Carlo method is used, for n1*n sample points, there are only
%points near the n sample points that were drawn, and their %
%KDE is not 0. The so-called nearby refers to its boundary bound_o                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,l_sort] = sort(t1);
% Select estimated points where KDE is not 0
% Background--Time
temp1_bt = zeros(n,1);
temp1_bt(bt_sort) = 1:length(bt_sort);
temp3_bt = find(temp1_bt~=0);
temp2_bt = single(temp1_bt(temp3_bt));
temp4_bt = [temp3_bt(1);temp3_bt(2:end-1)-temp3_bt(1:end-2);n-temp3_bt(end-1)];
temp_bt = repelem(temp2_bt,temp4_bt);
%Background--space
[~,y_sort] = sort(space_time_data(:,2));
temp1_bs = zeros(n,1);
temp1_bs(bs_sort) = 1:length(bs_sort);
temp2_bs = temp1_bs(y_sort);
temp3_bs = find(temp2_bs~=0);
temp_bs = 1:length(bs_sort);
temp4_bs = [temp3_bs(1);temp3_bs(2:end-1)-temp3_bs(1:end-2);n-temp3_bs(end-1)];
temp_bs = [repelem(single(temp_bs)',temp4_bs),y_sort];%Memory consumption here
temp_bs = sortrows(temp_bs,2);
temp_bs(:,2) = [];
%Descendants
temp1_o = zeros((n1-1)*(n-1),1);
temp1_o(o_sort) = 1:length(o_sort);
temp2_o = temp1_o(l_sort);
temp3_o = find(temp2_o~=0);
temp_o = 1:length(o_sort);
temp4_o = [temp3_o(1);temp3_o(2:end-1)-temp3_o(1:end-2);(n1-1)*(n-1)-temp3_o(end-1)];
temp_o = [repelem(single(temp_o)',temp4_o),l_sort];% Memory consumption here
temp_o = sortrows(temp_o,2);
temp_o(:,2) = [];
