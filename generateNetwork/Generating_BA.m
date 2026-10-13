clear;clc;close all;

N=400;
netType = 'BA';

% m: The number of added edges when one node is added
for m = 40:1:40
%     m=3;
    m0=m+1; %The initial number of nodes
    A=BAnetwork(N,m,m0);
    A = full(A);
    D = sum(A,2);
    avDegree = round(mean(D));
    % Save the adjacent and Laplacian matrix in degree descending order
    func_SaveSortNetData(netType,A);
    fprintf('Degree = %d, Done£¡\n',avDegree);
end
