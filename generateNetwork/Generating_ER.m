clear;clc;close all;

N=400; % The number of nodes in the network
netType = 'ER';

for k = 80:1:80
    p = k/N;
    [A,~]=ERnetwork(N,p);
    A = full(A);
    D=sum(A,2);
    avDegree = round(mean(D));
    
    % Save the adjacent and Laplacian matrix in degree descending order
    func_SaveSortNetData(netType,A);
    fprintf('Degree = %d, Done£¡\n',avDegree);
    
end

