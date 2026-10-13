
function A = BAnetwork(N,m,m0)
% This function is used to generate scale free network
% input: N - the (final) number of the nodes in the network
%        m - the number of edges added during a node being added
%        m0 - the initial number of nodes  

% Construct the initial adjacent matrix of the network 
A = sparse(N, N);
A(1:m0,1:m0) = round(rand(m0));
A = tril(A);
A = A+A';
A = A - diag(diag(A)); 
for i = m0+1 : N
    Degree = sum(A(1:i-1,1:i-1));

    
    % Use the roulette bet to calculate the probability of connecting each 
    % of the existing nodes
    for j=2:i-1
        Degree(j) = Degree(j) + Degree(j-1); % accumulate degree
    end
    LinksNum = 0;
    
    % The process of adding edges
    while LinksNum< m
        dice = rand * Degree(i-1); 
        dice = ceil(dice);   
        for j=1:i-2
            if dice<=Degree(1) && A(i,1)==0
                A(i,1) = 1;
                A(1,i) = 1;
                LinksNum = LinksNum+1;  
            elseif dice>Degree(j) && dice<=Degree(j+1) && A(i,j+1)==0
                A(i,j+1) = 1;
                A(j+1,i) = 1;
                LinksNum = LinksNum+1;
            end
        end
    end
    
    %         for j = 1 : i-1
    %             p = Degree(j)/Sum_Degree;
    %             if rand<=p && sum(A(i,:))< m
    %                 A(i,j)=1;
    %                 A(j,i)=1;
    %             end
    %         end
    
end

% Calculate the number of nodes for each degree
Degree = sum(A);
list = unique(Degree);
num = zeros(1,length(list));
for i=1:length(list)
    num(i) = length(find(list(i)==Degree));
end

% plot double-log figures of degree distribution
loglog(list,num ./ sum(num),'.','markersize',20)  %listÎªºá×ø±ê
xlabel('k'),ylabel('P(k)')