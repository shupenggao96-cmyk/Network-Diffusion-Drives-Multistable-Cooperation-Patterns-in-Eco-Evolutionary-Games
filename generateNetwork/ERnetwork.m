function [A,k]=ERnetwork(N,p)
% This function is used to generate random network
% input: N - the (final) number of the nodes in the network
%        p - the probability of rewiring edges
% N*p is approximately the mean average degree of the network

A=zeros(N,N);
for m=1:N
    for n=m+1:N
        if(rand(1,1)<p)
            A(m,n)=1;
            A(n,m)=1;
        end
    end
end
D=sum(A,2);
k=mean(D);

end