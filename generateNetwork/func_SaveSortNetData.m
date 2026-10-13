function func_SaveSortNetData(netType,A)
% input: adjacent matrix , network type
% output: sorted adjacent matrix and Laplacian matrix saved in .mat format
%         sorted Laplacian matrix  saved in .txt format

N = size(A,1);
D = sum(A,2);
avDegree = round(mean(D));
L = A - diag(D);

P = speye(N);
P = sortrows([diag(L),P]);
P = P(:,2:end);
L = P*L*P';
clear P;
A = L - diag(diag(L));

fileRoot = ['.\LaplacianMatrixSort_', num2str(N),'\'];
create_folder(fileRoot);

save([fileRoot,netType,num2str(avDegree),'.mat'],'A','L');

% take the up triangle of the adjacent matrix A (including the diagonal)
% L = triu(L);

% find the coordinate of nozero values
[x,y,~] = find(L~=0);

% Write the Laplacian matrix into a txt file, nodeIndex start at 0
% EXAMPLE
% nodIindex-1 nondeIndex-1 L(nondeIndex,nondeIndex)
fid = fopen([fileRoot,netType,num2str(avDegree),'.txt'],'w');
for j = 1:length(x)
    fprintf(fid,'%d %d %d\n',x(j)-1,y(j)-1, L(x(j),y(j)));
end
fclose(fid);

end

function create_folder(data_folder)
if ~isfolder(data_folder)
    mkdir(data_folder);
end
end