%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sum-product algorithm initialization  
%
% Assumption: size(node(i).CPT)=1 for all leaf nodes  
%
% Input arguments:
% - N: number of nodes in the factor graph
% - node[n]: node structure
% - Adj[nxn]: adjacency matrix of the factor graph
% - Message[nxn]: array of empty matrices representing the messages
%
% Output arguments:
% - Message[nxn]: array of matrices representing the messages after
%                   initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Message=InitMessages(N,node,Adj,Message)

% scan all nodes in the factor graph
Message=cell(N,N);
for X=1:N
    % find all neighboring nodes of X
    neighbors=find(Adj(X,:));
    % init only leaf node messages
    if length(neighbors)==1
        % init leaf variable nodes messages
        if node(X).type=='variable'
            % send  all one message to neighboring function nodes
            for n=1:length(neighbors)
                Message{X,neighbors(n)}=ones(size(node(X).values));
            end
        % init leaf function nodes messages   
        else    
            % send function node to neighboring variable nodes
            for n=1:length(neighbors)
                Message{X,neighbors(n)}=node(X).CPT;
            end
        end
    end
end