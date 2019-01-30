%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sum-product algorithm with flooding schedule: 
% update until each edge in factor graph has been traversed in each direction  
%
% Input arguments:
% - N: number of nodes in the factor graph
% - node[N]: node structure
% - Adj[NxN]: adjacency matrix of the factor graph
% - Message[NxN]: array of empty matrices representing the messages
% - MaxIterationNumber: max. number of iterations
%
% Output arguments:
% - Message[NxN]: array of matrices representing the messages after
%                 the sum-product update
% - Marginal[N]: array of matrices representing the marginals
% - Entropy[N,MaxIterationNumber]: Entropy for each variable node and each iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Message,Marginal,Entropy]=SumProductAlgorithmLoopy(N,node,Adj,Message,MaxIterationNumber)

% scan all nodes while a message to a neighbor is empty 
% MODIF loop=1;
%while loop
for iteration=1:MaxIterationNumber
   for X=1:N
        % find all neighboring nodes of X
        neighbors=find(Adj(X,:));
        % compute a message towards all neighbors of X - if possible
        for n=1:length(neighbors)
            % collect the set of neighbors\neighbor(n)
            delta=setdiff(neighbors,neighbors(n));
            % test if the message towards neighbors(n) is empty and computable
            %MODIF if isempty(Message{X,neighbors(n)}) && sum(cellfun(@isempty,Message(delta,X)))==0  
            if sum(cellfun(@isempty,Message(delta,X)))==0
                % update variable to function node messages
                if node(X).type=='variable'
                    % compute the product of all incoming messages from delta
                    Message{X,neighbors(n)}=ones(size(node(X).values));
                    for m=1:length(delta)
                        Message{X,neighbors(n)}=Message{X,neighbors(n)}.*Message{delta(m),X};
                    end
                % update function to variable node messages 
                elseif node(X).type=='function'
                    % - compute function node x product of all incoming
                    % messages from delta conserving the order of the variables
                    w=1;
                    for m=1:length(delta)
                        w=kron(Message{delta(m),X},w);
                    end
                    % - and marginalize out all variables except neighbors(n)
                    % for all values assumed by variable node neighbors(n)
                    for v=1:length(node(neighbors(n)).values)
                        % shift the dimension corresponding to neighbors(n) towards first dimension
                        [B,I]=find(neighbors(n)==neighbors);
                        prob=shiftdim(node(X).CPT,I-1);
                        Message{X,neighbors(n)}(v)=sum(prob(v,:).*w);
                    end
                end
            end
        end
   end
    
   % MODIF stopping rule: check if at least one message in the array is missing
   %loop=sum(cellfun(@isempty,Message(find(Adj))));
    
    % Termination: as soon as each edge in factor graph has been traversed in each direction
    % Collect set of variable nodes
    variable_nodes=[];
    for n=1:N
        if(node(n).type=='variable')
            variable_nodes=[variable_nodes,n];
        end
    end

    % MODIF: Compute marginal for each variable node for each iteration
    for X=1:length(variable_nodes)
         % find all neighboring nodes of X
        neighbors=find(Adj(X,:));
        % compute marginal of X
        Marginal{X}=ones(size(node(X).values));
        for m=1:length(neighbors)
            Marginal{X}=Marginal{X}.*Message{neighbors(m),X};
        end
        % normalize marginal of X to probability mass function
        Marginal{X}=Marginal{X}./sum(Marginal{X});
        % MODIF: Compute Entropy of variable node X for each iteration
        Entropy(X,iteration)=sum(-log2(Marginal{X}).*Marginal{X});
    end
end