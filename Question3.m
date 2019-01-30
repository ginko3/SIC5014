%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Model the Lung Cancer Diagnosis Problem, Ignoring Fatigue Feature
% - Compute Prior feature probabilities using the sum-product algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% Define Model1

% variables true and false
t=1; f=2;


% DEFINE STRUCTURE FOR THE NODES

% VARIABLE NODES
% DEFINE NODE1
node(1).name='S';                 % fill in name
node(1).type='variable';          % fill in type
node(1).values=[t f];             % assumed values
% DEFINE NODE2
node(2).name='B';                 % fill in name
node(2).type='variable';          % fill in type
node(2).values=[t f];             % assumed values
% DEFINE NODE3
node(3).name='L';                 % fill in name
node(3).type='variable';          % fill in type
node(3).values=[t f];             % assumed values
% DEFINE NODE4
node(4).name='X';                 % fill in name
node(4).type='variable';          % fill in type
node(4).values=[t f];             % assumed values

% FUNCTION NODES
% DEFINE NODE5
node(5).name='p(S)';              % fill in name
node(5).type='function';          % fill in type
node(5).CPT=[0.2 0.8];            % conditional probability table
% DEFINE NODE6
node(6).name='p(B|S)';            % fill in name
node(6).type='function';          % fill in type
node(6).CPT=[0.25 0.75;
             0.05 0.95];          % conditional probability table
% DEFINE NODE7
node(7).name='p(L|S)';            % fill in name
node(7).type='function';          % fill in type
node(7).CPT=[0.003   0.997;
             0.00005 0.99995];    % conditional probability table
% DEFINE NODE8
node(8).name='p(X|L)';            % fill in name
node(8).type='function';          % fill in type
node(8).CPT=[0.6   0.4;
             0.02 0.98];          % conditional probability table


% DEFINE ADJACENCY MATRIX OF THE FACTOR GRAPH
B=[1 1 1 0;0 1 0 0 ;0 0 1 1;0 0 0 1];
Adj=[zeros(4),B;B.',zeros(4)];

% draw factor graph using graphlayout package
figure(1);
A=sparse(Adj);addpath(genpath([pwd,'/graphlayout']));drawFG(A,node); title('Factor graph model without fatigue feature');
% NUMBER OF NODES IN THE FACTOR GRAPH
N=size(Adj,1);

% CREATE AN ARRAY [nxn] EMPTY MATRICES REPRESENTING THE MESSAGES
Message=cell(size(Adj));

% INIT MESSAGES
Message=InitMessages(N,node,Adj,Message);
%celldisp(Message)
% SUM-PRODUCT ALGORITHM ON A TREE IS EXACT
[Message,Marginal]=SumProductAlgorithm(N,node,Adj,Message);
%celldisp(Message)

% DISPLAY PRIOR NODE PROBABILITIES
disp(['Marginals obtained from sum-product algorithm']);
for X = 1:N
    if(node(X).type=='variable')
        disp(['P(' node(X).name, ') = ' num2str(Marginal{X})]);
    end
end

% PLOT PRIOR NODE PROBABILITIES
figure(2);
% Concatenate marginal probabilities for each variable node
x = cat(1,cell2mat(Marginal'));
% Collect variable node names
nodeNames=[];
for X=1:N
    if(node(X).type=='variable')
        nodeNames=[nodeNames;node(X).name];
    end
end

% Bar plot of the marginal probabilities for each variable node
bar(x, 'stacked'); set(gca, 'xticklabel', nodeNames);
ylabel('Probability');
title('Prior Marginal Probabilities');
legend('true', 'false');
hold on;

% loi de probabilité conjointe p(S,B,L,X)
% trouver la dimension de p
dim=[];
for X=1:N
    if(node(X).type=='variable')
        dim=[dim,length(node(X).values)];
    end
end

LoiConjointe=zeros(dim);
for valeurS=t:f
    for valeurB=t:f
        for valeurL=t:f
            for valeurX=t:f
                LoiConjointe(valeurS,valeurB,valeurL,valeurX)=node(5).CPT(valeurS)*node(6).CPT(valeurS,valeurB)*node(7).CPT(valeurS,valeurL)*node(8).CPT(valeurL,valeurX);
            end
        end
    end
end

% calcul des lois marginales par la méthode directe
for X=1:N
    if(node(X).type=='variable')
        % shifter la dimension sur laquelle on ne
        % marginalise pas en première position
        p=shiftdim(LoiConjointe,X-1);
        tmp=sum(sum(sum(p,2),3),4);
        Marginal{X}=reshape(tmp,1,length(node(X).values));
        Marginal{X}=Marginal{X}./sum(Marginal{X});
    end
end

% DISPLAY PRIOR NODE PROBABILITIES COMPUTED FROM DIRECT APPROACH
disp(['Marginals obtained from direct approach']);
for X = 1:N
    if(node(X).type=='variable')
        disp(['P(' node(X).name, ') = ' num2str(Marginal{X})]);
    end
end

