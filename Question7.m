%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Model the Lung Cancer Diagnosis Problem, Ignoring Fatigue Feature
% - Compute Prior feature probabilities using the sum-product algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% Define Model2

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
% DEFINE NODE5
node(5).name='F';                 % fill in name
node(5).type='variable';          % fill in type
node(5).values=[t f];             % assumed values

% FUNCTION NODES
% DEFINE NODE6
node(6).name='p(S)';              % fill in name
node(6).type='function';          % fill in type
node(6).CPT=[0.2 0.8];            % conditional probability table
% DEFINE NODE7
node(7).name='p(B|S)';            % fill in name
node(7).type='function';          % fill in type
node(7).CPT=[0.25 0.75;
             0.05 0.95];          % conditional probability table
% DEFINE NODE8
node(8).name='p(L|S)';            % fill in name
node(8).type='function';          % fill in type
node(8).CPT=[0.003   0.997;
             0.00005 0.99995];    % conditional probability table
% DEFINE NODE9
node(9).name='p(X|L)';            % fill in name
node(9).type='function';          % fill in type
node(9).CPT=[0.6   0.4;
             0.02 0.98];          % conditional probability table
% DEFINE NODE10
node(10).name='p(F|B,L)';            % fill in name
node(10).type='function';          % fill in type
node(10).CPT(:,:,t)=[0.75   0.1;
                    0.5 0.05];          % conditional probability table
node(10).CPT(:,:,f)=1-node(10).CPT(:,:,t);

% DEFINE ADJACENCY MATRIX OF THE FACTOR GRAPH
B=[1 1 1 0 0;0 1 0 0 1;0 0 1 1 1;0 0 0 1 0;0 0 0 0 1];
Adj=[zeros(5),B;B.',zeros(5)];

% draw factor graph using graphlayout package
figure(1);
A=sparse(Adj);addpath(genpath([pwd,'/graphlayout']));drawFG(A,node); title('Factor graph model without fatigue feature');
% NUMBER OF NODES IN THE FACTOR GRAPH
N=size(Adj,1);

% CREATE AN ARRAY [nxn] EMPTY MATRICES REPRESENTING THE MESSAGES
Message=cell(size(Adj));

% INIT MESSAGES
Message=InitMessagesLoopy(N,node,Adj,Message);
%celldisp(Message)
% SUM-PRODUCT ALGORITHM ON A TREE IS EXACT
MaxIterationNumber=100;
[Message,Marginal,Entropy]=SumProductAlgorithmLoopy(N,node,Adj,Message,MaxIterationNumber);
%celldisp(Message)

% DISPLAY PRIOR NODE PROBABILITIES
disp(['Marginals obtained from sum-product algorithm']);
for X = 1:N
    if(node(X).type=='variable')
        disp(['P(' node(X).name, ') = ' num2str(Marginal{X})]);
    end
end

% DISPLAY Entropy of each variable node
for X = 1:N
    if(node(X).type=='variable')
        subplot(1,5,X);
        plot((1:MaxIterationNumber),Entropy(X,:));
        xlabel('Iteration');
        ylabel(['Entropy of variable node' num2str(node(X).name) ' .']);
    end
end


