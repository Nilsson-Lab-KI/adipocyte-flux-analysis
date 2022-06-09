clear all
changeCobraSolver('gurobi');

%%
%Model Generation%
A = readtable("reaction-list.csv", 'Format', '%s%s%s'); %Make sure in the same directory
reactionNames = table2array(A(:,1));
reactionFormulae = table2array(A(:,2));

model = createModel(reactionNames, reactionNames, reactionFormulae);

B = readtable("met-formulae.xlsx"); %Make sure in the same directory
model.metFormulas = table2array(B(:,2));

%%
%Verification%

%Check if all reactions are mass balanced. 
%Note: exchange reaction will not be balanced
results = verifyModel(model, 'massBalance', true);

%Check for flux consistency
cnt=1;
tol = 1e-6;

TableChecks{cnt,1} = 'Check for flux consistency';
param.epsilon=1e-4;
param.modeFlag=0;
%param.method='null_fastcc';
param.method='fastcc';
printLevel = 1;
[fluxConsistentMetBool,fluxConsistentRxnBool,fluxInConsistentMetBool,fluxInConsistentRxnBool,model] = findFluxConsistentSubset(model,param,printLevel);
if isempty(find(fluxInConsistentRxnBool))
    TableChecks{cnt,2} = 'Model is flux consistent.';
else
    TableChecks{cnt,2} = 'Model is NOT flux consistent';
end
cnt = cnt + 1;

inconsistent_metabolites = [model.mets, num2cell(fluxInConsistentMetBool)];
inconsistent_reactions = [model.rxns, num2cell(fluxInConsistentRxnBool)];

%%
%Key metrics of model%
S = full(model.S); %Complete Stoichiometry Matrix

[rrfS, p] = rref(S);%Reduced row echelon form of S
rankS = rank(S);

%%
%Finding the Null Space of S
i = 1:size(reactionFormulae,1);
temp = ismember(i, p);
freerxns = [];
for i = 1:size(reactionFormulae,1)
    if temp(i) == 0
        freerxns = [freerxns i];
    end
end

temp = [];
for j = freerxns
    temp = [temp, rrfS(1:size(reactionFormulae,1),j)];
end

ident = eye(size(freerxns,2));

for i=1:size(freerxns,2)
    temp(freerxns(i), :) = ident(i,:);
end

nullS = temp;

%%
%Changing the model bounds to 10000
ub = 10000*ones(size(model.ub));
model.ub = ub;

lb = model.lb;
for i=1:size(model.lb)
    if model.lb(i) == -1000
        lb(i) = -10000;
    else
        lb(i) = lb(i);
    end
end
model.lb = lb;

%%
%Save Model
save('model.mat', "model");