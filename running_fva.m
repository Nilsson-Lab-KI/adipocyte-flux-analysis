clear all
changeCobraSolver('gurobi');

%%
%Create model again
A = readtable("reaction-list.csv", 'Format', '%s%s%s'); %Make sure in the same directory
reactionNames = table2array(A(:,1));
reactionFormulae = table2array(A(:,2));

model = createModel(reactionNames, reactionNames, reactionFormulae);

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
%Or load model from model.mat (make sure it in the directory)
model = load('model.mat');
model = model.model;
%%
%FVA with 95% Confidence Interval Original Control Experimental Data 

%Changing Measured Uptake/Release Fluxes
model = changeRxnBounds(model, 'EX_glc(e)', -319.565, 'l');
model = changeRxnBounds(model, 'EX_glc(e)', -214.967, 'u');

model = changeRxnBounds(model, 'EX_lac_L(e)', 33.31415, 'l');
model = changeRxnBounds(model, 'EX_lac_L(e)', 50.43576, 'u');

model = changeRxnBounds(model, 'EX_pyr(e)', -32.3, 'l');
model = changeRxnBounds(model, 'EX_pyr(e)', -30.8,'u');

model = changeRxnBounds(model, 'EX_glyc(e)', 1.853402, 'l');
model = changeRxnBounds(model, 'EX_glyc(e)', 3.33263, 'u');

%Changing TAG Prediction
model = changeRxnBounds(model, 'EX_tag_hs(e)', 1.63, 'l');
model = changeRxnBounds(model, 'EX_tag_hs(e)', 4.89, 'u');

%Adjusting Boundaries based on media composition
model = changeRxnBounds(model, 'EX_hdca(e)', 0, 'l');%Only FA release possible
model = changeRxnBounds(model, 'EX_o2(e)', 0, 'u');%Only uptake possible
model = changeRxnBounds(model, 'EX_co2(e)', 0, 'l');%Only release possible

[minFlux, maxFlux] = fluxVariability(model, 100,'printLevel', 2, 'allowLoops', 'original', 'threads', 1);

%%
%FVA with 95% Confidence Interval Original IL Experimental Data
%Create or load model^^

%Changing Measured Uptake/Release Fluxes
model = changeRxnBounds(model, 'EX_glc(e)', -255.02, 'l');
model = changeRxnBounds(model, 'EX_glc(e)', -216.63, 'u');

model = changeRxnBounds(model, 'EX_lac_L(e)', 72.69, 'l');
model = changeRxnBounds(model, 'EX_lac_L(e)', 88.30, 'u');

model = changeRxnBounds(model, 'EX_pyr(e)', -33.50, 'l');
model = changeRxnBounds(model, 'EX_pyr(e)', -29.64,'u');

model = changeRxnBounds(model, 'EX_glyc(e)', 0.60, 'l');
model = changeRxnBounds(model, 'EX_glyc(e)', 1.43, 'u');

%Changing TAG Prediction
model = changeRxnBounds(model, 'EX_tag_hs(e)', 2.69, 'l');
model = changeRxnBounds(model, 'EX_tag_hs(e)', 8.08, 'u');

%Adjusting Boundaries based on media composition
model = changeRxnBounds(model, 'EX_o2(e)', 0, 'u');
model = changeRxnBounds(model, 'EX_co2(e)', 0, 'l');

[minFlux, maxFlux] = fluxVariability(model, 100,'printLevel', 2, 'allowLoops', 'original', 'threads', 1);

%%
%FVA with 95% Confidence Interval Control Experimental Data Corrected for
%Matrigel uptake
%Create or load model^^

%Changing Measured Uptake/Release Fluxes
model = changeRxnBounds(model, 'EX_glc(e)', -234.05, 'l');
model = changeRxnBounds(model, 'EX_glc(e)', -129.45, 'u');

model = changeRxnBounds(model, 'EX_lac_L(e)', 75.44, 'l');
model = changeRxnBounds(model, 'EX_lac_L(e)', 92.56, 'u');

model = changeRxnBounds(model, 'EX_pyr(e)', -20.07, 'l');
model = changeRxnBounds(model, 'EX_pyr(e)', -18.61,'u');

model = changeRxnBounds(model, 'EX_glyc(e)', 1.85, 'l');
model = changeRxnBounds(model, 'EX_glyc(e)', 3.33, 'u');

%Changing TAG Prediction
model = changeRxnBounds(model, 'EX_tag_hs(e)', 1.63, 'l');
model = changeRxnBounds(model, 'EX_tag_hs(e)', 4.89, 'u');

%Adjusting Boundaries based on media composition
model = changeRxnBounds(model, 'EX_hdca(e)', 0, 'l');%Only FA release possible
model = changeRxnBounds(model, 'EX_o2(e)', 0, 'u');%Only uptake possible
model = changeRxnBounds(model, 'EX_co2(e)', 0, 'l');%Only release possible

[minFlux, maxFlux] = fluxVariability(model, 100,'printLevel', 2, 'allowLoops', 'original', 'threads', 1);

%%
%FVA with 95% Confidence Interval IL Experimental Data Corrected for
%Matrigel uptake
%Create or load model^^

%Changing Measured Uptake/Release Fluxes
model = changeRxnBounds(model, 'EX_glc(e)', -169.50, 'l');
model = changeRxnBounds(model, 'EX_glc(e)', -131.11, 'u');

model = changeRxnBounds(model, 'EX_lac_L(e)', 114.82, 'l');
model = changeRxnBounds(model, 'EX_lac_L(e)', 130.43, 'u');

model = changeRxnBounds(model, 'EX_pyr(e)', -20.88,'l');
model = changeRxnBounds(model, 'EX_pyr(e)', -17.36,'u');

model = changeRxnBounds(model, 'EX_glyc(e)', 0.60, 'l');
model = changeRxnBounds(model, 'EX_glyc(e)', 1.43, 'u');

%Changing TAG Prediction
model = changeRxnBounds(model, 'EX_tag_hs(e)', 2.69, 'l');
model = changeRxnBounds(model, 'EX_tag_hs(e)', 8.08, 'u');

%Adjusting Boundaries based on media composition
model = changeRxnBounds(model, 'EX_o2(e)', 0, 'u');
model = changeRxnBounds(model, 'EX_co2(e)', 0, 'l');

[minFlux, maxFlux] = fluxVariability(model, 100,'printLevel', 2, 'allowLoops', 'original', 'threads', 1);