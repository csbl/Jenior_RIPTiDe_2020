%% Comparing different integration algorithms using the same set of data
% 12/15/2019 BVD
clear all
close all
clc

% Load the original model and make the updates to the model
initCobraToolbox(false)
changeCobraSolver('gurobi6', 'all');

% Load the Ecoli model
ecoli_load = readSBML('iJO1366.xml', 1000);

% Change the objective function
% In base BIGG model, objective function = 0.9865
ecoli_load = changeObjective(ecoli_load, 'BIOMASS_Ec_iJO1366_WT_53p95M');

% Change all exchange reaction bounds that are not 1000, 0, or -1000
lbx = find(ecoli_load.lb ~= -1000 & ecoli_load.lb ~= 0);
ecoli_load.lb(lbx) = -1000;
ubx = find(ecoli_load.ub ~= 1000 & ecoli_load.ub ~= 0);
ecoli_load.ub(ubx) = 1000;

% Load in the transcript data
% m9 expression data
expressData = readtable('m9_aerobic_expression.csv');
% %Generate expressData variable
express_data.gene = expressData.gene;
express_data.value = expressData.reads;

% Map expression to reactions in the model
[expressionRxns parsedGPR] = mapExpressionToReactions(model, express_data);

%% Running algorithms using the createTissueSpecificModel function
% GIMME
clear options
model = ecoli_load;
options.solver = 'GIMME';
%   for GIMME
%       options.expressionRxns       reaction expression, expression data corresponding to model.rxns.
%                                    Note : If no gene-expression data are
%                                    available for the reactions, set the
%                                    value to -1
%       options.threshold            expression threshold, reactions below this are minimized
%       options.obj_frac*            minimum fraction of the model objective function
%                                    (default - 0.9)
options.expressionRxns = expressionRxns;
% use the median of the transcript data as the cutoff for transcript
% presence
options.threshold = 12;
options.obj_frac = 0.8;

[model_GIMME] = createTissueSpecificModel(model, options);

solution_GIMME = optimizeCbModel(model_GIMME)
%% iMAT or Shlomi et al algorithm
clear options
%    for iMAT
%       options.expressionRxns       reaction expression, expression data corresponding to model.rxns.
%                                    Note : If no gene-expression data are
%                                    available for the reactions, set the value to -1
%       options.threshold_lb         lower bound of expression threshold, reactions with
%                                    expression below this value are "non-expressed"
%       options.threshold_ub         upper bound of expression threshold, reactions with
%                                    expression above this value are
%                                    "expressed"
%       options.tol*                 minimum flux threshold for "expressed" reactions
%                                    (default 1e-8)
%       options.core*                cell with reaction names (strings) that are manually put in
%                                    the high confidence set (default - no core reactions)
%       options.logfile*             name of the file to save the MILP log (defaut - 'MILPlog')
%       options.runtime*             maximum solve time for the MILP (default - 7200s)
%       options.epsilon*             small value to consider when modeling
%                                    flux (default 1)
model = ecoli_load;
options.solver = 'iMAT';
options.expressionRxns = expressionRxns;
options.threshold_lb = 10;
options.threshold_ub = 900;

[model_iMAT] = createTissueSpecificModel(model, options);

% Model no longer has the objective function - need to add back in
add_rxns = {'BIOMASS_Ec_iJO1366_WT_53p95M'};
remove_rxns = {};
for rxn = 1:length(ecoli_load.rxns)
    if sum(strcmp(ecoli_load.rxns(rxn), model_iMAT.rxns)) == 1
    elseif sum(strcmp(ecoli_load.rxns(rxn), add_rxns)) == 1
    else
        remove_rxns{end+1,1} = ecoli_load.rxns{rxn};
    end
end

model_iMAT = removeRxns(ecoli_load, remove_rxns);
solution_iMAT = optimizeCbModel(model_iMAT);
%% fastcc to create flux consistent model
clear options

model = model;

% Create globally consistent network:
% checking for consistency of the network
% A - n x 1 boolean vector constaining consistent reactions
epsilon = 1e-4; % smallest flux considered non-zero
printLevel = 2; % summary print level
[A,modelFlipped,V] = fastcc(model, epsilon, printLevel);

% A - 5837 reactions/8336 reactions - number of reactions doesn't change
% based on epsilon
remove = setdiff(1:numel(model.rxns), A);
rxnRemoveList = model.rxns(remove);
base_model_fluxConsistent = removeRxns(model, rxnRemoveList);

%% MBA - using createTissueSpecificModel
% Use the FastCore model from above since MBA calls fastcc during algorithm
clear options
%   for MBA
%       options.medium_set           list of reaction names with medium confidence
%       options.high_set             list of reaction names with high confidence
%       options.tol*                 minimum flux threshold for "expressed" reactions
%                                    (default - 1e-8)

[A,modelFlipped,V] = fastcc(base_model_fluxConsistent, epsilon, printLevel);

% A - 5837 reactions/8336 reactions - number of reactions doesn't change
% based on epsilon
remove = setdiff(1:numel(base_model_fluxConsistent.rxns), A);
rxnRemoveList = base_model_fluxConsistent.rxns(remove);
base_model_MBA = removeRxns(base_model_fluxConsistent, rxnRemoveList);

model = base_model_MBA;

% Need to re-map expression to reactions
[expressionRxns parsedGPR] = mapExpressionToReactions(model, express_data);

options.solver = 'MBA';
options.medium_set = base_model_MBA.rxns(expressionRxns > 10 & expressionRxns < 900); % base_model_FastCore.rxns(expressionRxns_FastCore >= 0);
options.high_set = base_model_MBA.rxns(expressionRxns > 900);
options.tol = 1e-8;

[model_MBA] = createTissueSpecificModel(model, options);

% Model no longer has the objective function - need to add back in
add_rxns = {'BIOMASS_Ec_iJO1366_WT_53p95M'};
remove_rxns = {};
for rxn = 1:length(ecoli_load.rxns)
    if sum(strcmp(ecoli_load.rxns(rxn), model_MBA.rxns)) == 1
    elseif sum(strcmp(ecoli_load.rxns(rxn), add_rxns)) == 1
    else
        remove_rxns{end+1,1} = ecoli_load.rxns{rxn};
    end
end

model_MBA = removeRxns(ecoli_load, remove_rxns);
solution_MBA = optimizeCbModel(model_MBA);
%% CORDA 
% Function not in the new COBRA toolbox
%   metTests - metabolic tests to be included in the reconstruction. This
%       argument should be a cell array of strings of size nx2, where n is 
%       the number of metabolic tests to be performed. Column 1 should be 
%       the name of the reaction to be included and the corresponding column 
%       2 should be the reaction to be included. For example, to test for
%       the production of pep and pyruvate, metTests should be equal to
%       {'DM_pep[c]' 'pep[c] -> ';'DM_pyr[c]' 'pyr[c] -> '}. May be left
%       empty.
%   ES - High confidence reactions. Reactions to be included in the model. Cell
%       array of strings.
%   PR - Medium confidence reactions to be included in the model if they do 
%       not depend on too many NP reactions to carry a flux. Cell array of 
%       strings.
%   NP - Negatice confidence reactions not to be included in the model. 
%       These reactions will be included in the tissue model only if they 
%       are necessary for the flux of ES reactions or for the flux of PRtoNP 
%       or more PR reactions. Cell array of strings.

% Define medium reactions as those associated with medium-level proteins
% Classifications for reactions: -1, 1, 2
% ES: 2
% PR: 1
% NP: -1
% Map expression to reactions in the model

% Remove fields from structure, throws error with CORDA
%ecoli_load = rmfield(ecoli_load, 'csense');
%ecoli_load = rmfield(ecoli_load, 'osenseStr');

[expressionRxns parsedGPR] = mapExpressionToReactions(ecoli_load, express_data);

ES = ecoli_load.rxns(expressionRxns >= 900);
ES{end+1,1} = 'BIOMASS_Ec_iJO1366_WT_53p95M';
PR = ecoli_load.rxns(expressionRxns > 10 & expressionRxns < 900);
PR(650:675,:) = [];
NP = {};

constraint = 1e-4;
[model_CORDA, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA(ecoli_load, {}, ES,PR,NP, 2, constraint);
solution_CORDA = optimizeCbModel(model_CORDA);

%% Generate a model using RegrEx
base_model = base_model_fluxConsistent;
base_model = convertToIrreversible(base_model);

% Map expression to reactions in the model
[expressionRxns parsedGPR] = mapExpressionToReactions(base_model, express_data);

% add in reversible annotation
base_model.rev = 1;

% add in rxnGeneMat and grRules
base_model = buildRxnGeneMat(base_model);
base_model = creategrRulesField(base_model);

% Running the RegrEx method 
% GEM - COBRA model with no blocked reactions
% D - data mapped to reactions in the network which a value for each
% reaction
[Sol,ContextCOBRA] = RegrEx(base_model,expressionRxns);

% Model no longer has the objective function - need to add back in
RegrEx_rxns = ContextCOBRA.rxns;
RegrEx_rxns{end+1,1} = 'BIOMASS_Ec_iJO1366_WT_53p95M';
remove_rxns = setdiff(base_model.rxns, RegrEx_rxns);

model_RegrEx = removeRxns(base_model, remove_rxns);
solution_RegrEx = optimizeCbModel(model_RegrEx);