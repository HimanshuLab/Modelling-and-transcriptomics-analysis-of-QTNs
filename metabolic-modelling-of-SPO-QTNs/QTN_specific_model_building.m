clear
% initiating COBRA toolbox
initCobraToolbox(false)
changeCobraSolver('gurobi','all') % Changing the COBRA solver to gurobi
model = readCbModel('yeastGEM.mat'); % loading the Yeast 8 model

%% Constriaining the exchange reactions
model = changeRxnBounds(model,'r_1106',-1000,'l'); % acetate[c] to acetate [e]
model = changeRxnBounds(model,'r_1106',0,'u'); % acetate[c] to acetate [e]
model = changeRxnBounds(model,'r_1634',-1000,'l'); % Exchange reaction of acetate[e]
model = changeRxnBounds(model,'r_1634',0,'u'); % Exchange reaction of acetate[e]

model = changeRxnBounds(model,'r_1992',-1000,'l'); % Exchange reaction of oxygen
model = changeRxnBounds(model,'r_1979',-1000,'l'); % oxygen[c] to oxygen[e]

model = changeRxnBounds(model,'r_1714',0,'l'); % Glucose exchange
model = changeRxnBounds(model,'r_1115',0,'l'); % Ammonia exchange


%% To get IMAT models using localgini thresholding
geneExpData = readtable('tpm_counts_Average.csv'); % The gene expression data
expData.value = geneExpData(:,2:17);
expData.value = table2array(geneExpData(:,2:17));
expData.genes = table2cell(geneExpData(:,1));
context = geneExpData.Properties.VariableNames;
context(:,1) = [];
expData.context = context;
% Model extraction method
MeM = 'iMAT';
contexts = {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16'};
ut = 90; % Upper threshold
lt = 10; % Lower threshold
ThS = 1; % impliying at gene level
%Tolerance level above which reactions are considered as expressed
tol = 1e-4;
filename = './'; 
% Giving the biomass, nucleotide synthesis reaction and maintanence
% reaction a higher importance
biomass_id = find(strcmp(model.rxns,'r_2111'));
nucleotide_synthesis_id = find(strcmp(model.rxns,'r_0466'));
NGAM_id = find(strcmp(model.rxns,'r_4046'));

coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
[Models,RxnImp] = buildContextmodels(expData,model,MeM,contexts,ut,lt,ThS,coreRxn,filename,tol);

% Removing the unused genes in the built models
 for i=1:16
    temp_var = removeUnusedGenes(Models{i});
    Models{i} = temp_var;
end

C1_model = Models{1};
C2_model = Models{2};
C3_model = Models{3};
C4_model = Models{4};
C5_model = Models{5};
C6_model = Models{6};
C7_model = Models{7};
C8_model = Models{8};
C9_model = Models{9};
C10_model = Models{10};
C11_model = Models{11};
C12_model = Models{12};
C13_model = Models{13};
C14_model = Models{14};
C15_model = Models{15};
C16_model = Models{16};




