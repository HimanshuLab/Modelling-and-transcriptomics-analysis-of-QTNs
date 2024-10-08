function CreateSNPSpecificModel(MeM,threshold)
changeCobraSolverParams('LP', 'feasTol', 1e-9);
% loading the Yeast 8 model
model = readCbModel('yeastGEM.mat'); 
model = buildRxnGeneMat(model);
% Constriaining the exchange reactions
model = changeRxnBounds(model,'r_1106',-1000,'l'); % acetate[c] to acetate [e]
model = changeRxnBounds(model,'r_1106',0,'u'); % acetate[c] to acetate [e]
model = changeRxnBounds(model,'r_1634',-1000,'l'); % Exchange reaction of acetate[e]
model = changeRxnBounds(model,'r_1634',0,'u'); % Exchange reaction of acetate[e]

model = changeRxnBounds(model,'r_1992',-1000,'l'); % Exchange reaction of oxygen
model = changeRxnBounds(model,'r_1979',-1000,'l'); % oxygen[c] to oxygen[e]

model = changeRxnBounds(model,'r_1714',0,'l'); % Glucose exchange
model = changeRxnBounds(model,'r_1115',0,'l'); % Ammonia exchange


geneExpData = readtable('tpm_counts_Average.csv'); % The gene expression data
expData.value = geneExpData(:,2:17);
expData.value = table2array(geneExpData(:,2:17));
expData.genes = table2cell(geneExpData(:,1));
context = geneExpData.Properties.VariableNames;
context(:,1) = [];
expData.context = context;
% for LocalGini Thresholding
if strcmp(threshold, 'LocalGini')
    
    ut = 90; % Upper threshold
    lt = 10; % Lower threshold
    ThS = 1; % impliying at gene level
    %Tolerance level above which reactions are considered as expressed
    tol = 1e-8;
    mkdir(['./LocalGini/',MeM])
    filename = ['./LocalGini/',MeM,'/']; 
    % Giving the biomass, nucleotide synthesis reaction and maintanence
    % reaction a higher importance
    biomass_id = find(strcmp(model.rxns,'r_2111'));
    nucleotide_synthesis_id = find(strcmp(model.rxns,'r_0466'));
    NGAM_id = find(strcmp(model.rxns,'r_4046'));
    coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
    [Models,RxnImp] = buildContextmodels(expData,model,MeM,contexts,ut,lt,ThS,coreRxn,filename,fastcc(model,1e-8),tol);
elseif strcmp(threshold,'LocalT1')
    Gmodel = model;
    modelData.valuebyTissue = expData.value;
    modelData.gene = expData.genes;
    modelData.Tissue = expData.context;
    modelData = getModelData(modelData,Gmodel);

    [thrRxnData,rxnTisMat] = getLocalT1_case(modelData,Gmodel,90);
    % Giving the biomass, nucleotide synthesis reaction and maintanence
    % reaction a higher importance
    biomass_id = find(strcmp(Gmodel.rxns,'r_2111'));
    nucleotide_synthesis_id = find(strcmp(Gmodel.rxns,'r_0466'));
    NGAM_id = find(strcmp(Gmodel.rxns,'r_4046'));
    coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
    % maximum weights to manually added core reactions
    thrRxnData.value(coreRxn,:) = 10*log(2);
    % for those reaction for which there is no expression level in the data
    thrRxnData.value(isnan(thrRxnData.value)) =-2;
    
    if strcmp(MeM,'FASTCORE')
        core = thrRxnData.value>=5*log(2);
        a= fastcc(Gmodel,1e-8);
        ConsModel = removeRxns(Gmodel,Gmodel.rxns(setdiff(1:numel(Gmodel.rxns),a)));
        core = core(a,:);
        filename='./LocalT1/FASTCORE/';
        mkdir(filename)
        
        for j =1:size(core,2)
            model = fastcore(ConsModel,find(core(:,j)), 1e-8);
            save(['./LocalT1/FASTCORE/',contexts{j}],'model')
        end
    elseif strcmp(MeM,'iMAT')
        core = thrRxnData.value;
        filename='./LocalT1/iMAT/';
        mkdir(filename)
        for j =1:size(core,2)
            model = iMAT(Gmodel,core(:,j),5*log(2),5*log(2),1e-8,{});
            save(['./LocalT1/iMAT/',contexts{j}],'model')
        end
    elseif strcmp(MeM,'INIT')
        filename = './LocalT1/INIT/';
        mkdir(filename)
        core =thrRxnData.value;
        core(coreRxn,:) =max(max(core));
        core(find(sum(Gmodel.rxnGeneMat,2)==0),:)=0;
        for j =1:size(core,2)
            model = INIT(Gmodel,core(:,j),1e-8);
            save(['./LocalT1/INIT/',contexts{j}],'model')
        end
    end

elseif strcmp(threshold,'LocalT2')
    Gmodel = model;
    modelData.valuebyTissue = expData.value;
    modelData.gene = expData.genes;
    modelData.Tissue = expData.context;
    modelData = getModelData(modelData,Gmodel);
    [thrRxnData,rxnTisMat] = getLocalT2_case(modelData,Gmodel,10,90);
    % Giving the biomass, nucleotide synthesis reaction and maintanence
    % reaction a higher importance
    biomass_id = find(strcmp(Gmodel.rxns,'r_2111'));
    nucleotide_synthesis_id = find(strcmp(Gmodel.rxns,'r_0466'));
    NGAM_id = find(strcmp(Gmodel.rxns,'r_4046'));
    coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
    % maximum weights to manually added core reactions
    thrRxnData.value(coreRxn,:) = 10*log(2);
    % for those reaction for which there is no expression level in the data
    thrRxnData(isnan(thrRxnData.value)) =-2;
    if strcmp(MeM,'FASTCORE')
        core = thrRxnData.value>=5*log(2);
        a= fastcc(Gmodel,1e-8);
        ConsModel = removeRxns(Gmodel,Gmodel.rxns(setdiff(1:numel(Gmodel.rxns),a)));
        core = core(a,:);
        filename='./LocalT2/FASTCORE/';
        mkdir(filename)
        
        for j =1:size(core,2)
            model = fastcore(ConsModel,find(core(:,j)), 1e-8);
            save(['./LocalT2/Fastcore/',contexts{j}],'model')
        end
    elseif strcmp(MeM,'iMAT')
        core = thrRxnData.value;
        filename='./LocalT2/iMAT/';
        mkdir(filename)
        for j =1:size(core,2)
            model = iMAT(Gmodel,core(:,j),5*log(2),5*log(2),1e-8,{});
            save(['./LocalT2/iMAT/',contexts{j}],'model')
        end
    elseif strcmp(MeM,'INIT')
        filename = './LocalT2/INIT/';
        mkdir(filename)
        core =thrRxnData.value;
        core(coreRxn,:) =max(max(core));
        core(find(sum(Gmodel.rxnGeneMat,2)==0),:)=0;
        for j =1:size(core,2)
            model = INIT(Gmodel,core(:,j),1e-8);
            save(['./LocalT2/INIT/',contexts{j}],'model')
        end
    end

elseif strcmp(threshold,'StanDep')
    Gmodel = model;
    modelData.valuebyTissue = expData.value;
    modelData.gene = expData.genes;
    modelData.Tissue = expData.context;
    modelData = getModelData(modelData,Gmodel);

    spec = getSpecialistEnzymes(Gmodel);  
    prom = getPromEnzymes(Gmodel);
    enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

    edgeX = [-3:6]; % bins  
    distMethod = 'euclidean'; % distance method  
    linkageMethod = 'complete'; % linkage metric for hierarchical clustering

    clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,20,distMethod,linkageMethod);
    if strcmp(MeM,'FASTCORE')
        core= models4mClusters1(clustObj,enzymeData.Tissue,Gmodel,edgeX,[],[],false,0,[1 1]);
        biomass_id = find(strcmp(Gmodel.rxns,'r_2111'));
        nucleotide_synthesis_id = find(strcmp(Gmodel.rxns,'r_0466'));
        NGAM_id = find(strcmp(Gmodel.rxns,'r_4046'));
        coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
        core(coreRxn,:)=true;
        a = fastcc(Gmodel,1e-8);
        ConsModel = removeRxns(Gmodel,Gmodel.rxns(setdiff(1:numel(Gmodel.rxns),a)));
        core = core(a,:);
        filename='./StanDep/FASTCORE/';
        mkdir(filename)
        for j =1:size(core,2)
            model = fastcore(ConsModel,find(core(:,j)),1e-8);
            save(['./StanDep/FASTCORE/',contexts{j}],'model')
        end
    elseif strcmp(MeM,'iMAT')
        core= models4mClusters1(clustObj,enzymeData.Tissue,Gmodel,edgeX,[],[],false,0,[1 1]);
        biomass_id = find(strcmp(Gmodel.rxns,'r_2111'));
        nucleotide_synthesis_id = find(strcmp(Gmodel.rxns,'r_0466'));
        NGAM_id = find(strcmp(Gmodel.rxns,'r_4046'));
        coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
        core(coreRxn,:)=true;
        core = double(core);
        core(find(sum(Gmodel.rxnGeneMat,2)==0),:)=0.5;
        filename='./StanDep/iMAT/';
        mkdir(filename)
        for j =1:size(core,2)
           model = iMAT(Gmodel,core(:,j),0.2,0.7,1e-8,{});
           save(['./StanDep/iMAT/',contexts{j}],'model')
        end
    elseif strcmp(MeM,'INIT')
        core = getINITweights(clustObj,edgeX,Gmodel);
        core(core==-inf)=min(min(core(core~=-inf)));
        core(core==inf)=max(max(core(core~=inf)));
        core(isnan(core))=0;
        biomass_id = find(strcmp(Gmodel.rxns,'r_2111'));
        nucleotide_synthesis_id = find(strcmp(Gmodel.rxns,'r_0466'));
        NGAM_id = find(strcmp(Gmodel.rxns,'r_4046'));
        coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
        core(coreRxn,:)=1;
        core(find(sum(Gmodel.rxnGeneMat,2)==0),:)=0;
        filename='./StanDep/INIT/';
        mkdir(filename)
        for j =1:size(core,2)
           model = INIT(Gmodel,core(:,j),1e-8);
           save(['./StanDep/INIT/',contexts{j}],'model')
        end
    end

elseif strcmp(threshold,'Top10')
    Gmodel = model;
    modelData.valuebyTissue = expData.value;
    modelData.gene = expData.genes;
    modelData.Tissue = expData.context;
    modelData = getModelData(modelData,Gmodel);
    parsedGPR = GPRparser(Gmodel);
    for i=1:numel(contexts)
        rxn_exp(:,i) = selectGeneFromGPR(Gmodel,modelData.gene,modelData.value(:,i), parsedGPR); % gene to reaction mapping
    end
    biomass_id = find(strcmp(Gmodel.rxns,'r_2111'));
    nucleotide_synthesis_id = find(strcmp(Gmodel.rxns,'r_0466'));
    NGAM_id = find(strcmp(Gmodel.rxns,'r_4046'));
    coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
    rxn_exp(coreRxn,:) = max(max(rxn_exp)); 
    rxn_exp(find(sum(Gmodel.rxnGeneMat,2)==0),:)=0;
    if strcmp(MeM,'FASTCORE')
        core = rxn_exp;
        a= fastcc(Gmodel,1e-8);
        ConsModel = removeRxns(Gmodel,Gmodel.rxns(setdiff(1:numel(Gmodel.rxns),a)));
        filename='./Top10/FASTCORE/';
        mkdir(filename)
        for i=1:numel(contexts)
            model = fastcore(ConsModel,find(core(:,i)>prctile(core(:,i),90)),1e-8);
            save(['./Top10/FASTCORE/',contexts{i}],'model')
        end
    elseif strcmp(MeM,'INIT')
        core = rxn_exp;
        filename='./Top10/INIT/';
        mkdir(filename)
        for i=1:numel(contexts)
            temp =5*log2(core(:,i)/prctile(core(:,i),90));
            temp(temp==-inf)=min(temp(temp~=-inf));
            model = INIT(Gmodel,temp,1e-8);
            save(['./Top10/INIT/',contexts{i}],'model')
        end

    elseif strcmp(MeM,'iMAT')
        core = rxn_exp;
        filename='./Top10/iMAT/';
        mkdir(filename)
        for i=1:numel(contexts)
            temp = core(:,i);
            temp = min(temp(temp~=0));
            % lower threshold excludes 0 expression reactions alone and
            % upper threshold is min of top 10 percentile 
           model = iMAT(Gmodel,core(:,i),temp,prctile(core(:,i),90),1e-8,{});
           save(['./Top10/iMAT/',contexts{i}],'model')
        end
    end

end
end
    
