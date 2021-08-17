%% Community Modelling
%% Startup 
initCobraToolbox(false)

%Use cplex 12.9 but also have 12.8 installed 
solverOK=changeCobraSolver('ibm_cplex','all');

projectPath = 'C:\Users\maria_2kg7ogk\Documents\GitHub\community_modelling_';
cd(projectPath);

%% Download AGORA Models
%AGORA 1.03 - without mucin
system('curl -O https://www.vmh.life/files/reconstructions/AGORA/1.03/Agora-1.03.zip')
unzip('Agora-1.03.zip','AGORA103')

%% Creation of Pan Models
%Path to AGORA models 
agoraPath = [projectPath filesep 'AGORA-1.03-Mucin' filesep 'mat' filesep]; 

%Directory for PanModels
mkdir('panModels');
panPath = [projectPath filesep 'panModels' filesep]; 

%Taxon level of the Pan Models
taxonLevel = 'Genus'; 

%Create function 
createPanModels(agoraPath,panPath,taxonLevel);

%% Meta Data Sorting 

data = DEcancerhealthy; 
metadata = SOCSmetadata; 

cancer = string.empty; 
no_cancer = string.empty; 

m = 1;
n = 1; 

for i=1:length(metadata)
    for j=1:size(data,2) 
    k = contains(data(1,j),metadata(i,1));
        if k == 1 
            index = j;
            break
        end
    end
    if metadata(i,2) == "cancer"
        cancer(:,m) = data(:,index);
        m = m +1; 
    elseif metadata(i,2) == "no_cancer"
        no_cancer(:,n) = data(:,index);
        n = n +1;
    end
    disp(m)
    disp(n)
end 

writetable(array2table(cancer), 'DE_sorted_cancer.csv')
writetable(array2table(no_cancer), 'DE_sorted_healthy.csv')

%% Translate Metagenome to AGORA(pan) 

%Inputs
data = {'CRC', 'CRC_rep', 'DE', 'DE_rep'};
data_files = {'CRC_data_map2AGORA.csv', 'CRC_data_map2AGORA_cancer_healthy.csv', 'DE_cancer_healthy_2AGORA.csv', 'DE_cancer_healthy_2AGORA_rep.csv'};

sequencingDepth = 'g__';
i = 1;

%Function
MetagenomeAbundancePath = [projectPath filesep data_files{i}];
[translatedAbundances,normalizedAbundances,unmappedRows] = translateMetagenome2AGORA(MetagenomeAbundancePath, sequencingDepth);

%Save
writetable(cell2table(normalizedAbundances), strcat('normalizedAbundances_', data{i}, '.csv'), 'writevariablenames', false)

%% MgPipe

%Inputs - Consistent
objre={'EX_biomass(e)'}; %Objective function
figForm = '-dpng'; % the output is a vectorized picture, change to '-dpng' for .png
numWorkers = 3;% number of cores dedicated for parallelization
autoFix = true; % autofix for names mismatch
compMod = false; % if outputs in open formats should be produced for each section 
indInfoFilePath='none'; % if documentation (.csv) on stratification criteria is available
rDiet = false; % to enable also rich diet simulations 
extSolve = false; % if to use an external solver and save models with diet
fvaType = true; % the type of FVA function to use to solve
autorun = true; % to turn off the autorun to be able to manually execute each part of the pipeline
global CBTDIR
dietFilePath=[CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'resources' filesep 'AverageEuropeanDiet'];
modPath = [projectPath filesep 'panModels'];

%Inputs - Varied 
data = {'CRC', 'CRC_rep', 'DE', 'DE_rep'};
results_folders = {'results_CRC_individual', 'results_CRC_rep', 'results_DE_individual', 'results_DE_rep'};

%Change i in order to change data input
i = 1;

abunData = [projectPath filesep strcat('normalizedAbundances_', data{i}, '.csv')]; % char with path and name of file from which to retrieve information
abunFilePath = [projectPath filesep strcat('normalizedAbundances_', data{i}, '.csv')];

mkdir(string(results_folders{i}));
resPath = [projectPath filesep string(results_folders{i})];

%Function 
[init,modPath,~,resPath,dietFilePath,abunData,indInfoFilePath,objre,figForm,numWorkers,autoFix,compMod,rDiet,extSolve,fvaType,autorun]= initMgPipe(modPath, CBTDIR, resPath, dietFilePath, abunFilePath, indInfoFilePath, objre, figForm, numWorkers, autoFix, compMod, rDiet,extSolve,fvaType,autorun);

%% Simulations - MgPipe 

[ID,fvaCt,nsCt,presol,inFesMat] = microbiotaModelSimulator(resPath,setup,sampName,dietFilePath,rDiet,0,extSolve,patNumb,fvaType); 
[Fsp,Y]= mgSimResCollect(resPath,ID,sampName,rDiet,0,patNumb,indInfoFilePath,fvaCt,figForm);
[finRes] = extractFullRes(resPath, ID, 'sDiet', sampName, fvaCt, nsCt);

%% Pairwise Modelling

%Loading the models
healthy_colon = readCbModel('model_healthy_colon.mat');
cancer_colon = readCbModel('model_cancer_colon.mat'); 
healthy_microbiota_DE = readCbModel('microbiota_model_samp_healthy.mat');
cancer_microbiota_DE = readCbModel('microbiota_model_samp_cancer.mat');

%adding model ID for the name tag of the host 
healthy_colon.modelID = 'healthy'; 
cancer_colon.modelID = 'cancer'; 

%model arrays 
colon_models = {healthy_colon, cancer_colon};
microbiota_models = {healthy_microbiota_DE, cancer_microbiota_DE} ;

%Creation of Host-Microbiome Join Models (createMultipleSpeciesModel) 
models={};
nameTagsModels={};
bioID={};

for i=1:length(colon_models)
    display(i)
    modelHost= colon_models{i}; 
    nameTagHost = modelHost.modelID; 
    for j=1:length(microbiota_models)
        display(j)
        models{1,1}= microbiota_models{j};
        nameTagsModels{1,1} =  char(strcat(models{1,1}.name));
        bioID{1,1} = models{1,1}.rxns(find(strncmp(models{1,1}.rxns, 'biomass',7)));
        [modelJoint] = createMultipleSpeciesModel(models,'nameTagsModels',nameTagsModels,'modelHost',modelHost,'nameTagHost',nameTagHost,'mergeGenesFlag',false);
        = modelJoint
        model_name = strcat('modelJoint_DE_rep_', nameTagHost, '_', nameTagsModels{1,1}, '_microbiota');
        modelJoint.name = model_name;
        save([model_name '.mat'], 'modelJoint')
    end
end 

%% Simulations - Pairwise Models

models = {modelJoint_healthy_colon_healthy_micro, modelJoint_healthy_colon_cancer_micro, modelJoint_cancer_colon_healthy_micro, modelJoint_cancer_colon_cancer_micro}; 

results = cell(1,4);

for i=1:length(models)
    display(i)
    %checking for cancerous or healthy models in order to adjust objective approriately 
    if contains(models{i}.name,'cancer_colon')
        objective = 'cancerbiomass_reaction';
    else 
        objective = 'healthybiomass_reaction'; 
    end
    %change objective 
    model_temp = changeObjective(models{i},objective);
    %simulate 
    solution =solveCobraLP(model_temp);
    %optimizeCbModel(model_temp)
    %save result 
    results{i} = solution.obj;
end 

%compiling results 
results_final = table('Size', [4,2], 'VariableTypes', {'string', 'double'});
results_final.Properties.VariableNames = {'Model' 'Biomass'};

for j=1:length(models)
    results_final(j,1) = {models{j}.name};
    results_final(j,2) = results(j);
end 

disp(results_final)
writetable(cell2table(results_final), 'biomass_cancer_colon_dif_micro_results.csv')


%% Simulations - Pairwise Models

models = {modelJoint_cancer_colon_healthy_micro, modelJoint_cancer_colon_cancer_micro};

for i=1:length(models)
    model_temp = changeObjective(models{i},'cancerbiomass_reaction');
    %FVA
    [minFlux,maxFlux,optsol,ret]=fastFVA(model_temp,90,'max','ibm_cplex');
    %[minFlux, maxFlux, Vmin, Vmax] = fluxVariability(model_temp,90,'max',model_temp.rxns,0,true,'FBA');
    %Save results 
    writetable([minFlux, maxFlux, optsol, ret], stringcat('pairwise_fastFVA_', models{i}, '.csv')
end 
