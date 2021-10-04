import os

WORKDIR=os.getcwd()
region=config["region"]
species=config["species"]

rule all:
    input:
        'Results/'+species+'/'+region+'/Plots/lr/ROC.pdf',
        'Results/'+species+'/'+region+'/Plots/rf/ROC.pdf',
        'Results/'+species+'/'+region+'/Plots/bt/ROC.pdf',
        'Results/'+species+'/'+region+'/Plots/lr/Feature_importance.pdf',
        'Results/'+species+'/'+region+'/Plots/rf/Feature_importance.pdf',
        'Results/'+species+'/'+region+'/Plots/bt/Feature_importance.pdf',
        'Results/'+species+'/'+region+'/Plots/Model_comparison.pdf',
        'Results/'+species+'/'+region+'/Models/lr_best.json',
        'Results/'+species+'/'+region+'/Models/rf_best.json',
        'Results/'+species+'/'+region+'/Models/bt_best.json',
        'Results/'+species+'/'+region+'/Plots/features.pdf'

rule create_negative_set:
    params:
        script = WORKDIR+"/scripts/createNegativeSet.py",
        out = WORKDIR+'/Results/'+species+'/'+region+'/Data/',
        beds = config['beds'],
        p = config['pval'],
        regions = region,
        maskY = config['maskY'],
        nd = config['negativeDistance']
    output:
        'Results/'+species+'/'+region+'/Data/exampleSet.csv'
    threads: 1
    resources:
        memory = 10,
        time = 1
    shell:
        "python3 {params.script} -o {params.out} -b {params.beds} -p {params.p} -r {params.regions} --maskY {params.maskY} --negativeDistance {params.nd}"

rule extract_sequences:
    input:
        fa = config['fasta'],
        sets = 'Results/'+species+'/'+region+'/Data/exampleSet.csv'
    output:
        'Results/'+species+'/'+region+'/Data/exampleSequences.csv'
    params:
        script = WORKDIR+"/scripts/extractSequences.py",
        out = WORKDIR+'/Results/'+species+'/'+region+'/Data/',
        window = config['window'],
        hybrid = config['hybrid']
    resources:
        memory = 20,
        time = 1
    threads: 1
    shell:
        "python3 {params.script} -o {params.out} -f {input.fa} -w {params.window} --hybrid {params.hybrid}"

rule create_general_features:
    input:
        'Results/'+species+'/'+region+'/Data/exampleSequences.csv'
    output:
        'Results/'+species+'/'+region+'/Data/generalFeatures.csv'
    params:
        script = WORKDIR+"/scripts/createGeneralFeatures.py",
        out = WORKDIR+'/Results/'+species+'/'+region+'/Data/',
        window = config['window'],
        td = config['thermodynamics'],
        di = config['dinuc_stacking'],
        vienna = config['vienna']
    resources:
        memory = 20,
        time = 1
    threads: 1
    shell:
        "python3 {params.script} -o {params.out} -w {params.window} --td {params.td} --di {params.di} -v {params.vienna}"

rule create_all_features:
    input:
        'Results/'+species+'/'+region+'/Data/generalFeatures.csv'
    output:
        'Results/'+species+'/'+region+'/Data/allFeatures.csv'
    params:
        script = WORKDIR+"/scripts/createHumanFeatures.py",
        out = WORKDIR+'/Results/'+species+'/'+region+'/Data/',
        window = config['window'],
        tf = config['tfs'],
        rbp = config['rbps'],
        shape = config['DNAshape'],
        species = species
    resources:
        memory = 20,
        time = 1
    threads: 1
    run:
        if params.species=='human':
            shell("python3 {params.script} -o {params.out} -w {params.window} --tf {params.tf} -m {params.rbp} --shape {params.shape}")
        else:
            shell("cp {params.out}generalFeatures.csv {params.out}allFeatures.csv")

rule plot_features:
    input:
        'Results/'+species+'/'+region+'/Data/allFeatures.csv'
    output:
        'Results/'+species+'/'+region+'/Plots/features.pdf'
    params:
        script = WORKDIR+"/scripts/plotFeatures.R"
    resources:
        memory = 20,
        time = 1
    threads: 1
    run:
        shell("Rscript {params} {input} {output}")

rule normalize_features:
    input:
        'Results/'+species+'/'+region+'/Data/allFeatures.csv'
    output:
        'Results/'+species+'/'+region+'/Features/normalizedFeatures.npy',
        'Results/'+species+'/'+region+'/Features/normalizedLabels.npy',
        'Results/'+species+'/'+region+'/Features/featureList.txt'
    resources:
        memory = 20,
        time = 4
    threads: 1
    script:
        "scripts/normalizedFeatures.py"

rule create_ML_sets:
    input:
        'Results/'+species+'/'+region+'/Features/normalizedFeatures.npy',
        'Results/'+species+'/'+region+'/Features/normalizedLabels.npy'
    output:
        'Results/'+species+'/'+region+'/Features/validationX.npy',
        'Results/'+species+'/'+region+'/Features/trainX.npy',
        'Results/'+species+'/'+region+'/Features/testX.npy',
        'Results/'+species+'/'+region+'/Features/validationY.npy',
        'Results/'+species+'/'+region+'/Features/trainY.npy',
        'Results/'+species+'/'+region+'/Features/testY.npy'
    resources:
        memory = 20,
        time = 1
    threads: 1
    script:
        "scripts/createMLSets.py"

rule perform_gridSearch:
    input:
        'Results/'+species+'/'+region+'/Features/validationX.npy',
        'Results/'+species+'/'+region+'/Features/validationY.npy'
    output:
        'Results/'+species+'/'+region+'/Models/lr_random.pickle',
        'Results/'+species+'/'+region+'/Models/rf_random.pickle',
        'Results/'+species+'/'+region+'/Models/bt_random.pickle'
    params:
        threads = 80
    resources:
        memory = 900,
        time = 14
    threads: 80
    script:
        "scripts/gridSearch.py"

rule find_best_model:
    input:
        'Results/'+species+'/'+region+'/Models/lr_random.pickle',
        'Results/'+species+'/'+region+'/Models/rf_random.pickle',
        'Results/'+species+'/'+region+'/Models/bt_random.pickle'
    output:
        'Results/'+species+'/'+region+'/Validation/Validation_scores.csv',
        'Results/'+species+'/'+region+'/Plots/Validation.pdf'
    resources:
        memory = 20,
        time = 1
    threads: 10
    script:
        "scripts/findBestModel.py"

rule train_best_model:
    input:
        'Results/'+species+'/'+region+'/Models/{model}_random.pickle',
        'Results/'+species+'/'+region+'/Features/trainX.npy',
        'Results/'+species+'/'+region+'/Features/trainY.npy',
        'Results/'+species+'/'+region+'/Features/testX.npy',
        'Results/'+species+'/'+region+'/Features/testY.npy'
    output:
        'Results/'+species+'/'+region+'/Models/{model}_best.pickle',
        'Results/'+species+'/'+region+'/Plots/{model}/PR.pdf',
        'Results/'+species+'/'+region+'/Plots/{model}/ROC.pdf',
        'Results/'+species+'/'+region+'/Models/{model}_best.json'
    resources:
        memory = 20,
        time = 1
    threads: 10
    script:
        "scripts/trainBestModel.py"

rule permutation_importance:
    input:
        'Results/'+species+'/'+region+'/Models/{model}_best.pickle',
        'Results/'+species+'/'+region+'/Features/testX.npy',
        'Results/'+species+'/'+region+'/Features/testY.npy',
        'Results/'+species+'/'+region+'/Features/featureList.txt'
    output:
        'Results/'+species+'/'+region+'/FI/{model}_feature_imporance.npy',
        'Results/'+species+'/'+region+'/Plots/{model}/Feature_importance.pdf',
    params:
        jobs = 20
    resources:
        memory = 200,
        time = 3
    threads: 20
    script:
        "scripts/computeFeatureImportance.py"

rule compare_best_models:
    input:
        'Results/'+species+'/'+region+'/Features/trainX.npy',
        'Results/'+species+'/'+region+'/Features/trainY.npy',
        'Results/'+species+'/'+region+'/Features/testX.npy',
        'Results/'+species+'/'+region+'/Features/testY.npy',
        'Results/'+species+'/'+region+'/Models/lr_best.pickle',
        'Results/'+species+'/'+region+'/Models/rf_best.pickle',
        'Results/'+species+'/'+region+'/Models/bt_best.pickle'
    output:
        'Results/'+species+'/'+region+'/Plots/Model_comparison.pdf'
    resources:
        memory = 20,
        time = 1
    threads: 1
    script:
        "scripts/compareBestModels.py"
