# CoLiPE - Comprehensive Link Prediction Evaluation
A framework for bias-aware evaluation of link prediction algorithms focusing on under-studied proteins. 

### Running the evaluation
Please run the following script first to perform the evaluation of the link prediction algorithms:
```
demo_00_run_evaluation.m
```

### View 1: Bias in algorithms
This view investigates the disposition of the algorithms towards well-studied (high-degree) nodes. To perform this analysis, please run: 
```
demo_01_bias_in_algorithms
```

### View 2: Bias in benchmarking data
This view investigates (in imbalance analysis) how much the under-studied entities are under-represented in the benchmarking data (i.e., training/test splits) and (in separability analysis) how much an evaluation setting favors algorithms that bring forward high-degree nodes. To perform this analysis, please run: 
```
demo_02_imbalance_analysis
demo_02_separability_analysis
```

### View 3: Balanced evaluation to prevent under-studied entities from being under-represented
This view measures the prediction performance of the algorithms on standard (unweighted, edge-uniform) and balanced (weighted, node-uniform) settings. To perform this analysis, please run: 
```
demo_03_evaluation_barplots
demo_03_evaluation_pr_curves
```

### View 4: Stratified analysis to focus on under-studied entities
This view decomposes the prediction performance into different categories based on the degrees of the involved nodes for the predicted interactions. Thus, this analysis aims to give insight into the discovery chances of new interactions involving under-studied proteins with the use of a given algorithm. To perform this analysis, please run: 
```
demo_04_stratified_heatmaps
demo_04_stratified_pr_curves
```

### View 5: Simple and comprehensive summary for bias-aware evaluation
To make a comprehensive and bias-aware evaluation in a simple manner, this view considers five aspects: Early/late curve prediction performance for under-studied/well-studied nodes, and the disposition of an algorithm regarding degree bias. To perform this analysis, please run: 
```
demo_05_summary_spiderplot_5way
```



