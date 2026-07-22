# Drug Effects Predictive Model Machine Learning
**Method summary**

This R script implements the in-silico treatment predictive machine learning using a Random Forest (RF) multi-class classification approach (see Figure 1) described in Palm et,al (2026) (Publishing pending). 
In this study we modelled the effects of multiple perturbagens by training a Random Forest multi-class classification model with existing human bulk transcriptomic data (randomForest 4.7-1.1 package) to distinguish between the five distinct transcriptomic plaque subtypes using the 5000-plaque subtype defining genes recently identified by Mokry et al., (2023)4 (Figure 6). These genes were identified using shared nearest neighbour modularity optimization -based clustering algorithm implemented in the Seurat package 21. To develop this model in RStudio, we first identified differentially expressed genes (DEGs, FDR < 0.05) following in-vitro drug experiments or extracted signatures from public omics databases (Pvalue<0.05). Only genes that overlapped with both the 5,000-plaque subtype-defining genes, were present in the human Bulk RNA dataset, and were selected by LASSO were used for model training. For model training, the data was partitioned into 1,000 randomly generated subsets, each split into 80% training and 20% testing sets. Trained ML model performance was evaluated using multiclass area under the ROC curve (AUC) from the pROC 1.18.5 package, with an average AUC of 75% or higher considered indicative of a well performing model.

To simulate in-silico drug treatment, the expression changes of the same DESeq2-derived, LASSO-selected DEGs were projected onto human transcriptomic plaque data. This was achieved by multiplying the expression values by the antilogarithmic transformed fold changes, resulting in a “perturbed” dataset. The trained predictive model was then applied to classify the perturbed dataset, and these classifications were used to quantify the effects of the in-silico drug treatment. After 1,000 iterations, an averaged output was generated, revealing predicted changes in plaque classification after in-silico treatment. 


Libraries used:
progress 1.2.3,
pROC 1.18.5,
caret 6.0-94,
glmnet 4.1-8,
randomForest 4.7-1.1,
EnsDb.Hsapiens.v86 2.99.0,
biomaRT 2.60.1,
ggalluvial 0.12.5,
limma 3.60.4,
ggplot2 3.5.1
dplyr 1.1.4

Figure 1 Overview of model building 
<img width="4100" height="4584" alt="General Drug VAR flowchart" src="https://github.com/user-attachments/assets/c2a0e9b4-23d5-449a-9b6e-386a4fb06924" />

Numbered boxes indicate corresponding phase in step-by-step description. 

## Getting started: step by step description 
Before initiating the predictive Random Forest model loop, two main dataframes need to be prepared:

### 1. Before treatment data frame 
This data frame will be used to **train and test** the RF model. It contains the original normalized bulk RNA expression of the patient population, which will be used to train the mode. Additionally, it includes the known plaque subtype for each patient. 

### 2. After in-silico treatment data frame
This data frame represents data after simulated (in-silico) treatment. It will be **fed into the trained RF model** to make treatment predictions based on the learned patterns from the training phase. It **does not** contain any plaque subtype information. 


## Steps to prepare data frames (script 1)

1. **Load and normalize data**:  
   The raw bulk RNA dataset is loaded (with patients as rows and genes as columns), scaled, and quantile normalized using the 'limma' package and transposed 

2. **Load in-vitro cell drug treatment data**:  
  The in-vitro cell drug treatment expression data (for the drug of interest or controls) is loaded and saved in a table with the following variables: **`gene`** (ENSG), **`symbol`**, and **`FoldChange`**. Genes that are not differentially expressed after treatment (FRD>0.1) are filtered out. The variable **`FoldChange`** is created by transforming the existing log2 fold change values using the formula:
   
   $$
   \text{FoldChange} = 2^{\text{log fold change}}
   $$

3. **Filter normalized bulk RNA data frame**:  
   The normalized bulk RNA data frame is filtered for plaque subtype-defining genes (No. = 5,000) and genes that were differentially expressed after the in-vitro drug experiment (FDR < 0.1). Additionally, genes with low counts (count <= 10) are excluded.

4. **Create after in-silico treatment data frame**:  
   The in-silico treatment data frame is created by multiplying the expression values from the 'Before Treatment DataFrame' with the corresponding gene ‘FoldChange’ values from the in-vitro drug experiment table. 

## RF model building and performance evaluation (script 2)

5. **Initialize RF Prediction model loop**:  
   The RF model is trained and tested using the 'Before Treatment DataFrame' (See Figure 2). Model performance is evaluated using the Area Under the Curve (AUC), and the predicted plaque subtype is calculated by feeding the trained RF model with the **`After In-Silico Treatment DataFrame.`** To overcome overfitting, the RF prediction model is performed 1,000 times, with each iteration utilizing a different subset of the patient population created using createDataPartition(y=..., p=0.80, times = 1000, list = TRUE). This results in a final output that represents the average prediction over 1,000 separate models. 
 

Figure 2: Zoomed-in overview of model building and performance evaluation 
<img width="2563" height="2897" alt="Var model building and QC" src="https://github.com/user-attachments/assets/65e86357-949c-4565-b238-f384e836b8ee" />


## Model results & visualization (script 3)
6. Differences in the proportion of plaque subtype before and after in-silico treatment were assessed Stuart-Maxwell or McNemar’s Chi-squared test, which evaluates paired data to determine whether the observed shifts (e.g, from more vulnerable to less vulnerable plaques) are meaningful. To further interpret results, odds ratios (ORs) were calculated by dividing the number of plaques for a group of interest before in-silico treatment (e.g., plaque subtype 2) by the number of plaques predicted to shift after in-silico treatment. ORs were subsequently log-transformed, where values > 0 represent a higher log odd for that selected group of interest. Corresponding confidence intervals were calculated using standard formulas and p-values were derived from McNemar's Chi-squared test. To assess the magnitude of the drug effect on plaque subtype reclassification, net marginal change for each subtype was calculated as the difference between the number of patients transitioning into and out of that subtype. Then, net percentage change (NPC) was calculated by dividing the net transition count by the total study population. Plaque shifts before and after in-silico treatments were visualized using galluvial 0.12.5 package. 

Formulas used to calculate Odds ratios and corresponding upper and lower CIs:
   a= No. shifted less vulnerable plaques after in-silico treatment  
   b= No. shifted more vulnerable plaques after in-silico treatment
   
   $$
   \text{Odds Ratio} = \frac{a}{b}
   $$
 
   $$
   \text{SE}_{\log(\text{OR})} = \sqrt{\frac{1}{a} + \frac{1}{b}}
   $$

   $$
   \log(\text{OR}) = \log(\text{Odds Ratio})
   $$
   
   $$
   \text{Lower CI} = \exp(\log(\text{OR}) - 1.96 \times \text{SE}_{\log(\text{OR})})
   $$
  
   $$
   \text{Upper CI} = \exp(\log(\text{OR}) + 1.96 \times \text{SE}_{\log(\text{OR})})
   $$

# RF Example Output

![image](https://github.com/user-attachments/assets/3faa77ec-978b-4137-af16-4f63c2fe438e)

![image](https://github.com/user-attachments/assets/c488dafd-de8f-44c5-a66e-aa2bb9656cf1)

![image](https://github.com/user-attachments/assets/ca400bd0-5556-4d4b-91dc-add4431b6359)







