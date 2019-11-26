#run the following script after creating hypermutation files with the associated R file then go back to that file to plot the results
#FOR python2
import pandas as pd
import os
from collections import Counter

pathPrefix = '/Users/friedman/Desktop/mnt' #replace with '' if you are running this script on the cluster without a mount

#validates that the hypermutation cohort identified by the clustering fulfills two additional conditions:
#1. the majority signature in the hypermutated cases is distinct from the majority signature in the non-hypermutated cases
#2. the hypermutated cases represent <50% of the total cases

###NOTE: in the file /juno/work/taylorlab/friedman/myAdjustedDataFiles/impactSignatureCalls_Nov20_2019.tsv, I combined the following signatures together:
#6+15+20+21+26=MMR  #2+13=APOBEC  #4+24+29=SMOKING
###NOTE: in the file /juno/work/taylorlab/friedman/myAdjustedDataFiles/impactSignatureCalls_Nov20_2019.tsv, I combined the following signatures together:

def validate_hypermutation_decomposition(hypermutationFileDir = pathPrefix + '/juno/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/',
                                        signaturesFilePath = pathPrefix + '/juno/work/taylorlab/friedman/myAdjustedDataFiles/impactSignatureCalls_Nov20_2019.tsv',
                                        minNormalCasesToCompare = 25 #dont change anything with the normal cases if there are fewer than this many signatures
                                        ):
    
    signaturesData = pd.read_table(signaturesFilePath)
    signaturesData = signaturesData[signaturesData['dominantSignature'] != 'insufficientMutBurden']
    dominantSignatureDict = dict(zip(signaturesData['Tumor_Sample_Barcode'], signaturesData['dominantSignature']))
    for f in os.listdir(hypermutationFileDir):
        print f
        fullFilePath = os.path.join(hypermutationFileDir, f)
        cancerTypeDataFrame = pd.read_table(fullFilePath)
        cancerTypeDataFrame['dominantSignature'] = cancerTypeDataFrame['Tumor_Sample_Barcode'].apply(lambda x: dominantSignatureDict[x] if x in dominantSignatureDict else None)
        
        allCases = set(cancerTypeDataFrame['Tumor_Sample_Barcode'])
        hypermutatedCases = set(cancerTypeDataFrame[cancerTypeDataFrame['hypermutantClassification'] == 'Hypermutated']['Tumor_Sample_Barcode'])
        normalCases = set(cancerTypeDataFrame[cancerTypeDataFrame['hypermutantClassification'] == 'Normal']['Tumor_Sample_Barcode'])
        
        hypermutatedDominantSignatures = signaturesData[signaturesData['Tumor_Sample_Barcode'].isin(hypermutatedCases)]['dominantSignature']
        normalDominantSignatures = signaturesData[signaturesData['Tumor_Sample_Barcode'].isin(normalCases)]['dominantSignature']
        
        if len(hypermutatedDominantSignatures) > 0 and len(normalDominantSignatures) > minNormalCasesToCompare:
            mostCommonHyperSig = Counter(hypermutatedDominantSignatures).most_common()[0][0]
            mostCommonNormalSig = Counter(normalDominantSignatures).most_common()[0][0]
            if mostCommonNormalSig == 'mean_MMR':
                mostCommonNormalSig = 'mean_1' #fix mmr dominant signature calls in non hypermutated which is probably a mistake
            if mostCommonNormalSig == mostCommonHyperSig:
                cancerTypeDataFrame['hypermutantClassification'] = cancerTypeDataFrame.apply(lambda row:
                            'highMutationBurden' if (row['dominantSignature'] == mostCommonNormalSig) and (row['hypermutantClassification'] == 'Hypermutated')
                            else row['hypermutantClassification'], axis=1)
        
        #RESAVE THE NEW HYPERMUTATION INFO
        print 'saving new version'
        cancerTypeDataFrame.to_csv(fullFilePath, index=False, sep='\t')
