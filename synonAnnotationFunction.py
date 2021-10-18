#Function version of synonymousAnnotation.py
#Sidhant Puntambekar

import pandas as pd
import numpy as np

# Function takes in a curated synonymous mutation file,
# a genome in a bottle file describing hard to sequence sites.
# Make sure the column names in the GIAB and synon mutation file are the same as in this function.
def synonAnnotate(synonMutationFile, genomeBottleFile): 
        chromosomeArray = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                          'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                           'chr20', 'chr21', 'chr22', 'chrX', 'chrY']	

        dfSynonMutation = pd.read_csv(synonMutationFile, delimiter="\t", index_col=False)
        if dfSynonMutation.shape[0] == 0:
                return None
        print(dfSynonMutation.iloc[:,0])
        print(dfSynonMutation.iloc[:,1])
        print(dfSynonMutation.iloc[:,2])
        dfGenomeBottle = pd.read_csv(genomeBottleFile, delimiter="\t")
        
        chromDfUnique = dfGenomeBottle.chrom.unique()
        print("Unique Chromosomes GIAB: ", chromDfUnique)

        synonDfUnique = dfSynonMutation.iloc[:,2].unique()
        print("Unique Chromsomes Synon: ", synonDfUnique)
        
        dfChromosomes = dict()
        
        # Split up the GIAB dataframe by chromosome
        for i in range(len(chromDfUnique)):
                dfChromosomes[i] = dfGenomeBottle[dfGenomeBottle['chrom'] == chromDfUnique[i]]

        # print(dfChromosomes)

        # Get the first and last range starts from the starting chromosome GIAB list
        firstChrom = synonDfUnique[0]
        chromosomeInterestedIn = chromosomeArray[firstChrom-1]
        lb = dfChromosomes[firstChrom-1]['left_sites'].iloc[0]
        ub = dfChromosomes[firstChrom-1]['left_sites'].iloc[-1]
        
        #chromosomeInterestedIn = chromosome
        
        #print(chromosomeInterestedIn)
        #print(lb) 
        #print(ub)
        
        dfGenomeBottleChromBoundFiltered = dfGenomeBottle[dfGenomeBottle['chrom'] == chromosomeInterestedIn]

        lenSynonMutation = len(dfSynonMutation)
        dfSynonMutation.insert(dfSynonMutation.shape[1], 'site_in_genome_bottle', [False for x in range(lenSynonMutation)])
        
        chromCol = dfSynonMutation.iloc[:, 2]
        chromArray = chromCol.values
        
        siteCol = dfSynonMutation.iloc[:, 3]
        siteArray = siteCol.values
        
        leftSites = dfGenomeBottleChromBoundFiltered.loc[:, 'left_sites']
        leftSitesArray = leftSites.values
        
        rightSites = dfGenomeBottleChromBoundFiltered.loc[:, 'right_sites']
        rightSitesArray = rightSites.values
        
        inSiteBool = np.array([False] * lenSynonMutation)
        indexInSite = np.array([0] * lenSynonMutation)
        
        j = 0
        
        for i in range(lenSynonMutation):
                # Update the chromosome if necessary
                if (chromArray[i] != chromArray[i-1]):
                        newChrom = chromArray[i]
                        chromosomeInterestedIn = chromosomeArray[newChrom-1]
                        lb = dfChromosomes[newChrom-1]['left_sites'].iloc[0]
                        ub = dfChromosomes[newChrom-1]['left_sites'].iloc[-1]
                        
                        print("Chromosome: ", chromosomeInterestedIn)
                        print("Lower bound", lb) 
                        print("Upper bound", ub)
	                
                        dfGenomeBottleChromBoundFiltered = dfGenomeBottle[dfGenomeBottle['chrom'] == chromosomeInterestedIn]
			
                        leftSites = dfGenomeBottleChromBoundFiltered.loc[:, 'left_sites']
                        leftSitesArray = leftSites.values

                        rightSites = dfGenomeBottleChromBoundFiltered.loc[:, 'right_sites']
                        rightSitesArray = rightSites.values
                        
                        j = 0	
                        print(leftSitesArray)
                        print(rightSitesArray)

                while leftSitesArray[j] <= siteArray[i]:
                        if j == (len(leftSitesArray)-1):
                                break
                        j = j + 1

                
                if ((siteArray[i] >= leftSitesArray[j-1]) and (siteArray[i] < rightSitesArray[j-1])):
                        inSiteBool[i] = True
                        indexInSite[i] = i
                else:
                        inSiteBool[i] = False

        dfSynonMutation.loc[indexInSite, 'site_in_genome_bottle'] = True
        print(dfSynonMutation.head(100))
        print(dfSynonMutation.tail(100))
        
        return dfSynonMutation

def main():
        synonMutation1 = synonAnnotate('../data/synonMutations/fake_transcript_variants_sorted_v3_ad_syn_processed',
                                       '../data/GRCh38/GRCh38_alldifficultregions.bed')
        synonMutation2 = synonAnnotate('../data/synonMutations/fake_transcript_variants_sorted_v3_ad_syn_processed_chromChange',
                                       '../data/GRCh38/GRCh38_alldifficultregions.bed') 

if __name__ == "__main__":
        main()
