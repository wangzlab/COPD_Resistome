import csv
from ete3 import NCBITaxa

ncbi = NCBITaxa()

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}


FileID_byDate = "2021.4" 
taxidFile = "2.out_all_taxids\\\%s.txt" % (FileID_byDate)  #check the correct file name
inputObj  = open(taxidFile, "r")

taxids = []
for line in inputObj.readlines():
    taxids.append(int(line))
inputObj.close()

if __name__ == '__main__':
    #taxids = [1281578,  1469502, 920]
    desired_ranks = ['domain','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    results = list()
    absent_ids = []
    for taxid in taxids:
        try:
            get_desired_ranks(taxid, desired_ranks)
            results.append(list())
            results[-1].append(str(taxid))
            ranks = get_desired_ranks(taxid, desired_ranks)
            for key, rank in ranks.items():
                if rank != '<not present>':
                    results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
                else:
                    results[-1].append(rank)
        except:
            absent_ids.append(taxid)
            pass
            
    
    outFile = "3.out_taxid_mapping\\\%s_out_taxid_mapping.txt" % (FileID_byDate)
    outputObj = open(outFile, "w")
    
    #generate the header
    header = ['Original_query_taxid']
    header.extend(desired_ranks)
    #print('\t'.join(header))
    outputObj.write("%s\n" % ('\t'.join(header)))
    
    
    #print the results
    for result in results:
        #print('\t'.join(result))
        outputObj.write("%s\n" % ('\t'.join(result)))
    
    outputObj.close()