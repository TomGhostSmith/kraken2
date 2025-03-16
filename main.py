import os
import sys
import json
import multiprocessing
from tqdm import tqdm
sys.path.append('..')

from config import config


def getDomain():
    taxID2seqID = dict()
    seqID2domain = dict()
    taxID2domain = dict()
    with open("/Data/ICTVData/model/kraken2/seqid2taxid.map") as fp:
        for line in fp:
            seqID, taxID = line.strip().split('|')[-1].split('\t')
            taxID2seqID[taxID] = seqID
    
    with open("/Data/ICTVData/model/kraken2/library_report.tsv") as fp:
        for line in fp:
            domain, seqID = line.strip().split('\t')[:2]
            seqID = seqID.split(' ')[0][1:]
            seqID2domain[seqID] = domain
    
    for taxID, seqID in taxID2seqID.items():
        if seqID in seqID2domain:
            taxID2domain[taxID] = seqID2domain[seqID]

    with open('/Data/ICTVData/model/kraken2/taxID2domain.json', 'wt') as fp:
        json.dump(taxID2domain, fp, indent=2)
    with open('/Data/ICTVData/model/kraken2/taxID2seqID.json', 'wt') as fp:
        json.dump(taxID2seqID, fp, indent=2)
    with open('/Data/ICTVData/model/kraken2/seqID2domain.json', 'wt') as fp:
        json.dump(seqID2domain, fp, indent=2)
    
def getViralSequences():
    with open("/Data/ICTVData/model/kraken2/library_viral.tsv") as fp:
        lines = fp.readlines()
    
    params = list()
    for line in lines:
        terms = line.strip().split('\t')
        URL = terms[2]
        URL = "https" + URL[3:]
        seqID = terms[1].split(' ')[0][1:]
        params.append((seqID, URL))
    with tqdm(total=len(params)) as progress:
        with multiprocessing.Pool(processes=64) as pool:
            for param in params:
                pool.apply_async(getVirusSeq, param, callback=lambda _:updateProgress(progress))
            
            pool.close()
            pool.join()

def updateProgress(progress):
    progress.update(1)

def getVirusSeq(seqID, URL):
    os.system(f'wget -qO- "{URL}" | gunzip > /Data/ICTVData/dataset/ncbi/fasta/{seqID}.fasta')

# getDomain()
getViralSequences()
# getVirusSeq("test", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/574/705/GCF_018574705.1_ASM1857470v1/GCF_018574705.1_ASM1857470v1_genomic.fna.gz")
