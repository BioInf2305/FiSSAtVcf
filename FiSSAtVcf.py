##load the config file 

configfile: "config.json"

##path of the programs to be used in the pipeline

plink = config["tools"]["PLINK"]
king = config["tools"]["KING"]
fastindep = config["tools"]["fastIndep"]

##path of the input and output directory

outDir = config["base_out"]["outputDir"]
inDir = config["base_out"]["inputDir"]

##thresholds for the various tools

ch_list = config["params"]["chrmList"]
kinThreshold = float(config["params"]["kinThreshold"])
runs = int(config["params"]["runs"])
seed = int(config["params"]["seeds"])
threads = int(config["params"]["threads"])
extensionList = ["bim", "fam", "bed", "nosex", "log"]

##output of the final steps

rule all:
    input:
        expand(outDir+"chr{chrm}_unrelatedSamplesWithTags.vcf.gz",chrm=ch_list)

##combined vcf files using bcftools

rule vcfConcat:
    input:
        expand(inDir+"{chrm}.vcf.gz",chrm=ch_list)
    output:
        expand(outDir+"Chrms_Merged.vcf.gz")
    shell:
        "bcftools concat -o {output} -O z {input}"

##run plink to generate binary file

rule runPlink:
    input:
        expand(outDir+"Chrms_Merged.vcf.gz")
    output:
        expand(outDir+"Chrms_Merged.{extension}",extension=extensionList)
    params:
        out1 = expand(outDir+"Chrms_Merged")
    shell:
        plink +" --vcf {input} --make-bed --out {params.out1}"

##run KING to estimate kinship coefficients

rule runKing:
    input:
        expand(outDir+"Chrms_Merged.bed")
    output:
        expand(outDir+"Chrms_Merged.kin0")
    params:
        out1 = expand(outDir+"Chrms_Merged")
    shell:
        king + " -b {input} --kinship --degree 0 --prefix {params.out1}"

##python codes to prepare input for Fastindep

def extractTotalNumber(inFile):

    totalSamples = []
    header = 0
    with open(inFile) as source:
        for line in source:
            if header == 0:
                header += 1
            else:
                a = line.split()
                if a[1] not in totalSamples:
                    totalSamples.append(a[0])
                if a[3] not in totalSamples:
                    totalSamples.append(a[3])
    return len(totalSamples)


def prepareInput(kingMatrixIn,nSample,fastIndepIn):

    nSample = int(nSample)
    totalList = []
    diagonal = [10000]
    sampleList = []
    orderList = []
    dest = open(fastIndepIn,"w")
    header = 0
    for i in range(nSample):
        tmpList = (nSample)*diagonal
        totalList.append(tmpList)
    with open(kingMatrixIn) as source:
        for line in source:
            if header == 0:
                header += 1
            else:
                line = line.split()
                if line[1] not in sampleList:
                    start = len(sampleList[:])
                    populate = start+1
                    sampleList.append(line[1])
                    if start == len(totalList)-2:
                        sampleList.append(line[3])
                    if len(orderList) == 0:
                        orderList.append(line[1])
                totalList[start][populate] = line[7]
                totalList[populate][start] = line[7]
                populate = populate+1
                if start == 0:
                    orderList.append(line[3])
                else:
                    if orderList[populate-1] != line[3]:
                        print("Order does not match",populate,orderList[populate-1],line[3])
                        sys.exit(0)
    dest.write("\t"+"\t".join(sampleList)+"\n")
    for i,v in enumerate(totalList):
        dest.write(sampleList[i]+"\t")
        pHetList = list(map(str,v))
        dest.write("\t".join(pHetList)+"\n")
    dest.close()

## prepare input for Fastindep

rule convertInput:
    input:
        expand(outDir+"Chrms_Merged.kin0")
    output:
        expand(outDir+"Chrms_Merged.kin0.sorted.matrix"),
    params:
        out1 = expand(outDir+"Chrms_Merged.kin0.sorted"),
        out2 = expand(outDir+"Chrms_Merged.kin0.sorted.matrix"),
    run:
        shell("sort -k1,1 {input} > {params.out1}")
        nSamples = extractTotalNumber(input[0])
        prepareInput(params.out1[0],nSamples,params.out2[0])

##run Fastindep to extract unrelated individuals

rule runFastindep:
    input:
        file1 = expand(outDir+"Chrms_Merged.kin0.sorted.matrix"),
    output:
        expand(outDir+"outfile.txt")
    run:
        shell(fastindep+" -t {kinThreshold} -n {runs} -s {seed} -i {input} -o {output}",kinThreshold=kinThreshold,runs=runs,seed=seed)

def extractUnrelatedIndi(inFile,outFile):

    dest = open(outFile,"w")
    count = 0
    with open(inFile) as source:
        for line in source:
            if "Set" in line:
                line = line.rstrip().split()
                if line[0] == "Set" and line[1] == "Size" and line[2] == "=" and count == 0:
                    count += 1
                    indiv = line[4:]
                    dest.write("\n".join(indiv))
    dest.close()

## prepare vcf file containing only unrelated individuals extracted in the previous step, note that 
## bcftools, by default, will add "AC" and "AN" tag at this stage. 

rule extractUnrelatedSamples:
    input:
        expand(outDir+"outfile.txt")
    output:
        expand(outDir+"Chrms_Merged_UnrelatedSamples.vcf.gz")
    params:
        inp1 = expand(outDir+"UnrelatedSamples.txt"),
        inp2 = expand(outDir+"Chrms_Merged.vcf.gz")
    run:
        extractUnrelatedIndi(input[0],params.inp1[0])
        shell("bcftools view -S {params.inp1} -o {output} -O z {params.inp2}")


##python codes to remove monomorphic SNPs, add tags in parallel, and delete temporary files

import sys
import decimal
from pysam import VariantFile
from multiprocessing import Pool
from collections import OrderedDict
import os,glob
import os.path

class AddVcfTag:


    def __init__(self,vcfIn,chrmList,nCores,outDir):

        self.vcfIn = vcfIn
        self.chrmList = chrmList.split(",")
        self.chrmSizeDict = OrderedDict()
        self.nCores = int(nCores)
        self.chrmOutList = OrderedDict()
        self.globalDpPDict = OrderedDict()
        self.outDir = outDir

    ##method to extract chromosome size from the vcf file, helps in running "vcfAddTags" function in parallel

    def populatateChrmSize(self):

        tmpVcfIn = VariantFile(self.vcfIn)
        for chrom in tmpVcfIn.header.records:
            if chrom.key == "contig":
                if chrom.values()[0] in self.chrmList:
                    self.chrmSizeDict[chrom.values()[0]] = int(chrom.values()[1])

    ##method to add tags

    def vcfAddTags(self,inputList):
        localDpPList=[]
        vcfIn=VariantFile(self.vcfIn)
        vcfIn.header.info.add("AF","1","Float","alternate allele frequency")
        vcfIn.header.info.add("MAF","1","Float","minor allele frequency")
        vcfIn.header.info.add("AVG_DP","1","Float","average depth across samples")
        vcfIn.header.info.add("N_ALT_HET","1","Integer","number of heterozygous alternate carriers")
        vcfIn.header.info.add("N_ALT_HOM","1","Float","number of homozygous alternate carriers")
        vcfIn.header.info.add("AVG_DP_P","1","Float","average depth's percentile values")
        vcfOutP = os.path.join(self.outDir,inputList[0])
        vcfOut = VariantFile(vcfOutP,'w',header=vcfIn.header)
        sampleList = list(vcfIn.header.samples)
        genotypeDict = {(None,None):0, (0, 1):1, (1, 0):1, (1, 1):2, (0, 0):0}
        hetCarriers = {(0, 1):1, (1, 0):1, (0, 0):0, (1, 1):0}
        homCarriers = {(0, 1):0, (1, 0):0, (0, 0):0, (1, 1):1}
        for rec in vcfIn.fetch(inputList[1],inputList[2],inputList[3]):
            ac = rec.info["AC"][0]
            an = rec.info["AN"]
            avgDp = 0
            nAltHet = 0
            nAltHom = 0
            sampleCount = 0
            for sample in sampleList:
                if rec.samples[sample]["GT"] != (None,None):
                    sampleCount += 1
                    avgDp += rec.samples[sample]["DP"]
                    nAltHet += hetCarriers[rec.samples[sample]["GT"]]
                    nAltHom += homCarriers[rec.samples[sample]["GT"]]
            if ac != sampleCount and ac > 0:
                af = ac/an
                rec.info["AF"] = af
                rec.info["MAF"] = min(1-af,af)
                avgDp = str(avgDp/sampleCount)
                avgDp = float(decimal.Decimal(avgDp).quantize(decimal.Decimal('1e-3')))
                rec.info["AVG_DP"] = avgDp
                rec.info["N_ALT_HET"] = nAltHet
                rec.info["N_ALT_HOM"] = nAltHom
                localDpPList.append(avgDp)
                vcfOut.write(rec)
        return localDpPList

    ##method to run "vcfAddTags" code in parallel

    def runParallel(self):

        self.populatateChrmSize()
        sortedChrmSizeDict = dict(sorted(self.chrmSizeDict.items(), key=lambda item: item[1]))
        listRecords = []
        for chrom in sortedChrmSizeDict:
            self.chrmOutList[chrom] = []
        if self.nCores >= len(sortedChrmSizeDict):
            divideChrm = self.nCores//len(list(sortedChrmSizeDict.keys()))
            remainChrm = self.nCores%len(list(sortedChrmSizeDict.keys()))
            for chrom in sortedChrmSizeDict:
                if chrom == list(sortedChrmSizeDict.keys())[-1]:
                    chromDivision = divideChrm+remainChrm
                else:chromDivision = divideChrm
                equalDiv = sortedChrmSizeDict[chrom]//chromDivision
                start = 0
                vcfOut = 0
                for i in range(chromDivision):
                    vcfOut += 1
                    vcfOutName = chrom+"_"+str(vcfOut)+".tmp.vcf.gz"
                    self.chrmOutList[chrom].append(vcfOutName)
                    tmpList = []
                    tmpList.append(vcfOutName)
                    tmpList.append(chrom)
                    if i == chromDivision-1:
                        tmpList.append(start)
                        tmpList.append(sortedChrmSizeDict[chrom])
                    else:
                        tmpList.append(start)
                        tmpList.append(start+equalDiv)
                    listRecords.append(tmpList[:])
                    start = start+equalDiv
            pool = Pool(self.nCores)
            localDpPListAppended = pool.map(self.vcfAddTags,listRecords)
        else:
            chromList = list(sortedChrmSizeDict.keys())
            pool = Pool(self.nCores)
            chromCoresList = [chromList[i:i+self.nCores] for i in range(0,len(chromList),self.nCores)]
            listRecords = []
            for i,v in enumerate(chromCoresList):
                for i in range(len(v)):
                    tmpList = []
                    vcfOutName = v[i]+".tmp.vcf.gz"
                    self.chrmOutList[v[i]].append(vcfOutName)
                    tmpList.append(vcfOutName)
                    tmpList.append(v[i])
                    tmpList.append(0)
                    tmpList.append(sortedChrmSizeDict[v[i]])
                    listRecords.append(tmpList[:])
            localDpPListAppended = pool.map(self.vcfAddTags,listRecords)
        localDpPListFlat = [DP_P for DP_PList in localDpPListAppended for DP_P in DP_PList]
        localDpPListFlat = sorted(localDpPListFlat)
        for i,v in enumerate(localDpPListFlat):
            if str(v) not in self.globalDpPDict:
                self.globalDpPDict[str(v)] = (100*i)/len(localDpPListFlat)
        self.combinedVcf()
        self.removeTmpFiles()

    def combinedVcf(self):

        for chrm in self.chrmOutList:
            vcfOutList = self.chrmOutList[chrm]
            firstVcfFile = os.path.join(self.outDir,vcfOutList[0])
            firstVcfFile = VariantFile(firstVcfFile)
            outFileP = os.path.join(self.outDir,str(chrm)+"_unrelatedSamplesWithTags.vcf.gz")
            outFile = VariantFile(outFileP,"w",header=firstVcfFile.header)
            for rec in firstVcfFile.fetch():
                tmpDpValue = float(decimal.Decimal(rec.info["AVG_DP"]).quantize(decimal.Decimal('1e-3')))
                tmpDpValue = str(tmpDpValue)
                rec.info["AVG_DP_P"] = float(self.globalDpPDict[tmpDpValue])
                outFile.write(rec)
            for vcfFile in vcfOutList[1:]:
                otherVcfFileP = os.path.join(self.outDir,vcfFile)
                otherVcfFile = VariantFile(otherVcfFileP)
                for rec in otherVcfFile:
                    tmpDpValue = float(decimal.Decimal(rec.info["AVG_DP"]).quantize(decimal.Decimal('1e-3')))
                    tmpDpValue = str(tmpDpValue)
                    rec.info["AVG_DP_P"] = float(self.globalDpPDict[tmpDpValue])
                    outFile.write(rec)

    ##method to remove temporary files

    def removeTmpFiles(self):
        os.chdir(self.outDir)
        for fil in glob.glob("Chrms_Merged*"):
            os.remove(fil)
        for fil in glob.glob("*tmp*"):
            os.remove(fil)
        os.remove("outfile.txt")

##add vcf tags using the above python codes

rule AddTagVcf:
    input:
        expand(outDir+"Chrms_Merged_UnrelatedSamples.vcf.gz")
    output:
        expand(outDir+"chr{chrm}_unrelatedSamplesWithTags.vcf.gz",chrm=ch_list)
    params:
        inp1=",".join(list(map(lambda x:"chr"+x,list(map(str,ch_list))))),
        inp2=threads
    run:
        shell("bcftools index {input}")
        runAddTags = AddVcfTag(input[0],params.inp1,params.inp2,outDir)
        runAddTags.runParallel()
