import sys, os, shutil

from Bio import SeqIO
from Bio import Phylo

INFILE = os.environ.get("INFILE",None).split()

rule all:
    input: expand("{infile}.aln.skyride.pdf",infile=INFILE)

rule aln:
    input: "{infile}.fas"
    output: "{infile}.aln"
    #log: "{infile}.mafft.log"
    message: "Aligning sequences using MAFFT"
    shell: "mafft --auto {input} > {output}"  

rule iqtree:
    input: "{infile}.aln"
    output: "{infile}.aln.iqtree","{infile}.aln.treefile","{infile}.aln.mldist","{infile}.aln.log" 
    message: "Building phylogeny using IQ-TREE"
    shell: "iqtree-omp -s {input} -m TEST -nt 2"

rule rtt:
    input: "{infile}.aln.treefile"
    output: "{infile}.aln.treefile.rtt.tre", "{infile}.aln.treefile.rtt.txt"
    message: "Generating Root-to-Tip regression tree. Calculating substitution rate and tMRCA"
    shell: "Rscript ~/bin/rtt2.R {input} {input}.rtt"  

rule treble:
    input: "{infile}.aln.treefile.rtt.tre", "{infile}.aln.treefile"   
    output: "{infile}.aln.treefile.treble.txt"
    message: "TREBLE estimates substitution rate and tMRCA from rooted tree"
    shell: "Rscript ~/bin/treble2.R {input[0]} {input[1]}.treble"

rule getDates:
    input: "{infile}.aln.treefile.treble.txt", "{infile}.aln.treefile"
    output: temp("{infile}.aln.treefile.temp"), temp("{infile}.dates"), temp("{infile}.map")
    message: "Preparing tree file for LSD run"
    run:
        # read the tree
        tr = Phylo.read(input[1],'newick')
        
        # get the name of the tips
        clades = [term.name for term in tr.get_terminals()]
        
        # get number of taxa in the tree
        numTaxa = len(clades)
        
        # read the tree as a string for renaming taxa names
        treeStr = [line.strip() for line in open(input[1])][0]
        
        ## get the tip.labels from the tree and assign new names
        i  = 1
        
        # create the str object to write dates for LSD
        dStr = str(numTaxa)
        
        # create the mapping string
        mapStr = str(numTaxa)
        
        # get dates for all the taxa
        for clade in clades:
            cld = str(clade)
            ind = cld.rfind('_')
            cName = cld[:ind]
            cDate = cld[ind+1:]
            nName = 'X' + str(i)
            i += 1
            dStr += '\n' + nName + '\t' + cDate
            
            # replace old names with simplified names
            treeStr = treeStr.replace(cld,nName)
            
            # add mapping to the string
            mapStr += '\n' + nName + '\t' + cld
            
        # write the new tree
        with open(output[0],'w') as fh:
            fh.write(treeStr)
        
        # write the dates file
        with open(output[1],'w') as fh:
            fh.write(dStr)
        
        # write the mapping file
        with open(output[2],'w') as fh:
            fh.write(mapStr)

rule runLSD:
    input: name="{infile}.aln.treefile.temp", dates="{infile}.dates", map="{infile}.map" 
    output: temp("{infile}.aln.treefile.temp.result"), temp("{infile}.aln.treefile.temp.result.nexus"), temp("{infile}.aln.treefile.temp.result.newick"), temp("{infile}.aln.treefile.temp.result.date.newick") 
    message: "LSD roots the tree, generates subsRate and tMRCA"
    shell: "lsd -i {input.name} -d {input.dates} -c -r a"

 
rule processLSDoutput:
    input: "{infile}.map", "{infile}.aln.treefile.temp.result","{infile}.aln.treefile.temp.result.nexus","{infile}.aln.treefile.temp.result.newick", "{infile}.aln.treefile.temp.result.date.newick"    
    output: "{infile}.aln.treefile.result","{infile}.aln.treefile.result.nexus","{infile}.aln.treefile.result.newick", "{infile}.aln.treefile.result.date.newick"
    message: "Creating LSD output files"
    run:
        # read in the mapping file
        mapLines = [line.strip() for line in open(input[0])]
        
        # get number of taxa
        nTaxa = int(mapLines[0])
        
        # create the result file
        shutil.copy(input[1],output[0])
        
        # read the LSD output files as strings
        newickStr = [line.strip() for line in open(input[3])][0]
        newickDateStr = [line.strip() for line in open(input[4])][0]
        
        # replace simplified taxa names with actual names
        for i in range(nTaxa,0,-1):
            token = mapLines[i].split('\t')
            #print(token[0],'\t',token[1])
            newickStr = newickStr.replace(token[0],token[1])
            newickDateStr = newickDateStr.replace(token[0],token[1])
        
        # write the new tree files with original taxa names
        with open(output[2],'w') as fh:
            fh.write(newickStr)
        
        with open(output[3],'w') as fh:
            fh.write(newickDateStr)
        
        Phylo.convert(output[2],'newick',output[1],'nexus') 
        
        
rule rateLSD:
  input: "{infile}.aln.treefile.result"
  output: "{infile}.lsd.rates.txt"
  message: "Extracts the substitution rate and tMRCA estimated by LSD"
  run:
    lines = [line.strip('\n') for line in open(input[0],'r')]
    for line in lines:
      if 'tMRCA' in line:
        words = line.split(' ')
        rate = words[3]
        tmrca = words[5]
        
        fh = open(output[0],'w')
        fh.write('TMRCA\tSubstRate\n%s\t%s\n' %(tmrca,rate))
        fh.close()
        
rule treedater:
    input: "{infile}.lsd.rates.txt", "{infile}.aln.treefile", "{infile}.aln"
    output: "{infile}.treedater.tre","{infile}.treedater.png", "{infile}.treedater.result"
    message: "Running treeDater to get time-stamped tree, tMRCA and subtitution rates"
    shell: "Rscript ~/bin/treedater.R {input[1]} {input[2]} {output[0]} {output[1]} {output[2]}"        
        
rule genieR:
    input: "{infile}.treedater.tre"
    output: "{infile}.genieR.txt"    
    message: "Fitting demographic models CONST, EXPO and LOG on the time-stamped tree using GenieR package"
    shell: "Rscript ~/bin/genieR.R {input[0]} {output[0]}"    
    
    
rule skyRidesPhyloDyn:
  input: "{infile}.aln.treefile.result.date.newick", "{infile}.fas", "{infile}.genieR.txt"
  output: "{infile}.aln.skyride.pdf"
  message: "Construct skyride using Phylodyn"
  shell: "Rscript ~/bin/skyridePhylodyn.R {input[0]} {output}"    
