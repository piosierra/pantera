from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import subprocess
import pandas as pd
import os
import re


cores = 4

input_file = sys.argv[1]
    
def write_sequences_file(sequences, filename):

    try:
        SeqIO.write(sequences, filename, "fasta")
    except FileNotFoundError:
        print("FATAL ERROR: I couldn't find the file, please check: '" + filename + "'. Path not found")
        sys.exit(0)
    except PermissionError:
        print("FATAL ERROR: I couldn't access the files, please check: '" + filename + "'. I don't have permissions.")
        sys.exit(0)
    except Exception as exp :
        print("FATAL ERROR: There is a unknown problem writing sequences in : '" + filename + "'.")
        print(exp)
        sys.exit(0)



def find_TRs2(te, minLTR, minTIR, minpolyA):
    lenLTR = 0
    lenTIR = 0
    lenPolyA = 0
    te = te.upper()
    mismatches = minpolyA // 8
    if te.seq[:minpolyA + mismatches].count("A") >= minpolyA:
        lenPolyA = max(lenPolyA, te.seq[:minpolyA + mismatches].count("A") + len(re.match("(.*?)[^A]",str(te.seq)[minpolyA:]).group())-2)
    if te.seq[:minpolyA + mismatches].count("T") >= minpolyA:
        lenPolyA = max(lenPolyA, te.seq[:minpolyA + mismatches].count("T") + len(re.match("(.*?)[^T]",str(te.seq)[minpolyA:]).group())-2)
    if te.seq[-(minpolyA + mismatches):].count("A") >= minpolyA:
        lenPolyA = max(lenPolyA, te.seq[-(minpolyA + mismatches):].count("A") + len(re.match("(.*?)[^A]",str(te.seq)[:len(te.seq)-minpolyA][::-1]).group())-2)
    if te.seq[-(minpolyA + mismatches):].count("T") >= minpolyA:       
        lenPolyA = max(lenPolyA, te.seq[-(minpolyA + mismatches):].count("T") + len(re.match("(.*?)[^T]",str(te.seq)[:len(te.seq)-minpolyA][::-1]).group())-2)
    try: 
 #       os.mkdir("temp" + str(te.id[:14]))
        os.mkdir("temp") 
    except FileExistsError:
        pass   
 #   os.chdir("temp" + str(te.id[:14]))
    os.chdir("temp")
    seq_name = str(te.id).split("#")[0]
    write_sequences_file(te, seq_name + ".fa")
    output = subprocess.run(['makeblastdb', '-in',  seq_name + ".fa", '-dbtype', 'nucl'], stdout=subprocess.PIPE, text=True)
    output = subprocess.run(["blastn -query " +  seq_name + ".fa -db " + seq_name + ".fa -num_threads " + str(cores) + " -outfmt 6 -word_size 11 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 | cut -f 1,7-10 | sed 's/#/-/g' > " + str(seq_name) + ".blast"], shell=True)
    blastresult = pd.read_table( str(seq_name) + ".blast", sep='\t', names=['qseqid', 'qstart', 'qend', 'sstart', 'send'],  dtype = {'qseqid':str, 'qstart':int, 'qend':int, 'sstart':int, 'send':int} )
    # print(blastresult)
    blastHits = blastresult.shape[0]
    if blastresult.shape[0] > 1:
        blastresult = blastresult[1:]
        blastresult["len"] = abs(blastresult["qend"]- blastresult["qstart"])
    #    blastresult.drop(blastresult[blastresult.len > len(te.seq)/2].index, inplace=True)
        ltr_check = blastresult[(((blastresult.qend-blastresult.qstart)*(blastresult.send-blastresult.sstart))>0) & (blastresult[['qstart','qend','sstart','send']].min(axis=1) < len(te)*0.1) & (blastresult[['qstart','qend','sstart','send']].max(axis=1) > len(te)*0.9) ]
        try: 
            lenLTR = max(ltr_check[ltr_check.len>=minLTR].len) 
        except: 
            lenLTR = 0
        tir_check = blastresult[(((blastresult.qend-blastresult.qstart)*(blastresult.send-blastresult.sstart))<0) & (blastresult[['qstart','qend','sstart','send']].min(axis=1) < len(te)*0.1) & (blastresult[['qstart','qend','sstart','send']].max(axis=1) > len(te)*0.9) ]
        try: 
            lenTIR = max(tir_check[tir_check.len>=minTIR].len) 
        except: 
            lenTIR = 0
    os.chdir("..")
    os.system("rm -r temp")
#    os.system("rm -r temp" + str(te.id[:14]))

    return lenLTR, lenTIR, lenPolyA, blastHits

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
with open(input_file + ".benchmark2", 'a') as out_file:
    for fasta in fasta_sequences:
        name = fasta.id
        out_file.write(name + " " + str(len(fasta.seq)) + " ")
        LTR, TIR, PA, BH =find_TRs2(fasta,10,6,8)
        out_file.write(str(LTR)+" "+ str(TIR)+ " "+str(PA)+" "+str(BH)+"\n")








