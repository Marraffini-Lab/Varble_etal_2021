from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import csv
import pickle

handle1 = open("testR1.fastq", "rU")
records1 = SeqIO.parse(handle1, 'fastq')


outReport = open("parser_report.txt", "wb")
out1 = csv.writer(open("output1.csv", "wb"))


totalrecords = 0
phagereads = 0
repeatreads = 0
spacer1 = 0
spacer2 = 0
spacer3 = 0
spacer4 = 0
spacer5 = 0
spacer6 = 0

spacer = ''
ribobc = ''
crispbc = ''


phage_start = "TTTAGTTTAC"
DR_5prime = "TATGCTGTTTTG"
riboB = 'CCAGCTGTACAGTCT'
crispB= 'CATTGGGATATATCAA'

spacers_list = ["CTTCTTGCGCTTTTT", "TCAATTTGTAAAAAA",
            "AATTAATTGCGCTCT", "TTAGGTGCGCTTGGC", "GGTAAACCGTGCTTT",
            "CTTGTTGAGTTCCAT"]

spacer_listcount = [spacer1, spacer2, spacer3, spacer4, spacer5, spacer6]

for r1 in records1:
    totalrecords += 1
    seq1 = str(r1.seq)
    seq1rc = str(r1.seq.reverse_complement())

    target = ""

    pos1 = nt_search(seq1, phage_start)
    if len(pos1) > 1:
        target2 = seq1
    else:
        pos2 = nt_search(seq1rc, phage_start)
        if len(pos2) > 1:
            target2 = seq1rc

    if (target2 != ""):
        #print("Found barcode: " + str(barcodes[x]))
        phagereads += 1
        pos1e = nt_search(target2, DR_5prime)

        if len(pos1e) > 1:
            repeatreads += 1
            for x in range(6):
                pos1b = nt_search(target2, spacers_list[x])

                if len(pos1b) > 1:
                    spacer_listcount[x] += 1
                    spacer_pos = pos1b[1]


            if len(pos1e) > 2:
                repeatreads += 1
                print target2




    if ((totalrecords % 25) == 0):
        print("Total records: " + str(totalrecords))
        print("My phage reads: " + str(phagereads))
        print("My repeat reads: " + str(repeatreads))
        print ('Spacer1: ' + str(spacer1))
        print ('Spacer2: ' + str(spacer2))
        print ('Spacer3: ' + str(spacer3))
        print ('Spacer4: ' + str(spacer4))
        print ('Spacer5: ' + str(spacer5))
        print ('Spacer6: ' + str(spacer6))



outReport.write("Total records: " + str(totalrecords)+'\n')
outReport.write("My phage reads: " + str(phagereads)+'\n')
outReport.write("My repeat reads: " + str(repeatreads) + '\n')
outReport.write('Spacer 1: ' + str(spacer1) + '\n')
outReport.write('Spacer 2: ' + str(spacer2) + '\n')
outReport.write('Spacer 3: ' + str(spacer3) + '\n')
outReport.write('Spacer 4: ' + str(spacer4) + '\n')
outReport.write('Spacer 5: ' + str(spacer5) + '\n')
outReport.write('Spacer 6: ' + str(spacer6))

