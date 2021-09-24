import csv

#Genome Size
gen_size = 56025
#Normalize ratio
ratio = 1

#open a wig file where the output of this code will be saved
output_file = open ('CEM_1.wig', 'w')

#create a header needed  in the wig file
output_file.write("track"+"\t"+"type=wiggle_0"+"\n"+"color"+"\t"+"5:150:55"+"\t"+"225:0:0"+"\n"+"variableStep"+"\t"+"chrom=Phiap1"+"\t"+"span=1"+"\n")


#create empty directionary to save mapping data with keys for every bp position and values for the number of reads mapped at that bp
for_strand = {}
rev_strand = {}


                 #fill the dictionary with keys each bp and items to 0
key_bp= 1
for_strand[key_bp] = 0
rev_strand[key_bp] = 0
#here change the value to the genome size (e.g SA2 = 138341)
while key_bp != gen_size :
    key_bp = int(key_bp +1)
    for_strand[key_bp] = 0
    rev_strand[key_bp] = 0



# open the SAM file and specify that the file needs to be sparce based on tab delimitation
with open('CEM_1.sam') as sam_file:
    sam_reader = csv.reader(sam_file, delimiter='\t')
    	
    #my SAM files a weird header made of 3 lines, the code below skip them
    skip_firstline = next(sam_reader)
    skip_secondline = next(sam_reader)
    skip_thirdline = next(sam_reader)
        
    #this loop through each line of the SAM file starting at line 4
    for line in sam_reader :
            
        #if the line is labeled 163 which are reads aligned to top strand of SA2.
        #the code will retrived the starting and end position and will add 1 to each bp location associated with that.
        if int(line[1])==163 :
            for bp in range(int(line[3]), (int(line[7])+1)):
                for_strand[bp] +=1
                
        #same thing as above but for line labeled 99 which are the bottom strand.
        elif int(line[1]) == 99:
            for bp in range(int(line[3]), (int(line[7])+1)):
                rev_strand[bp] +=1
                    
        #if the line is not labeled 163 or 99 that means it was not aligned to the viral genome so skip        
        else:
            pass


for bp in range(1, (gen_size+1)):
    #get forward and reverse read value from each bp location from forward and reverse dictionary and multiply it by the normalization factor
    norm_value_for = int(for_strand.get(bp))*(ratio)
    norm_value_rev = int(rev_strand.get(bp))*(-ratio)
    
    #write a line in the outpu file containing the bp location and then the forward and reverse read value
    output_file.write(str(bp) + "\t" + str(norm_value_for) + "\t" + str(norm_value_rev) + "\n")

output_file.close()

    
    
#163 is top strand
#99 is bottom strannd

