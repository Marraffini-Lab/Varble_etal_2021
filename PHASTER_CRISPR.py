import requests
import json
import csv

#write output files
output = open ('results round 3 - 1.txt', 'w')
output.write( "accession_number"+"\t"+"CRISPR_start"+"\t"+"CRISPR_end"+"\t"+"Phage_start"+"\t"+"Phage_end"+"\t""Phage_score"+"\n")
output_2 = open('accession number with no summarry round 3 - 1.txt', 'w')

output_3 = open ('accession with no json round 3 - 1.txt', 'w')

output_4 = open ('accession that timed out round 3 - 1.txt', 'w')

counter= 0
with open('accession with no json round 3.csv','rU') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        counter = counter + 1
        print (counter)
        accession = str(row["name"])
        CRISPR_start = int(row["CRISPR_start"])
        CRISPR_end_top = int(row["CRISPR_end_top"])
        CRISPR_end_bottom = int(row["CRISPR_end_bottom"])
        CRISPR_length = int(row["CRISPR_length"])
        CRISPR_class=str(row["CRISPR_class"])
        url = "http://phaster.ca/phaster_api?acc="+accession 


        try:
            x = requests.get(url, timeout=10)
        except:
            output_4.write(str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(CRISPR_end_bottom)+"\t"+str(CRISPR_length)+"\t"+str(CRISPR_class)+"\t"+ str(accession)+"\n")
            continue
    
        # decode the json response
        try:
            json_response = x.json()
        
        except ValueError:
            output_3.write(str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(CRISPR_end_bottom)+"\t"+str(CRISPR_length)+"\t"+str(CRISPR_class)+"\t"+ str(accession)+"\n")
            continue
            
            
        #it sees the file as being a dictionary, one entry of the dict would be "summary" whioch contains all the infos we need from the results, if there is no summary that means it was not in the database so the accession number was written to output_2
        try:
            useful_results = json_response['summary']
        
        except KeyError:
            output_2.write(str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(CRISPR_end_bottom)+"\t"+str(CRISPR_length)+"\t"+str(CRISPR_class)+"\t"+ str(accession)+"\n")
            continue



        #split the results into a list for ecah line 
        split_results = useful_results.split("\n")

        #find the position of this line in the results
        try: 
            position = split_results.index("                                 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

        except ValueError:
            continue
            
        #remove the first lines from the results that we dont't need
        trimmed_split_results = split_results[position+1:-1]


        # go through the trimmed list of results
        for item in trimmed_split_results:
            split_item=item.split(" ")
            while("" in split_item) : 
                split_item.remove("")
            
            if split_item[2] != "":
                phage_location_range = split_item[4]
                position_of_dash = phage_location_range.find("-")
                begining_phage = int(phage_location_range [:position_of_dash])
                end_phage = int(phage_location_range [position_of_dash+1:])

        
                if abs(CRISPR_start - begining_phage) <= 2000:
                
                    output.write( str(accession)+"\t"+str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(begining_phage)+"\t"+str(end_phage)+"\t"+str(split_item[2])+"\t"+str(CRISPR_class)+"\n")
            
                elif abs(CRISPR_start - end_phage) <= 2000:
                    output.write( str(accession)+"\t"+str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(begining_phage)+"\t"+str(end_phage)+"\t"+str(split_item[2])+"\t"+str(CRISPR_class)+"\n")
                    
                elif abs(CRISPR_end_top - begining_phage) <= 2000: 
                    output.write( str(accession)+"\t"+str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(begining_phage)+"\t"+str(end_phage)+"\t"+str(split_item[2])+"\t"+str(CRISPR_class)+"\n")
            
                elif abs(CRISPR_end_top - end_phage) <= 2000:
                    output.write( str(accession)+"\t"+str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(begining_phage)+"\t"+str(end_phage)+"\t"+str(split_item[2])+"\t"+str(CRISPR_class)+"\n")
                    
                elif abs(CRISPR_end_bottom - begining_phage) <= 2000: 
                    output.write( str(accession)+"\t"+str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(begining_phage)+"\t"+str(end_phage)+"\t"+str(split_item[2])+"\t"+"*"+str(CRISPR_class)+"\n")
            
                elif abs(CRISPR_end_bottom - end_phage) <= 2000:
                    output.write( str(accession)+"\t"+str(CRISPR_start)+"\t"+str(CRISPR_end_top)+"\t"+str(begining_phage)+"\t"+str(end_phage)+"\t"+str(split_item[2])+"\t"+"*"+str(CRISPR_class)+"\n")
                
                else:
                    pass
        
            
        
            else:
                pass
            
         
output.close()
output_2.close()
output_3.close()
output_4.close()