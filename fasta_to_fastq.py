read_list=[]#make fake reads data


with open ('hiv.fasta', 'r') as fasta_file:#open orgrinal HIV virus fasta reference data
    for line in fasta_file:
        print (line)
        read_list.append(line)

print(read_list)
del read_list[0] #delete the first which contains description info
print(read_list)

read_list_withoutLinefeed=[] #create empty list to put reads without linefeed
for read in read_list:
    read=read.strip("\n") #remove the linefeed
    print(read)
    read_list_withoutLinefeed.append(read)


print(read_list_withoutLinefeed)

fastq_listfile=[]
read_id_number=1 #make read id number
read_id=""#make read id string as empty
for read in read_list_withoutLinefeed:
    fastq_list=[]
    read_id="@"+str(read_id_number)
    quality_score_per_read = len(read) * "F"
    fastq_string=read_id+","+read+","+"+"+","+quality_score_per_read  #get 4 line info into a string to save into list
    fastq_list=fastq_string.split(",") #the list contains 4 line info for one read
    read_id_number=read_id_number+1
    print(fastq_list)
    fastq_listfile.append(fastq_list)#get the total reads info into a list, each element contains one read info

print(fastq_listfile)
fastq_listfile.pop()
print(fastq_listfile)

fastq_file=open('/home/weiboyong/HIV/hiv.fastq', 'w')#make a new txt document to store fastq file

for read in fastq_listfile:
    for per_read_info in read:
        fastq_file.write(per_read_info)
        fastq_file.write('\n')

fastq_file.close()


#to check if the fq document is right
z=open("/home/weiboyong/HIV/hiv.fastq")
print(z.read())

z.close()
















