import os
with open ('/home/weiboyong/HIV/hiv.fasta', 'r') as fasta_file:#open orgrinal HIV virus fasta reference data
 raw_string=fasta_file.read()
 print(raw_string)
 print(type(raw_string))




print(len(raw_string))
print(len(">NC_001802.1 Human immunodeficiency virus 1, complete genome"))
sequence_string=raw_string.replace(">NC_001802.1 Human immunodeficiency virus 1, complete genome", '')
print(sequence_string)
#sequence_string is the pure atcg sequence

print(len(sequence_string))
total_base_number=len(sequence_string)
print(total_base_number//20)
print(465*20)
start_point=0
end_point=start_point+20
#while end_point != total_base_number:
#twenty_bps_variable=sequence_string[start_point:end_point]
new_ref_string=""
with open('/home/weiboyong/practice/GCF_000001405.39_GRCh38.p13_genomic.fna','r')as reference_file:
    raw_ref_string=reference_file.read()
    for string in raw_ref_string:#t_string means transfer the bases from lower letter ino capital letter
        if string=='a':
            t_string='A'
        elif string=='t':
            t_string='T'
        elif string=='c':
            t_string='C'
        elif string=='g':
            t_string='G'
        else:
            t_string=string
        new_ref_string=new_ref_string+t_string
test_file=open('/home/weiboyong/HIV/test_ref.txt','w')
test_file.write(new_ref_string)
test_file.close()