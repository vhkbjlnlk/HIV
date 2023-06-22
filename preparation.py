'''chromosome_list_with_row_number = []
chromosome_list=[]
with open('/home/weiboyong/HIV/chromosome.txt', 'r') as chromosome_info:
    for line in chromosome_info:
        line = line.strip()

        print(line, type(line))
        chromosome_list_with_row_number.append(line)
print(chromosome_list_with_row_number)
print(chromosome_info)
for info in chromosome_list_with_row_number:
    x=info.find(":")
    info=info[(x+1):]
    print(info)
    chromosome_list.append(info)


with open ('/home/weiboyong/HIV/test_ref.txt','r') as ref:
    reference_sequence=ref.read()
print(reference_sequence)
'''
with open ('/home/weiboyong/HIV/hiv.fasta','r') as HIV:
    HIV_sequence=HIV.read()
x=HIV_sequence.find(">NC_001802.1 Human immunodeficiency virus 1, complete genome")
print(len(">NC_001802.1 Human immunodeficiency virus 1, complete genome"))
HIV_sequence=HIV_sequence[x+60:]
print(HIV_sequence)