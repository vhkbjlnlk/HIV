import os

with open('/mnt/ddngs/errand/HIV/hiv.fasta', 'r') as HIV:
    HIV_sequence = HIV.read()
x = HIV_sequence.find(">NC_001802.1 Human immunodeficiency virus 1, complete genome")
print(len(">NC_001802.1 Human immunodeficiency virus 1, complete genome"))
HIV_sequence = HIV_sequence[x + 60:]
startpoint = 1
endpoint = startpoint + 20
cut_hiv_20bps_sequence = HIV_sequence[startpoint:endpoint]
HIV_sequence = HIV_sequence.replace('\n', '')
print(HIV_sequence)
print(HIV_sequence[0])
chromosome_list_with_row_number = []  # waiting to fill,use grep in shell to get
chromosome_list = []
with open('/mnt/ddngs/errand/HIV/chromosome.txt', 'r') as chromosome_info:
    for line in chromosome_info:
        line = line.strip()

        print(line, type(line))
        chromosome_list_with_row_number.append(line)
print(chromosome_list_with_row_number)
print(chromosome_info)
for info in chromosome_list_with_row_number:
    x = info.find(":")
    info = info[(x + 1):]
    print(info)
    chromosome_list.append(info)
with open('/mnt/ddngs/errand/HIV/test_ref.txt', 'r') as ref:
    reference_sequence = ref.read()
sgRNA_seq_length = 35
print("sadfwe", chromosome_list[-1])
del (chromosome_list[-1])
del (chromosome_list[-1])
print(chromosome_list)


def mkdir(dict_path):
    exist_situation = os.path.exists(dict_path)
    if exist_situation:
        print("folder exists")
    else:
        os.mkdir(dict_path)


overall_result = open("/mnt/ddngs/errand/HIV/mapping_result/overall_result.txt", "w")

for compatible_number in range(10, 21):
    for chromosome_number in chromosome_list:
        compatible_result = [
            ("compatible_number" + "\t" + str(
                compatible_number) + "\t" + "chromosome_number" + "\t" + chromosome_number + '\n')]
        print("compatible list", compatible_result)
        if chromosome_list.index(chromosome_number) == (len(chromosome_list) - 1):
            chromosome_sequence = reference_sequence[
                                  ((reference_sequence.find(chromosome_number)) + len(chromosome_number)):(
                                      reference_sequence.find(
                                          ">NC_012920.1 Homo sApiens miToChondrion, CompleTe Genome"))]
        else:
            chromosome_sequence = reference_sequence[
                                  ((reference_sequence.find(chromosome_number)) + len(chromosome_number)):(
                                      reference_sequence.find(
                                          chromosome_list[(chromosome_list.index(chromosome_number) + 1)]))]

        chromosome_sequence = chromosome_sequence.replace('\n', '')

        set_result_list = []
        startpoint = 0
        endpoint = 0

        while endpoint <= len(HIV_sequence):
            endpoint = startpoint + compatible_number + 1
            set_seq = HIV_sequence[startpoint:endpoint]
            d = chromosome_sequence.find(set_seq)

            set_result_list.append([set_seq, d, startpoint, endpoint])
            if len(set_result_list) >= (sgRNA_seq_length - compatible_number):
                if len(set_result_list) > (sgRNA_seq_length - compatible_number):
                    del set_result_list[0]
            match_number_calculate = 0
            for find_return_value in set_result_list:
                if find_return_value[1] == -1:
                    match_number_calculate = match_number_calculate + 1
                else:
                    break
            if match_number_calculate == (sgRNA_seq_length - compatible_number):
                coordinate = (set_result_list[0])[2]
                sgRNA_seq = HIV_sequence[coordinate:(coordinate + sgRNA_seq_length)]

                compatible_result.append(sgRNA_seq + '\t'+str(coordinate) + '\n')
            startpoint = startpoint + 1
        with open("/mnt/ddngs/errand/HIV/mapping_result/overall_result.txt", "a") as result_file:
            print(compatible_result)
            for x in compatible_result:
                result_file.write(x)
        overall_result.close()
