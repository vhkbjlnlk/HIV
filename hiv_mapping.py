with open ('/home/weiboyong/HIV/hiv.fasta','r') as HIV:
    HIV_sequence=HIV.read()
x=HIV_sequence.find(">NC_001802.1 Human immunodeficiency virus 1, complete genome")
print(len(">NC_001802.1 Human immunodeficiency virus 1, complete genome"))
HIV_sequence=HIV_sequence[x+60:]
startpoint = 1
endpoint = startpoint + 20
cut_hiv_20bps_sequence = HIV_sequence[startpoint:endpoint]
chromosome_list_with_row_number = []# waiting to fill,use grep in shell to get
chromosome_list = []
with open('/home/weiboyong/HIV/chromosome.txt', 'r') as chromosome_info:
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
with open ('/home/weiboyong/HIV/test_ref.txt','r') as ref:
    reference_sequence=ref.read()


all_chromosome_ATCG_index_dict = {}  # to put every chromosome atcg location info in dictionary as : 'chromosome name':[A base string,T base string]
for chromosome_number in chromosome_list:
    if chromosome_list.index(chromosome_number) == (len(chromosome_list) - 1):
        chromosome_sequence = reference_sequence[
                              ((reference_sequence.find(chromosome_number)) + len(chromosome_number)):]
    else:
        chromosome_sequence = reference_sequence[
                              ((reference_sequence.find(chromosome_number)) + len(chromosome_number)):(
                                  reference_sequence.find(
                                      chromosome_list[(chromosome_list.index(chromosome_number) + 1)]))]

    A_base_single_location = 0
    A_base_list = []
    end_number = 0
    while end_number != 1:
        A_base_single_location = chromosome_sequence.find("A", A_base_single_location) + 1
        end_number = chromosome_sequence.find("A", A_base_single_location) + 2
        A_base_list.append(A_base_single_location)

    T_base_single_location = 0
    T_base_list = []
    end_number = 0
    while end_number != 1:
        T_base_single_location = chromosome_sequence.find("T", T_base_single_location) + 1
        end_number = chromosome_sequence.find("T", T_base_single_location) + 2
        T_base_list.append(T_base_single_location)

    C_base_single_location = 0
    C_base_list = []
    end_number = 0
    while end_number != 1:
        C_base_single_location = chromosome_sequence.find("C", C_base_single_location) + 1
        end_number = chromosome_sequence.find("C", C_base_single_location) + 2
        C_base_list.append(C_base_single_location)

    G_base_single_location = 0
    G_base_list = []
    end_number = 0
    while end_number != 1:
        G_base_single_location = chromosome_sequence.find("G", G_base_single_location) + 1
        end_number = chromosome_sequence.find("G", G_base_single_location) + 2
        G_base_list.append(G_base_single_location)

    all_chromosome_ATCG_index_dict[chromosome_number] = [A_base_list, T_base_list, C_base_list, G_base_list,
                                                         "in order of ATCG"]
sequence_result_list = []
list_to_compare = []
statistic_summary_list = []
alignment_sequence = ""
sequence_length = 20
result_file=open('/home/weiboyong/HIV/result.txt','w')

for compatible_number in range(2, 20):
    final_result = []


    for chromosome_number in chromosome_list:
        startpoint = 0
        endpoint=0
        sequence_result_list=[]
        while endpoint < len(HIV_sequence):
            endpoint = startpoint + compatible_number + 1
            alignment_sequence = HIV_sequence[startpoint:endpoint]
            placement_number = 0
            list_to_compare = []
            for base in alignment_sequence:
                if base == "A":
                    middle_list = (all_chromosome_ATCG_index_dict.get(chromosome_number))[0]
                    middle_list = [i - placement_number for i in middle_list]

                elif base == "T":
                    middle_list = (all_chromosome_ATCG_index_dict.get(chromosome_number))[1]
                    middle_list = [i - placement_number for i in middle_list]

                elif base == "C":
                    middle_list = (all_chromosome_ATCG_index_dict.get(chromosome_number))[2]
                    middle_list = [i - placement_number for i in middle_list]

                else:
                    middle_list = (all_chromosome_ATCG_index_dict.get(chromosome_number))[3]
                    middle_list = [i - placement_number for i in middle_list]
                placement_number = placement_number + 1

                list_to_compare.append(middle_list)
            result = set(list_to_compare[0]).intersection(*list_to_compare[1:(compatible_number + 1)])
            result_list = [result, alignment_sequence[0], startpoint, endpoint,chromosome_number]
            sequence_result_list.append(result_list)
            if len(sequence_result_list) >= (sequence_length-compatible_number):
                if len(sequence_result_list) > (sequence_length-compatible_number):
                    del sequence_result_list[0]
                a = 0
                for result_list in sequence_result_list:
                    if not result_list[0]:
                        a = a + 1
                    else:
                        break
                if a == (sequence_length-compatible_number):
                    final_result.append([HIV_sequence[(sequence_result_list[0])[2]:(sequence_result_list[0])[3]],
                                         (sequence_result_list[0])[2], (sequence_result_list[0])[3],
                                         (sequence_result_list[0])[1],(sequence_result_list[0])[4]])
            startpoint = startpoint + 1

    print(final_result)
    result_file.write(final_result)
    statistic_list = [compatible_number, len(final_result)]
    statistic_summary_list.append(statistic_list)

print(statistic_summary_list)
result_file.write(statistic_summary_list)
result_file.close()