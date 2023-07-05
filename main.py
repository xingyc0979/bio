from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import os
import re
import csv
import time

def main():
    csv_title = ['sequence', 'protein', 'has_sth']
    record_file = './record.csv'
    # with open('./record.csv', mode='w', newline='') as f:
    #     writer = csv.writer(f)  # 创建csv写入对象
    #     writer.writerow(csv_title)  # 将列表写入csv文件
    #     f.close()

    # result_handle = open("/home/xyc/biopython/results/103/WP_223065108.1.xml")
    # blast_record = NCBIXML.parse(result_handle)
    # for record_seq in blast_record:
    #     for alig in record_seq.alignments:
    #         print(alig.length)
    #         print(alig.accession)
    #         print(alig.title)
    #         print(alig.hit_id)
    #         print(alig.hit_def)
    #     break
    record = SeqIO.parse("test.fasta", format="fasta")
    sequences = []
    for iter in record:
        sequences.append(iter)
        if sequences.__len__() == 1:
            flag=False
            while not flag:
                flag=NCBIWWW.qblast_seqs("blastp", "nr", sequences, format_type="XML", hitlist_size=1000, expect=20000,alignments=1000,descriptions=1000,record_index=1)
                if not flag:
                    print("retry with after 30 secs:",sequences[0].id)
                    time.sleep(3)
            print("iter for 1 sequences :", iter.id)
            sequences.clear()
    NCBIWWW.qblast_seqs("blastp", "nr", sequences, format_type="XML", hitlist_size=1000, expect=20000,alignments=1000,descriptions=1000,record_index=1)
    print("iter for last sequences :", iter.id)
    sequences.clear()
    print("YES!")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
