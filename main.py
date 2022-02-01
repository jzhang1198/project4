# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
from os import listdir

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    path2data = '../project4/data/'
    sequences = [path2data + sequence for sequence in listdir('../project4/data') if 'test' not in sequence]
    sequences.remove(path2data + 'Homo_sapiens_BRD2.fa')
    reference, _ = read_fasta(path2data + 'Homo_sapiens_BRD2.fa')

    sub_matrix_file = '../project4/substitution_matrices/BLOSUM62.mat'
    gap_open = -10
    gap_extend = -1
    f = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)

    scores_alignments = []
    for sequence in sequences:
        splits = sequence.split('/')[-1].split('_')
        species = splits[0] + ' ' + splits[1]
        seq,_ = read_fasta(sequence)
        score,reference_align,seq_align = f.align(reference,seq)
        scores_alignments.append((score,reference_align,seq_align,species))

    scores_alignments.sort(reverse=True)

    for i in range(0,len(scores_alignments)):
        print(f'Alignment score between Homo sapiens BRD2 and {scores_alignments[i][-1]} BRD2 is: {scores_alignments[i][0]}')
        print('')

    return

if __name__ == "__main__":
    main()
