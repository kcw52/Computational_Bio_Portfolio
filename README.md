# Computational_Bio_Portfolio
This is the portfolio of python code that I learned and created during BISC 5163.

About the author:

Kaylie Canerday is an undergraduate student at Louisiana Tech University. She is double majored in Biology and Environmental Science. 


## Sequence Objects
```python
from Bio.Seq import Seq
```

```python
my_seq = Seq("GATCG")
```

```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G

```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5

```python
print(my_seq[0])
```

    G

```python
print(my_seq[4])
```

    G

```python
print(my_seq[2])
```

    T

```python
# We can also do a . count function
Seq("AAAA").count("AA")
```

    2

```python
# We can also print the length of a new sequence 
my_seq = Seq("GATCGGCTAATCGGTAGCAAATGAAAGTTAG")
```

```python
len(my_seq)
```

    31

```python
# We can also have it count specific letters in the new sequence
my_seq.count("G")
```

    9

```python
# We can measure the GC Count which is good for designing primers
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```

    41.935483870967744

```python
# Common things are already built into Biopython

from Bio.SeqUtils import gc_fraction
```

```python
my_seq = Seq("GATCGATGGGCCTATATGGATCGAAAATCGC")
```

```python
gc_fraction(my_seq)
```

    0.4838709677419355

```python
# We can slice sequences into mulitple parts 
my_seq[0::3]
```

    Seq('GCTGTAGCATC')

```python
my_seq[4:12]
```

    Seq('GATGGGCC')

```python
my_seq[1::3]
```

    Seq('AGGCATAGAC')

```python
my_seq[2:3] 
```

```python
# We can also choose our start position
my_seq[2:3]
```
 
    Seq('T')

```python
# We can use negatives to start from the back end 
my_seq[::-1]
```

    Seq('CGCTAAAAGCTAGGTATATCCGGGTAGCTAG')

```python
# We can turn a seq back into a string
str(my_seq)
```

    'GATCGATGGGCCTATATGGATCGAAAATCGC'

```python
# A verb placeholder can be used
fasta_format_string = ">Name\n%s\n" % my_seq
```

```python
print(fasta_format_string)
```

    >Name
    GATCGATGGGCCTATATGGATCGAAAATCGC
    
```python
# We can also cocantenate 
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```

```python
seq1 + seq2
```

    Seq('ACGTAACCGG')

```python
seq2 + seq1
```

    Seq('AACCGGACGT')

```python
contigs = {Seq("ATG"), Seq("ATCCCG"), Seq("TTTGCA")}
```

```python
from Bio.Seq import Seq
```

```python
contigs = {Seq("ATG"), Seq("ATCCCG"), Seq("TTTGCA")}
```

```python
spacer = Seq("N" *10)
```

```python
spacer.join(contigs)
```

    Seq('ATCCCGNNNNNNNNNNATGNNNNNNNNNNTTTGCA')

```python
dna_seq = Seq("acgtACGT")
```

```python
dna_seq
```

    Seq('acgtACGT')

```python
# We can make sure everything is uppercase
dna_seq.upper()
```

    Seq('ACGTACGT')

```python
# We can make sure everything is lowercase
dna_seq.lower()
```

    Seq('acgtacgt')

```python
dna_seq.upper()
```

    Seq('ACGTACGT')

```python
# We can search for specific codons
"gtac" in dna_seq
```

    False

```python
"GTAC" in dna_seq
```
    False

```python
dna_seq
```

    Seq('acgtACGT')

```python
dna_seq = dna_seq.upper()
```

```python
"GTAC" in dna_seq
```

    True

```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAATCGC")
```

```python
# We can get the complement
my_seq.complement()
```

    Seq('CTAGCTACCCGGATATATCCTAGCTTTAGCG')

```python
# We can get the reverse complement
my_seq.reverse_complement()
```

    Seq('GCGATTTCGATCCTATATAGGCCCATCGATC')

```python
# We can get protein seq.
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```

    Seq('EBYNTM')

```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```

```python
coding_dna
```

    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')

```python
template_dna = coding_dna.reverse_complement()
```

```python
template_dna
```

    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')

```python
coding_dna
```

    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')

```python
messenger_rna = coding_dna.transcribe()
```

```python
messenger_rna
```

    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')

```python
template_dna.reverse_complement().transcribe()
```

    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')

```python
messenger_rna.back_transcribe()
```

    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')

```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# Do transcription
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')

```python
coding_dna.translate(table ="Vertebrate Mitochondrial")
```

    Seq('MAIVMGRWKGAR*')

```python
coding_dna.translate(table = 2)
```

    Seq('MAIVMGRWKGAR*')

```python
coding_dna.translate(to_stop = True)
```

    Seq('MAIVMGR')

```python
coding_dna.translate(table = 2, to_stop=True)
```

    Seq('MAIVMGRWKGAR')

```python
coding_dna.translate(table = 2, stop_symbol = "!")
```

    Seq('MAIVMGRWKGAR!')

```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCTGATGGAGCACAGGCTGCGGAAATTACGTTGTCCGTCAGGTTACGTAAACGGGTTCACTGGTACCCATTGACCTTGCCAGGGTAA")
```

```python
gene.translate(table = "Bacterial")
```

    Seq('VKKMQSIVLALSLVLVAPDGAQAAEITLSVRLRKRVHWYPLTLPG*')

```python
gene.translate(table = "Bacterial", to_stop = True)
```

    Seq('VKKMQSIVLALSLVLVAPDGAQAAEITLSVRLRKRVHWYPLTLPG')

```python
gene.translate(table = "Bacterial", cds = True)
```

    Seq('MKKMQSIVLALSLVLVAPDGAQAAEITLSVRLRKRVHWYPLTLPG')

```python
from Bio.Data import CodonTable
```

```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```

```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```

```python
print (standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--

```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--

```python
mito_table.stop_codons
```

    ['TAA', 'TAG', 'AGA', 'AGG']

```python
mito_table.start_codons
```

    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']

```python
seq = Seq("ACGT")
```

```python
"ACGT" == seq1
```

    True

```python
seq1 == "ACGT"
```

    True

```python
unknown_seq = Seq(None, 10)
```

```python
unknown_seq
```

    Seq(None, length=10)

```python
len(unknown_seq)
```

    10

```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```

```python
seq[1000:1020]
```

    Seq(None, length=20)

```python
seq[117512690:117512700]
```

    Seq('CCTGAATGTG')

```python
seq[117512670:]
```

    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)

```python
seq = Seq("ACGT")
```

```python
undefined_seq = Seq(None, length = 10)
```


```python
seq + undefined_seq +seq
```

    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)

```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```

```python
from Bio.Seq import MutableSeq
```

```python
mutable_seq = MutableSeq(my_seq)
```

```python
mutable_seq
```

    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')

```python
mutable_seq[5] = "C"
```

```python
mutable_seq
```

    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')

```python
mutable_seq.remove("T")
```

```python
mutable_seq
```

    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')

```python
mutable_seq.reverse()
```

```python
mutable_seq
```

    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')

```python
new_seq = Seq(mutable_seq)
```

```python
new_seq
```

    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')

```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```

```python
my_string = "GCTTATGGGTCGTTGGAAGGGTGGTCGCTGCTGGTTAG"
```

```python
reverse_complement(my_string)
```

    'CTAACCAGCAGCGACCACCCTTCCAACGACCCATAAGC'

```python
back_transcribe(my_string)
```

    'GCTTATGGGTCGTTGGAAGGGTGGTCGCTGCTGGTTAG'

```python
translate(my_string)
```
    'AYGSLEGWSLLV'


## Sequence Annotations

