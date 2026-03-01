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


## Sequence I/O

## Multiple Sequence Alignment

```python
# https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/PF05371_seed.sth
```


```python
from Bio import AlignIO
```


```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print("Alignment length %i" % alignment.get_alignment_length())
```

    Alignment length 52



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
align1 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTGCTAGCTAG"), id = "Alpha"),
        SeqRecord(Seq("ACT-CTAGCTAG"), id = "Beta"),
        SeqRecord(Seq("ACTGCTAGDTAG"), id = "Gamma"),
    ]
)
align2 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("GTCAGC-AG"), id = "Delta"),
        SeqRecord(Seq("GACAGCTAG"), id = "Epsilon"),
        SeqRecord(Seq("GTCAGCTAG"), id = "Zeta"),
    ]
)
align3 = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTAGTACAGCTG"), id = "Eta"),
        SeqRecord(Seq("ACTAGTACAGCT-"), id = "Theta"),
        SeqRecord(Seq("-CTACTACAGGTG"), id = "Iota"),
    ]
)
```


```python
my_alignments = [align1, align2, align3]
```


```python
my_alignments
```




    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7fac001d6210>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7fac001d6110>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7fac001d6310>]




```python
print(my_alignments)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7fac001d6210>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7fac001d6110>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7fac001d6310>]



```python
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
alignments = list(AlignIO.parse("my_example.phy", "phylip" ))
```


```python
last_align = alignments[-1]
```


```python
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
first_align = alignments[0]
```


```python
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
from Bio import AlignIO
```


```python
count = AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
alignments = AlignIO.parse("PF05371_seed.sth.txt", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i 
```


```python
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
AlignIO.write([alignment], "PF05371_seed.phy.txt", "phylip")
```




    1




```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print("Number of rows: %i" % len(alignment))
```

    Number of rows: 7



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
print(alignment[2, 6])
```

    T



```python
print(alignment[2].seq[6])
```

    T



```python
print(alignment[:, 6])
```

    TTT---T



```python
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
print(alignment[:, :6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited.sort()
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Align import MultipleSeqAlignment
```


```python
alignment = MultipleSeqAlignment(
[
    SeqRecord(Seq("ACTCCTA"), id = "seq1"),
    SeqRecord(Seq("AAT-CTA"), id = "seq2"), 
    SeqRecord(Seq("CCTACT-"), id = "seq3"),
    SeqRecord(Seq("TCTCCTC"), id = "seq4"),

])
```


```python
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
substitutions = alignment.substitutions
```


```python
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
m = substitutions.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
me = substitutions.select("ACTG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
import Bio.Align.Applications
```


```python
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
# https://github.com/biopython/biopython/blob/master/Doc/examples/opuntia.aln
```


```python
from Bio import AlignIO
```


```python
align = AlignIO.read("opuntia.aln", "clustal")
```


```python
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
from Bio import Phylo
```


```python
tree = Phylo.read("opuntia.dnd", "newick")
```


```python
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    



```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
alligner = Align.PairwiseAligner(match_score = 1.0)
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    3.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
aligner.mode = "local"
```


```python
target = "AGAACTC"
```


```python
query = "GAACT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    5.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
aligner.algorithm
```




    'Smith-Waterman'




```python
aligner.epsilon
```




    1e-06




```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
alignments = aligner.align(target, query)
```


```python
alignment = alignments[0]
```


```python
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7fecdd11c610>




```python
alignment.score
```




    3.0




```python
alignment.target
```




    'GAACT'




```python
alignment.query
```




    'GAT'




```python
print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    



```python
alignment.coordinates
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
len(alignment)
```




    2




```python
alignment.shape
```




    (2, 5)




```python
aligner.mode = "local"
```


```python
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
local_alignment = local_alignments[0]
```


```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
local_alignment.shape
```




    (2, 4)




```python
aligner.mode = "global"
```


```python
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
```


```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
local_alignment.sort()
```


```python
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    



```python
from Bio import Align
```


```python
from Bio.Seq import reverse_complement
```


```python
target = "AAACCC"
```


```python
query = "AACC"
```


```python
aligner = Align.PairwiseAligner(mismatch_score = -1, internal_gap_score = -1)
```


```python
aligner.score(target, query)
```




    4.0




```python
aligner.score(target, reverse_complement(query))
```




    0.0




```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
aligner.score(target, query, strand = "-")
```




    0.0




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, query, strand = "-")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
print(alignments[1])
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
aligner.left_gap_score = -0.5
```


```python
aligner.right_gap_score = -0.2
```


```python
aligner.score(target, query)
```




    3.3




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments)
```

    <Bio.Align.PairwiseAlignments object at 0x7fecdd11cd90>



```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
aligner.score(target, query, strand = "+")
```




    3.3

```python

```

## Blast

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import os


record = SeqIO.read("m_cold.fasta", format="fasta")


os.makedirs("blast_results", exist_ok=True)
output_file = "blast_results/my_blast.xml"


NCBIWWW.email = "kcc025@email.latech.edu"


result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

blast_xml = result_handle.read()   # read **once**
result_handle.close()               # close the handle
with open(output_file, "w") as out_handle:
    out_handle.write(blast_xml)

print("BLAST results saved to:", output_file)

with open(output_file) as result_handle:
    blast_record = NCBIXML.read(result_handle)

print("Query:", blast_record.query)
```

    BLAST results saved to: blast_results/my_blast.xml
    Query: No definition line



```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


record = SeqIO.read("m_cold.fasta", "fasta")


clean_record = SeqRecord(
    seq=record.seq,
    id="MP14H09",
    description="Mesembryanthemum crystallinum cold acclimation cDNA"
)


SeqIO.write(clean_record, "m_cold_clean.fasta", "fasta")
```


    1


```python
from Bio.Blast import NCBIWWW, NCBIXML
import os


os.makedirs("blast_results", exist_ok=True)
output_file = "blast_results/my_blast.xml"


NCBIWWW.email = "kcc025@email.latech.edu"
record = SeqIO.read("m_cold_clean.fasta", "fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)


blast_xml = result_handle.read()
result_handle.close()
with open(output_file, "w") as out_handle:
    out_handle.write(blast_xml)


with open(output_file) as result_handle:
    blast_record = NCBIXML.read(result_handle)

print("Query:", blast_record.query)  # Now shows your clean ID
```

    Query: No definition line



```python
s
import os
print(os.path.exists("m_cold_clean.fasta"))  # Should print True


with open("m_cold_clean.fasta") as f:
    for _ in range(3):
        print(f.readline().strip())
```

    True
    >MP14H09 Mesembryanthemum crystallinum cold acclimation cDNA
    CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNTGTGAAC
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTA



```python
s
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import os
```


```python

record = SeqIO.read("m_cold.fasta", "fasta")

clean_record = SeqRecord(
    seq=record.seq,
    id="MP14H09",
    description="Mesembryanthemum crystallinum cold acclimation cDNA"
)


clean_fasta_file = "m_cold_clean.fasta"
SeqIO.write(clean_record, clean_fasta_file, "fasta")
print("Cleaned FASTA saved:", clean_fasta_file)
```

    Cleaned FASTA saved: m_cold_clean.fasta



```python

os.makedirs("blast_results", exist_ok=True)
blast_xml_file = "blast_results/my_blast.xml"
print("BLAST results will be saved to:", blast_xml_file)
```

    BLAST results will be saved to: blast_results/my_blast.xml



```python
# Run remote BLAST 
NCBIWWW.email = "kcc025@email.latech.edu"  # Replace with your email

# Run BLAST and read results 
record = SeqIO.read(clean_fasta_file, "fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

# Save BLAST results to XML 
blast_xml = result_handle.read()
result_handle.close()

with open(blast_xml_file, "w") as out_handle:
    out_handle.write(blast_xml)

print("BLAST results saved successfully!")
```

    BLAST results saved successfully!



```python
# Parse the saved BLAST XML
with open(blast_xml_file) as result_handle:
    blast_record = NCBIXML.read(result_handle)

print("Query sequence:", blast_record.query)
print("Number of hits:", len(blast_record.alignments))
```

    Query sequence: No definition line
    Number of hits: 50



```python
#  Print top 5 hits and alignment info
for alignment in blast_record.alignments[:5]:  # top 5 hits
    print("Hit:", alignment.hit_def)
    for hsp in alignment.hsps[:1]:  # top HSP per hit
        print("Score:", hsp.score, "E-value:", hsp.expect)
        print("Query snippet:", hsp.query[:50], "...")
    print("-"*50)
```

    Hit: PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
    Score: 482.0 E-value: 6.88175e-117
    Query snippet: ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAAT ...
    --------------------------------------------------
    Hit: PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
    Score: 463.0 E-value: 1.84663e-111
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------
    Hit: PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
    Score: 443.0 E-value: 4.9552e-106
    Query snippet: TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA ...
    --------------------------------------------------
    Hit: PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
    Score: 441.0 E-value: 1.72953e-105
    Query snippet: AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCG ...
    --------------------------------------------------
    Hit: PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    Score: 439.0 E-value: 6.03667e-105
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------



```python
print("Number of hits:", len(blast_record.alignments))

for alignment in blast_record.alignments[:5]:  # top 5 hits
    print("Hit:", alignment.hit_def)
    for hsp in alignment.hsps[:1]:  # top HSP per hit
        print("Score:", hsp.score, "E-value:", hsp.expect)
        print("Query snippet:", hsp.query[:50], "...")
    print("-"*50)
```

    Number of hits: 50
    Hit: PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
    Score: 482.0 E-value: 6.88175e-117
    Query snippet: ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAAT ...
    --------------------------------------------------
    Hit: PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
    Score: 463.0 E-value: 1.84663e-111
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------
    Hit: PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
    Score: 443.0 E-value: 4.9552e-106
    Query snippet: TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA ...
    --------------------------------------------------
    Hit: PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
    Score: 441.0 E-value: 1.72953e-105
    Query snippet: AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCG ...
    --------------------------------------------------
    Hit: PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    Score: 439.0 E-value: 6.03667e-105
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------



```python
blast_record.query = clean_record.id
print("Query sequence:", blast_record.query)
```

    Query sequence: MP14H09



```python
for alignment in blast_record.alignments[:5]:
    print("Hit:", alignment.hit_def)
    for hsp in alignment.hsps[:1]:
        print("Score:", hsp.score, "E-value:", hsp.expect)
        print("Query snippet:", hsp.query[:50], "...")
    print("-"*50)
```

    Hit: PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
    Score: 482.0 E-value: 6.88175e-117
    Query snippet: ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAAT ...
    --------------------------------------------------
    Hit: PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
    Score: 463.0 E-value: 1.84663e-111
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------
    Hit: PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
    Score: 443.0 E-value: 4.9552e-106
    Query snippet: TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA ...
    --------------------------------------------------
    Hit: PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
    Score: 441.0 E-value: 1.72953e-105
    Query snippet: AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCG ...
    --------------------------------------------------
    Hit: PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    Score: 439.0 E-value: 6.03667e-105
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------



```python
# Optional – Print top 5 hits and alignment info
for alignment in blast_record.alignments[:5]:  # top 5 hits
    print("Hit:", alignment.hit_def)
    for hsp in alignment.hsps[:1]:  # top HSP per hit
        print("Score:", hsp.score, "E-value:", hsp.expect)
        print("Query snippet:", hsp.query[:50], "...")
    print("-"*50)
```

    Hit: PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
    Score: 482.0 E-value: 6.88175e-117
    Query snippet: ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAAT ...
    --------------------------------------------------
    Hit: PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
    Score: 463.0 E-value: 1.84663e-111
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------
    Hit: PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
    Score: 443.0 E-value: 4.9552e-106
    Query snippet: TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA ...
    --------------------------------------------------
    Hit: PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
    Score: 441.0 E-value: 1.72953e-105
    Query snippet: AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCG ...
    --------------------------------------------------
    Hit: PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    Score: 439.0 E-value: 6.03667e-105
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------



```python
# Inspect top BLAST hits
for alignment in blast_record.alignments[:5]:  # top 5 hits
    print("Hit:", alignment.hit_def)
    for hsp in alignment.hsps[:1]:  # top HSP per hit
        print("Score:", hsp.score, "E-value:", hsp.expect)
        print("Query snippet:", hsp.query[:50], "...")
    print("-"*50)
```

    Hit: PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
    Score: 482.0 E-value: 6.88175e-117
    Query snippet: ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAAT ...
    --------------------------------------------------
    Hit: PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
    Score: 463.0 E-value: 1.84663e-111
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------
    Hit: PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
    Score: 443.0 E-value: 4.9552e-106
    Query snippet: TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA ...
    --------------------------------------------------
    Hit: PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
    Score: 441.0 E-value: 1.72953e-105
    Query snippet: AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCG ...
    --------------------------------------------------
    Hit: PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    Score: 439.0 E-value: 6.03667e-105
    Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC ...
    --------------------------------------------------



```python
# Inspect top 10 BLAST hits
top_hits = 10  # Number of hits to display

if len(blast_record.alignments) == 0:
    print("No BLAST hits found.")
else:
    print(f"Showing top {min(top_hits, len(blast_record.alignments))} hits:\n")
    for i, alignment in enumerate(blast_record.alignments[:top_hits], start=1):
        print(f"Hit #{i}: {alignment.hit_def}")  # Description of the hit
        for hsp in alignment.hsps[:1]:  # Only show the top HSP per hit
            print(f"  Score: {hsp.score}, E-value: {hsp.expect}")
            print(f"  Query snippet: {hsp.query[:50]}...")  # First 50 bases
        print("-"*60)
```

    Showing top 10 hits:
    
    Hit #1: PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
      Score: 482.0, E-value: 6.88175e-117
      Query snippet: ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAAT...
    ------------------------------------------------------------
    Hit #2: PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
      Score: 463.0, E-value: 1.84663e-111
      Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC...
    ------------------------------------------------------------
    Hit #3: PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
      Score: 443.0, E-value: 4.9552e-106
      Query snippet: TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA...
    ------------------------------------------------------------
    Hit #4: PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
      Score: 441.0, E-value: 1.72953e-105
      Query snippet: AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCG...
    ------------------------------------------------------------
    Hit #5: PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
      Score: 439.0, E-value: 6.03667e-105
      Query snippet: AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGC...
    ------------------------------------------------------------
    Hit #6: PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122306609), transcript variant X1, mRNA
      Score: 434.0, E-value: 7.35417e-104
      Query snippet: ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGT...
    ------------------------------------------------------------
    Hit #7: PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122306609), transcript variant X2, mRNA
      Score: 434.0, E-value: 7.35417e-104
      Query snippet: ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGT...
    ------------------------------------------------------------
    Hit #8: PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2 (LOC108995251), transcript variant X2, mRNA
      Score: 430.0, E-value: 8.95921e-103
      Query snippet: AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT-----...
    ------------------------------------------------------------
    Hit #9: PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2 (LOC108995251), transcript variant X1, mRNA
      Score: 430.0, E-value: 8.95921e-103
      Query snippet: AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT-----...
    ------------------------------------------------------------
    Hit #10: PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X3, mRNA
      Score: 413.0, E-value: 6.88783e-98
      Query snippet: AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCC...
    ------------------------------------------------------------



```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|1219041180|ref|XM_021875076.1| PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
    length: 1173
    e value: 6.88175e-117
    ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTC...
    || ||||||||| |||| | |||| ||  |||| |||| | |||| ||| | |||| ||| ||| ||||| | ||...
    ACCGAAAATGGGCAGAGGAGTGAATTATATGGCAATGACACCTGAGCAACTAGCCGCGGCCAATTTGATCAACTC...
    ****ALIGNMENT****
    sequence: gi|2514617377|ref|XM_021992092.2| PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
    length: 752
    e value: 1.84663e-111
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||| |||  |||| | || ||||| |||||||| || ||||| |||| ||| ||| ||||||||||||||...
    AAAATGGGTAGACGAATGGATTATTTGGCGATGAAAACCGAGCAATTAGCCGCGGCCAATTTGATCGATTCCGAT...
    ****ALIGNMENT****
    sequence: gi|2518612504|ref|XM_010682658.3| PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
    length: 621
    e value: 4.9552e-106
    TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACA...
    ||||||||||||||||| ||| ||||  |||||||| |||| ||||  ||||| ||||| ||||| || ||    ...
    TTGGCCATGAAAACTGAGCAAATGGCGTTGGCTAATTTGATAGATTATGATATGAATGAACTTAAGATCGCTTTG...
    ****ALIGNMENT****
    sequence: gi|2031543140|ref|XM_041168865.1| PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
    length: 1020
    e value: 1.72953e-105
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    ||||||||| |||  | |  | |||||||||||||||||||    ||||  |||  || ||||||| || |||| ...
    AATGGGGAG-GAA--GGATAATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATAA...
    ****ALIGNMENT****
    sequence: gi|2618480339|ref|XM_048479995.2| PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    length: 1028
    e value: 6.03667e-105
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  ||||| |||| |||||||| |   |||  |||| |  ||||  |||| |||...
    AAAATGGGGAGG---ATGGAGTTTTTGGCTATGAGAACTGATCCA---GCCACGGCTGACTTGATAAATTCTGAT...
    ****ALIGNMENT****
    sequence: gi|2082357253|ref|XM_043119041.1| PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122306609), transcript variant X1, mRNA
    length: 1020
    e value: 7.35417e-104
    ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATC...
    |||||||| |||  | | || |||||||||||||||||||    ||||  |||  || ||||||| || ||||||...
    ATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATATC...
    ****ALIGNMENT****
    sequence: gi|2082357255|ref|XM_043119049.1| PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122306609), transcript variant X2, mRNA
    length: 1036
    e value: 7.35417e-104
    ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATC...
    |||||||| |||  | | || |||||||||||||||||||    ||||  |||  || ||||||| || ||||||...
    ATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATATC...
    ****ALIGNMENT****
    sequence: gi|1882610310|ref|XM_035691634.1| PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2 (LOC108995251), transcript variant X2, mRNA
    length: 909
    e value: 8.95921e-103
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT---------GGCCGTGGCTAATATGATCGA...
    ||||||||| |||  | | || |||||||||||||||||||             ||||  |||  || |||||||...
    AATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCCGGCCACGGCCACGGCCACGGCGGATTTGATCGA...
    ****ALIGNMENT****
    sequence: gi|1882610309|ref|XM_018970776.2| PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2 (LOC108995251), transcript variant X1, mRNA
    length: 1025
    e value: 8.95921e-103
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT---------GGCCGTGGCTAATATGATCGA...
    ||||||||| |||  | | || |||||||||||||||||||             ||||  |||  || |||||||...
    AATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCCGGCCACGGCCACGGCCACGGCGGATTTGATCGA...
    ****ALIGNMENT****
    sequence: gi|2395983798|ref|XM_006466623.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X3, mRNA
    length: 1052
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1350315641|ref|XM_024180293.1| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X4, mRNA
    length: 868
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1350315634|ref|XM_006425717.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X1, mRNA
    length: 952
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1350315636|ref|XM_006425716.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X2, mRNA
    length: 881
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2395983796|ref|XM_025094967.2| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X1, mRNA
    length: 980
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2395983800|ref|XM_006466626.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X5, mRNA
    length: 913
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2395983799|ref|XM_006466625.3| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X4, mRNA
    length: 978
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1350315638|ref|XM_006425719.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X3, mRNA
    length: 893
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2395983797|ref|XM_006466624.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X2, mRNA
    length: 968
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1204884098|ref|XM_021445554.1| PREDICTED: Herrania umbratica cold-regulated 413 plasma membrane protein 2-like (LOC110429488), mRNA
    length: 905
    e value: 6.88783e-98
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || |||||  ||| ||||...
    AAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGATCAGTTCTGATA...
    ****ALIGNMENT****
    sequence: gi|1227938481|ref|XM_022049453.1| PREDICTED: Carica papaya cold-regulated 413 plasma membrane protein 2-like (LOC110820077), mRNA
    length: 1009
    e value: 2.92878e-96
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    |||||||||||||    ||| | || ||||| ||||| ||||||||   ||||   ||| || | |||  ||| |...
    AGAAAATGGGGAGG---ATGGAATATTTGGCTATGAAGACTGATCA---GGCCACTGCTGATCTCATCACTTCTG...
    ****ALIGNMENT****
    sequence: gi|1063463253|ref|XM_007047033.2| PREDICTED: Theobroma cacao cold-regulated 413 plasma membrane protein 2 (LOC18611025), transcript variant X2, mRNA
    length: 1071
    e value: 1.24535e-94
    TGTGAACAGA-AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGAT...
    || |||| || |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || ||||...
    TGAGAACTGAGAAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGAT...
    ****ALIGNMENT****
    sequence: gi|1063463252|ref|XM_007047032.2| PREDICTED: Theobroma cacao cold-regulated 413 plasma membrane protein 2 (LOC18611025), transcript variant X1, mRNA
    length: 1065
    e value: 1.24535e-94
    TGTGAACAGA-AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGAT...
    || |||| || |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || ||||...
    TGAGAACTGAGAAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGAT...
    ****ALIGNMENT****
    sequence: gi|1269881407|ref|XM_022895605.1| PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X3, mRNA
    length: 1069
    e value: 4.3467e-94
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ****ALIGNMENT****
    sequence: gi|1269881403|ref|XM_022895603.1| PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X1, mRNA
    length: 1072
    e value: 4.3467e-94
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ****ALIGNMENT****
    sequence: gi|1269881405|ref|XM_022895604.1| PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X2, mRNA
    length: 1091
    e value: 4.3467e-94
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ****ALIGNMENT****
    sequence: gi|2082386146|ref|XM_043113302.1| PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122301958), transcript variant X2, mRNA
    length: 824
    e value: 1.51715e-93
    ATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAA...
    ||||| |||||||| |||||||||||||    |||  |||   || ||||||  || ||||||||||| || || ...
    ATGAATTACTTGGCTATGAAAACTGATCC---GGCAATGGAGGATTTGATCGGCTCTGATATCAATGACCTCAAG...
    ****ALIGNMENT****
    sequence: gi|2082386143|ref|XM_043113301.1| PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122301958), transcript variant X1, mRNA
    length: 844
    e value: 1.51715e-93
    ATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAA...
    ||||| |||||||| |||||||||||||    |||  |||   || ||||||  || ||||||||||| || || ...
    ATGAATTACTTGGCTATGAAAACTGATCC---GGCAATGGAGGATTTGATCGGCTCTGATATCAATGACCTCAAG...
    ****ALIGNMENT****
    sequence: gi|1954740698|ref|XM_038867092.1| PREDICTED: Tripterygium wilfordii cold-regulated 413 plasma membrane protein 2 (LOC120014952), mRNA
    length: 999
    e value: 5.29536e-93
    GAACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGAT...
    ||| ||||||||||||||   | | | || ||||| ||||| |||||||    ||  ||||   || |||||   ...
    GAAAAGAAAATGGGGAGA---ACGGATTATTTGGCGATGAAGACTGATCC---GGTTGTGGACGATTTGATCAGC...
    ****ALIGNMENT****
    sequence: gi|1882636119|ref|XM_018974650.2| PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC108998174), mRNA
    length: 1015
    e value: 6.45107e-92
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    |||||||||    ||||| || ||||| |||||||||||||    |||  ||| | || ||||||  || |||||...
    AATGGGGAGG---ATGAATTATTTGGCTATGAAAACTGATCC---GGCAATGGATGATTTGATCGGCTCTGATAT...
    ****ALIGNMENT****
    sequence: gi|2526866810|ref|XM_057645500.1| PREDICTED: Actinidia eriantha cold-regulated 413 plasma membrane protein 2-like (LOC130785340), mRNA
    length: 1152
    e value: 6.45107e-92
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    ||||||||||||   ||| | || ||||| ||||| || |||| |  |||  || | ||| ||||| |||| || ...
    AAAATGGGGAGA---ATGGATTATTTGGCGATGAAGACCGATCCAGCGGC--TGCCGAAT-TGATCAATTCGGAC...
    ****ALIGNMENT****
    sequence: gi|1187397285|gb|KX009413.1| Santalum album COR413-PM2 mRNA, complete cds
    length: 837
    e value: 2.25165e-91
    AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    |||||||||    ||| | | ||||||||||||||| ||||    |||||  |   || ||||| ||||||| ||...
    AATGGGGAGG---ATGGATTTCTTGGCCATGAAAACAGATCCCGCGGCCGCCG---ATTTGATCAATTCCGACAT...
    ****ALIGNMENT****
    sequence: gi|2550782781|ref|XM_058372567.1| PREDICTED: Rhododendron vialii cold-regulated 413 plasma membrane protein 2 (LOC131336659), mRNA
    length: 1110
    e value: 2.74307e-90
    GCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGT...
    ||||  ||| | |||||||| || |||||||| ||||| || || ||  | | | | || || |   | || |  ...
    GCCGATGCTGAAATGATCGACTCGGATATCAACGAGCTGAAGATCGCGGCCAAGCGACTGATTAGCCACGCCACC...
    ****ALIGNMENT****
    sequence: gi|2806124758|ref|XM_068481225.1| PREDICTED: Pyrus communis cold-regulated 413 plasma membrane protein 2-like (LOC137741519), mRNA
    length: 850
    e value: 9.57424e-90
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    |||| ||||| |||||||| ||||| || || ||| | | ||| ||||||| |||||| |  | ||| |||  ||...
    TGATAGATTCAGATATCAAAGAGCTCAAGATTGCAGCCAAGAGACTCATCAGTGATGCCACCAAGCTTGGTGGTT...
    ****ALIGNMENT****
    sequence: gi|2532162279|ref|XM_058104265.1| PREDICTED: Malania oleifera cold-regulated 413 plasma membrane protein 2-like (LOC131152402), mRNA
    length: 2364
    e value: 1.16638e-88
    GAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA...
    ||||||||||||      | |||||||||||||||||||||||| |  ||| |  |   || |||||  ||| ||...
    GAAAATGGGGAGGTC---GGAGTACTTGGCCATGAAAACTGATCCAGCGGCTGCCG---ATTTGATCAGTTCGGA...
    ****ALIGNMENT****
    sequence: gi|2250518185|ref|XM_009343631.3| PREDICTED: Pyrus x bretschneideri cold-regulated 413 plasma membrane protein 2 (LOC103933927), mRNA
    length: 787
    e value: 1.16638e-88
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    |||| ||||| |||||||| ||||| || || ||| | | ||| ||||||| |||||| |  | ||| |||  ||...
    TGATAGATTCAGATATCAAAGAGCTCAAGATTGCAGCCAAGAGACTCATCAGTGATGCCACCAAGCTTGGTGGTT...
    ****ALIGNMENT****
    sequence: gi|1350280614|ref|XM_024170292.1| PREDICTED: Morus notabilis cold-regulated 413 plasma membrane protein 2 (LOC21394987), mRNA
    length: 1020
    e value: 4.07107e-88
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||||||||| ||       || |||||||||||||| || | |   |||  |||| || ||||  |||| ||||...
    AAATGGGGAGGGAT------TATTTGGCCATGAAAACGGACCCA---GCCACGGCTGATTTGATAAATTCTGATA...
    ****ALIGNMENT****
    sequence: gi|743838297|ref|XM_011027373.1| PREDICTED: Populus euphratica cold-regulated 413 plasma membrane protein 2 (LOC105126500), transcript variant X2, mRNA
    length: 1132
    e value: 4.07107e-88
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || |||||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGATTCCGAT...
    ****ALIGNMENT****
    sequence: gi|743838293|ref|XM_011027372.1| PREDICTED: Populus euphratica cold-regulated 413 plasma membrane protein 2 (LOC105126500), transcript variant X1, mRNA
    length: 980
    e value: 4.07107e-88
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || |||||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGATTCCGAT...
    ****ALIGNMENT****
    sequence: gi|1768569081|ref|XM_031406607.1| PREDICTED: Pistacia vera cold-regulated 413 plasma membrane protein 2-like (LOC116120644), mRNA
    length: 982
    e value: 4.95958e-87
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT-GGCCGTGGCTAATATGATCGATTCCGA...
    |||||||||||    ||| | ||  |||  ||||||||||| ||||  ||     ||| |  ||||  | || ||...
    AAAATGGGGAGG---ATGGATTATCTGGGAATGAAAACTGA-CAATCAGGTTACTGCTGAGGTGATTAACTCTGA...
    ****ALIGNMENT****
    sequence: gi|2396494060|ref|XM_052454347.1| PREDICTED: Populus trichocarpa cold-regulated 413 plasma membrane protein 2 (LOC18101203), transcript variant X1, mRNA
    length: 1018
    e value: 1.73106e-86
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ****ALIGNMENT****
    sequence: gi|2396494064|ref|XM_024605027.2| PREDICTED: Populus trichocarpa cold-regulated 413 plasma membrane protein 2 (LOC18101203), transcript variant X2, mRNA
    length: 1178
    e value: 1.73106e-86
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ****ALIGNMENT****
    sequence: gi|1585724761|ref|XM_028202722.1| PREDICTED: Camellia sinensis cold-regulated 413 plasma membrane protein 2-like (LOC114262355), mRNA
    length: 910
    e value: 6.042e-86
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    |||||||||||||  ||||| |||| ||||| ||||| || |||||    |||    |  |   |||  ||||||...
    AGAAAATGGGGAGGAAAATGGAGTATTTGGCAATGAAGACCGATCATCCAGCCCCAACCCAATCGATGAATTCCG...
    ****ALIGNMENT****
    sequence: gi|2537663858|ref|XM_021815584.2| PREDICTED: Hevea brasiliensis cold-regulated 413 plasma membrane protein 2 (LOC110658100), mRNA
    length: 945
    e value: 2.10887e-85
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    ||||||||||    ||| ||||||||   ||||  ||||||||| |  |   |||  || ||||| | || ||| ...
    AAATGGGGAGG---ATGGAGTACTTGAAAATGAGTACTGATCAAGTACC---GGCCGATTTGATCAAGTCTGATC...
    ****ALIGNMENT****
    sequence: gi|2960782598|ref|XM_035077206.2| PREDICTED: Populus alba cold-regulated 413 plasma membrane protein 2 (LOC118063227), mRNA
    length: 937
    e value: 2.10887e-85
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   | |||| | || || ||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGGTAATTTAATTGAGTCCGAT...
    ****ALIGNMENT****
    sequence: gi|2645357626|ref|XM_062094449.1| PREDICTED: Populus nigra cold-regulated 413 plasma membrane protein 2-like (LOC133673573), mRNA
    length: 1175
    e value: 2.10887e-85
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ****ALIGNMENT****
    sequence: gi|1162571918|ref|XM_007202530.2| PREDICTED: Prunus persica cold-regulated 413 plasma membrane protein 2 (LOC18770198), transcript variant X1, mRNA
    length: 811
    e value: 2.56913e-84
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    ||||  |||| || |||||||| || || || ||| | | ||  |||||||||||||| | || ||| |||   |...
    TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAATGATGCCACTAAGCTTGGTGGGT...
    ****ALIGNMENT****
    sequence: gi|1162571919|ref|XM_020568695.1| PREDICTED: Prunus persica cold-regulated 413 plasma membrane protein 2 (LOC18770198), transcript variant X2, mRNA
    length: 929
    e value: 2.56913e-84
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    ||||  |||| || |||||||| || || || ||| | | ||  |||||||||||||| | || ||| |||   |...
    TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAATGATGCCACTAAGCTTGGTGGGT...
    ****ALIGNMENT****
    sequence: gi|2583747300|ref|XM_059787294.1| PREDICTED: Cornus florida cold-regulated 413 plasma membrane protein 2-like (LOC132285128), mRNA
    length: 1126
    e value: 2.56913e-84
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    |||||||||||||| |   | |||| ||||| |||||||||||||    ||||   ||  |  ||||| ||||||...
    AGAAAATGGGGAGAAA---GGAGTATTTGGCTATGAAAACTGATCC---GGCCACAGCCGAATTGATCAATTCCG...
    ****ALIGNMENT****
    sequence: gi|1229761331|ref|XM_022277554.1| PREDICTED: Momordica charantia cold-regulated 413 plasma membrane protein 2-like (LOC111005887), mRNA
    length: 850
    e value: 8.96713e-84
    ATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATTACGGGT...
    |||| |||||||| ||||||||||| ||| | | ||||||  |  |  |||| |  | |||||||    | ||  ...
    ATTCTGATATCAACGAGCTTAAAATTGCAGCCACGAGGCTTCTTGAACATGCCACCAAGCTCGGTGGAAAGGGCC...
    ****ALIGNMENT****
    sequence: gi|2118882425|ref|XM_044613294.1| PREDICTED: Mangifera indica cold-regulated 413 plasma membrane protein 2-like (LOC123198583), mRNA
    length: 1083
    e value: 3.12984e-83
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    || ||||||||    ||| | ||  |||| |||||||||||  |   ||  |  ||| |  ||||| | || |||...
    AAGATGGGGAGG---ATGGATTATCTGGCAATGAAAACTGACGATCAGGTTGCTGCTGACTTGATCAACTCTGAT...



```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("Sequence:", alignment.hit_def)
            print("Length:", alignment.length)
            print("E-value:", hsp.expect)
            print("Query  :", hsp.query[0:75] + "...")
            print("Match  :", hsp.match[0:75] + "...")
            print("Subject:", hsp.sbjct[0:75] + "...")
            print("-"*60)
```

    ****ALIGNMENT****
    Sequence: PREDICTED: Chenopodium quinoa cold-regulated 413 plasma membrane protein 2-like (LOC110697660), mRNA
    Length: 1173
    E-value: 6.88175e-117
    Query  : ACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTC...
    Match  : || ||||||||| |||| | |||| ||  |||| |||| | |||| ||| | |||| ||| ||| ||||| | ||...
    Subject: ACCGAAAATGGGCAGAGGAGTGAATTATATGGCAATGACACCTGAGCAACTAGCCGCGGCCAATTTGATCAACTC...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Spinacia oleracea cold-regulated 413 plasma membrane protein 2-like (LOC110787470), mRNA
    Length: 752
    E-value: 1.84663e-111
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : |||||||| |||  |||| | || ||||| |||||||| || ||||| |||| ||| ||| ||||||||||||||...
    Subject: AAAATGGGTAGACGAATGGATTATTTGGCGATGAAAACCGAGCAATTAGCCGCGGCCAATTTGATCGATTCCGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Beta vulgaris subsp. vulgaris cold-regulated 413 plasma membrane protein 2 (LOC104895996), mRNA
    Length: 621
    E-value: 4.9552e-106
    Query  : TTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACA...
    Match  : ||||||||||||||||| ||| ||||  |||||||| |||| ||||  ||||| ||||| ||||| || ||    ...
    Subject: TTGGCCATGAAAACTGAGCAAATGGCGTTGGCTAATTTGATAGATTATGATATGAATGAACTTAAGATCGCTTTG...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Juglans microcarpa x Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC121265293), mRNA
    Length: 1020
    E-value: 1.72953e-105
    Query  : AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    Match  : ||||||||| |||  | |  | |||||||||||||||||||    ||||  |||  || ||||||| || |||| ...
    Subject: AATGGGGAG-GAA--GGATAATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATAA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    Length: 1028
    E-value: 6.03667e-105
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : |||||||||||    ||| |||  ||||| |||| |||||||| |   |||  |||| |  ||||  |||| |||...
    Subject: AAAATGGGGAGG---ATGGAGTTTTTGGCTATGAGAACTGATCCA---GCCACGGCTGACTTGATAAATTCTGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122306609), transcript variant X1, mRNA
    Length: 1020
    E-value: 7.35417e-104
    Query  : ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATC...
    Match  : |||||||| |||  | | || |||||||||||||||||||    ||||  |||  || ||||||| || ||||||...
    Subject: ATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATATC...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122306609), transcript variant X2, mRNA
    Length: 1036
    E-value: 7.35417e-104
    Query  : ATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATC...
    Match  : |||||||| |||  | | || |||||||||||||||||||    ||||  |||  || ||||||| || ||||||...
    Subject: ATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCC---GGCCACGGCGGATTTGATCGACTCTGATATC...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2 (LOC108995251), transcript variant X2, mRNA
    Length: 909
    E-value: 8.95921e-103
    Query  : AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT---------GGCCGTGGCTAATATGATCGA...
    Match  : ||||||||| |||  | | || |||||||||||||||||||             ||||  |||  || |||||||...
    Subject: AATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCCGGCCACGGCCACGGCCACGGCGGATTTGATCGA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2 (LOC108995251), transcript variant X1, mRNA
    Length: 1025
    E-value: 8.95921e-103
    Query  : AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT---------GGCCGTGGCTAATATGATCGA...
    Match  : ||||||||| |||  | | || |||||||||||||||||||             ||||  |||  || |||||||...
    Subject: AATGGGGAG-GAA--GGATTATTTGGCCATGAAAACTGATCCGGCCACGGCCACGGCCACGGCGGATTTGATCGA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X3, mRNA
    Length: 1052
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X4, mRNA
    Length: 868
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X1, mRNA
    Length: 952
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X2, mRNA
    Length: 881
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X1, mRNA
    Length: 980
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X5, mRNA
    Length: 913
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X4, mRNA
    Length: 978
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X3, mRNA
    Length: 893
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X2, mRNA
    Length: 968
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||    |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    Subject: AAATGGGGAGAT---TGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Herrania umbratica cold-regulated 413 plasma membrane protein 2-like (LOC110429488), mRNA
    Length: 905
    E-value: 6.88783e-98
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || |||||  ||| ||||...
    Subject: AAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGATCAGTTCTGATA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Carica papaya cold-regulated 413 plasma membrane protein 2-like (LOC110820077), mRNA
    Length: 1009
    E-value: 2.92878e-96
    Query  : AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    Match  : |||||||||||||    ||| | || ||||| ||||| ||||||||   ||||   ||| || | |||  ||| |...
    Subject: AGAAAATGGGGAGG---ATGGAATATTTGGCTATGAAGACTGATCA---GGCCACTGCTGATCTCATCACTTCTG...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Theobroma cacao cold-regulated 413 plasma membrane protein 2 (LOC18611025), transcript variant X2, mRNA
    Length: 1071
    E-value: 1.24535e-94
    Query  : TGTGAACAGA-AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGAT...
    Match  : || |||| || |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || ||||...
    Subject: TGAGAACTGAGAAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Theobroma cacao cold-regulated 413 plasma membrane protein 2 (LOC18611025), transcript variant X1, mRNA
    Length: 1065
    E-value: 1.24535e-94
    Query  : TGTGAACAGA-AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGAT...
    Match  : || |||| || |||||||||||   ||| | || ||||| |||||||| ||||   | || |  |   || ||||...
    Subject: TGAGAACTGAGAAATGGGGAGA---ATGGACTATTTGGCTATGAAAACAGATCCTGTAGCAGAAG---ATTTGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X3, mRNA
    Length: 1069
    E-value: 4.3467e-94
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    Subject: AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X1, mRNA
    Length: 1072
    E-value: 4.3467e-94
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    Subject: AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Durio zibethinus cold-regulated 413 plasma membrane protein 2 (LOC111300020), transcript variant X2, mRNA
    Length: 1091
    E-value: 4.3467e-94
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||||   ||| |||| ||||| |||||||||||||   | || |  |  ||| |||||  ||| ||||...
    Subject: AAATGGGGAGA---ATGGAGTATTTGGCTATGAAAACTGATCCTGTAGCTGAAG--AAT-TGATCAGTTCTGATA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122301958), transcript variant X2, mRNA
    Length: 824
    E-value: 1.51715e-93
    Query  : ATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAA...
    Match  : ||||| |||||||| |||||||||||||    |||  |||   || ||||||  || ||||||||||| || || ...
    Subject: ATGAATTACTTGGCTATGAAAACTGATCC---GGCAATGGAGGATTTGATCGGCTCTGATATCAATGACCTCAAG...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Carya illinoinensis cold-regulated 413 plasma membrane protein 2-like (LOC122301958), transcript variant X1, mRNA
    Length: 844
    E-value: 1.51715e-93
    Query  : ATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAA...
    Match  : ||||| |||||||| |||||||||||||    |||  |||   || ||||||  || ||||||||||| || || ...
    Subject: ATGAATTACTTGGCTATGAAAACTGATCC---GGCAATGGAGGATTTGATCGGCTCTGATATCAATGACCTCAAG...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Tripterygium wilfordii cold-regulated 413 plasma membrane protein 2 (LOC120014952), mRNA
    Length: 999
    E-value: 5.29536e-93
    Query  : GAACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGAT...
    Match  : ||| ||||||||||||||   | | | || ||||| ||||| |||||||    ||  ||||   || |||||   ...
    Subject: GAAAAGAAAATGGGGAGA---ACGGATTATTTGGCGATGAAGACTGATCC---GGTTGTGGACGATTTGATCAGC...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Juglans regia cold-regulated 413 plasma membrane protein 2-like (LOC108998174), mRNA
    Length: 1015
    E-value: 6.45107e-92
    Query  : AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    Match  : |||||||||    ||||| || ||||| |||||||||||||    |||  ||| | || ||||||  || |||||...
    Subject: AATGGGGAGG---ATGAATTATTTGGCTATGAAAACTGATCC---GGCAATGGATGATTTGATCGGCTCTGATAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Actinidia eriantha cold-regulated 413 plasma membrane protein 2-like (LOC130785340), mRNA
    Length: 1152
    E-value: 6.45107e-92
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : ||||||||||||   ||| | || ||||| ||||| || |||| |  |||  || | ||| ||||| |||| || ...
    Subject: AAAATGGGGAGA---ATGGATTATTTGGCGATGAAGACCGATCCAGCGGC--TGCCGAAT-TGATCAATTCGGAC...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: Santalum album COR413-PM2 mRNA, complete cds
    Length: 837
    E-value: 2.25165e-91
    Query  : AATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATAT...
    Match  : |||||||||    ||| | | ||||||||||||||| ||||    |||||  |   || ||||| ||||||| ||...
    Subject: AATGGGGAGG---ATGGATTTCTTGGCCATGAAAACAGATCCCGCGGCCGCCG---ATTTGATCAATTCCGACAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Rhododendron vialii cold-regulated 413 plasma membrane protein 2 (LOC131336659), mRNA
    Length: 1110
    E-value: 2.74307e-90
    Query  : GCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGT...
    Match  : ||||  ||| | |||||||| || |||||||| ||||| || || ||  | | | | || || |   | || |  ...
    Subject: GCCGATGCTGAAATGATCGACTCGGATATCAACGAGCTGAAGATCGCGGCCAAGCGACTGATTAGCCACGCCACC...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Pyrus communis cold-regulated 413 plasma membrane protein 2-like (LOC137741519), mRNA
    Length: 850
    E-value: 9.57424e-90
    Query  : TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    Match  : |||| ||||| |||||||| ||||| || || ||| | | ||| ||||||| |||||| |  | ||| |||  ||...
    Subject: TGATAGATTCAGATATCAAAGAGCTCAAGATTGCAGCCAAGAGACTCATCAGTGATGCCACCAAGCTTGGTGGTT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Malania oleifera cold-regulated 413 plasma membrane protein 2-like (LOC131152402), mRNA
    Length: 2364
    E-value: 1.16638e-88
    Query  : GAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA...
    Match  : ||||||||||||      | |||||||||||||||||||||||| |  ||| |  |   || |||||  ||| ||...
    Subject: GAAAATGGGGAGGTC---GGAGTACTTGGCCATGAAAACTGATCCAGCGGCTGCCG---ATTTGATCAGTTCGGA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Pyrus x bretschneideri cold-regulated 413 plasma membrane protein 2 (LOC103933927), mRNA
    Length: 787
    E-value: 1.16638e-88
    Query  : TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    Match  : |||| ||||| |||||||| ||||| || || ||| | | ||| ||||||| |||||| |  | ||| |||  ||...
    Subject: TGATAGATTCAGATATCAAAGAGCTCAAGATTGCAGCCAAGAGACTCATCAGTGATGCCACCAAGCTTGGTGGTT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Morus notabilis cold-regulated 413 plasma membrane protein 2 (LOC21394987), mRNA
    Length: 1020
    E-value: 4.07107e-88
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : |||||||||| ||       || |||||||||||||| || | |   |||  |||| || ||||  |||| ||||...
    Subject: AAATGGGGAGGGAT------TATTTGGCCATGAAAACGGACCCA---GCCACGGCTGATTTGATAAATTCTGATA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Populus euphratica cold-regulated 413 plasma membrane protein 2 (LOC105126500), transcript variant X2, mRNA
    Length: 1132
    E-value: 4.07107e-88
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || |||||||||...
    Subject: AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGATTCCGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Populus euphratica cold-regulated 413 plasma membrane protein 2 (LOC105126500), transcript variant X1, mRNA
    Length: 980
    E-value: 4.07107e-88
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || |||||||||...
    Subject: AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGATTCCGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Pistacia vera cold-regulated 413 plasma membrane protein 2-like (LOC116120644), mRNA
    Length: 982
    E-value: 4.95958e-87
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATT-GGCCGTGGCTAATATGATCGATTCCGA...
    Match  : |||||||||||    ||| | ||  |||  ||||||||||| ||||  ||     ||| |  ||||  | || ||...
    Subject: AAAATGGGGAGG---ATGGATTATCTGGGAATGAAAACTGA-CAATCAGGTTACTGCTGAGGTGATTAACTCTGA...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Populus trichocarpa cold-regulated 413 plasma membrane protein 2 (LOC18101203), transcript variant X1, mRNA
    Length: 1018
    E-value: 1.73106e-86
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    Subject: AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Populus trichocarpa cold-regulated 413 plasma membrane protein 2 (LOC18101203), transcript variant X2, mRNA
    Length: 1178
    E-value: 1.73106e-86
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    Subject: AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Camellia sinensis cold-regulated 413 plasma membrane protein 2-like (LOC114262355), mRNA
    Length: 910
    E-value: 6.042e-86
    Query  : AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    Match  : |||||||||||||  ||||| |||| ||||| ||||| || |||||    |||    |  |   |||  ||||||...
    Subject: AGAAAATGGGGAGGAAAATGGAGTATTTGGCAATGAAGACCGATCATCCAGCCCCAACCCAATCGATGAATTCCG...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Hevea brasiliensis cold-regulated 413 plasma membrane protein 2 (LOC110658100), mRNA
    Length: 945
    E-value: 2.10887e-85
    Query  : AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    Match  : ||||||||||    ||| ||||||||   ||||  ||||||||| |  |   |||  || ||||| | || ||| ...
    Subject: AAATGGGGAGG---ATGGAGTACTTGAAAATGAGTACTGATCAAGTACC---GGCCGATTTGATCAAGTCTGATC...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Populus alba cold-regulated 413 plasma membrane protein 2 (LOC118063227), mRNA
    Length: 937
    E-value: 2.10887e-85
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : |||||||||||    ||| |||  |||   ||||| |||||| |    | |   | |||| | || || ||||||...
    Subject: AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGGTAATTTAATTGAGTCCGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Populus nigra cold-regulated 413 plasma membrane protein 2-like (LOC133673573), mRNA
    Length: 1175
    E-value: 2.10887e-85
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : |||||||||||    ||| |||  |||   ||||| |||||| |    | |   |||||| | || || ||||||...
    Subject: AAAATGGGGAGG---ATGGAGTTTTTGAAGATGAAGACTGATGATGAAGTCAGCGCTAATTTAATTGAGTCCGAT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Prunus persica cold-regulated 413 plasma membrane protein 2 (LOC18770198), transcript variant X1, mRNA
    Length: 811
    E-value: 2.56913e-84
    Query  : TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    Match  : ||||  |||| || |||||||| || || || ||| | | ||  |||||||||||||| | || ||| |||   |...
    Subject: TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAATGATGCCACTAAGCTTGGTGGGT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Prunus persica cold-regulated 413 plasma membrane protein 2 (LOC18770198), transcript variant X2, mRNA
    Length: 929
    E-value: 2.56913e-84
    Query  : TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    Match  : ||||  |||| || |||||||| || || || ||| | | ||  |||||||||||||| | || ||| |||   |...
    Subject: TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAATGATGCCACTAAGCTTGGTGGGT...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Cornus florida cold-regulated 413 plasma membrane protein 2-like (LOC132285128), mRNA
    Length: 1126
    E-value: 2.56913e-84
    Query  : AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    Match  : |||||||||||||| |   | |||| ||||| |||||||||||||    ||||   ||  |  ||||| ||||||...
    Subject: AGAAAATGGGGAGAAA---GGAGTATTTGGCTATGAAAACTGATCC---GGCCACAGCCGAATTGATCAATTCCG...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Momordica charantia cold-regulated 413 plasma membrane protein 2-like (LOC111005887), mRNA
    Length: 850
    E-value: 8.96713e-84
    Query  : ATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATTACGGGT...
    Match  : |||| |||||||| ||||||||||| ||| | | ||||||  |  |  |||| |  | |||||||    | ||  ...
    Subject: ATTCTGATATCAACGAGCTTAAAATTGCAGCCACGAGGCTTCTTGAACATGCCACCAAGCTCGGTGGAAAGGGCC...
    ------------------------------------------------------------
    ****ALIGNMENT****
    Sequence: PREDICTED: Mangifera indica cold-regulated 413 plasma membrane protein 2-like (LOC123198583), mRNA
    Length: 1083
    E-value: 3.12984e-83
    Query  : AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    Match  : || ||||||||    ||| | ||  |||| |||||||||||  |   ||  |  ||| |  ||||| | || |||...
    Subject: AAGATGGGGAGG---ATGGATTATCTGGCAATGAAAACTGACGATCAGGTTGCTGCTGACTTGATCAACTCTGAT...
    ------------------------------------------------------------

```python

```

## Challenge 1 - Blast of KRAS Gene

```python
from Bio import Entrez

Entrez.email = "kcc025@email.latech.edu"

handle = Entrez.efetch(
    db="nucleotide",
    id="NM_004985.5",  # KRAS mRNA
    rettype="fasta",
    retmode="text"
)

sequence = handle.read()

with open("KRAS.fasta", "w") as f:
    f.write(sequence)

print("KRAS.fasta downloaded")
```

    KRAS.fasta downloaded



```python
from Bio import Entrez

Entrez.email = "kcc025@email.latech.edu"  # REQUIRED by NCBI

handle = Entrez.efetch(
    db="nucleotide",
    id="NM_004985.5",  # KRAS mRNA accession
    rettype="fasta",
    retmode="text"
)

sequence = handle.read()

with open("KRAS.fasta", "w") as f:
    f.write(sequence)

print("KRAS.fasta downloaded successfully")
```

    KRAS.fasta downloaded successfully



```python
import os
print(os.path.exists("KRAS.fasta"))
```

    True



```python
from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(
    query="KRAS.fasta",
    db="your_database_prefix",
    out="kras_blast.xml",
    outfmt=5
)

stdout, stderr = blastn_cline()
print(stderr)
```


```python
from Bio.Blast.Applications import NcbimakeblastdbCommandline

make_db = NcbimakeblastdbCommandline(
    input_file="m_cold.fasta",  # <-- database FASTA
    dbtype="nucl",
    out="localdb"               # <-- database prefix name
)

stdout, stderr = make_db()
print("makeblastdb stderr:\n", stderr)
```

    makeblastdb stderr:
     



```python
import glob
print(glob.glob("localdb.n*"))  # should show localdb.nsq etc
```

    ['localdb.nin', 'localdb.nsq', 'localdb.nhr']



```python
from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(
    query="KRAS.fasta",
    db="localdb",
    out="kras_blast.xml",
    outfmt=5
)

stdout, stderr = blastn_cline()
print("blastn stderr:\n", stderr)
```

    blastn stderr:
     



```python
from Bio.Blast import NCBIXML
with open("kras_blast.xml") as handle:
    record = NCBIXML.read(handle)
print("Alignments:", len(record.alignments))
```

    Alignments: 0



```python
from Bio import SeqIO
rec = next(SeqIO.parse("KRAS.fasta", "fasta"))
print(rec.id)
print("Length:", len(rec.seq))
print("First 60 bp:", rec.seq[:60])
```

    NM_004985.5
    Length: 5306
    First 60 bp: CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCGGCGGCAGTGGCGGCGGCGAAGGT



```python
import os
print("XML size:", os.path.getsize("kras_blast.xml"))
```

    XML size: 2046



```python
from Bio.Blast.Applications import NcbimakeblastdbCommandline

make_db = NcbimakeblastdbCommandline(
    input_file="KRAS.fasta",
    dbtype="nucl",
    out="kras_db"
)
stdout, stderr = make_db()
print(stderr)
```

    



```python
import glob, os
print("DB files:", glob.glob("kras_db.n*"))
for f in glob.glob("kras_db.n*"):
    print(f, os.path.getsize(f))
```

    DB files: ['kras_db.nsq', 'kras_db.nin', 'kras_db.nhr']
    kras_db.nsq 1328
    kras_db.nin 88
    kras_db.nhr 150



```python
from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(
    query="KRAS.fasta",
    db="kras_db",
    out="kras_vs_kras.xml",
    outfmt=5
)
stdout, stderr = blastn_cline()
print("blastn stderr:\n", stderr)
```

    blastn stderr:
     



```python
from Bio.Blast import NCBIXML

with open("kras_vs_kras.xml") as handle:
    record = NCBIXML.read(handle)

print("Alignments:", len(record.alignments))
if record.alignments:
    print("Top hit:", record.alignments[0].title)
```

    Alignments: 1
    Top hit: gnl|BL_ORD_ID|0 NM_004985.5 Homo sapiens KRAS proto-oncogene, GTPase (KRAS), transcript variant b, mRNA



```python
for alignment in record.alignments[:5]:
    print("Hit:", alignment.hit_def)
    for hsp in alignment.hsps[:1]:
        print("Score:", hsp.score, "E-value:", hsp.expect)
        print("Query snippet:", hsp.query[:50], "...")
    print("-"*50)
```

    Hit: NM_004985.5 Homo sapiens KRAS proto-oncogene, GTPase (KRAS), transcript variant b, mRNA
    Score: 5306.0 E-value: 0.0
    Query snippet: CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCGGCGGCAGTGGCGG ...
    --------------------------------------------------



```python

for alignment in record.alignments[:5]:  # top 5 hits
    print("Hit:", alignment.hit_def)
    for hsp in alignment.hsps[:1]:  # top HSP per hit
        print("Score:", hsp.score, "E-value:", hsp.expect)
        print("Query snippet:", hsp.query[:50], "...")
    print("-"*50)
```

    Hit: NM_004985.5 Homo sapiens KRAS proto-oncogene, GTPase (KRAS), transcript variant b, mRNA
    Score: 5306.0 E-value: 0.0
    Query snippet: CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCGGCGGCAGTGGCGG ...
    --------------------------------------------------



```python

top_hits = 10  # Number of hits to display

if len(record.alignments) == 0:
    print("No BLAST hits found.")
else:
    print(f"Showing top {min(top_hits, len(record.alignments))} hits:\n")
    for i, alignment in enumerate(record.alignments[:top_hits], start=1):
        print(f"Hit #{i}: {alignment.hit_def}")  # Description of the hit
        for hsp in alignment.hsps[:1]:  # Only show the top HSP per hit
            print(f"  Score: {hsp.score}, E-value: {hsp.expect}")
            print(f"  Query snippet: {hsp.query[:50]}...")  # First 50 bases
        print("-"*60)
```

    Showing top 1 hits:
    
    Hit #1: NM_004985.5 Homo sapiens KRAS proto-oncogene, GTPase (KRAS), transcript variant b, mRNA
      Score: 5306.0, E-value: 0.0
      Query snippet: CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCGGCGGCAGTGGCGG...
    ------------------------------------------------------------



```python
E_VALUE_THRESH = 0.04

for alignment in record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("Hit:", alignment.hit_def)
            print("Score:", hsp.score)
            print("E-value:", hsp.expect)
            print("-"*50)
```

    Hit: NM_004985.5 Homo sapiens KRAS proto-oncogene, GTPase (KRAS), transcript variant b, mRNA
    Score: 5306.0
    E-value: 0.0
    --------------------------------------------------



```python
E_VALUE_THRESH = 0.04

for alignment in record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gnl|BL_ORD_ID|0 NM_004985.5 Homo sapiens KRAS proto-oncogene, GTPase (KRAS), transcript variant b, mRNA
    length: 5306
    e value: 0.0
    CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCGGCGGCAGTGGCGGCGGCGAAGGTGGCGGCGGCTCGGCC...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCGGCGGCAGTGGCGGCGGCGAAGGTGGCGGCGGCTCGGCC...



```python
# The evalue when compared to chimpanzee is 0.0
```
