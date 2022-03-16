####Benchmark4####
FASTADIR=./dockground/fasta/
META=./dockground/dockground.csv
MODE='bench4'
OUTDIR=./dockground/fasta/merged/
#python3 ./merge_fasta.py --fastadir $FASTADIR --meta $META --mode $MODE --outdir $OUTDIR

#####Marks#######
#Positive set
FASTADIR=./marks/fasta/
META=./marks/final_ids_lens.csv
MODE='marks'
OUTDIR=./marks/fasta/merged/
#python3 ./merge_fasta.py --fastadir $FASTADIR --meta $META --mode $MODE --outdir $OUTDIR

#Negative set
FASTADIR=./marks/marks_negative_set/fasta/
META=./marks/marks_negative_set/neg_set_ids.csv
MODE='marks_neg'
OUTDIR=./marks/marks_negative_set/fasta/merged/
#python3 ./merge_fasta.py --fastadir $FASTADIR --meta $META --mode $MODE --outdir $OUTDIR

#####CASP14#######
FASTADIR=./CASP14/fasta/
META=./CASP14/casp14_multimeric_targets.csv
MODE='CASP14'
OUTDIR=./CASP14/fasta/merged/
#python3 ./merge_fasta.py --fastadir $FASTADIR --meta $META --mode $MODE --outdir $OUTDIR


######Negatome#######
FASTADIR=./negatome/fasta/
META=./negatome/negatome_ids.csv
MODE='negatome'
OUTDIR=./negatome/fasta/merged/
python3 ./merge_fasta.py --fastadir $FASTADIR --meta $META --mode $MODE --outdir $OUTDIR
