EVAL=../eval_docking.py

#Evaluate the Docking
##########new_dimers#######
##AF2
DOCKQFILES=./dockqstats_newdimers_af2   #Path to file location
OUTFILE=./dockqstats_newdimers_af.csv
#python3 $EVAL --dockqfiles $DOCKQFILES --outfile $OUTFILE

##########Marks###########

##AF2
DOCKQFILES=./dockqstats_marks_af2   #Path to file location
OUTFILE=./dockqstats_marks_af.csv
python3 $EVAL --dockqfiles $DOCKQFILES --outfile $OUTFILE

#########Benchmark 4########
#Model 1
#Evaluate the docking of all different runs
DOCKQFILES=./dockqstats_bench4_af2   #Path to file location with dockq stats
OUTFILE=./bench4_dockqstats.csv #Where to write the merged csv with extracted dockq scores
python3 $EVAL --dockqfiles $DOCKQFILES --outfile $OUTFILE
