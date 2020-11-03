ssh -t ks2232@compute17.fractus.cs.cornell.edu python3 partitionspread_withchildren_optimized.py;
scp ks2232@compute17.fractus.cs.cornell.edu:/home/ks2232/output.pickle output.pickle;
python3 loadpltfigure.py
