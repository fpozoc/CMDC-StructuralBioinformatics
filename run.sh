echo "Welcome to CMDC.py You are predicting correlating mutations and distance correlations."

mkdir -p CMDCresults && cd CMDCresults/ && mkdir -p CMDC_{01..03}Example_results && cd ../ &&

for I in CMDC_{01..03}Example_results; do
    cp -f CMDC.py -r CMDCresults/${I}/
    cp -f distances.py -r CMDCresults/${I}/
    cp -f mi.py -r CMDCresults/${I}/
    cp -f extract_sequences.py -r CMDCresults/${I}/
done

cd CMDCresults/CMDC_01Example_results/ && python3 CMDC.py 5cyt -atom CB -seqs 600 -msa clustalo  && rm -r obsolete/ && rm -r __pycache__/ && mkdir -p outputs && mv -f my_predicted_residuecontacts* outputs/ && mv -f *_blast* outputs/ && mkdir -p std.sys/ && mv -f stderr.txt std.sys/ && mv -f stdout.txt std.sys/ && mv -f  *.fa outputs/ && mkdir -p my_scripts && mv -f CMDC.py my_scripts/ && mv -f mi.py my_scripts/ && mv -f extract_sequences.py my_scripts/ && mv -f distances.py my_scripts/ && cd ../ && echo "First example has already finished" &&
cd CMDC_02Example_results/ && python3 CMDC.py 1rbb -gaps -msa muscle  && rm -r obsolete/ && rm -r __pycache__/ && mkdir -p outputs && mv -f my_predicted_residuecontacts* outputs/ && mv -f *_blast* outputs/ && mkdir -p std.sys/ && mv -f stderr.txt std.sys/ && mv -f stdout.txt std.sys/ && mv -f  *.fa outputs/ && mkdir -p my_scripts && mv -f CMDC.py my_scripts/ && mv -f mi.py my_scripts/ && mv -f extract_sequences.py my_scripts/ && mv -f distances.py my_scripts/ && cd ../ && echo "Second example has already finished" &&
cd CMDC_03Example_results/ && python3 CMDC.py 3bp2 -a CA -seqs 50 -msa clustalw  && rm -r obsolete/ && rm -r __pycache__/ && mkdir -p outputs && mv -f my_predicted_residuecontacts* outputs/ && mv -f *_blast* outputs/ && mkdir -p std.sys/ && mv -f stderr.txt std.sys/ && mv -f stdout.txt std.sys/ && mv -f  *.fa outputs/ && mkdir -p my_scripts && mv -f CMDC.py my_scripts/ && mv -f mi.py my_scripts/ && mv -f extract_sequences.py my_scripts/ && mv -f distances.py my_scripts/ && cd ../ && echo "Third example has already finished" &&

echo "All of your examples executed have finished correctly. Please, check that your results in corresponding directories. Also check if all of outputs are totally complete. If not, may be could be a problem of your Internet Conexion."
