
if [ $# -le 1 ]
then
	echo "Usage: $0 name_of_directory <tool>"
	echo "tool can be 'scarpa', 'encore' or 'hapsembler'"
	exit
fi

DNAME=$1

mkdir ./$DNAME
mkdir ./$DNAME/bin

if [ "$2" == scarpa  ]
then
	mkdir ./$DNAME/test
	cp BiEdge.h ContigEdge.h DirNode.h DnaRead.h HapSuite.h HapUtils.h Naive.h ReadAln.h Scaffer.h SmithWaterman.h UEdge.h UGraph.h UNode.h ./$DNAME
	cp BiEdge.cpp ContigEdge.cpp DirNode.cpp DnaRead.cpp HapUtils.cpp Naive.cpp ReadAln.cpp Scaffer.cpp SmithWaterman.cpp UEdge.cpp UGraph.cpp UNode.cpp ./$DNAME
	cp BiNode.h ContigNode.h GraphDef.h LinkedIter.h LinkedList.h LinkedNode.h ScafferMain.cpp ./$DNAME
	cp lp_Hash.h lp_lib.h lp_matrix.h lp_mipbb.h lp_SOS.h lp_types.h lp_utils.h liblpsolve55.a liblpsolve55.so ./$DNAME
	cp makefile_scarpa ./$DNAME/makefile
	cp scarpa_parser ./$DNAME/bin
	cp scarpa_process ./$DNAME/bin
	cp ./documentation/TEST.README ./$DNAME/test
	cat ./documentation/SCARPA.README | awk '{sub(/scarpa-X.X/, "'$DNAME'"); print}' > ./$DNAME/SCARPA.README
	cp ./documentation/LICENSE ./$DNAME
	cp -r ./test ./$DNAME/
	chmod 644 ./$DNAME/*.h
	chmod 644 ./$DNAME/*.cpp
	chmod 644 ./$DNAME/test/*.*
	chmod 744 ./$DNAME/bin/scarpa_parser
	chmod 744 ./$DNAME/bin/scarpa_process
fi

if [ "$2" == encore  ]
then
	cp Dna.h DnaRead.h HapSuite.h HapUtils.h KmerHash.h KmerIter.h KmerList.h Naive.h ReadAln.h ReadMatch.h SmithWaterman.h ./$DNAME
	cp Dna.cpp DnaRead.cpp HapUtils.cpp KmerHash.cpp KmerIter.cpp KmerList.cpp Naive.cpp ReadAln.cpp ReadMatch.cpp SmithWaterman.cpp ./$DNAME
	cp GraphDef.h KmerNode.h DnaAlignment.cpp ./$DNAME
	cp makefile_encore ./$DNAME/makefile
	cat ./documentation/ENCORE.README | awk '{sub(/encore-X.X/, "'$DNAME'"); print}' > ./$DNAME/ENCORE.README
	cp ./documentation/LICENSE ./$DNAME
	chmod 644 ./$DNAME/*.h
	chmod 644 ./$DNAME/*.cpp
fi

if [ "$2" == hapsembler  ]
then
	mkdir ./$DNAME/doc
	mkdir ./$DNAME/sample
	cp BiEdge.h BiNode.h Containee.h ContigEdge.h ContigNode.h DirNode.h Dna.h DnaRead.h DnaTempl.h EdgeTuple.h GraphDef.h HapSuite.h HapUtils.h KmerHash.h KmerIter.h KmerList.h KmerNode.h ./$DNAME
	cp LinkedIter.h LinkedList.h LinkedNode.h MateEdge.h MatePair.h Naive.h OvlEdge.h OvlGraph.h OvlGraphTempl.h PathEdge.h PathNode.h ReadAln.h ReadMatch.h ReadNode.h ./$DNAME
	cp Scaffer.h SmithWaterman.h UEdge.h UGraph.h UNode.h ./$DNAME
	cp BiEdge.cpp Containee.cpp ContigEdge.cpp DirNode.cpp Dna.cpp DnaAlignment.cpp DnaRead.cpp EdgeTuple.cpp HapsemblerCore.cpp HapUtils.cpp KmerHash.cpp KmerIter.cpp KmerList.cpp ./$DNAME
	cp MateEdge.cpp Naive.cpp OvlEdge.cpp OvlGraph.cpp PathEdge.cpp ReadAln.cpp ReadMatch.cpp Scaffer.cpp ScafferMain.cpp SmithWaterman.cpp UEdge.cpp UGraph.cpp UNode.cpp ./$DNAME
	cp lp_Hash.h lp_lib.h lp_matrix.h lp_mipbb.h lp_SOS.h lp_types.h lp_utils.h liblpsolve55.a liblpsolve55.so ./$DNAME
	cp hapsemble ./$DNAME/bin/hapsemble
	cp makefile ./$DNAME
	cp ./sample/sample_illu.fastq ./$DNAME/sample/
	cp ./sample/sample_illu.info ./$DNAME/sample/
	cat ./documentation/SAMPLE.README | awk '{sub(/hapsembler-X.X/, "'$DNAME'"); print}' > ./$DNAME/sample/SAMPLE.README
	cat ./documentation/README | awk '{sub(/hapsembler-X.X/, "'$DNAME'"); print}' > ./$DNAME/README
	cp ./documentation/INSTALL ./$DNAME
	cp ./documentation/LICENSE ./$DNAME
	cp ./documentation/HapsemblerDocumentation.pdf ./$DNAME/doc/${DNAME}_manual.pdf
	cp scarpa_parser ./$DNAME/bin
	cp scarpa_process ./$DNAME/bin
	chmod 644 ./$DNAME/*.h
	chmod 644 ./$DNAME/*.cpp
	chmod 644 ./$DNAME/doc/*.pdf
	chmod 644 ./$DNAME/sample/*.*
	chmod 744 ./$DNAME/bin/hapsemble
	chmod 744 ./$DNAME/bin/scarpa_parser
	chmod 744 ./$DNAME/bin/scarpa_process
fi

tar -zcvf $DNAME".tar.gz" ./$DNAME
rm -r ./$DNAME
