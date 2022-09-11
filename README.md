

VIMTM processing steps

    python3 miRExplore/python/vimtm/create_virus_synonyms.py --obo miRExplore_PMID_PMC/obodir/ncit.obo --syn vimtm/syns/virus.syn

    python3 miRExplore/python/vimtm/create_vim_synonyms.py --virus-syn vimtm/syns/virus.syn --syn vimtm/syns/target_word.syn --out vimtm/syns/vimtm.syn

    python3 miRExplore/python/vimtm/splitXMLIntoChunks.py vimtm/litcovid/litcovid2pubtator.xml vimtm/litcovid/litcovid2_chunks

    python3 miRExplore/python/vimtm/extract_litcovid.py INT 32 vimtm/litcovid/ &> nohup_extract_litcovid2


    bash vimtm/runSyngrep.sh miRExplore/python/textmining/textmineDocument.py "./vimtm/results.litcovid/virus_terms/" ./vimtm/litcovid/ "" ./vimtm/syns/virus.syn


    bash vimtm/runSyngrep.sh miRExplore/python/textmining/textmineDocument.py "./vimtm/results.litcovid/vim_targets/" ./vimtm/litcovid/ "" ./vimtm/syns/vimtm.syn

    cat vimtm/results.litcovid/vim_targets/*.index | cut -f 1 | cut -f 1 -d. | sort | uniq > vimtm/tm.pmc.pmid
    sort vimtm/covid_gold.pmc.pmid vimtm/tm.pmc.pmid vimtm/tm.pmc.pmid | uniq -u

    python3 miRExplore/python/vimtm/download_europepmc_xmls.py vimtm/tm.pmc.pmid vimtm/reltexts/
    python3 miRExplore/python/vimtm/europepmc_supplement_discovery.py vimtm/reltexts/ vimtm/relsents/

    python3 miRExplore/python/vimtm/download_europepmc_xmls.py vimtm/covid_gold.pmc.pmid vimtm/covidtexts/
    python3 miRExplore/python/vimtm/europepmc_supplement_discovery.py vimtm/covidtexts/ vimtm/covidsents/


    cd vimtm/blastdb

    /mnt/raidtmp/joppich/pubmed_pmc/pmc/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/update_blastdb.pl --decompress ref_viruses_rep_genomes


    export PATH=/mnt/raidtmp/joppich/pubmed_pmc/pmc/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/:$PATH

    /mnt/raidtmp/joppich/pubmed_pmc/pmc/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/update_blastdb.pl --decompress human_genome
    /mnt/raidtmp/joppich/pubmed_pmc/pmc/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/update_blastdb.pl --decompress Betacoronavirus

    
