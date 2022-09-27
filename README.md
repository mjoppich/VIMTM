

# VIMTM processing steps

## Creating synonym dictionaries
    python3 miRExplore/python/vimtm/create_virus_synonyms.py --obo miRExplore_PMID_PMC/obodir/ncit.obo --syn vimtm/syns/virus.syn

    python3 miRExplore/python/vimtm/create_vim_synonyms.py --virus-syn vimtm/syns/virus.syn --syn vimtm/syns/target_word.syn --out vimtm/syns/vimtm.syn


## Extracting Sentences
    python3 miRExplore/python/vimtm/splitXMLIntoChunks.py vimtm/litcovid/litcovid2pubtator.xml vimtm/litcovid/litcovid2_chunks
    python3 miRExplore/python/vimtm/extract_litcovid.py INT 32 vimtm/litcovid/ &> nohup_extract_litcovid2


## NER
    bash vimtm/runSyngrep.sh "miRExplore/python/textmining/textmineDocument.py --nosuffix --only-longest-match" "./vimtm/results.litcovid/virus_terms/" ./vimtm/litcovid/ "" ./vimtm/syns/virus.syn

    bash vimtm/runSyngrep.sh "miRExplore/python/textmining/textmineDocument.py --only-longest-match" "./vimtm/results.litcovid/vim_targets/" ./vimtm/litcovid/ "" ./vimtm/syns/vimtm.syn


## Acquiring relevant literature
    cat vimtm/results.litcovid/vim_targets/*.index | cut -f 1 | cut -f 1 -d. | sort | uniq > vimtm/tm.pmc.pmid
    #sort vimtm/covid_gold.pmc.pmid vimtm/tm.pmc.pmid vimtm/tm.pmc.pmid | uniq -u   


## Downloading full texts and supplements
    python3 miRExplore/python/vimtm/download_europepmc_xmls.py vimtm/tm.pmc.pmid vimtm/reltexts/ > vimtm/reltexts/conversion.tsv
    python3 miRExplore/python/vimtm/europepmc_supplement_discovery.py vimtm/reltexts/ vimtm/relsents/

    python3 miRExplore/python/vimtm/download_europepmc_xmls.py vimtm/covid_gold.pmc.pmid vimtm/covidtexts/
    python3 miRExplore/python/vimtm/europepmc_supplement_discovery.py vimtm/covidtexts/ vimtm/covidsents/


## Fetching TaxID for documents
    python3 miRExplore/python/vimtm/get_major_document_taxid.py vimtm/results.litcovid/vim_targets/ vimtm/syns/virus.syn_ncit2tax.obo vimtm/syns/virus.syn vimtm/reltexts/conversion.tsv vimtm/doc2tax

    python3 miRExplore/python/vimtm/get_major_document_taxid.py vimtm/results.litcovid/virus_terms/ vimtm/syns/virus.syn_ncit2tax.obo vimtm/syns/virus.syn vimtm/reltexts/conversion.tsv vimtm/doc2tax


## Preparing BLAST database for checking miRNA sequences

    cd vimtm/blastdb

    # export PATH=/mnt/raidtmp/joppich/pubmed_pmc/pmc/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/:$PATH

    /mnt/raidtmp/joppich/pubmed_pmc/pmc/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/update_blastdb.pl --decompress ref_viruses_rep_genomes

    /mnt/raidtmp/joppich/pubmed_pmc/pmc/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/update_blastdb.pl --decompress human_genome
    /mnt/raidtmp/joppich/pubmed_pmc/pmc/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/update_blastdb.pl --decompress Betacoronavirus

    
    # run extract_mirnas.ipynb