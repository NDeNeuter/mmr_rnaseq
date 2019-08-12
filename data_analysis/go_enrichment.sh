source activate mmr_rna

POPFILE="../data/go_enrichment/uniprot_ids/measured_genes_protein_ids.txt"
GOATOOL="../code/find_enrichment_goatools.py"
OBOFILE="../original_data/go.obo"
ASSFILE="../data/go_enrichment/parsed_goa_human.gaf"

#DIRPATH="../data/go_enrichment/uniprot_ids"
#for FILE in "${DIRPATH}/*protein_ids.txt"
#do
#    echo $FILE
#    OUTFILE=${FILE/protein_ids.txt/out.txt}
#    OUTFILE=${OUTFILE/uniprot_ids/go_results}
#    echo $OUTFILE
#    if [ -e "${OUTFILE}" ]
#    then
#        echo "${OUTFILE} already exists"
#        continue
#    else
#        python "${GOATOOL}" "${FILE}" "${POPFILE}" "${ASSFILE}" --outfile "$OUTFILE" --obo "${OBOFILE}"
#    fi
#    echo
#done

DIRPATH="../data/deseq2"
for FILE in $(find $DIRPATH -name '*_uniprot_ids.txt')
do
    echo $FILE
    OUTFILE=${FILE/uniprot_ids.txt/go.txt}
    echo $OUTFILE
    if [ -e "${OUTFILE}" ]
    then
        echo "${OUTFILE} already exists"
        continue
    else
        python "${GOATOOL}" "${FILE}" "${POPFILE}" "${ASSFILE}" --outfile "$OUTFILE" --obo "${OBOFILE}"
    fi
done
    