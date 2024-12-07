awk -v gene_name=$1 '
BEGIN { FS="\t"}
{
    # If the line is an operon number, we store it
    if ($1 ~ /^[0-9]+$/) {
        operon_number = $1
    }
    # If the line contains the target gene, print the operon number and gene details
    else if ($2 == gene_name) {
        print "Operon: " operon_number "\n"  $0
    }
    # Print the details of genes when we have a match; useful for printing rest of genes in the operon
    else if (operon_number && $2 ~ /gene_name/) {
        print $0
    }

}' $2
