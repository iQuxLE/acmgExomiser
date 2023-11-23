package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

public class BS1Assigner {

    /*
    //BS1 - "Allele frequency is greater than expected for disorder" - requires disease incidence info from Orphanet and
    // implementation of allele frequency cut-off using CardioDB algorithm. Varsome uses a more tractable implementation
    // which we could do using gnomAD and ClinVar data.
    //
    //Varsome:
    //BS1
    //Allele frequency is greater than expected for disorder. (Benign, Strong)
    //
    //Here we find the highest GnomAD allele frequency for the variant across the main population ethnicities and compare
    // this to the benign cut-off frequency derived from the gene statistics. If there are too few known variants
    // (fewer than 4), we use a much higher default threshold, 0.015, for rare diseases.
    //
    //For mitochondrial variants, a single frequency of 0.005 is used per ClinGen guidelines.
    //
    //In order to avoid double-counting, rule BS1 is not evaluated if either rules BA1 or PM2 were triggered first.
    //
    //Rule BS1 will not trigger if there is strong pathogenic clinical evidence.
     */
}
