package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;

public class BS1BS2Assigner {

    private final VariantDataService variantDataService;

    public BS1BS2Assigner(VariantDataService variantDataService){

        this.variantDataService = variantDataService;
    }

    //BS1 - "Allele frequency is greater than expected for disorder" - requires disease incidence info from Orphanet and
    // implementation of allele frequency cut-off using CardioDB algorithm. Varsome uses a more tractable implementation
    // which we could do using gnomAD and ClinVar data.

    //BS2 - "Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked
    // (hemizygous) disorder, with full penetrance expected at an early age" - needs het/hom counts from gnomAD
    // (these should be in the newer v14 variant store, but I need to double-check)

    public void assign(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation){
        var frequencydata = variantEvaluation.getFrequencyData();
        System.out.println(frequencydata);
    }

}
