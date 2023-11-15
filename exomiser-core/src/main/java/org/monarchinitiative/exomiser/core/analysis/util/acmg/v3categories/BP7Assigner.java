package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import htsjdk.variant.variantcontext.VariantContext;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;

public class BP7Assigner {

    private final VariantDataService variantDataService;

    public BP7Assigner(VariantDataService variantDataService) {
        this.variantDataService = variantDataService;
    }

    //"A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus
    // sequence nor the creation of a new splice site AND the nucleotide is not highly conserved" - Should get this
    // already from Jannovar? Need PhyloP or other conservation score.

//
// Varsome
//For non-mitochondrial variants, splicing is checked as follows:
//
//    the variant is found more than 2 bases away from the next splice site,
//    it isn't predicted splicing using splice-site prediction.

    //Splice-Site Prediction
    //We use the scSNV database and MaxEntScan for splice-site prediction in rules for rules BP7 and PP3.
    // This score is limited to single-nucleotide variants only. A variant will be predicted splicing if
    // 'ADA Boost Splicing' threshold is greater than 0.958.
    public void assign(AcmgEvidence.Builder acmgEvidence, VariantEvaluation variantEvaluation){
        // Jannovar annotates with variantEffect.splice_Acceptor or splice_donor, so when its synonymous, cant be something else
        VariantEffect variantEffect = variantEvaluation.getVariantEffect();
        if (variantEffect == VariantEffect.SYNONYMOUS_VARIANT && nucleotideNotHighlyConserved(variantEvaluation)) {
            acmgEvidence.add(AcmgCriterion.BP7);
        }
    }

    public boolean nucleotideNotHighlyConserved(VariantEvaluation variantEvaluation){
        // check for phylopScore
        // add to vcf when processing ?
        variantDataService.getPhyloScore(variantEvaluation);
//        int phyloScore = variantEvaluation.getPhyloScore(); //SNPiff

        //VARSOME:
        // Conservation
        //We use PhyloP100Way for conservation tests, this is available for nearly all positions in both genomes, and proves to be a useful indication of whether a variant may be benign or pathogenic.
        //
        //Conservation is used:
        //
        //To exclude highly-conserved variants from rules BP4 and BP7.
        //As a last-resort fallback in rules PP3 and BP4 if no other in-silico predictions are available.
        //Thresholds have been carefully calibrated to maximise accuracy whilst not over-calling, see in-silico predictions
        //
        //Supporting Benign: if the score is less than 3.58
        //Supporting Pathogenic: if the score is greater than 7.52,
        //Moderate Pathogenic: if the score is greater than 9.88.
    }
}
