package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;

import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.BS1;

public class BS1BS2Assigner {

    private final VariantDataService variantDataService;

    public BS1BS2Assigner(VariantDataService variantDataService){

        this.variantDataService = variantDataService;
    }

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

    //BS2 - "Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked
    // (hemizygous) disorder, with full penetrance expected at an early age" - needs het/hom counts from gnomAD
    // (these should be in the newer v14 variant store, but I need to double-check)\
    //Varsome:
    //We first determine the mode of inheritance of the gene, then compares the allele count (see allele frequency for quality checks) to the corresponding threshold:
    //
    //recessive or X-linked genes: allele count greater than 2,
    //dominant genes: allele count greater than 5.

    public void assign(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation, ModeOfInheritance modeOfInheritance) {
        // Check if BA1 or PM2 have been triggered already, and if strong pathogenic evidence exists
//        if (variantEvaluation.isBA1Triggered() || variantEvaluation.isPM2Triggered() || variantEvaluation.hasStrongPathogenicEvidence()) {
//            return;
//        }

        FrequencyData frequencyData = variantEvaluation.getFrequencyData();
        float benignFrequencyThreshold = getBenignFrequencyThreshold(variantDataService, variantEvaluation);

        // For mitochondrial variants, use a specific threshold
        if (modeOfInheritance == ModeOfInheritance.MITOCHONDRIAL) {
            benignFrequencyThreshold = CLINGEN_MITO_THRESHOLD;
        }

        // Check if the frequency is greater than the benign threshold
        if (frequencyData.hasFrequencyOverPercentageValue(benignFrequencyThreshold)) {
            acmgEvidenceBuilder.add(BS1); // Evidence.STRONG
        }
    }

    private float getBenignFrequencyThreshold(VariantDataService variantDataService, VariantEvaluation variantEvaluation) {
        // Get the count of known variants for the gene (use geneSymbol)
        int knownVariantsCount = variantDataService.getKnownVariantsCount(variantEvaluation.getGeneSymbol());

        // Use default threshold if too few known variants
        if (knownVariantsCount < 4) {
            return DEFAULT_THRESHOLD_RARE_DISEASE;
        } else {
            // Here we would normally calculate a gene-specific benign threshold
            // This is a placeholder for the actual implementation which would require additional data
            return BENIGN_THRESHOLD;
        }
    }
}

}
