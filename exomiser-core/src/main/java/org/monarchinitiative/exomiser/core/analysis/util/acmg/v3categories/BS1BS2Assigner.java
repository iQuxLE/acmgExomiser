package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;

import java.util.Arrays;
import java.util.List;

import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.BS1;

public class BS1BS2Assigner {

    private final VariantDataService variantDataService;

    private final float CLINGEN_MITO_THRESHOLD = 0.005f;
    private final float DEFAULT_THRESHOLD_RARE_DISEASE = 0.015f;

    private final float BENIGN_CUT_OFF_FREQUENCY;

    public BS1BS2Assigner(VariantDataService variantDataService){
        this.BENIGN_CUT_OFF_FREQUENCY = 1;
//                variantDataService.calculateGeneSpecificBenignCutOffFrequency();
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
        List<Double> benignFrequencies = Arrays.asList(0.0, 0.2);
        List<Double> pathogenicFrequencies = Arrays.asList(0.01, 0.02);
        BenignCutOffFrequency cutoffCalculator = new BenignCutOffFrequency(benignFrequencies, pathogenicFrequencies);
        double cutoff = cutoffCalculator.calculateCutoffFrequency();

//        int knownVariantsCount = variantDataService.getKnownVariantsCount(variantEvaluation.getGeneSymbol());
        int knownVariantsCount = 10;

        // Use default threshold if too few known variants
        if (knownVariantsCount < 4) {
            return DEFAULT_THRESHOLD_RARE_DISEASE;
        } else {
            // Here we would normally calculate a gene-specific benign threshold --> cutOFF freqency Gene specific
            // This is a placeholder for the actual implementation which would require additional data
//            return new BENIGN_CUT_OFF_FREQUENCY; // gene specific
             return (float) cutoff;
        }
    }
}


