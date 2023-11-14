package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;

import java.util.Map;

import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.BP1;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PP2;

public class BP1PP2Assigner {

    private final VariantDataService variantDataService;

    public BP1PP2Assigner(VariantDataService variantDataService){

        this.variantDataService = variantDataService;
    }

    public void assign(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation){
        String geneSymbol = variantEvaluation.getGeneSymbol();
        ClinVarGeneStats clinVarGeneStats = variantDataService.getClinVarGeneStats(geneSymbol);
        if (clinVarGeneStats == null){
            return;
        }
        if (missenseVariantsAreCommonMechanismOfDisease(clinVarGeneStats)) {
            acmgEvidenceBuilder.add(PP2);
        }
    }

    private boolean missenseVariantsAreCommonMechanismOfDisease(ClinVarGeneStats clinVarGeneStats) {
        /*
        algorithm for determining if missense variants are a common mechanism of disease
        --> see percentage of missense variants that are pathogenic in comparison to truncating in this gene
        --> Varsome Germline Classification
        PP2 checks that the ratio of pathogenic missense variants over all non-VUS missense variants is greater than 0.808

         */
        // define Thresholds
        double pathogenicMissenseThreshold = 0.808;
        double ratioPathMissenseVariantsOverNonVusMissense = calculateRatioPathogenicMissenseToTruncatingVariants(clinVarGeneStats);
        return ratioPathMissenseVariantsOverNonVusMissense > pathogenicMissenseThreshold;
    }

    private double calculateRatioPathogenicMissenseToTruncatingVariants(ClinVarGeneStats clinVarGeneStatsMap){
        ClinVarGeneStats geneStats = clinVarGeneStatsMap.getVariantEffects().contains(VariantEffect.MISSENSE_VARIANT) ? clinVarGeneStatsMap : null;
        if ( geneStats == null) return 0;
        // only be inside Missense
        Map<ClinVarData.ClinSig, Integer> clinSigMap = geneStats.getClinSigMap(VariantEffect.MISSENSE_VARIANT);
        int pathogenicityCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.PATHOGENIC, 0);
        int likelyPathogenicCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.LIKELY_PATHOGENIC, 0);
        int totalMissense = clinSigMap.values().stream().mapToInt(Integer::intValue).sum();
        int missenseVusCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);
        double nonVusMissenseVariants = totalMissense - missenseVusCount;
        return (pathogenicityCount + likelyPathogenicCount) / nonVusMissenseVariants;
    }

    /*
BP1 "Missense variant in a gene for which primarily truncating variants are known to cause disease"
 */
    public void assignBP1(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation) {
        String geneSymbol = variantEvaluation.getGeneSymbol();
        ClinVarGeneStats clinVarGeneStats = variantDataService.getClinVarGeneStats(geneSymbol);
        if (truncatingVariantsAreKnownToCauseDisease(clinVarGeneStats) && calculateRatioBenignMissenseVariantsOverAllNonVus(clinVarGeneStats) > 0.569) {
            acmgEvidenceBuilder.add(BP1);
        }
    }

    private boolean truncatingVariantsAreKnownToCauseDisease(ClinVarGeneStats clinVarGeneStats) {
        /*
        algorithm for determining if truncating variants are known to cause disease
        --> BP1 conversely checks that the ratio of benign missense variants over all non-VUS missense variants is greater than 0.569.
        --> see percentage of truncating variants that are pathogenic in comparison to missense in this gene
        -->> supporting BENIGN, so if we are in a GENE for which primarily truncating Variants are know to cause disease and are pathogenic, finding benign Missense
        variants here is a supporting evidence for a benign variant
         */

        // currently not checking if gene has primarily truncating variants
        // would be a function like primarly truncating ? and do ratio of pathogenic truncating to pathogenic (all other) and then a default threshold of like 5x .

        double benignMissenseThreshold = 0.569;
        double ratioBenignMissenseVariantsOverAllNonVusMissense = calculateRatioBenignMissenseVariantsOverAllNonVus(clinVarGeneStats);
        return ratioBenignMissenseVariantsOverAllNonVusMissense > benignMissenseThreshold ;
    }

    private double calculateRatioBenignMissenseVariantsOverAllNonVus(ClinVarGeneStats clinVarGeneStatsMap){
        ClinVarGeneStats geneStats = clinVarGeneStatsMap.getVariantEffects().contains(VariantEffect.MISSENSE_VARIANT) ? clinVarGeneStatsMap : null;
        if ( geneStats == null) return 0;
        Map<ClinVarData.ClinSig, Integer> clinSigMap = geneStats.getClinSigMap(VariantEffect.MISSENSE_VARIANT);
        int benignCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.BENIGN, 0);
        // over all non vus missense bitte beachten
        int likelyBenignCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.LIKELY_BENIGN, 0);
        int totalMissense = clinSigMap.values().stream().mapToInt(Integer::intValue).sum();
        int vusCountMissense = clinSigMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);
        double nonVusMissenseVariants = totalMissense - vusCountMissense;
        return (benignCount + likelyBenignCount) / nonVusMissenseVariants;
    }
}
