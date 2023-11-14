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

    private final double truncatingVariantsArePrimarilyDiseaseCausing = 0.9; //MOCKED NOT REAL

    private final double missenseIsAcommonMechanismOfDisease = 0.66; // MOCKED NOT REAL (NEED STATS FOR THIS) --> LIKELYHOOD

    public BP1PP2Assigner(VariantDataService variantDataService){

        this.variantDataService = variantDataService;
    }

    // PP2: "Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common mechanism of disease."
    public void assign(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation){
        double pathogenicMissenseThreshold = 0.808;
        String geneSymbol = variantEvaluation.getGeneSymbol();
        ClinVarGeneStats clinVarGeneStats = variantDataService.getClinVarGeneStats(geneSymbol);
        if (clinVarGeneStats == null){
            return;
        }
        if (variantEffectIsPrimarilyDiseaseCausing(clinVarGeneStats, variantEvaluation.getVariantEffect()) && calculateRatioPathogenicMissenseOverAllNonVusMissense(clinVarGeneStats) > pathogenicMissenseThreshold) {
            acmgEvidenceBuilder.add(PP2);
        }
    }

    private double calculateRatioPathogenicMissenseOverAllNonVusMissense(ClinVarGeneStats clinVarGeneStatsMap){
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
        double benignMissenseThreshold = 0.569;
        String geneSymbol = variantEvaluation.getGeneSymbol();
        ClinVarGeneStats clinVarGeneStats = variantDataService.getClinVarGeneStats(geneSymbol);
        if (variantEffectIsPrimarilyDiseaseCausing(clinVarGeneStats, variantEvaluation.getVariantEffect()) && calculateRatioBenignMissenseVariantsOverAllNonVusMissense(clinVarGeneStats) > benignMissenseThreshold) {
            acmgEvidenceBuilder.add(BP1);
        }
    }

    private double calculateRatioBenignMissenseVariantsOverAllNonVusMissense(ClinVarGeneStats clinVarGeneStatsMap){
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

    public boolean variantEffectIsPrimarilyDiseaseCausing(ClinVarGeneStats geneStats, VariantEffect variantEffect) {
        Map<ClinVarData.ClinSig, Integer> dataMap = variantEffect == VariantEffect.MISSENSE_VARIANT ? geneStats.getMissenseData() : geneStats.getTruncatingData();
        int totalPathogenic = dataMap.getOrDefault(ClinVarData.ClinSig.PATHOGENIC, 0) + dataMap.getOrDefault(ClinVarData.ClinSig.LIKELY_PATHOGENIC, 0);
        int totalNonVus = dataMap.values().stream().mapToInt(Integer::intValue).sum() - dataMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);

        double ratio = totalNonVus == 0 ? 0 : (double) totalPathogenic / totalNonVus;

        if (variantEffect == VariantEffect.MISSENSE_VARIANT) {
            return ratio > missenseIsAcommonMechanismOfDisease;
        } else {
            return ratio > truncatingVariantsArePrimarilyDiseaseCausing;
        }
    }


            /*BP1 BP1 BP1
        algorithm for determining if truncating variants are known to cause disease
        --> BP1 conversely checks that the ratio of benign missense variants over all non-VUS missense variants is greater than 0.569.
        --> see percentage of truncating variants that are pathogenic in comparison to missense in this gene
        -->> supporting BENIGN, so if we are in a GENE for which primarily truncating Variants are know to cause disease and are pathogenic, finding benign Missense
        variants here is a supporting evidence for a benign variant
         */

    // currently not checking if gene has primarily truncating variants
    // would be a function like primarly truncating ? and do ratio of pathogenic truncating to pathogenic (all other) and then a default threshold of like 5x .

            /*
        algorithm for determining if missense variants are a common mechanism of disease
        --> see percentage of missense variants that are pathogenic in comparison to truncating in this gene
        --> Varsome Germline Classification
        PP2 checks that the ratio of pathogenic missense variants over all non-VUS missense variants is greater than 0.808

         */
    // define Thresholds

    private double calculateMissenseVariantRatio(ClinVarGeneStats clinVarGeneStatsMap, ClinVarData.ClinSig clinSig) {
        if (!clinVarGeneStatsMap.getVariantEffects().contains(VariantEffect.MISSENSE_VARIANT)) {
            return 0;
        }

        Map<ClinVarData.ClinSig, Integer> clinSigMap = clinVarGeneStatsMap.getClinSigMap(VariantEffect.MISSENSE_VARIANT);
        int countOfInterest = 0;
        if (clinSig == ClinVarData.ClinSig.PATHOGENIC || clinSig == ClinVarData.ClinSig.LIKELY_PATHOGENIC) {
            countOfInterest = clinSigMap.getOrDefault(ClinVarData.ClinSig.PATHOGENIC, 0) + clinSigMap.getOrDefault(ClinVarData.ClinSig.LIKELY_PATHOGENIC, 0);
        } else if (clinSig == ClinVarData.ClinSig.BENIGN || clinSig == ClinVarData.ClinSig.LIKELY_BENIGN) {
            countOfInterest = clinSigMap.getOrDefault(ClinVarData.ClinSig.BENIGN, 0) + clinSigMap.getOrDefault(ClinVarData.ClinSig.LIKELY_BENIGN, 0);
        }

        int totalMissense = clinSigMap.values().stream().mapToInt(Integer::intValue).sum();
        int missenseVusCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);
        double nonVusMissenseVariants = totalMissense - missenseVusCount;

        return nonVusMissenseVariants == 0 ? 0 : (double) countOfInterest / nonVusMissenseVariants;
    }


    // unified function
    public void assignBP1orPP2(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation) {
        String geneSymbol = variantEvaluation.getGeneSymbol();
        ClinVarGeneStats clinVarGeneStats = variantDataService.getClinVarGeneStats(geneSymbol);
        if (clinVarGeneStats == null || variantEvaluation.getVariantEffect() != VariantEffect.MISSENSE_VARIANT) {
            return;
        }

        // Calculate the pathogenic ratio for missense variants
        double pathogenicRatio = calculateMissenseVariantRatio(clinVarGeneStats, ClinVarData.ClinSig.PATHOGENIC);
        // Calculate the benign ratio for missense variants
        double benignRatio = calculateMissenseVariantRatio(clinVarGeneStats, ClinVarData.ClinSig.BENIGN);

        // Check for PP2: Pathogenic missense variant ratio
        if (pathogenicRatio > missenseIsAcommonMechanismOfDisease) {
            acmgEvidenceBuilder.add(PP2);
        }

        // Check for BP1: Benign missense variant ratio
        if (benignRatio > truncatingVariantsArePrimarilyDiseaseCausing) {
            acmgEvidenceBuilder.add(BP1);
        }
    }




}
