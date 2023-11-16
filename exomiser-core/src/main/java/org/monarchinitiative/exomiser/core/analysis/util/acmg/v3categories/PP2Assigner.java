package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;

import java.util.Map;

import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PP2;

public class PP2Assigner {

    private final VariantDataService variantDataService;

    private final double missenseIsAcommonMechanismOfDisease = 0.66; // MOCKED -> NEED STATS FOR THIS -> LIKELYHOOD

    public PP2Assigner(VariantDataService variantDataService){

        this.variantDataService = variantDataService;
    }
    /*
        TODO:Algorithm for determining if missense variants are a common mechanism of disease
        --> see percentage of missense variants that are pathogenic in comparison to truncating in this gene

        Varsome Germline Classification:
        --> PP2 checks that the ratio of pathogenic missense variants over all non-VUS missense variants is greater than 0.808

    */

    // PP2: "Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common mechanism of disease."
    public void assign(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation){
        double pathogenicMissenseThreshold = 0.808;
        String geneSymbol = variantEvaluation.getGeneSymbol();
        ClinVarGeneStats clinVarGeneStats = variantDataService.getClinVarGeneStats(geneSymbol);
        if (clinVarGeneStats == null){
            System.out.println("geneStats is null");
            return;
        }
        if (missenseIsAcommonMechanismOfDisease(clinVarGeneStats, variantEvaluation.getVariantEffect()) && calculateRatioPathogenicMissenseOverAllNonVusMissense(clinVarGeneStats) > pathogenicMissenseThreshold) {
            acmgEvidenceBuilder.add(PP2);
        }
    }

    private double calculateRatioPathogenicMissenseOverAllNonVusMissense(ClinVarGeneStats clinVarGeneStatsMap){
        Map<ClinVarData.ClinSig, Integer> clinSigMap = clinVarGeneStatsMap.getClinSigMap(VariantEffect.MISSENSE_VARIANT);
        int pathogenicityCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.PATHOGENIC, 0);
        int likelyPathogenicCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.LIKELY_PATHOGENIC, 0);
        int totalMissense = clinSigMap.values().stream().mapToInt(Integer::intValue).sum();
        int missenseVusCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);
        double nonVusMissenseVariants = totalMissense - missenseVusCount;
        System.out.println(""+(pathogenicityCount + likelyPathogenicCount) / nonVusMissenseVariants);
        return (pathogenicityCount + likelyPathogenicCount) / nonVusMissenseVariants;
    }

    public boolean missenseIsAcommonMechanismOfDisease(ClinVarGeneStats geneStats, VariantEffect variantEffect) {
        if (variantEffect != VariantEffect.MISSENSE_VARIANT) {
            return false;
        }
        Map<ClinVarData.ClinSig, Integer> dataMap = geneStats.getMissenseData();
        int totalPathogenic = dataMap.getOrDefault(ClinVarData.ClinSig.PATHOGENIC, 0) + dataMap.getOrDefault(ClinVarData.ClinSig.LIKELY_PATHOGENIC, 0);
        int totalNonVus = dataMap.values().stream().mapToInt(Integer::intValue).sum() - dataMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);
        double ratio = totalNonVus == 0 ? 0 : (double) totalPathogenic / totalNonVus;
        System.out.println(ratio > missenseIsAcommonMechanismOfDisease);
        return ratio > missenseIsAcommonMechanismOfDisease;
    }

}
