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

public class BP1Assigner {

    private final VariantDataService variantDataService;

    private final double truncatingVariantsArePrimarilyDiseaseCausing = 0.9; //MOCKED

    public BP1Assigner(VariantDataService variantDataService){

        this.variantDataService = variantDataService;
    }

            /*
        TODO: Algorithm for determining if truncating variants are known to cause disease
        --> "BP1 conversely checks that the ratio of benign missense variants over all non-VUS missense variants is greater than 0.569."
        --> see percentage of truncating variants that are pathogenic in comparison to missense in this gene
        -->> supporting BENIGN, so if we are in a GENE for which primarily truncating Variants are know to cause disease and are pathogenic, finding benign Missense
        variants here is a supporting evidence for a benign variant
         */


    //BP1 "Missense variant in a gene for which primarily truncating variants are known to cause disease"
    public void assignBP1(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation) {
        double benignMissenseThreshold = 0.569;
        VariantEffect variantEffect = variantEvaluation.getVariantEffect();
        String geneSymbol = variantEvaluation.getGeneSymbol();
        if (variantEffect == null || geneSymbol == null ){
            return;
        }
        ClinVarGeneStats clinVarGeneStats = variantDataService.getClinVarGeneStats(geneSymbol);
        if(clinVarGeneStats == null){
            return;
        }
        if (truncatingVariantsArePrimarilyDiseaseCausing(clinVarGeneStats) && calculateRatioBenignMissenseVariantsOverAllNonVusMissense(clinVarGeneStats, variantEffect) > benignMissenseThreshold) {
            acmgEvidenceBuilder.add(BP1);
        }
    }

    private double calculateRatioBenignMissenseVariantsOverAllNonVusMissense(ClinVarGeneStats clinVarGeneStatsMap, VariantEffect variantEffect){
        if (variantEffect != VariantEffect.MISSENSE_VARIANT){
            System.out.println("missense = null");
            return 0;
        }
        Map<ClinVarData.ClinSig, Integer> clinSigMap = clinVarGeneStatsMap.getClinSigMap(VariantEffect.MISSENSE_VARIANT);
        int benignCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.BENIGN, 0);
        int likelyBenignCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.LIKELY_BENIGN, 0);
        int totalMissense = clinSigMap.values().stream().mapToInt(Integer::intValue).sum();
        int vusCountMissense = clinSigMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);
        double nonVusMissenseVariants = totalMissense - vusCountMissense;
        System.out.println((benignCount + likelyBenignCount) / nonVusMissenseVariants);
        return (benignCount + likelyBenignCount) / nonVusMissenseVariants;
    }

    public boolean truncatingVariantsArePrimarilyDiseaseCausing(ClinVarGeneStats geneStats) {
        System.out.println(geneStats.getVariantEffects());
        Map<ClinVarData.ClinSig, Integer> dataMap = geneStats.getTruncatingData();
        System.out.println(dataMap);
        int totalPathogenic = dataMap.getOrDefault(ClinVarData.ClinSig.PATHOGENIC, 0) + dataMap.getOrDefault(ClinVarData.ClinSig.LIKELY_PATHOGENIC, 0);
        int totalNonVus = dataMap.values().stream().mapToInt(Integer::intValue).sum() - dataMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);
        double ratio = totalNonVus == 0 ? 0 : (double) totalPathogenic / totalNonVus;
        System.out.println(ratio);
        System.out.println(ratio > truncatingVariantsArePrimarilyDiseaseCausing);
        return ratio > truncatingVariantsArePrimarilyDiseaseCausing;

    }
}
