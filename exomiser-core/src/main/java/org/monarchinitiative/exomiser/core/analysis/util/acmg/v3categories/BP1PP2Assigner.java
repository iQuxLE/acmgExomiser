package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;

import java.util.Collections;
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
        double ratioPathMissenseVariantsOverNonVus = calculatePercentagePathogenicVariantsForEffect(clinVarGeneStats, VariantEffect.MISSENSE_VARIANT);
        return ratioPathMissenseVariantsOverNonVus > pathogenicMissenseThreshold;
    }

    private double calculatePercentagePathogenicVariantsForEffect(ClinVarGeneStats clinVarGeneStatsMap, VariantEffect variantEffect){
        ClinVarGeneStats geneStats = clinVarGeneStatsMap.getVariantEffects().contains(variantEffect) ? clinVarGeneStatsMap : null;
        if ( geneStats == null) return 0;
        Map<ClinVarData.ClinSig, Integer> clinSigMap = geneStats.getClinSigMap(variantEffect);
        int pathogenicityCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.PATHOGENIC, 0);

        // PP2 checks that the ratio of pathogenic missense variants over all non-VUS missense variants is greater than 0.808
        int likelyPathogenicCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.LIKELY_PATHOGENIC, 0);
        int total = clinSigMap.values().stream().mapToInt(Integer::intValue).sum();
        int vusCount = clinSigMap.getOrDefault(ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE, 0);

        double everythingBesidesVus = total - vusCount;

        return (pathogenicityCount + likelyPathogenicCount) / everythingBesidesVus;
    }

    /*
BP1 "Missense variant in a gene for which primarily truncating variants are known to cause disease"
 */
    public void assignBP1(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation) {
        String geneSymbol = variantEvaluation.getGeneSymbol();
        ClinVarGeneStats clinVarGeneStats = variantDataService.getClinVarGeneStats(geneSymbol);
        if (truncatingVariantsAreKnownToCauseDisease(clinVarGeneStats)) {
            acmgEvidenceBuilder.add(BP1);
        }
    }

    private boolean truncatingVariantsAreKnownToCauseDisease(ClinVarGeneStats clinVarGeneStats) {
        /*
        algorithm for determining if truncating variants are known to cause disease
        --> see percentage of truncating variants that are pathogenic in comparison to missense in this gene
        BP1 conversely checks that the ratio of benign missense variants over all non-VUS missense variants is greater than 0.569.

         */
        double benignMissenseThreshold = 0.569;
        double missenseVariantsPathogenicity = calculatePercentagePathogenicVariantsForEffect(clinVarGeneStats, VariantEffect.MISSENSE_VARIANT);
        return true;


    }
}
