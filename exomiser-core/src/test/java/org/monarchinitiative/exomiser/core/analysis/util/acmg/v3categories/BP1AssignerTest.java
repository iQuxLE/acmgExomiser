package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.genome.TestVariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;
import org.monarchinitiative.svart.Coordinates;
import org.monarchinitiative.svart.Strand;

import java.util.Map;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.is;
import static org.junit.jupiter.api.Assertions.*;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.BP1;

class BP1AssignerTest {

    @Test
    void assignBP1() {
    }

    @Test
    void assignWithBP1ThresholdBelowEdgeCases() {
        ClinVarGeneStats expectedMapBP1Below = new ClinVarGeneStats("TestBP1Below", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.BENIGN, 80, ClinVarData.ClinSig.PATHOGENIC, 10)));
        TestVariantDataService serviceBP1Below = TestVariantDataService.builder()
                .put("TestBP1Below", expectedMapBP1Below)
                .build();
        var assignerBP1Below = new BP1Assigner(serviceBP1Below);
        VariantEvaluation variantEvaluationBP1Below = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestBP1Below")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBP1Below = new AcmgEvidence.Builder();

        assignerBP1Below.assignBP1(acmgBuilderBP1Below, variantEvaluationBP1Below);
        assertThat(acmgBuilderBP1Below.contains(BP1), is(false));
    }

    @Test
    void assignWithBP1ThresholdAboveEdgeCases() {
        ClinVarGeneStats expectedMapBP1Above = new ClinVarGeneStats("TestBP1Above",
                Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 1, ClinVarData.ClinSig.BENIGN, 10), VariantEffect.STOP_GAINED, Map.of(ClinVarData.ClinSig.PATHOGENIC, 125, ClinVarData.ClinSig.BENIGN, 10)));
        TestVariantDataService serviceBP1Above = TestVariantDataService.builder()
                .put("TestBP1Above", expectedMapBP1Above)
                .build();
        var assignerBP1Above = new BP1Assigner(serviceBP1Above);
        VariantEvaluation variantEvaluationBP1Above = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestBP1Above")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBP1Above = new AcmgEvidence.Builder();

        assignerBP1Above.assignBP1(acmgBuilderBP1Above, variantEvaluationBP1Above);
        assertThat(acmgBuilderBP1Above.contains(BP1), is(true));
    }
}