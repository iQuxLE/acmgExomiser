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
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PP2;

class PP2AssignerTest {

    @Test
    void assignWithThresholdBelowEdgeCases() {
        // Test case for pathogenicMissenseThreshold just below the threshold
        ClinVarGeneStats expectedMapBelow = new ClinVarGeneStats("TestBelow",
                Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 10, ClinVarData.ClinSig.BENIGN, 1), VariantEffect.FRAMESHIFT_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 10, ClinVarData.ClinSig.BENIGN, 5)));

        TestVariantDataService serviceBelow = TestVariantDataService.builder()
                .put("TestBelow", expectedMapBelow)
                .build();
        var assignerBelow = new PP2Assigner(serviceBelow);
        VariantEvaluation variantEvaluationBelow = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestBelow")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBelow = new AcmgEvidence.Builder();

        assignerBelow.assign(acmgBuilderBelow, variantEvaluationBelow);
        assertThat(acmgBuilderBelow.contains(PP2), is(true));
    }

    @Test
    void dontAssignNoMissense() {
        // Test case for pathogenicMissenseThreshold just below the threshold
        ClinVarGeneStats expectedMapBelow = new ClinVarGeneStats("TestBelow",
                Map.of(VariantEffect.FRAMESHIFT_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 10, ClinVarData.ClinSig.BENIGN, 5)));

        TestVariantDataService serviceBelow = TestVariantDataService.builder()
                .put("TestBelow", expectedMapBelow)
                .build();
        var assignerBelow = new PP2Assigner(serviceBelow);
        VariantEvaluation variantEvaluationBelow = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestBelow")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBelow = new AcmgEvidence.Builder();

        assignerBelow.assign(acmgBuilderBelow, variantEvaluationBelow);
        assertThat(acmgBuilderBelow.contains(PP2), is(false));
    }

    @Test
    void assignWithThresholdAboveEdgeCases() {
        // Test case for pathogenicMissenseThreshold  below the threshold
        ClinVarGeneStats expectedMapBelow = new ClinVarGeneStats("TestBelow", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 100, ClinVarData.ClinSig.BENIGN, 95)));
        TestVariantDataService serviceBelow = TestVariantDataService.builder()
                .put("TestBelow", expectedMapBelow)
                .build();
        var assignerBelow = new PP2Assigner(serviceBelow);
        VariantEvaluation variantEvaluationBelow = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestBelow")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBelow = new AcmgEvidence.Builder();

        assignerBelow.assign(acmgBuilderBelow, variantEvaluationBelow);
        assertThat(acmgBuilderBelow.contains(PP2), is(false));
    }

    @Test
    void assignWithNullGeneSymbolAndClinVarGeneStats() {
        TestVariantDataService service = TestVariantDataService.builder()
                .put("Test", null)
                .build();
        var assigner = new BP1PP2Assigner(service);
        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("Test")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilder = new AcmgEvidence.Builder();

        assigner.assign(acmgBuilder, variantEvaluation);
        assertThat(acmgBuilder.contains(PP2), is(false));
    }
}