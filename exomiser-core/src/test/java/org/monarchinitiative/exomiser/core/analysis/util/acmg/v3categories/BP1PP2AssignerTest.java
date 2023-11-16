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

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.is;
import static org.junit.jupiter.api.Assertions.*;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.BP1;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PP2;

class BP1PP2AssignerTest {


    @Test
    void assignMapOfMissenseShouldBeTrue() {
        ClinVarGeneStats expectedMap = new ClinVarGeneStats("Test", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 12, ClinVarData.ClinSig.BENIGN, 2)));
        TestVariantDataService service = TestVariantDataService.builder()
                .put("Test", expectedMap)
                .build();
        var assigner = new BP1PP2Assigner(service);
        VariantEvaluation variantEvaluation = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("Test")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilder = new AcmgEvidence.Builder();

        assigner.assign(acmgBuilder, variantEvaluation);
        assertThat(acmgBuilder.contains(PP2), is(true));
    }

    void assignMapOfFrameshiftShouldBeFalse() {
        ClinVarGeneStats expectedMapNonsense = new ClinVarGeneStats("TestNonsense", Map.of(VariantEffect.FRAMESHIFT_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 10, ClinVarData.ClinSig.BENIGN, 1)));
        TestVariantDataService serviceNonsense = TestVariantDataService.builder()
                .put("TestNonsense", expectedMapNonsense)
                .build();
        var assignerNonsense = new BP1PP2Assigner(serviceNonsense);
        VariantEvaluation variantEvaluationNonsense = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestNonsense")
                .variantEffect(VariantEffect.FRAMESHIFT_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderNonsense = new AcmgEvidence.Builder();

        assignerNonsense.assign(acmgBuilderNonsense, variantEvaluationNonsense);
        assertThat(acmgBuilderNonsense.contains(PP2), is(false));
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

    @Test
    void assignWithThresholdBelowEdgeCases() {
        // Test case for pathogenicMissenseThreshold just below the threshold
        ClinVarGeneStats expectedMapBelow = new ClinVarGeneStats("TestBelow", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 10, ClinVarData.ClinSig.BENIGN, 5)));
        TestVariantDataService serviceBelow = TestVariantDataService.builder()
                .put("TestBelow", expectedMapBelow)
                .build();
        var assignerBelow = new BP1PP2Assigner(serviceBelow);
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
        // Test case for pathogenicMissenseThreshold just below the threshold
        ClinVarGeneStats expectedMapBelow = new ClinVarGeneStats("TestBelow", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 100, ClinVarData.ClinSig.BENIGN, 95)));
        TestVariantDataService serviceBelow = TestVariantDataService.builder()
                .put("TestBelow", expectedMapBelow)
                .build();
        var assignerBelow = new BP1PP2Assigner(serviceBelow);
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
    void assignWithBP1ThresholdBelowEdgeCases() {
        ClinVarGeneStats expectedMapBP1Below = new ClinVarGeneStats("TestBP1Below", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.BENIGN, 80, ClinVarData.ClinSig.PATHOGENIC, 10)));
        TestVariantDataService serviceBP1Below = TestVariantDataService.builder()
                .put("TestBP1Below", expectedMapBP1Below)
                .build();
        var assignerBP1Below = new BP1PP2Assigner(serviceBP1Below);
        VariantEvaluation variantEvaluationBP1Below = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestBP1Below")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBP1Below = new AcmgEvidence.Builder();

        assignerBP1Below.assign(acmgBuilderBP1Below, variantEvaluationBP1Below);
        assertThat(acmgBuilderBP1Below.contains(BP1), is(false));
    }

    @Test
    void assignWithBP1ThresholdAboveEdgeCases() {
        ClinVarGeneStats expectedMapBP1Above = new ClinVarGeneStats("TestBP1Above", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.BENIGN, 95, ClinVarData.ClinSig.PATHOGENIC, 5)));
        TestVariantDataService serviceBP1Above = TestVariantDataService.builder()
                .put("TestBP1Above", expectedMapBP1Above)
                .build();
        var assignerBP1Above = new BP1PP2Assigner(serviceBP1Above);
        VariantEvaluation variantEvaluationBP1Above = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestBP1Above")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBP1Above = new AcmgEvidence.Builder();

        assignerBP1Above.assign(acmgBuilderBP1Above, variantEvaluationBP1Above);
        assertThat(acmgBuilderBP1Above.contains(BP1), is(true));
    }

    @Test
    void assignWithPP2ThresholdBelowEdgeCases() {
        ClinVarGeneStats expectedMapPP2Below = new ClinVarGeneStats("TestPP2Below", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 65, ClinVarData.ClinSig.BENIGN, 100)));
        TestVariantDataService servicePP2Below = TestVariantDataService.builder()
                .put("TestPP2Below", expectedMapPP2Below)
                .build();
        var assignerPP2Below = new BP1PP2Assigner(servicePP2Below);
        VariantEvaluation variantEvaluationPP2Below = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestPP2Below")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderPP2Below = new AcmgEvidence.Builder();

        assignerPP2Below.assign(acmgBuilderPP2Below, variantEvaluationPP2Below);
        assertThat(acmgBuilderPP2Below.contains(PP2), is(false));
    }

    @Test
    void assignWithPP2ThresholdAboveEdgeCases() {
        ClinVarGeneStats expectedMapPP2Above = new ClinVarGeneStats("TestPP2Above", Map.of(VariantEffect.MISSENSE_VARIANT, Map.of(ClinVarData.ClinSig.PATHOGENIC, 67, ClinVarData.ClinSig.BENIGN, 100)));
        TestVariantDataService servicePP2Above = TestVariantDataService.builder()
                .put("TestPP2Above", expectedMapPP2Above)
                .build();
        var assignerPP2Above = new BP1PP2Assigner(servicePP2Above);
        VariantEvaluation variantEvaluationPP2Above = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("TestPP2Above")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderPP2Above = new AcmgEvidence.Builder();

        assignerPP2Above.assign(acmgBuilderPP2Above, variantEvaluationPP2Above);
        assertThat(acmgBuilderPP2Above.contains(PP2), is(true));
    }
}