package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.genome.TestVariantDataService;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;
import org.monarchinitiative.svart.Coordinates;
import org.monarchinitiative.svart.Strand;

import java.util.HashMap;
import java.util.Map;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.is;
import static org.junit.jupiter.api.Assertions.*;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.BP1;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PP2;

class BP1AssignerTest {

    private ClinVarGeneStats parseGeneStats(String string) {
        if (string == null || string.trim().isEmpty() || string.trim().equalsIgnoreCase("null")) {
            return null;
        }
        String[] parts = string.split("-");
        String geneSymbol = parts[0].trim();
        String[] effectData = parts[1].split(";");

        Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>> variantEffectMap = new HashMap<>();

        for (String data : effectData) {
            String[] cvData = data.trim().split(",");
            VariantEffect variantEffect = VariantEffect.valueOf(cvData[0].trim());

            for (int i = 1; i < cvData.length; i += 2) {
                ClinVarData.ClinSig clinSig = ClinVarData.ClinSig.valueOf(cvData[i].trim());
                int count = Integer.parseInt(cvData[i + 1].trim());
                variantEffectMap.computeIfAbsent(variantEffect, k -> new HashMap<>()).put(clinSig, count);
            }
        }
        return variantEffectMap.isEmpty() ? null : new ClinVarGeneStats(geneSymbol, variantEffectMap);
    }

    @ParameterizedTest
    @CsvSource(value = {
            "TestBP1BelowThreshold-STOP_GAINED, BENIGN, 80, PATHOGENIC, 10 | false", // is truncating variant but below Threshold + no missense to check above non vus missense
            "TestBP1AboveThreshold-MISSENSE_VARIANT, PATHOGENIC, 1, BENIGN, 10; STOP_GAINED, PATHOGENIC, 125, BENIGN, 10 | true", // stop gained is truncating and above threshold
            "TestNonsense-FRAMESHIFT_VARIANT, PATHOGENIC, 10, BENIGN, 1 | false",
            "null | false",
    }, delimiter = '|')
    void assignBP1(String testData, boolean expected) {
        String[] parts = testData.split("\\|");
        String geneStatsString = parts[0].trim();
        ClinVarGeneStats expectedMap = geneStatsString.isEmpty() ? null : parseGeneStats(geneStatsString);
        TestVariantDataService serviceBP1Below = TestVariantDataService.builder()
                .put("Test", expectedMap)
                .build();
        var assignerBP1Below = new BP1Assigner(serviceBP1Below);
        VariantEvaluation variantEvaluationBP1Below = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("Test")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBP1Below = new AcmgEvidence.Builder();

        assignerBP1Below.assignBP1(acmgBuilderBP1Below, variantEvaluationBP1Below);
        assertThat(acmgBuilderBP1Below.contains(BP1), is(expected));
        // might want to adjust geneSymbol in Map and geneSymbol from put in TestVariantDataService to be the same, by now could be different
        // or just only allow one assignment
    }
}