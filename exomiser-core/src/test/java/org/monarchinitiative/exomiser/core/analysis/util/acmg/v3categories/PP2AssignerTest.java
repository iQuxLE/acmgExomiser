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
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PP2;

class PP2AssignerTest {

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
            "TestPP2AboveTheshold-MISSENSE_VARIANT, PATHOGENIC, 10, BENIGN, 1; FRAMESHIFT_VARIANT, PATHOGENIC, 10, BENIGN, 5 | true",
            "DontAssignNoMissense-FRAMESHIFT_VARIANT, PATHOGENIC, 10, BENIGN, 5 | false",
            "TestPP2BelowThreshold-MISSENSE_VARIANT, PATHOGENIC, 100, BENIGN, 95 | false",
            "null | false",
    }, delimiter = '|')
    void assignPP2(String testData , boolean expected) {
        String[] geneStats = testData.split("\\|");
        String geneStatsString = geneStats[0].trim();
        ClinVarGeneStats expectedMap = parseGeneStats(geneStatsString);
        TestVariantDataService serviceBelow = TestVariantDataService.builder()
                .put("Test", expectedMap)
                .build();
        var assignerBelow = new PP2Assigner(serviceBelow);
        VariantEvaluation variantEvaluationBelow = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(1), Strand.POSITIVE, Coordinates.oneBased(123, 123), "T", "A")
                .geneSymbol("Test")
                .variantEffect(VariantEffect.MISSENSE_VARIANT)
                .build();
        AcmgEvidence.Builder acmgBuilderBelow = new AcmgEvidence.Builder();

        assignerBelow.assign(acmgBuilderBelow, variantEvaluationBelow);
        assertThat(acmgBuilderBelow.contains(PP2), is(expected));
    }
}