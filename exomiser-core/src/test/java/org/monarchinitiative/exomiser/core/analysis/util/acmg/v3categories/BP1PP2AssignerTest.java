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
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PP2;

class BP1PP2AssignerTest {


    @Test
    void assign() {
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
}