package org.monarchinitiative.exomiser.core.model.pathogenicity;

import static org.junit.jupiter.api.Assertions.*;


import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.junit.jupiter.api.Test;

import java.util.Map;


public class ClinVarGeneStatisticsTest {

    @Test
    public void testClinVarGeneStatistics() {
        ClinVarGeneStatistics.Builder builder = new ClinVarGeneStatistics.Builder();
        builder.addEntry("BRCA1", VariantEffect.MISSENSE_VARIANT, ClinVarData.ClinSig.PATHOGENIC)
                .addEntry("BRCA1", VariantEffect.MISSENSE_VARIANT, ClinVarData.ClinSig.BENIGN)
                .addEntry("BRCA1", VariantEffect.MISSENSE_VARIANT, ClinVarData.ClinSig.LIKELY_BENIGN)
                .addEntry("BRCA1", VariantEffect.START_LOST, ClinVarData.ClinSig.PATHOGENIC);

        ClinVarGeneStatistics clinVarGeneStatistics = builder.build();

        // Retrieve the map
        Map<String, Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>>> map = clinVarGeneStatistics.getMap();

        // Verify counts for BRCA1 and MISSENSE
        assertEquals(1, (int) map.get("BRCA1").get(VariantEffect.MISSENSE_VARIANT).get(ClinVarData.ClinSig.PATHOGENIC));
        assertEquals(1, (int) map.get("BRCA1").get(VariantEffect.MISSENSE_VARIANT).get(ClinVarData.ClinSig.BENIGN));
        assertEquals(1, (int) map.get("BRCA1").get(VariantEffect.MISSENSE_VARIANT).get(ClinVarData.ClinSig.LIKELY_BENIGN));

        // Verify counts for BRCA1 and START_LOST
        assertEquals(1, (int) map.get("BRCA1").get(VariantEffect.START_LOST).get(ClinVarData.ClinSig.PATHOGENIC));
    }
}

