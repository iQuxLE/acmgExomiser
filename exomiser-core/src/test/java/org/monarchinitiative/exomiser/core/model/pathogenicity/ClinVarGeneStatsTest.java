package org.monarchinitiative.exomiser.core.model.pathogenicity;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.junit.jupiter.api.Test;

import java.util.*;

import static org.hamcrest.CoreMatchers.*;
import static org.hamcrest.MatcherAssert.assertThat;

public class ClinVarGeneStatsTest {

//    @Test
//    public void testEmptyBuilder() {
//        ClinVarGeneStats stats = ClinVarGeneStats.builder().geneSymbol("BRCA1").build();
//
//        assertThat(stats.getInnerMap(), equalTo(Collections.emptyMap()));
//    }

//    @Test
//    public void testBuilderWithValues() {
//        String geneSymbol = "BRCA1";
//        VariantEffect variantEffect = VariantEffect.MISSENSE_VARIANT;
//        ClinVarData.ClinSig primaryInterpretation = ClinVarData.ClinSig.PATHOGENIC;
//        ClinVarData clinVarData = ClinVarData.builder()
//                .geneSymbol(geneSymbol)
//                .variantEffect(variantEffect)
//                .primaryInterpretation(primaryInterpretation)
//                .build();
//        ClinVarGeneStats stats = ClinVarGeneStats.builder()
//                .geneSymbol(geneSymbol)
//                .addClinVarData(clinVarData)
//                .build();
//
//        Map<ClinVarData.ClinSig, Integer> expectedMap = new HashMap<>();
//        expectedMap.put(primaryInterpretation, 1);
//        assertThat(stats.getInnerMap(), equalTo(Map.of(variantEffect, expectedMap)));
//    }
//
//    @Test
//    public void testBuilderWithMultipleClinVarData() {
//        String geneSymbol = "BRCA1";
//        VariantEffect variantEffect1 = VariantEffect.MISSENSE_VARIANT;
//        ClinVarData.ClinSig pathogenic = ClinVarData.ClinSig.PATHOGENIC;
//        ClinVarData.ClinSig likelyPathogenic = ClinVarData.ClinSig.LIKELY_PATHOGENIC;
//
//        ClinVarData clinVarData1 = ClinVarData.builder()
//                .geneSymbol(geneSymbol)
//                .variantEffect(variantEffect1)
//                .primaryInterpretation(pathogenic)
//                .build();
//        ClinVarData clinVarData2 = ClinVarData.builder()
//                .geneSymbol(geneSymbol)
//                .variantEffect(variantEffect1)
//                .primaryInterpretation(likelyPathogenic)
//                .build();
//        ClinVarData clinVarData3 = ClinVarData.builder()
//                .geneSymbol(geneSymbol)
//                .variantEffect(variantEffect1)
//                .primaryInterpretation(likelyPathogenic)
//                .build();
//
//        ClinVarGeneStats.Builder builder = ClinVarGeneStats.builder()
//                .geneSymbol(geneSymbol);
//
//        List<ClinVarData> clinVarDataList = Arrays.asList(clinVarData1, clinVarData2, clinVarData3);
//        for (ClinVarData clinVarData : clinVarDataList) {
//            builder.addClinVarData(clinVarData);
//        }
//
//        ClinVarGeneStats stats = builder.build();
//
//
//        Map<ClinVarData.ClinSig, Integer> expectedMap1 = new HashMap<>();
//        expectedMap1.put(pathogenic, 1);
//        expectedMap1.put(likelyPathogenic, 2);
//        System.out.println(stats.getInnerMap());
//        assertThat(stats.getInnerMap(), equalTo(Map.of(variantEffect1, expectedMap1)));
//    }

}
