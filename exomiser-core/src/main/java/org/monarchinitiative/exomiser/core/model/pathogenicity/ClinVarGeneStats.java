package org.monarchinitiative.exomiser.core.model.pathogenicity;

import de.charite.compbio.jannovar.annotation.VariantEffect;

import java.util.*;

public record ClinVarGeneStats(String geneSymbol, Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>> variantEffectCounts){

    private static final ClinVarGeneStats EMPTY_DATA = new ClinVarGeneStats.Builder().build();

    public ClinVarGeneStats(String geneSymbol, Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>> variantEffectCounts) {
        this.geneSymbol = geneSymbol;
        this.variantEffectCounts = variantEffectCounts;
    }

    public static Builder builder() {
        return new Builder();
    }

    public Set<VariantEffect> getVariantEffects() {
        return variantEffectCounts.keySet();
    }

    public Map<ClinVarData.ClinSig, Integer> getClinSigMap(VariantEffect variantEffect) {
        return variantEffectCounts.getOrDefault(variantEffect, Map.of());
    }


    public static ClinVarGeneStats empty(){
        return EMPTY_DATA;
    }

    @Override
    public String toString() {
        return "ClinVarGeneStats{" +
                "geneSymbol='" + geneSymbol + '\'' +
                ", innerMap=" + variantEffectCounts +
                '}';
    }

    public static class Builder {
        private String geneSymbol = "";
        private Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>> innerMap = new HashMap<>();

        public Builder geneSymbol(String geneSymbol) {
            Objects.requireNonNull(geneSymbol);
            this.geneSymbol = geneSymbol;
            return this;
        }

        public Builder addClinVarData(ClinVarData clinVarData) {
            Objects.requireNonNull(clinVarData);
            if (!clinVarData.getGeneSymbol().equals(this.geneSymbol)) {
                throw new IllegalArgumentException("Gene symbol mismatch");
            }
            VariantEffect variantEffect = clinVarData.getVariantEffect();
            ClinVarData.ClinSig primaryInterpretation = clinVarData.getPrimaryInterpretation();
            Map<ClinVarData.ClinSig, Integer> clinSigMap = innerMap.computeIfAbsent(variantEffect, k -> new HashMap<>());
            clinSigMap.merge(primaryInterpretation, 1, Integer::sum);
            return this;
        }

        public Builder addEntry(VariantEffect variantEffect, ClinVarData.ClinSig primaryInterpretation) {
            innerMap.computeIfAbsent(variantEffect, k -> new HashMap<>())
                    .merge(primaryInterpretation, 1, Integer::sum);
            return this;
        }

        public ClinVarGeneStats build() {
            return new ClinVarGeneStats(geneSymbol, innerMap);
        }
    }
}

