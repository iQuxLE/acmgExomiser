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

    public Map<ClinVarData.ClinSig, Integer> getMissenseData() {
        return getClinSigMap(VariantEffect.MISSENSE_VARIANT);
    }

    public Map<ClinVarData.ClinSig, Integer> getTruncatingData() {
        Set<VariantEffect> truncatingEffects = Set.of(VariantEffect.START_LOST, VariantEffect.STOP_GAINED, VariantEffect.FRAMESHIFT_TRUNCATION,  VariantEffect.SPLICE_ACCEPTOR_VARIANT, VariantEffect.SPLICE_DONOR_VARIANT, VariantEffect.FRAMESHIFT_VARIANT, VariantEffect.DISRUPTIVE_INFRAME_DELETION, VariantEffect.DISRUPTIVE_INFRAME_INSERTION, VariantEffect.INITIATOR_CODON_VARIANT, VariantEffect.INFRAME_INSERTION);
        return getClinSigMapForSet(truncatingEffects);
    }

    public Map<ClinVarData.ClinSig, Integer> getClinSigMapForSet(Set<VariantEffect> variantEffectSet) {
        Map<ClinVarData.ClinSig, Integer> aggregatedData = new HashMap<>();
        for (VariantEffect variantEffect : variantEffectSet) {
            getClinSigMap(variantEffect).forEach((clinSig, count) ->
                    aggregatedData.merge(clinSig, count, Integer::sum)
            );
        }
        return aggregatedData;
    }


    public static ClinVarGeneStats empty(){
        return EMPTY_DATA;
    }

    @Override
    public String toString() {
        return "ClinVarGeneStats{" +
                "geneSymbol='" + geneSymbol + '\'' +
                ", effectMap=" + variantEffectCounts +
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

        public ClinVarGeneStats build() {
            return new ClinVarGeneStats(geneSymbol, innerMap);
        }
    }
}

