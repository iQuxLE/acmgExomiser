package org.monarchinitiative.exomiser.core.model.pathogenicity;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Sets;
import de.charite.compbio.jannovar.annotation.VariantEffect;
import com.google.common.collect.ImmutableMap;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

//public class ClinVarGeneStatistics {
//    private final Map<String, Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>>> map;
//
//    private ClinVarGeneStatistics(Builder builder) {
//        this.map = ImmutableMap.copyOf(builder.map);
//    }
//
//    public static class Builder {
//
//        private String geneSymbol;
//        private final Map<String, Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>>> map;
//
//        public Builder() {
//            this.map = new HashMap<>();
//        }
//
//        public Builder addEntry(String geneSymbol, VariantEffect variantEffect, ClinVarData.ClinSig clinicalSignificance) {
//            map.computeIfAbsent(geneSymbol, k -> new HashMap<>())
//                    .computeIfAbsent(variantEffect, k -> new HashMap<>())
//                    .merge(clinicalSignificance, 1, Integer::sum);
//            return this;
//        }
//
//        public int getNumberOfEntries() {
//            return map.values().stream()
//                    .mapToInt(m -> m.values().stream()
//                            .mapToInt(Map::size)
//                            .sum())
//                    .sum();
//        }
//
//        public int getOccurrencesOfClinSigForVariantEffect(String geneSymbol, VariantEffect variantEffect, ClinVarData.ClinSig clinSig) {
//            return map.getOrDefault(geneSymbol, Collections.emptyMap())
//                    .getOrDefault(variantEffect, Collections.emptyMap())
//                    .getOrDefault(clinSig, 0);
//        }
//
//        public ClinVarGeneStatistics build() {
//            return new ClinVarGeneStatistics(this);
//        }
//    }
//
//    public Map<String, Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>>> getMap() {
//        return map;
//    }
//}

import java.util.*;

    public class ClinVarGeneStatistics {
        private final List<ClinVarData.ClinSig> clinSigList;
        private final List<VariantEffect> molecularConsequenceList;
        private final String geneSymbol;

        public ClinVarGeneStatistics(String geneSymbol, List<ClinVarData.ClinSig> clinSigList, List<VariantEffect> molecularConsequenceList) {
            this.geneSymbol = geneSymbol;
            this.clinSigList = clinSigList;
            this.molecularConsequenceList = molecularConsequenceList;
        }
    }

     class GeneStatsBuilder {
        private final Map<String, ClinVarGeneStatistics> geneStatsMap = new HashMap<>();
        private List<ClinVarData.ClinSig> clinSigList = new ArrayList<>();
        private List<VariantEffect> molecularConsequenceList = new ArrayList<>();
        private String currentGeneSymbol;

        public void addData(String geneSymbol, ClinVarData.ClinSig clinSig, VariantEffect molecularConsequence) {
            if (currentGeneSymbol == null) {
                currentGeneSymbol = geneSymbol;
            }

            if (!currentGeneSymbol.equals(geneSymbol)) {
                build();
                currentGeneSymbol = geneSymbol;
                clinSigList = new ArrayList<>();
                molecularConsequenceList = new ArrayList<>();
            }

            clinSigList.add(clinSig);
            molecularConsequenceList.add(molecularConsequence);
        }

        public void build() {
            ClinVarGeneStatistics stats = new ClinVarGeneStatistics(currentGeneSymbol, clinSigList, molecularConsequenceList);
            geneStatsMap.put(currentGeneSymbol, stats);
        }

        public Map<String, ClinVarGeneStatistics> getGeneStatsMap() {
            return geneStatsMap;
        }
    }



