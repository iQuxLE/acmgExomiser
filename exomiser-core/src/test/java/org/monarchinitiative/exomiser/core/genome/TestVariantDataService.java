/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2019 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.monarchinitiative.exomiser.core.genome;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.genome.dao.ClinVarDao;
import org.monarchinitiative.exomiser.core.genome.dao.ClinVarDaoMvStore;
import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.frequency.Frequency;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.*;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.svart.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Mock of VariantDataService to provide canned responses for variants. Enables
 * testing of the service without requiring any database back-end. This is
 * backed by maps mapping the variants to their respective frequency/pathogenicity/whatever data.
 *
 * @author Jules Jacobsen<jules.jacobsen@sanger.ac.uk>
 */
public class TestVariantDataService implements VariantDataService {
    private static final Logger logger = LoggerFactory.getLogger(TestVariantDataService.class);
    private static final VariantDataService STUB_SERVICE = new VariantDataServiceStub();
    private final GenomeAssembly genomeAssembly;

    private final Set<Variant> expectedWhiteList;
    private final Map<Variant, FrequencyData> expectedFrequencyData;
    private final Map<Variant, PathogenicityData> expectedPathogenicityData;
    private final MVStore mvStore;
    private final ClinVarDao clinVarDao;
    private final MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap;
    private final MVMap<String, ClinVarGeneStats> clinVarGeneStatsMap;



    private TestVariantDataService(Builder builder) {
        this.expectedWhiteList = ImmutableSet.copyOf(builder.expectedWhiteList);
        this.expectedFrequencyData = ImmutableMap.copyOf(builder.expectedFrequencyData);
        this.expectedPathogenicityData = ImmutableMap.copyOf(builder.expectedPathogenicityData);
        this.mvStore = builder.mvStore;
        this.genomeAssembly = builder.genomeAssembly;
        this.clinVarDao = new ClinVarDaoMvStore(builder.mvStore, builder.genomeAssembly);
        this.clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);
        this.clinVarGeneStatsMap = MvStoreUtil.openGeneStatsMVMap(mvStore);

        }

    /**
     * Returns a stub service containing no data.
     *
     * @return a stub service instance
     */
    public static VariantDataService stub() {
        return STUB_SERVICE;
    }

    /**
     * Adds the expected FrequencyData for the given Variant to the
     * VariantDataService.
     *
     * @param variant
     * @param frequencyData
     */
    public void put(Variant variant, FrequencyData frequencyData) {
        expectedFrequencyData.put(variant, frequencyData);
    }

    /**
     * Adds the expected PathogenicityData for the given Variant to the
     * VariantDataService.
     *
     * @param variant
     * @param pathogenicityData
     */
    public void put(Variant variant, PathogenicityData pathogenicityData) {
        expectedPathogenicityData.put(variant, pathogenicityData);
    }

    @Override
    public boolean variantIsWhiteListed(Variant variant) {
        return expectedWhiteList.contains(variant);
    }

    @Override
    public FrequencyData getVariantFrequencyData(Variant variant, Set<FrequencySource> frequencySources) {
        FrequencyData allFrequencyData = expectedFrequencyData.getOrDefault(variant, FrequencyData.empty());

        List<Frequency> wanted = allFrequencyData.getKnownFrequencies()
                .stream()
                .filter(frequency -> frequencySources.contains(frequency.getSource()))
                .collect(Collectors.toList());

        return FrequencyData.of(allFrequencyData.getRsId(), wanted);
    }

    @Override
    public PathogenicityData getVariantPathogenicityData(Variant variant, Set<PathogenicitySource> pathogenicitySources) {
        PathogenicityData pathData = expectedPathogenicityData.getOrDefault(variant, PathogenicityData.empty());

        List<PathogenicityScore> wanted = pathData.getPredictedPathogenicityScores()
                .stream()
                .filter(pathogenicity -> pathogenicitySources.contains(pathogenicity.getSource()))
                .collect(Collectors.toList());

        return PathogenicityData.of(pathData.getClinVarData(), wanted);
    }

    @Override
    public Map<GenomicVariant, ClinVarData> findClinVarDataOverlappingGenomicInterval(GenomicInterval genomicInterval) {
        return clinVarDao.findClinVarDataOverlappingGenomicInterval(genomicInterval);
    }

    @Override
    public ClinVarGeneStats getClinVarGeneStats(String geneSymbol) {
        return clinVarDao.getClinVarGeneStats(geneSymbol);
    }

    public static Builder builder() {
        return new Builder();
    }

    public static class Builder {
        private Set<Variant> expectedWhiteList = new HashSet<>();
        private Map<Variant, FrequencyData> expectedFrequencyData = new HashMap<>();
        private Map<Variant, PathogenicityData> expectedPathogenicityData = new HashMap<>();

        private GenomeAssembly genomeAssembly = GenomeAssembly.HG19;
        private MVStore mvStore = new MVStore.Builder().open();
        private MVMap<String, ClinVarGeneStats> clinVarGeneStatsMap = MvStoreUtil.openGeneStatsMVMap((mvStore));

        private MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);

        public Builder setGenomeAssembly(GenomeAssembly genomeAssembly){
            this.genomeAssembly = genomeAssembly;
            return this;
        }

        public Builder expectedWhiteList(Set<Variant> expectedWhiteList) {
            this.expectedWhiteList = expectedWhiteList;
            return this;
        }

        public Builder put(AlleleProto.AlleleKey AlleleKey, AlleleProto.ClinVar ClinVar){
            this.clinVarMap.put(AlleleKey, ClinVar);
            return this;
        }

        public Builder put(String geneSymbol, ClinVarGeneStats clinVarGeneStats) {
            this.clinVarGeneStatsMap.put(geneSymbol, clinVarGeneStats);
            return this;
        }

        public Builder addToWhitelist(Variant... variants) {
            this.expectedWhiteList.addAll(Arrays.asList(variants));
            return this;
        }

        public Builder expectedFrequencyData(Map<Variant, FrequencyData> expectedFrequencyData) {
            this.expectedFrequencyData = expectedFrequencyData;
            return this;
        }

        public Builder put(Variant variant, FrequencyData expectedFrequencyData) {
            this.expectedFrequencyData.put(variant, expectedFrequencyData);
            return this;
        }

        public Builder expectedPathogenicityData(Map<Variant, PathogenicityData> expectedPathogenicityData) {
            this.expectedPathogenicityData = expectedPathogenicityData;
            return this;
        }

        public Builder put(Variant variant, PathogenicityData expectedPathogenicityData) {
            this.expectedPathogenicityData.put(variant, expectedPathogenicityData);
            return this;
        }

        public TestVariantDataService build() {
            return new TestVariantDataService(this);
        }
    }


    private static class VariantDataServiceStub implements VariantDataService {

        @Override
        public Map<GenomicVariant, ClinVarData> findClinVarDataOverlappingGenomicInterval(GenomicInterval genomicInterval) {
            return Collections.emptyMap();
        }

        @Override
        public ClinVarGeneStats getClinVarGeneStats(String geneSymbol){
            return ClinVarGeneStats.empty();
        }

        @Override
        public boolean variantIsWhiteListed(Variant variant) {
            return false;
        }

        @Override
        public FrequencyData getVariantFrequencyData(Variant variant, Set<FrequencySource> frequencySources) {
            return FrequencyData.empty();
        }

        @Override
        public PathogenicityData getVariantPathogenicityData(Variant variant, Set<PathogenicitySource> pathogenicitySources) {
            return PathogenicityData.empty();
        }

    }

}
