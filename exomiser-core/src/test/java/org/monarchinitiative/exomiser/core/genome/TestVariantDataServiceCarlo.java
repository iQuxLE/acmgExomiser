package org.monarchinitiative.exomiser.core.genome;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.genome.dao.ClinVarDao;
import org.monarchinitiative.exomiser.core.genome.dao.ClinVarDaoMvStore;
import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
import org.monarchinitiative.exomiser.core.model.AlleleProtoAdaptor;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.frequency.Frequency;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.svart.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.stream.Collectors;

public class TestVariantDataServiceCarlo implements VariantDataService {

    private static final Logger logger = LoggerFactory.getLogger(TestVariantDataService.class);
    private final GenomeAssembly genomeAssembly;
    private final MVStore mvStore;
    private final ClinVarDao clinVarDao;

    private MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap;


    private TestVariantDataServiceCarlo(Builder builder) {
        this.mvStore = builder.mvStore;
        this.genomeAssembly = builder.genomeAssembly;
        this.clinVarDao = new ClinVarDaoMvStore(builder.mvStore, builder.genomeAssembly);
        this.clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);
    }

    @Override
    public boolean variantIsWhiteListed(Variant variant) {
        return false;
    }

    @Override
    public FrequencyData getVariantFrequencyData(Variant variant, Set<FrequencySource> frequencySources) {
        return null;
    }

    @Override
    public PathogenicityData getVariantPathogenicityData(Variant variant, Set<PathogenicitySource> pathogenicitySources) {
        return null;
    }

    public ClinVarData getClinVarDataFromGenomicVariant(@Nonnull GenomicVariant variant) {
        AlleleProto.AlleleKey alleleKey = AlleleProtoAdaptor.toAlleleKey(variant);
        AlleleProto.ClinVar clinVar = clinVarMap.get(alleleKey);
        return clinVar == null ? ClinVarData.empty() : AlleleProtoAdaptor.toClinVarData(clinVar);
    }

    public GenomicVariant alleleKeyToGenomicVariant(AlleleProto.AlleleKey alleleKey, Contig contig) {
        return GenomicVariant.builder()
                .variant(contig, Strand.POSITIVE, Coordinates.oneBased(alleleKey.getPosition(), alleleKey.getPosition()), alleleKey.getRef(), alleleKey.getAlt()).build();
    }

    @Override
    public Map<GenomicVariant, ClinVarData> findClinVarDataOverlappingGenomicInterval(GenomicInterval genomicInterval) {
        return clinVarDao.findClinVarDataOverlappingGenomicInterval(genomicInterval);
    }

    public static Builder builder() {
        return new Builder();
    }

    public static class Builder {

        private GenomeAssembly genomeAssembly;
        private MVStore mvStore;
        private MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap;

        public Builder setMVStore(MVStore store) {
            this.mvStore = store;
            if(mvStore == null) {
                throw new IllegalArgumentException("MVStore is not initialized!");
            }
            this.clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);
            return this;
        }

        public Builder putAK(AlleleProto.AlleleKey AlleleKey, AlleleProto.ClinVar ClinVar){
            if (clinVarMap == null) {
                throw new IllegalStateException("MVMap is not initialized. Ensure you have called setMVStore() before using put().");
            }
            this.clinVarMap.put(AlleleKey, ClinVar);
            System.out.println("Map size after put: " + clinVarMap.size());
            System.out.println("Contents: " + clinVarMap.entrySet());
            return this;
        }


        public Builder setGenomeAssembly(GenomeAssembly genomeAssembly){
            this.genomeAssembly = genomeAssembly;
            return this;
        }

        public TestVariantDataServiceCarlo build() {
            System.out.println("Map size before build: " + clinVarMap.size());
            return new TestVariantDataServiceCarlo(this);
        }
    }

}

