package org.monarchinitiative.exomiser.core.genome.dao.serialisers;
import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.genome.dao.ClinVarDaoMvStore;
import org.monarchinitiative.exomiser.core.model.AlleleProtoAdaptor;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.exomiser.core.proto.JannovarProto;
import org.monarchinitiative.svart.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.annotation.Nonnull;
import java.util.*;

public class VariantStoreDao {
    private static final Logger logger = LoggerFactory.getLogger(VariantStoreDao.class);

    private final MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap;

    public VariantStoreDao(MVStore mvStore) {
        clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);
    }

    class PeekingIterator<T> implements Iterator<T> {
        private final Iterator<T> iterator;
        private T nextItem;

        public PeekingIterator(Iterator<T> iterator) {
            this.iterator = iterator;
            advance();
        }

        public T peek() {
            return nextItem;
        }

        @Override
        public boolean hasNext() {
            return nextItem != null;
        }

        @Override
        public T next() {
            if (nextItem == null) throw new NoSuchElementException();
            T toReturn = nextItem;
            advance();
            return toReturn;
        }

        private void advance() {
            nextItem = iterator.hasNext() ? iterator.next() : null;
        }
    }

    public ClinVarData getClinVarData(@Nonnull Variant variant) {
        AlleleProto.AlleleKey alleleKey = AlleleProtoAdaptor.toAlleleKey(variant);
        AlleleProto.ClinVar clinVar = clinVarMap.get(alleleKey);
        return clinVar == null ? ClinVarData.empty() : AlleleProtoAdaptor.toClinVarData(clinVar);
    }

    public ClinVarData getClinVarDataFromGenomicVariant(@Nonnull GenomicVariant variant) {
        AlleleProto.AlleleKey alleleKey = AlleleProtoAdaptor.toAlleleKey(variant);
        AlleleProto.ClinVar clinVar = clinVarMap.get(alleleKey);
        return clinVar == null ? ClinVarData.empty() : AlleleProtoAdaptor.toClinVarData(clinVar);
    }

    public GenomicVariant alleleKeyToGenomicVariant(AlleleProto.AlleleKey alleleKey) {

        //problem is the genomeAssemby part
        // if HG 19 /./ and if HG38 ??

        return GenomicVariant.builder()
                .variant(GenomeAssembly.HG19.getContigById(alleleKey.getChr()), Strand.POSITIVE, Coordinates.oneBased(alleleKey.getPosition(), alleleKey.getPosition()), alleleKey.getRef(), alleleKey.getAlt()).build();
//                .variant(genomicInterval.contig(), Strand.POSITIVE, genomicInterval.coordinates(), alleleKey.getRef(), alleleKey.getAlt()).build();
    }

    public Map<GenomicVariant, ClinVarData> findEntriesInRange(GenomicInterval genomicInterval) {
        Map<GenomicVariant, ClinVarData> results = new LinkedHashMap<>();
        logger.info("open ClinVarMap");
//        MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> map = MvStoreUtil.openClinVarMVMap(mvStore);

        var genomicI = genomicInterval.contigId();

        //

        int chr = genomicInterval.contigId();
        Contig contig = genomicInterval.contig();
        int start = genomicInterval.start();
        int end = genomicInterval.end();
        Coordinates coordinates = genomicInterval.coordinates();

        GenomicVariant gv = GenomicVariant.builder()
                .variant(contig, Strand.POSITIVE, coordinates, "","" )
                .build();

        // make this to AlleleKey or map through this dunno

        // build Allele keys for map and define bounds
        AlleleProto.AlleleKey lowerBound = AlleleProto.AlleleKey.newBuilder()
                .setChr(chr)
                .setPosition(start - 2)
                .build();
        AlleleProto.AlleleKey upperBound = AlleleProto.AlleleKey.newBuilder()
                .setChr(chr)
                .setPosition(end + 2)
                .build();

        // add to map
        AlleleProto.AlleleKey floorKey = clinVarMap.floorKey(lowerBound);
        if (floorKey != null && floorKey.getPosition() < lowerBound.getPosition()) {
            lowerBound = floorKey;
        }

        //not really needed if compare fails with upperBound
        AlleleProto.AlleleKey ceilingKey = clinVarMap.ceilingKey(upperBound);
        if (ceilingKey != null && ceilingKey.getPosition() > upperBound.getPosition()) {
            upperBound = ceilingKey;
        }
        logger.info("Iterate through entries");
        PeekingIterator<AlleleProto.AlleleKey> keyIterator = new PeekingIterator<>(clinVarMap.keyIterator(lowerBound));

        while (keyIterator.hasNext()) {
            AlleleProto.AlleleKey ak = keyIterator.next();
            GenomicVariant gvFromAk = alleleKeyToGenomicVariant(ak);
            // don't process keys out of the initial boundaries
            if (ak.getPosition() < (start - 2) || ak.getPosition() > (end + 2)) {
                continue;
            }
//            AlleleProto.ClinVar cv = clinVarMap.get(ak);
            ClinVarData cvData = getClinVarDataFromGenomicVariant(gvFromAk);
            results.put(gvFromAk, cvData);
            if (keyIterator.peek() == null) // if (ak == null/false){ continue }
                continue;
            if (keyIterator.peek().getPosition() > upperBound.getPosition()){
                break;
            }
        }
        return results;

    }
//    @Override
//    public void close () throws Exception {
//        mvStore.close();
//    }
}







