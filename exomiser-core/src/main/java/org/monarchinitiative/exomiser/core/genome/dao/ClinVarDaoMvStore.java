package org.monarchinitiative.exomiser.core.genome.dao;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
import org.monarchinitiative.exomiser.core.model.AlleleProtoAdaptor;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.svart.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.annotation.Nonnull;
import javax.swing.text.Position;
import java.util.*;

/**
 * @since 14.0.0
 */
public class ClinVarDaoMvStore implements ClinVarDao {

    private static final Logger logger = LoggerFactory.getLogger(ClinVarDaoMvStore.class);

    private final GenomeAssembly genomeAssembly;
    private final MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> clinVarMap;
    private final MVMap<String, ClinVarGeneStats> clinVarGeneStatsMap;
    public ClinVarDaoMvStore(MVStore mvStore, GenomeAssembly genomeAssembly) {
        clinVarMap = MvStoreUtil.openClinVarMVMap(mvStore);
        clinVarGeneStatsMap = MvStoreUtil.openGeneStatsMVMap(mvStore);
        this.genomeAssembly = genomeAssembly;
    }

    @Override
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

    public GenomicVariant alleleKeyToGenomicVariant(AlleleProto.AlleleKey alleleKey, Contig contig) {
        return GenomicVariant.builder()
                .variant(contig, Strand.POSITIVE, Coordinates.ofAllele(CoordinateSystem.oneBased(), alleleKey.getPosition(), alleleKey.getRef()), alleleKey.getRef(), alleleKey.getAlt()).build();
    }

    private ClinVarData getClinVarData(AlleleProto.AlleleKey alleleKey) {
        AlleleProto.ClinVar clinVar = clinVarMap.get(alleleKey);
        return clinVar == null ? ClinVarData.empty() : AlleleProtoAdaptor.toClinVarData(clinVar);
    }

    @Override
    public Map<GenomicVariant, ClinVarData> findClinVarDataOverlappingGenomicInterval(GenomicInterval genomicInterval) {
        Map<GenomicVariant, ClinVarData> results = new LinkedHashMap<>();
        logger.debug("open ClinVarMap and show content: " + clinVarMap.size());
        Contig contig = genomicInterval.contig();
        if (!genomeAssembly.containsContig(contig)) {
            return Collections.emptyMap();
        }

        int chr = genomicInterval.contigId();
        int start = genomicInterval.start();
        int end = genomicInterval.end();

        // build Allele keys for map and define bounds
        AlleleProto.AlleleKey lowerBound = AlleleProto.AlleleKey.newBuilder()
                .setChr(chr)
                .setPosition(start)
                .build();
        AlleleProto.AlleleKey upperBound = AlleleProto.AlleleKey.newBuilder()
                .setChr(chr)
                .setPosition(end)
                .build();

        AlleleProto.AlleleKey floorKey = clinVarMap.floorKey(lowerBound);
        if (floorKey != null && floorKey.getPosition() < lowerBound.getPosition()) {
            lowerBound = floorKey;
        }

        logger.debug("Iterate through entries");
        Iterator<AlleleProto.AlleleKey> keyIterator = clinVarMap.keyIterator(lowerBound);

        while (keyIterator.hasNext() ) {
            AlleleProto.AlleleKey ak = keyIterator.next();
            // don't process keys out of the initial boundaries
            if (ak.getPosition() >= start && ak.getPosition() <= end ){

                GenomicVariant gvFromAk = alleleKeyToGenomicVariant(ak, contig);
                ClinVarData cvData = getClinVarDataFromGenomicVariant(gvFromAk);
                results.put(gvFromAk, cvData);
            }
            if (ak.getPosition() > upperBound.getPosition()){
                break;
            }
        }
        return results;

    }

    @Override
    public ClinVarGeneStats getClinVarGeneStats(String geneSymbol) {
        return clinVarGeneStatsMap.get(geneSymbol);


//        ClinVarGeneStats.Builder builder = ClinVarGeneStats.builder().geneSymbol(geneSymbol);
//        logger.info("clinVarMap size: " + clinVarMap.size());
//        for (AlleleProto.AlleleKey alleleKey : clinVarMap.keySet()) {
//            ClinVarData clinVarData = getClinVarData(alleleKey);
//            logger.info("clinVarData geneSymbol: " + clinVarData.getGeneSymbol());
//            if (clinVarData.getGeneSymbol().equals(geneSymbol)) {
//                VariantEffect variantEffect = clinVarData.getVariantEffect();
//                ClinVarData.ClinSig clinSig = clinVarData.getPrimaryInterpretation();
//                builder.addEntry(variantEffect, clinSig);
//            }
//        }
//        ClinVarGeneStats clinVarGeneStats = builder.build();
//        logger.info("Finished building ClinVarGeneStats for geneSymbol: " + geneSymbol);
//        return clinVarGeneStats;
//    }


    }

}
