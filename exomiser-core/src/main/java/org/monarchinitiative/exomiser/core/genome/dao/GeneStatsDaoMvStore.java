//package org.monarchinitiative.exomiser.core.genome.dao;
//
//import org.h2.mvstore.MVMap;
//import org.h2.mvstore.MVStore;
//import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
//import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
//import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStatistics;
//import org.monarchinitiative.exomiser.core.proto.AlleleProto;
//import org.monarchinitiative.svart.Contig;
//import org.monarchinitiative.svart.GenomicInterval;
//import org.monarchinitiative.svart.GenomicVariant;
//
//import java.util.Collections;
//import java.util.Iterator;
//import java.util.LinkedHashMap;
//import java.util.Map;
//
//public class GeneStatsDaoMvStore implements GeneStatsDao{
//
//    private final ClinVarGeneStatistics geneStats = ClinVarGeneStatistics.newBuilder();
//
//    private final MVMap<String, ClinVarGeneStatistics> geneStatsMap;
//
//    public GeneStatsDaoMvStore(MVStore mvStore) {
//        this.geneStatsMap = mvStore.openMap("ClinVarGeneStatistics"); // openGeneStatsMap/ add method in MvStoreUtil for more modular open and build
//    }
//
//    @Override
//    public ClinVarGeneStatistics getGeneStats(VariantEvaluation variantEvaluation) {
//        return geneStatsMap.get(variantEvaluation.getGeneSymbol());
//    }
//}
