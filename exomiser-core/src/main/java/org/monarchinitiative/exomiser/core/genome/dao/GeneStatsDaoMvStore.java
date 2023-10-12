package org.monarchinitiative.exomiser.core.genome.dao;

import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStatistics;

public class GeneStatsDaoMvStore implements GeneStatsDao{

    private final MVMap<String, ClinVarGeneStatistics> geneStatsMap;

    public GeneStatsDaoMvStore(MVStore mvStore) {
        this.geneStatsMap = mvStore.openMap("ClinVarGeneStatistics"); // add method in MvStoreUtil for more modular open and build
    }

    @Override
    public ClinVarGeneStatistics getGeneStats(VariantEvaluation variantEvaluation) {
        return geneStatsMap.get(variantEvaluation.getGeneSymbol());
    }
}
