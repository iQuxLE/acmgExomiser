package org.monarchinitiative.exomiser.core.genome.dao;

import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStatistics;

public interface GeneStatsDao {
    ClinVarGeneStatistics getGeneStats(VariantEvaluation variantEvaluation);

}
