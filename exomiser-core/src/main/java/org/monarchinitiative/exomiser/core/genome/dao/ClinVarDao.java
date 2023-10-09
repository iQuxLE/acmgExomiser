package org.monarchinitiative.exomiser.core.genome.dao;

import org.monarchinitiative.exomiser.core.analysis.util.acmg.GeneStatistics;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.svart.GenomicInterval;
import org.monarchinitiative.svart.GenomicVariant;

import java.util.List;
import java.util.Map;


/**
 * Interface for providing {@link ClinVarData} about a {@link Variant}.
 *
 * @since 14.0.0
 */
public interface ClinVarDao {

    ClinVarData getClinVarData(Variant variant);

//    GeneStatistics getCLinVarDataForGene(Variant variant);

//    List<VariantEvaluation> getClinVarDataForGene(GenomicInterval genomicInterval);
//    ClinVarData getClinVarData(Gene gene);

    Map<GenomicVariant, ClinVarData> findClinVarDataOverlappingGenomicInterval(GenomicInterval genomicInterval);

}
