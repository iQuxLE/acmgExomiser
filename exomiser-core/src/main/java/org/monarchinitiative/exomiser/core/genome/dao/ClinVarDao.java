package org.monarchinitiative.exomiser.core.genome.dao;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStats;
import org.monarchinitiative.svart.GenomicInterval;
import org.monarchinitiative.svart.GenomicVariant;

import java.util.Map;


/**
 * Interface for providing {@link ClinVarData} about a {@link Variant}.
 *
 * @since 14.0.0
 */
public interface ClinVarDao {

    ClinVarData getClinVarData(Variant variant);

    Map<GenomicVariant, ClinVarData> findClinVarDataOverlappingGenomicInterval(GenomicInterval genomicInterval);

     ClinVarGeneStats getClinVarGeneStats(String geneSymbol);

}
