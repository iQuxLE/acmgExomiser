package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.genome.VariantAnnotator;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.model.TranscriptAnnotation;
import org.monarchinitiative.exomiser.core.model.VariantAnnotation;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.svart.GenomicVariant;
import org.monarchinitiative.svart.VariantType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.Map;

import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PM5;
import static org.monarchinitiative.exomiser.core.analysis.util.acmg.AcmgCriterion.PS1;

public class PS1PM5Assigner {
    private static final Logger logger = LoggerFactory.getLogger(PS1PM5Assigner.class);
    private final VariantDataService variantDataService;
    private final VariantAnnotator variantAnnotator;


    public PS1PM5Assigner(VariantDataService variantDataService, VariantAnnotator variantAnnotator) {
        this.variantDataService = variantDataService;
        this.variantAnnotator = variantAnnotator;
    }

    /**
     * PS1 "Same amino acid change as a previously established pathogenic variant regardless of nucleotide change"
     * PM5 "Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before"
     */
    public void assignPS1orPM5(AcmgEvidence.Builder acmgEvidenceBuilder, VariantEvaluation variantEvaluation) {
        if (variantEvaluation.getVariantEffect() == VariantEffect.MISSENSE_VARIANT){
            List<TranscriptAnnotation> annotations = variantEvaluation.getTranscriptAnnotations();
            if (annotations == null || annotations.isEmpty()){
                logger.warn("TranscriptAnnotation is empty for variantEvaluation: {}", variantEvaluation);
                return;
            }
            TranscriptAnnotation transcriptAnnotation = annotations.get(0);
            String proteinChangeFromInput = transcriptAnnotation.getHgvsProtein();
            String cdnaChangeFromInput = transcriptAnnotation.getHgvsCdna();
            Map<GenomicVariant, ClinVarData> cvData = variantDataService.findClinVarDataOverlappingGenomicInterval(variantEvaluation.withPadding(2, 2));
            logger.debug("" + cvData);

            for (Map.Entry<GenomicVariant, ClinVarData> entry : cvData.entrySet()) {
                GenomicVariant clinVarVariant = entry.getKey();
                logger.debug("" + clinVarVariant);
                if (GenomicVariant.compare(clinVarVariant, variantEvaluation) == 0){
                    // should be assigning PP5 or BP6
                    continue;
                }
                ClinVarData clinVarData = entry.getValue();

                VariantType variantType = clinVarVariant.variantType();
                if (clinVarData.isPathOrLikelyPath() && clinVarData.starRating() >= 2 && (variantType == VariantType.SNV || variantType == VariantType.MNV)) {

                    List<VariantAnnotation> annotatedVariantList = variantAnnotator.annotate(clinVarVariant);
                    if (!annotatedVariantList.isEmpty()) {

                        VariantAnnotation variantAnnotation = annotatedVariantList.get(0);
                        logger.debug("AnnotatedClinVarData: " + variantAnnotation);

                        if (variantAnnotation.hasTranscriptAnnotations()) {
                            VariantEffect variantEffectFromVariantStore = variantAnnotation.getVariantEffect();
                            if (variantEffectFromVariantStore != VariantEffect.MISSENSE_VARIANT){
                                return;
                            }
                            TranscriptAnnotation transcriptAnnotationFromAnnotatedClinVarData = variantAnnotation.getTranscriptAnnotations().get(0);
                            String proteinChangeFromClinVar = transcriptAnnotationFromAnnotatedClinVarData.getHgvsProtein();
                            String cdnaChangeFromClinVar = transcriptAnnotationFromAnnotatedClinVarData.getHgvsCdna();

                            if (proteinChangeFromInput.equals(proteinChangeFromClinVar) && !cdnaChangeFromInput.equals(cdnaChangeFromClinVar)) {
                                acmgEvidenceBuilder.add(PS1);
                            }
                            if (!proteinChangeFromInput.equals(proteinChangeFromClinVar)) {
                                acmgEvidenceBuilder.add(PM5);
                            }
                        }
                    }
                }
            }
        }
    }
}
