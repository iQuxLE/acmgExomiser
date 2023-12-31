/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2021 Queen Mary University of London.
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

package org.monarchinitiative.exomiser.core.genome;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.monarchinitiative.exomiser.core.model.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Utility class for converting an HTSJDK {@link VariantContext} into a {@link Map} of {@link SampleGenotype} indexed by
 * sample name.
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 * @since 10.0.0
 */
public class VariantContextSampleGenotypeConverter {

    private static final Logger logger = LoggerFactory.getLogger(VariantContextSampleGenotypeConverter.class);

    private VariantContextSampleGenotypeConverter() {
    }

    /**
     * Function for converting the genotypes displayed in the {@link VariantContext} into the genotypes for the
     * {@param altAlleleId}.
     * <p>
     * For more complex cases such as a 1/2 genotype, because of the way exomiser evaluates individual alleles these will
     * be reported as -/1 in two separate alleles.
     *
     * @param variantContext
     * @param altAlleleId
     * @return
     */
    public static SampleGenotypes createAlleleSampleGenotypes(VariantContext variantContext, int altAlleleId) {
        Allele refAllele = variantContext.getReference();
        Allele altAllele = variantContext.getAlternateAllele(altAlleleId);
        logger.debug("Making sample genotypes for altAllele: {} {} {} {}", altAlleleId, refAllele, altAllele, variantContext);
        List<SampleData> samples = new ArrayList<>(variantContext.getNSamples());
        for (Genotype genotype : variantContext.getGenotypes()) {
            logger.debug("Building sample genotype for {}", genotype);
            String sampleName = genotype.getSampleName();
            SampleGenotype sampleGenotype = buildSampleGenotype(refAllele, altAllele, genotype);
            CopyNumber copyNumber = buildCopyNumber(genotype);
            if (sampleGenotype.isNoCall() && !copyNumber.isEmpty()) {
                // Canvas hack to
                logger.debug("Building sample genotype from CN {}", genotype);
                int CN = copyNumber.copies();
                // BUT chrX 140205371 Canvas:REF:chrX:140205371-140208082 N . 7.53 PASS DQ=31.0549859513643;dq20;END=140208082;CIPOS=-221,221;CIEND=-291,221 RC:BC:CN 56:5:1
                // MCC is a Canvas-specific major chromosome count - WT is 1 (1 maternal, 1 paternal)
                int MCC = parseIntAttribute(genotype, "MCC", -1);
                sampleGenotype = (CN == 0 || CN == MCC) ? SampleGenotype.homAlt() : SampleGenotype.het();
            }
            logger.debug("Variant [{} {}] sample {} {} has genotype {}", variantContext.getReference(), altAllele, genotype, genotype.getType(), sampleGenotype);
            SampleData sampleData = SampleData.of(sampleName, sampleGenotype, copyNumber);
            samples.add(sampleData);
        }
        return SampleGenotypes.of(samples);
    }

    private static CopyNumber buildCopyNumber(Genotype genotype) {
        if (genotype.hasExtendedAttribute("CN")) {
            int copyNumber = parseIntAttribute(genotype, "CN", -1);
            return CopyNumber.of(copyNumber);
        }
        return CopyNumber.empty();
    }

    private static int parseIntAttribute(Genotype genotype, String key, int defaultValue) {
        try {
            return Integer.parseInt((String) genotype.getExtendedAttribute(key));
        } catch (NumberFormatException e) {
            // swallow
        }
        return defaultValue;
    }

    private static SampleGenotype buildSampleGenotype(Allele refAllele, Allele altAllele, Genotype genotype) {
        if (genotype.hasAnyAttribute("GT")) {
            AlleleCall[] alleleCalls = buildAlleleCalls(refAllele, altAllele, genotype.getAlleles());
            return genotype.isPhased() ? SampleGenotype.phased(alleleCalls) : SampleGenotype.of(alleleCalls);
        }
        return SampleGenotype.empty();
    }

    private static AlleleCall[] buildAlleleCalls(Allele refAllele, Allele altAllele, List<Allele> genotypeAlleles) {
        AlleleCall[] alleleCalls = new AlleleCall[genotypeAlleles.size()];
        logger.trace("Checking genotype {} against {} {}", genotypeAlleles, refAllele, altAllele);
        for (int currentAlleleId = 0; currentAlleleId < genotypeAlleles.size(); currentAlleleId++) {
            Allele currentAllele = genotypeAlleles.get(currentAlleleId);
            alleleCalls[currentAlleleId] = determineAlleleCall(altAllele, currentAllele);
        }
        logger.trace("Assigned: {}", Arrays.asList(alleleCalls));
        return alleleCalls;
    }

    private static AlleleCall determineAlleleCall(Allele altAllele, Allele allele) {
        if (allele.isNoCall()) {
            return AlleleCall.NO_CALL;
        }
        if (allele.isReference()) {
           return AlleleCall.REF;
        }
        // this is the key bit - we've originally split the VCF into single ref/alt alleles so that:
        //     CHR POS REF ALT GT
        // VCF: 1 12345 A T,C 1/2
        // becomes two variant alleles:
        // 1 12345 A T  -/1 (altAlleleId = 0)
        // 1 12345 A C  -/1 (altAlleleId = 1)
        // so a heterozygous non-ref genotype 1/2 should become two hets - one for each alt allele *where the alt allele matches*
        //so the -/1 genotype represents a split genotype where '-' means 'present, but in another allele'
        if (allele.isNonReference() && allele.equals(altAllele)) {
            return AlleleCall.ALT;
        }
        //does this make sense for symbolic?
        return AlleleCall.OTHER_ALT;
    }
}
