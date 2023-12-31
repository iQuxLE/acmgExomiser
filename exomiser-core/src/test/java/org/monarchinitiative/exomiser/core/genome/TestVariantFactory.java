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
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.monarchinitiative.exomiser.core.model.SampleGenotype;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.List;

/**
 * Helper class for constructing {@link Variant} objects for tests.
 *
 * The construction of {@link Variant} objects is quite complex but for tests,
 * we would ideally have them for testing our data sets. This class helps us
 * with the construction.
 */
public class TestVariantFactory {

    private static final Logger logger = LoggerFactory.getLogger(TestVariantFactory.class);

    /**
     * Construct a new {@link Variant} object with the given values. n.b. this follows the VCF standard of being 1-based.
     *
     * @param chrom numeric chromosome id
     * @param pos one-based position of the variant
     * @param ref reference string
     * @param alt alt string
     * @param gt the Genotype to use
     * @param readDepth depth the read depth to use
     * @param qual phred-scale quality
     * @return {@link Variant} with the setting
     */
    public static VariantEvaluation buildVariant(int chrom, int pos, String ref, String alt, SampleGenotype gt, int readDepth, double qual) {
        if (alt.contains(",")) {
            throw new IllegalArgumentException("Only single ALT alleles supported.");
        }
        String vcfLine = String.format("%d %d . %s %s %f PASS . GT:RD %s:%d", chrom, pos, ref, alt, qual, gt, readDepth);
        VcfReader vcfReader = TestVcfReader.builder()
                .samples("sample")
                .vcfLines(vcfLine)
                .build();
        VariantFactory variantFactory = TestFactory.buildDefaultVariantFactory(vcfReader);
        return variantFactory.createVariantEvaluations().findFirst().get();
    }

    private VariantContext buildVariantContext(int chrom, int pos, String ref, String alt, SampleGenotype genotype, int readDepth, double qual) {
        Allele refAllele = Allele.create(ref, true);
        Allele altAllele = Allele.create(alt);

        // build Genotype
        GenotypeBuilder genotypeBuilder = buildGenotype(genotype, readDepth, refAllele, altAllele);

        // build VariantContext
        VariantContextBuilder vcBuilder = new VariantContextBuilder();
        vcBuilder.loc("chr" + chrom, pos, pos - 1L + ref.length());
//        vcBuilder.loc(Integer.toString(chrom), pos, pos - 1L + ref.length());
        vcBuilder.alleles(Arrays.asList(refAllele, altAllele));
        vcBuilder.genotypes(genotypeBuilder.make());
        vcBuilder.attribute("RD", readDepth);
        vcBuilder.log10PError(-0.1 * qual);

        return vcBuilder.make();
    }

    private GenotypeBuilder buildGenotype(SampleGenotype gt, int readDepth, Allele refAllele, Allele altAllele) {
        GenotypeBuilder gtBuilder = new GenotypeBuilder("sample");
        gtBuilder.alleles(genotype(refAllele, altAllele, gt));
        gtBuilder.attribute("RD", readDepth);
        return gtBuilder;
    }

    private List<Allele> genotype(Allele refAllele, Allele altAllele, SampleGenotype genotype) {
        if (genotype.isHomRef()) {
            return List.of(altAllele, altAllele);
        } else if (genotype.isHomRef()) {
            return List.of(refAllele, refAllele);
        } else if (genotype.isHet()) {
            return List.of(refAllele, altAllele);
        }
        return List.of();
    }

}
