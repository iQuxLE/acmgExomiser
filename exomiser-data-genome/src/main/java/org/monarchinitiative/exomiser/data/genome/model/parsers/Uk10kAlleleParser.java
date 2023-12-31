/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2019 Queen Mary University of London.
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

package org.monarchinitiative.exomiser.data.genome.model.parsers;

import org.monarchinitiative.exomiser.data.genome.model.Allele;
import org.monarchinitiative.exomiser.data.genome.model.AlleleProperty;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class Uk10kAlleleParser extends VcfAlleleParser {

    private static final Logger logger = LoggerFactory.getLogger(Uk10kAlleleParser.class);

    /**
     * Parses the AF value form the INFO line and adds it to the {@link Allele}. The AF is pre-calculated as the AC/AN.
     * AC and AN are the aggregated values of the ALSPAC and TWINSUK cohorts. These are non-rare disease cohorts so should
     * be suitable for rare disease analysis.
     * <p>
     * http://www.uk10k.org/studies/cohorts.html
     * <p>
     * ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
     * ##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false under the trained gaussian mixture model (GATK)">
     * ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in called genotypes">
     * ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
     * ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency in called genotypes">
     * ##INFO=<ID=AC_TWINSUK,Number=A,Type=Integer,Description="Allele count in called genotypes in TWINSUK cohort">
     * ##INFO=<ID=AN_TWINSUK,Number=1,Type=Integer,Description="Total number of alleles in called genotypes in TWINSUK cohort">
     * ##INFO=<ID=AF_TWINSUK,Number=A,Type=Float,Description="Allele frequency in called genotypes in TWINSUK cohort">
     * ##INFO=<ID=AC_ALSPAC,Number=A,Type=Integer,Description="Allele count in called genotypes in ALSPAC cohort">
     * ##INFO=<ID=AF_ALSPAC,Number=A,Type=Float,Description="Allele frequency in called genotypes in ASLPAC cohort">
     * ##INFO=<ID=AF_AFR,Number=1,Type=Float,Description="1000 Genomes Phase 1 Allele Frequency in African population (YRI,LWK,ASW)">
     * ##INFO=<ID=AN_ALSPAC,Number=1,Type=Integer,Description="Total number of alleles in called genotypes in ALSPAC cohort">
     * ##INFO=<ID=AF_AMR,Number=1,Type=Float,Description="1000 Genomes Phase 1 Allele Frequency in American population (MXL,CLM,PUR)">
     * ##INFO=<ID=AF_ASN,Number=1,Type=Float,Description="1000 Genomes Phase 1 Allele Frequency in Asian population (CHB,CHS,JPT)">
     * ##INFO=<ID=AF_EUR,Number=1,Type=Float,Description="1000 Genomes Phase 1 Allele Frequency in European population (CEU,TSI,FIN,GBR,IBS)">
     * ##INFO=<ID=AF_MAX,Number=1,Type=Float,Description="1000 Genomes Phase 1 Maximum Allele Frequency">
     * ##INFO=<ID=ESP_MAF,Number=3,Type=Float,Description="Minor allele frequecy in percent for European American, African American and All populations in the NHLBI Exome Sequencing Project (ESP)">
     * ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence of the ALT alleles from Ensembl 75 VEP v75, format transcriptId:geneName:consequence[:codingSeqPosition:proteinPosition:proteinAlleles:proteinPredictions]+...[+gerpScore]">
     * ##INFO=<ID=AC_TWINSUK_NODUP,Number=A,Type=Integer,Description="Allele count in called genotypes in TWINSUK cohort excluding 67 samples where a monozygotic or dyzygotic twin was included in the release">
     * ##INFO=<ID=AN_TWINSUK_NODUP,Number=1,Type=Integer,Description="Total number of alleles in called genotypes in TWINSUK cohort excluding 67 samples where a monozygotic or dyzygotic twin was included in the release">
     * ##INFO=<ID=AF_TWINSUK_NODUP,Number=A,Type=Float,Description="Allele frequency in called genotypes in TWINSUK cohort excluding 67 samples where a monozygotic or dyzygotic twin was included in the release">
     *
     * @param alleles
     * @param info
     * @return
     */
    @Override
    List<Allele> parseInfoField(List<Allele> alleles, String info) {
        List<String> alleleFrequencies = parseAlleleFrequencies(info);

        for (int i = 0; i < alleles.size(); i++) {
            Allele allele = alleles.get(i);
            if (!alleleFrequencies.isEmpty()) {
                String af = alleleFrequencies.get(i);
                if (!af.isEmpty() && !af.equals(".")) {
                    Float freq = 100f * Float.valueOf(af);
                    allele.addValue(AlleleProperty.UK10K, freq);
                }
            }
        }

        return alleles;
    }

    private List<String> parseAlleleFrequencies(String info) {
        //##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency in called genotypes">
        String[] infoFields = info.split(";");
        for (String infoField : infoFields) {
            if (infoField.startsWith("AF=")) {
                return parseFrequencyField(3, infoField);
            }
        }
        return Collections.emptyList();
    }

    private List<String> parseFrequencyField(int keyLength, String infoField) {
        String[] freqs = infoField.substring(keyLength).split(",");
        return Arrays.asList(freqs);
    }
}
